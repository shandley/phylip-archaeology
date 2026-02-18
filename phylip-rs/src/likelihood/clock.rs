//! Clock-constrained maximum likelihood (dnamlk).
//!
//! Implements PHYLIP's `dnamlk` approach: the same pruning algorithm as
//! unconstrained ML, but with the tree constrained to be ultrametric
//! (molecular clock). All tips are equidistant from the root, meaning
//! all lineages evolved at the same rate.
//!
//! The approach:
//! 1. Build an initial ultrametric tree using UPGMA on ML distances.
//! 2. Parameterize by node heights instead of branch lengths.
//! 3. Branch length from node to parent = height(parent) - height(node).
//! 4. Optimize node heights to maximize log-likelihood.
//! 5. Constraint: parent height > child height, all heights >= 0.
//!
//! # Example
//!
//! ```
//! use phylip_rs::likelihood::clock::{clock_ml_search, clock_log_likelihood};
//! use phylip_rs::likelihood::models::Jc69Model;
//! use phylip_rs::tree::{Alignment, Base, Sequence};
//!
//! let alignment = Alignment::new(vec![
//!     Sequence::new("A", vec![Base::A, Base::C, Base::G, Base::T]),
//!     Sequence::new("B", vec![Base::A, Base::C, Base::G, Base::T]),
//!     Sequence::new("C", vec![Base::A, Base::C, Base::A, Base::T]),
//!     Sequence::new("D", vec![Base::A, Base::T, Base::A, Base::C]),
//! ]).unwrap();
//!
//! let model = Jc69Model;
//! let result = clock_ml_search(&alignment, &model, None).unwrap();
//! assert!(result.lnl < 0.0);
//! assert!(result.root_height > 0.0);
//! ```
//!
//! # References
//!
//! - Felsenstein, J. (1981). Evolutionary trees from DNA sequences:
//!   a maximum likelihood approach. *J. Mol. Evol.*, 17, 368-376.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//!   Chapter 20: Likelihood methods for testing hypotheses.

use crate::distance::ml_distances::ml_distance_matrix;
use crate::distance::upgma::upgma;
use crate::likelihood::models::SubstitutionModel;
use crate::likelihood::pruning::log_likelihood;
use crate::tree::types::{Alignment, Tree};

/// Minimum node height for numerical stability.
const MIN_HEIGHT: f64 = 1e-8;
/// Default branch length (minimum) for stability.
const MIN_BRANCH_LENGTH: f64 = 1e-8;
/// Convergence tolerance for height optimization.
const HEIGHT_TOLERANCE: f64 = 1e-6;
/// Maximum rounds of height optimization.
const MAX_OPT_ROUNDS: usize = 50;

/// Result of clock-constrained ML analysis.
#[derive(Debug, Clone)]
pub struct ClockMLResult {
    /// The ultrametric tree with optimized heights.
    pub tree: Tree,
    /// Log-likelihood of the tree.
    pub lnl: f64,
    /// Height of the root (total evolutionary depth).
    pub root_height: f64,
    /// Name of the substitution model used.
    pub model_name: String,
}

/// Compute node heights from an ultrametric tree.
///
/// Height of a leaf is 0. Height of an internal node is the distance
/// from that node to any of its descendant leaves (they should all be
/// equal in an ultrametric tree).
fn compute_heights(tree: &Tree) -> Vec<f64> {
    let n = tree.num_nodes();
    let mut heights = vec![0.0; n];

    // For an ultrametric tree built by UPGMA, we can compute heights
    // by summing branch lengths from leaves upward.
    // Start from leaves (height 0) and go up.
    let postorder = tree.postorder();
    for &node_id in &postorder {
        if tree.nodes[node_id].is_leaf() {
            heights[node_id] = 0.0;
        } else {
            // Height = child height + branch length to child
            // (should be the same for all children in an ultrametric tree)
            let children = &tree.nodes[node_id].children;
            if let Some(&first_child) = children.first() {
                let child_bl = tree.nodes[first_child]
                    .branch_length
                    .unwrap_or(0.1);
                heights[node_id] = heights[first_child] + child_bl;
            }
        }
    }

    heights
}

/// Set branch lengths from heights.
///
/// For each non-root node: branch_length = height(parent) - height(node).
/// Ensures branch lengths are non-negative (at least MIN_BRANCH_LENGTH).
fn set_branch_lengths_from_heights(tree: &mut Tree, heights: &[f64]) {
    let n = tree.num_nodes();
    for i in 0..n {
        if let Some(parent_id) = tree.nodes[i].parent {
            let bl = (heights[parent_id] - heights[i]).max(MIN_BRANCH_LENGTH);
            tree.nodes[i].branch_length = Some(bl);
        }
    }
}

/// Evaluate likelihood of an ultrametric tree.
///
/// Sets branch lengths from the tree's current structure (which should
/// be ultrametric) and evaluates the log-likelihood using the standard
/// pruning algorithm.
pub fn clock_log_likelihood(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
) -> Result<f64, String> {
    log_likelihood(tree, alignment, model).map_err(|e| format!("{}", e))
}

/// Optimize node heights of an ultrametric tree to maximize likelihood.
///
/// Uses golden section search on each internal node's height iteratively:
/// 1. For each internal node, fix all other heights.
/// 2. Find the height that maximizes likelihood within
///    [max_child_height + epsilon, parent_height - epsilon].
/// 3. Update branch lengths accordingly.
/// 4. Iterate until convergence.
///
/// Returns the final log-likelihood.
pub fn optimize_heights(
    tree: &mut Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
) -> Result<f64, String> {
    let mut heights = compute_heights(tree);
    set_branch_lengths_from_heights(tree, &heights);

    let mut prev_lnl = clock_log_likelihood(tree, alignment, model)?;

    for _ in 0..MAX_OPT_ROUNDS {
        // Get internal nodes in postorder (children before parents)
        let postorder = tree.postorder();
        let internal_nodes: Vec<usize> = postorder
            .iter()
            .copied()
            .filter(|&id| !tree.nodes[id].is_leaf())
            .collect();

        for &node_id in &internal_nodes {
            let children = &tree.nodes[node_id].children;
            if children.is_empty() {
                continue;
            }

            // Lower bound: maximum height among children + epsilon
            let max_child_height = children
                .iter()
                .map(|&c| heights[c])
                .fold(0.0f64, f64::max);
            let lower = max_child_height + MIN_HEIGHT;

            // Upper bound: parent height - epsilon, or unbounded if root
            let upper = if let Some(parent_id) = tree.nodes[node_id].parent {
                heights[parent_id] - MIN_HEIGHT
            } else {
                // Root: upper bound is current height * 2 or a reasonable max
                (heights[node_id] * 3.0).max(lower + 1.0)
            };

            if lower >= upper {
                continue; // No room to optimize
            }

            // Golden section search for optimal height
            let golden_ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
            let mut a = lower;
            let mut b = upper;

            let mut x1 = b - golden_ratio * (b - a);
            let mut x2 = a + golden_ratio * (b - a);

            // Evaluate at x1
            heights[node_id] = x1;
            set_branch_lengths_from_heights(tree, &heights);
            let mut f1 = clock_log_likelihood(tree, alignment, model)
                .unwrap_or(f64::NEG_INFINITY);

            // Evaluate at x2
            heights[node_id] = x2;
            set_branch_lengths_from_heights(tree, &heights);
            let mut f2 = clock_log_likelihood(tree, alignment, model)
                .unwrap_or(f64::NEG_INFINITY);

            for _ in 0..30 {
                if (b - a) < HEIGHT_TOLERANCE {
                    break;
                }

                if f1 < f2 {
                    // Maximum is in [x1, b]
                    a = x1;
                    x1 = x2;
                    f1 = f2;
                    x2 = a + golden_ratio * (b - a);
                    heights[node_id] = x2;
                    set_branch_lengths_from_heights(tree, &heights);
                    f2 = clock_log_likelihood(tree, alignment, model)
                        .unwrap_or(f64::NEG_INFINITY);
                } else {
                    // Maximum is in [a, x2]
                    b = x2;
                    x2 = x1;
                    f2 = f1;
                    x1 = b - golden_ratio * (b - a);
                    heights[node_id] = x1;
                    set_branch_lengths_from_heights(tree, &heights);
                    f1 = clock_log_likelihood(tree, alignment, model)
                        .unwrap_or(f64::NEG_INFINITY);
                }
            }

            // Set to the best height found
            let best_height = (a + b) / 2.0;
            heights[node_id] = best_height;
            set_branch_lengths_from_heights(tree, &heights);
        }

        let lnl = clock_log_likelihood(tree, alignment, model)?;

        if (lnl - prev_lnl).abs() < HEIGHT_TOLERANCE {
            return Ok(lnl);
        }
        prev_lnl = lnl;
    }

    Ok(prev_lnl)
}

/// Check whether a tree is ultrametric (all root-to-tip distances are
/// equal within the given tolerance).
pub fn is_ultrametric(tree: &Tree, tolerance: f64) -> bool {
    let root = tree.root;
    let leaves: Vec<usize> = tree
        .nodes
        .iter()
        .filter(|n| n.is_leaf())
        .map(|n| n.id)
        .collect();

    if leaves.is_empty() {
        return true;
    }

    let distances: Vec<f64> = leaves
        .iter()
        .map(|&leaf| distance_to_root(tree, leaf, root))
        .collect();

    let first = distances[0];
    distances.iter().all(|&d| (d - first).abs() < tolerance)
}

/// Compute the sum of branch lengths from `node` up to `ancestor`.
pub fn distance_to_root(tree: &Tree, node: usize, ancestor: usize) -> f64 {
    let mut current = node;
    let mut total = 0.0;
    while current != ancestor {
        total += tree.nodes[current].branch_length.unwrap_or(0.0);
        match tree.nodes[current].parent {
            Some(p) => current = p,
            None => break,
        }
    }
    total
}

/// Clock-constrained ML search.
///
/// 1. Compute ML pairwise distances from the alignment.
/// 2. Build an initial ultrametric tree using UPGMA.
/// 3. Optimize node heights to maximize log-likelihood.
///
/// # Arguments
///
/// * `alignment` - the DNA alignment to analyze
/// * `model` - the substitution model (e.g., JC69 or F84)
/// * `seed` - optional random seed (currently unused; UPGMA is deterministic)
///
/// # Returns
///
/// A `ClockMLResult` with the optimized ultrametric tree and its
/// log-likelihood.
pub fn clock_ml_search(
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    _seed: Option<u64>,
) -> Result<ClockMLResult, String> {
    let ntaxa = alignment.ntaxa();
    if ntaxa < 2 {
        return Err("need at least 2 taxa for clock ML".to_string());
    }

    // Step 1: Compute ML pairwise distances
    let dist_matrix = ml_distance_matrix(alignment, model);

    // Step 2: Build initial ultrametric tree using UPGMA
    let mut tree = upgma(&dist_matrix);

    // Step 3: Optimize node heights
    let lnl = optimize_heights(&mut tree, alignment, model)?;

    // Compute root height from final tree
    let heights = compute_heights(&tree);
    let root_height = heights[tree.root];

    Ok(ClockMLResult {
        tree,
        lnl,
        root_height,
        model_name: model.name().to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::likelihood::models::{F84Model, Jc69Model};
    use crate::likelihood::pruning::{log_likelihood as free_log_likelihood, optimize_branch_lengths};
    use crate::tree::newick::parse_newick;
    use crate::tree::types::{Base, Sequence};

    fn make_alignment(data: &[(&str, &str)]) -> Alignment {
        let seqs: Vec<Sequence> = data
            .iter()
            .map(|(name, seq)| {
                let bases: Vec<Base> = seq
                    .chars()
                    .map(|c| Base::from_char(c).unwrap())
                    .collect();
                Sequence::new(*name, bases)
            })
            .collect();
        Alignment::new(seqs).unwrap()
    }

    // --- Basic clock ML tests ---

    #[test]
    fn test_clock_ml_basic() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGTACGT"),
            ("C", "TGCATGCA"),
        ]);
        let model = Jc69Model;
        let result = clock_ml_search(&aln, &model, None).unwrap();
        assert!(result.lnl < 0.0, "log-likelihood should be negative");
        assert!(
            result.root_height > 0.0,
            "root height should be positive"
        );
        assert_eq!(result.model_name, "JC69");
    }

    #[test]
    fn test_clock_ml_four_taxa() {
        let aln = make_alignment(&[
            ("A", "AAAAACCCCC"),
            ("B", "AAAAACCCCC"),
            ("C", "CCCCCAAAAA"),
            ("D", "CCCCCAAAAA"),
        ]);
        let model = Jc69Model;
        let result = clock_ml_search(&aln, &model, None).unwrap();
        assert_eq!(result.tree.num_leaves(), 4);
        assert!(result.lnl < 0.0);
    }

    #[test]
    fn test_clock_tree_is_ultrametric() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGG"),
        ]);
        let model = Jc69Model;
        let result = clock_ml_search(&aln, &model, None).unwrap();

        assert!(
            is_ultrametric(&result.tree, 1e-4),
            "clock tree should be ultrametric"
        );
    }

    #[test]
    fn test_clock_nonnegative_branch_lengths() {
        let aln = make_alignment(&[
            ("A", "ACGTACGTACGTACGT"),
            ("B", "ACGTACGTACGTACGA"),
            ("C", "TGCATGCATGCATGCA"),
            ("D", "TGCATGCATGCATGCC"),
        ]);
        let model = Jc69Model;
        let result = clock_ml_search(&aln, &model, None).unwrap();

        for node in &result.tree.nodes {
            if let Some(bl) = node.branch_length {
                assert!(
                    bl >= 0.0,
                    "branch length should be non-negative, got {} for node {:?}",
                    bl,
                    node.name
                );
            }
        }
    }

    #[test]
    fn test_clock_lnl_leq_free_lnl() {
        // Clock model is a restricted version of the free model,
        // so clock lnL should be <= free lnL.
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTTAACCGGTA"),
            ("C", "CCAATTGGCCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGGCCAATTGA"),
        ]);
        let model = Jc69Model;

        // Clock ML
        let clock_result = clock_ml_search(&aln, &model, None).unwrap();

        // Free ML (using the same tree topology but optimizing freely)
        let mut free_tree = clock_result.tree.clone();
        let free_lnl = optimize_branch_lengths(&mut free_tree, &aln, &model)
            .map_err(|e| format!("{}", e))
            .unwrap();

        assert!(
            clock_result.lnl <= free_lnl + 1e-4,
            "clock lnL ({}) should be <= free lnL ({})",
            clock_result.lnl,
            free_lnl
        );
    }

    #[test]
    fn test_clock_optimization_improves_over_upgma() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTA"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGA"),
        ]);
        let model = Jc69Model;

        // Build UPGMA tree without optimization
        let dist_matrix = ml_distance_matrix(&aln, &model);
        let upgma_tree = upgma(&dist_matrix);
        let upgma_lnl = free_log_likelihood(&upgma_tree, &aln, &model)
            .map_err(|e| format!("{}", e))
            .unwrap();

        // Clock ML with optimization
        let result = clock_ml_search(&aln, &model, None).unwrap();

        assert!(
            result.lnl >= upgma_lnl - 1e-4,
            "optimized clock lnL ({}) should be >= UPGMA lnL ({})",
            result.lnl,
            upgma_lnl
        );
    }

    #[test]
    fn test_clock_identical_sequences() {
        let aln = make_alignment(&[
            ("A", "ACGTACGTACGT"),
            ("B", "ACGTACGTACGT"),
            ("C", "ACGTACGTACGT"),
        ]);
        let model = Jc69Model;
        let result = clock_ml_search(&aln, &model, None).unwrap();
        assert!(result.lnl < 0.0);
        // Per-site should be reasonable
        let per_site = result.lnl / 12.0;
        assert!(
            per_site > -3.0,
            "identical sequences should have reasonable per-site lnl, got {}",
            per_site
        );
    }

    #[test]
    fn test_clock_log_likelihood_function() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGTACGA"),
            ("C", "TGCATGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.1,C:0.2);").unwrap();
        let model = Jc69Model;
        let lnl = clock_log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    #[test]
    fn test_clock_with_f84() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGG"),
        ]);
        let model = F84Model::equal_freqs(2.0);
        let result = clock_ml_search(&aln, &model, None).unwrap();
        assert!(result.lnl < 0.0);
        assert_eq!(result.model_name, "F84");
    }

    #[test]
    fn test_clock_five_taxa() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAA"),
            ("B", "AACCGGTTAA"),
            ("C", "CCAATTGGCC"),
            ("D", "CCAATTGGCC"),
            ("E", "GGTTAACCGG"),
        ]);
        let model = Jc69Model;
        let result = clock_ml_search(&aln, &model, None).unwrap();
        assert_eq!(result.tree.num_leaves(), 5);
        assert!(result.lnl < 0.0);
    }

    #[test]
    fn test_clock_like_data_small_lrt() {
        // When data is clock-like (equal divergences), the likelihood
        // ratio test statistic should be small.
        // Use identical sequences as the extreme case of clock-like data.
        let aln = make_alignment(&[
            ("A", "ACGTACGTACGTACGT"),
            ("B", "ACGTACGTACGTACGT"),
            ("C", "ACGTACGTACGTACGT"),
            ("D", "ACGTACGTACGTACGT"),
        ]);
        let model = Jc69Model;

        let clock_result = clock_ml_search(&aln, &model, None).unwrap();

        let mut free_tree = clock_result.tree.clone();
        let free_lnl = optimize_branch_lengths(&mut free_tree, &aln, &model)
            .map_err(|e| format!("{}", e))
            .unwrap();

        let lrt_stat = 2.0 * (free_lnl - clock_result.lnl);
        assert!(
            lrt_stat < 10.0,
            "LRT statistic should be small for clock-like data, got {}",
            lrt_stat
        );
    }

    #[test]
    fn test_compute_heights() {
        // Build a simple ultrametric tree: ((A:0.5,B:0.5):0.5,C:1.0)
        let tree = parse_newick("((A:0.5,B:0.5):0.5,C:1.0);").unwrap();
        let heights = compute_heights(&tree);

        // Leaves should have height 0
        for node in &tree.nodes {
            if node.is_leaf() {
                assert!(
                    heights[node.id].abs() < 1e-10,
                    "leaf '{}' should have height 0, got {}",
                    node.name.as_deref().unwrap_or("?"),
                    heights[node.id]
                );
            }
        }

        // Root should have height 1.0
        assert!(
            (heights[tree.root] - 1.0).abs() < 1e-10,
            "root height should be 1.0, got {}",
            heights[tree.root]
        );
    }

    #[test]
    fn test_is_ultrametric_true() {
        let tree = parse_newick("((A:0.5,B:0.5):0.5,C:1.0);").unwrap();
        assert!(is_ultrametric(&tree, 1e-8));
    }

    #[test]
    fn test_is_ultrametric_false() {
        let tree = parse_newick("((A:0.3,B:0.5):0.5,C:1.0);").unwrap();
        assert!(!is_ultrametric(&tree, 0.01));
    }

    #[test]
    fn test_set_branch_lengths_from_heights() {
        let mut tree = parse_newick("((A:0.5,B:0.5):0.5,C:1.0);").unwrap();
        let heights = compute_heights(&tree);

        // Verify heights are sensible
        assert!((heights[tree.root] - 1.0).abs() < 1e-10);

        // Set branch lengths and verify they match
        set_branch_lengths_from_heights(&mut tree, &heights);

        for node in &tree.nodes {
            if let Some(parent_id) = node.parent {
                let expected_bl = (heights[parent_id] - heights[node.id]).max(MIN_BRANCH_LENGTH);
                let actual_bl = node.branch_length.unwrap();
                assert!(
                    (actual_bl - expected_bl).abs() < 1e-10,
                    "branch length mismatch for node {:?}: expected {}, got {}",
                    node.name,
                    expected_bl,
                    actual_bl
                );
            }
        }
    }

    #[test]
    fn test_clock_ml_two_taxa() {
        let aln = make_alignment(&[
            ("A", "ACGTACGTACGT"),
            ("B", "ACGAACGAACGA"),
        ]);
        let model = Jc69Model;
        let result = clock_ml_search(&aln, &model, None).unwrap();
        assert_eq!(result.tree.num_leaves(), 2);
        assert!(result.lnl < 0.0);
    }
}

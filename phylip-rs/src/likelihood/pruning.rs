//! Felsenstein's pruning algorithm (1981).
//!
//! The foundational algorithm for computing the likelihood of sequence
//! data on a phylogenetic tree. This dynamic programming approach
//! traverses the tree from tips to root, computing conditional
//! likelihoods at each internal node.
//!
//! At each leaf node, the conditional likelihood vector is set from
//! the observed data (1.0 for consistent states, 0.0 otherwise).
//! At each internal node, the conditional likelihood for ancestral
//! state `b` is:
//!
//! ```text
//! L_node(b) = PRODUCT over children c of
//!               SUM over states j of P_c(b, j) * L_c(j)
//! ```
//!
//! where `P_c(b, j)` is the transition probability from state `b`
//! to state `j` along the branch to child `c`.
//!
//! The site log-likelihood at the root is:
//!
//! ```text
//! ln L_site = ln( SUM over b of pi_b * L_root(b) )
//! ```
//!
//! # Example
//!
//! ```
//! use phylip_rs::likelihood::models::{SubstitutionModel, Jc69Model};
//! use phylip_rs::likelihood::pruning::log_likelihood;
//! use phylip_rs::tree::{Alignment, Base, Sequence};
//! use phylip_rs::tree::newick::parse_newick;
//!
//! let alignment = Alignment::new(vec![
//!     Sequence::new("A", vec![Base::A, Base::C, Base::G, Base::T]),
//!     Sequence::new("B", vec![Base::A, Base::C, Base::A, Base::T]),
//!     Sequence::new("C", vec![Base::T, Base::G, Base::G, Base::A]),
//! ]).unwrap();
//!
//! let tree = parse_newick("((A:0.1,B:0.1):0.1,C:0.2);").unwrap();
//! let model = Jc69Model;
//! let lnl = log_likelihood(&tree, &alignment, &model).unwrap();
//! assert!(lnl < 0.0); // log-likelihood is always negative
//! ```
//!
//! # References
//!
//! - Felsenstein, J. (1981). Evolutionary trees from DNA sequences:
//!   a maximum likelihood approach. *Journal of Molecular Evolution*,
//!   17, 368-376.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//!   Chapter 16: Models of DNA evolution.

use std::fmt;

use crate::tree::types::{Alignment, Base, Tree};

use super::models::SubstitutionModel;

/// Minimum branch length for numerical stability.
const MIN_BRANCH_LENGTH: f64 = 1e-8;
/// Maximum branch length for optimization.
const MAX_BRANCH_LENGTH: f64 = 10.0;
/// Default branch length when none is specified.
const DEFAULT_BRANCH_LENGTH: f64 = 0.1;
/// Underflow threshold: rescale when max conditional likelihood drops below this.
const UNDERFLOW_THRESHOLD: f64 = 1e-100;

/// Errors in likelihood computation.
#[derive(Debug, Clone)]
pub enum LikelihoodError {
    /// Tree leaf names don't match alignment taxa.
    TaxaMismatch(String),
    /// Tree has no leaves.
    EmptyTree,
    /// All branch lengths are missing.
    NoBranchLengths,
    /// Numerical issue (e.g., all likelihoods are zero at a site).
    NumericalError(String),
}

impl fmt::Display for LikelihoodError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LikelihoodError::TaxaMismatch(name) => {
                write!(f, "leaf '{}' not found in alignment", name)
            }
            LikelihoodError::EmptyTree => write!(f, "tree has no leaves"),
            LikelihoodError::NoBranchLengths => {
                write!(f, "tree has no branch lengths")
            }
            LikelihoodError::NumericalError(msg) => {
                write!(f, "numerical error: {}", msg)
            }
        }
    }
}

impl std::error::Error for LikelihoodError {}

/// Conditional likelihood vector at a node for one site.
/// `[L(A), L(C), L(G), L(T)]`
type CondLike = [f64; 4];

/// Per-node data for the pruning algorithm.
#[derive(Debug, Clone)]
struct NodeLikelihoods {
    /// Conditional likelihood vectors, one per site.
    cond: Vec<CondLike>,
    /// Log-scaling factor per site (for underflow prevention).
    log_scale: Vec<f64>,
}

impl NodeLikelihoods {
    fn new(nsites: usize) -> Self {
        Self {
            cond: vec![[0.0; 4]; nsites],
            log_scale: vec![0.0; nsites],
        }
    }
}

/// Initialize conditional likelihood at a leaf node from observed data.
fn init_leaf(base: Base) -> CondLike {
    let bs = base.base_set(); // [A, C, G, T]
    [
        if bs[0] { 1.0 } else { 0.0 },
        if bs[1] { 1.0 } else { 0.0 },
        if bs[2] { 1.0 } else { 0.0 },
        if bs[3] { 1.0 } else { 0.0 },
    ]
}

/// Get the effective branch length for a node, using default if missing.
fn effective_branch_length(tree: &Tree, node_id: usize) -> f64 {
    tree.nodes[node_id]
        .branch_length
        .unwrap_or(DEFAULT_BRANCH_LENGTH)
        .max(MIN_BRANCH_LENGTH)
}

/// Run the pruning algorithm: compute conditional likelihoods at all nodes.
///
/// Returns per-node likelihood data, with the root containing the
/// final conditional likelihoods needed for the site log-likelihood.
fn pruning_pass(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
) -> Result<Vec<NodeLikelihoods>, LikelihoodError> {
    let nsites = alignment.nsites();
    let nnodes = tree.num_nodes();

    if tree.num_leaves() == 0 {
        return Err(LikelihoodError::EmptyTree);
    }

    // Build name -> alignment index map
    let name_to_idx: std::collections::HashMap<&str, usize> = alignment
        .sequences
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.as_str(), i))
        .collect();

    // Initialize node data
    let mut data = vec![NodeLikelihoods::new(nsites); nnodes];

    // Set leaf conditional likelihoods from observed data
    for node in &tree.nodes {
        if node.is_leaf() {
            if let Some(name) = &node.name {
                let taxon_idx = name_to_idx.get(name.as_str()).ok_or_else(|| {
                    LikelihoodError::TaxaMismatch(name.clone())
                })?;
                for site in 0..nsites {
                    data[node.id].cond[site] =
                        init_leaf(alignment.base(*taxon_idx, site));
                }
            }
        }
    }

    // Postorder traversal: pruning algorithm
    let postorder = tree.postorder();
    for &node_id in &postorder {
        let children = &tree.nodes[node_id].children;
        if children.is_empty() {
            continue; // leaf: already initialized
        }

        // For each child, precompute the transition probability matrix
        let child_data: Vec<(usize, [[f64; 4]; 4])> = children
            .iter()
            .map(|&child_id| {
                let t = effective_branch_length(tree, child_id);
                let p = model.transition_probs(t);
                (child_id, p)
            })
            .collect();

        for site in 0..nsites {
            let mut node_cond = [1.0f64; 4];
            let mut node_log_scale = 0.0f64;

            for &(child_id, ref p_matrix) in &child_data {
                // For each ancestral state b, compute:
                // partial_c[b] = SUM_j P(b->j; t_c) * L_c[j]
                let child_cond = &data[child_id].cond[site];
                node_log_scale += data[child_id].log_scale[site];

                let mut partial = [0.0f64; 4];
                for b in 0..4 {
                    let mut sum = 0.0;
                    for j in 0..4 {
                        sum += p_matrix[b][j] * child_cond[j];
                    }
                    partial[b] = sum;
                }

                // Multiply into the running product
                for b in 0..4 {
                    node_cond[b] *= partial[b];
                }
            }

            // Underflow prevention: rescale if values are too small
            let max_val = node_cond
                .iter()
                .copied()
                .fold(0.0f64, f64::max);

            if max_val > 0.0 && max_val < UNDERFLOW_THRESHOLD {
                for b in 0..4 {
                    node_cond[b] /= max_val;
                }
                node_log_scale += max_val.ln();
            }

            data[node_id].cond[site] = node_cond;
            data[node_id].log_scale[site] = node_log_scale;
        }
    }

    Ok(data)
}

/// Compute per-site log-likelihoods for a tree and alignment.
///
/// Returns a vector of log-likelihood values, one per alignment site.
pub fn site_log_likelihoods(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
) -> Result<Vec<f64>, LikelihoodError> {
    let data = pruning_pass(tree, alignment, model)?;
    let nsites = alignment.nsites();
    let root = tree.root;
    let pi = model.frequencies();

    let mut site_lnl = Vec::with_capacity(nsites);

    for site in 0..nsites {
        let root_cond = &data[root].cond[site];
        let log_scale = data[root].log_scale[site];

        // L_site = SUM_b pi[b] * L_root[b]
        let mut site_like = 0.0;
        for b in 0..4 {
            site_like += pi[b] * root_cond[b];
        }

        if site_like <= 0.0 {
            return Err(LikelihoodError::NumericalError(format!(
                "site {} has zero likelihood",
                site
            )));
        }

        site_lnl.push(site_like.ln() + log_scale);
    }

    Ok(site_lnl)
}

/// Compute the total log-likelihood of the alignment on the tree.
///
/// This is the sum of per-site log-likelihoods.
pub fn log_likelihood(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
) -> Result<f64, LikelihoodError> {
    let site_lnl = site_log_likelihoods(tree, alignment, model)?;
    Ok(site_lnl.iter().sum())
}

/// Optimize a single branch length to maximize the log-likelihood.
///
/// Uses modified Newton-Raphson with numerical derivatives, following
/// PHYLIP's `makenewv()` approach.
///
/// Returns the optimized branch length.
fn optimize_one_branch(
    tree: &mut Tree,
    node_id: usize,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    max_iter: usize,
) -> Result<f64, LikelihoodError> {
    let epsilon = 1e-6;

    let mut t = effective_branch_length(tree, node_id);

    for _ in 0..max_iter {
        let h = (t * 0.001).max(1e-8);

        // Evaluate likelihood at t, t+h, t-h
        tree.nodes[node_id].branch_length = Some(t);
        let lnl_center = log_likelihood(tree, alignment, model)?;

        tree.nodes[node_id].branch_length = Some((t + h).min(MAX_BRANCH_LENGTH));
        let lnl_plus = log_likelihood(tree, alignment, model)?;

        tree.nodes[node_id].branch_length = Some((t - h).max(MIN_BRANCH_LENGTH));
        let lnl_minus = log_likelihood(tree, alignment, model)?;

        // Numerical derivatives
        let slope = (lnl_plus - lnl_minus) / (2.0 * h);
        let curve = (lnl_plus - 2.0 * lnl_center + lnl_minus) / (h * h);

        // Newton step (maximizing, so we want to go uphill)
        let delta = if curve < -1e-12 {
            // Concave: standard Newton step
            -slope / curve
        } else if slope.abs() > 1e-12 {
            // Not concave: step in direction of slope
            slope.signum() * h * 10.0
        } else {
            break; // flat: converged
        };

        let t_new = (t + delta).clamp(MIN_BRANCH_LENGTH, MAX_BRANCH_LENGTH);

        if (t_new - t).abs() < epsilon {
            t = t_new;
            break;
        }
        t = t_new;
    }

    tree.nodes[node_id].branch_length = Some(t);
    Ok(t)
}

/// Optimize all branch lengths in the tree to maximize the log-likelihood.
///
/// Iterates over all branches multiple times until convergence, following
/// the approach in PHYLIP's `smooth()` function.
///
/// Returns the final log-likelihood.
pub fn optimize_branch_lengths(
    tree: &mut Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
) -> Result<f64, LikelihoodError> {
    let max_rounds = 20;
    let max_iter_per_branch = 10;
    let tolerance = 1e-4;

    // Ensure all branches have initial lengths
    for node in &mut tree.nodes {
        if node.parent.is_some() && node.branch_length.is_none() {
            node.branch_length = Some(DEFAULT_BRANCH_LENGTH);
        }
    }

    let mut prev_lnl = log_likelihood(tree, alignment, model)?;

    for _ in 0..max_rounds {
        // Optimize each branch
        let node_ids: Vec<usize> = tree
            .nodes
            .iter()
            .filter(|n| n.parent.is_some())
            .map(|n| n.id)
            .collect();

        for node_id in node_ids {
            optimize_one_branch(tree, node_id, alignment, model, max_iter_per_branch)?;
        }

        let lnl = log_likelihood(tree, alignment, model)?;

        if (lnl - prev_lnl).abs() < tolerance {
            return Ok(lnl);
        }
        prev_lnl = lnl;
    }

    Ok(prev_lnl)
}

/// Result of a likelihood evaluation.
#[derive(Debug, Clone)]
pub struct LikelihoodResult {
    /// The tree (possibly with optimized branch lengths).
    pub tree: Tree,
    /// Total log-likelihood.
    pub lnl: f64,
    /// Per-site log-likelihoods.
    pub site_lnl: Vec<f64>,
}

/// Evaluate the likelihood of a tree given an alignment and model.
///
/// If `optimize` is true, branch lengths will be optimized.
pub fn evaluate(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    optimize: bool,
) -> Result<LikelihoodResult, LikelihoodError> {
    let mut tree = tree.clone();

    if optimize {
        optimize_branch_lengths(&mut tree, alignment, model)?;
    }

    let site_lnl = site_log_likelihoods(&tree, alignment, model)?;
    let lnl = site_lnl.iter().sum();

    Ok(LikelihoodResult {
        tree,
        lnl,
        site_lnl,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::likelihood::models::{F84Model, Jc69Model};
    use crate::tree::newick::parse_newick;
    use crate::tree::types::{Base, BaseFrequencies, Sequence};

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

    // --- Basic likelihood tests ---

    #[test]
    fn test_log_likelihood_is_negative() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "TGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = Jc69Model;
        let lnl = log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0, "log-likelihood should be negative, got {}", lnl);
    }

    #[test]
    fn test_identical_sequences_high_likelihood() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGTACGT"),
            ("C", "ACGTACGT"),
        ]);
        let tree = parse_newick("((A:0.001,B:0.001):0.001,C:0.001);").unwrap();
        let model = Jc69Model;
        let lnl = log_likelihood(&tree, &aln, &model).unwrap();
        // With very short branches and identical sequences, likelihood should be high
        // (close to log(0.25) * nsites for each site having one state)
        assert!(lnl < 0.0);
        // Per-site should be close to ln(0.25) ≈ -1.386
        let per_site = lnl / 8.0;
        assert!(
            per_site > -2.0,
            "per-site lnl should be > -2.0 for identical seqs, got {}",
            per_site
        );
    }

    #[test]
    fn test_divergent_sequences_low_likelihood() {
        let aln = make_alignment(&[
            ("A", "AAAAAAAA"),
            ("B", "CCCCCCCC"),
            ("C", "GGGGGGGG"),
        ]);
        let tree = parse_newick("((A:0.01,B:0.01):0.01,C:0.01);").unwrap();
        let model = Jc69Model;
        let lnl = log_likelihood(&tree, &aln, &model).unwrap();
        // Very divergent on very short branches: low likelihood
        assert!(lnl < -20.0);
    }

    #[test]
    fn test_site_log_likelihoods_sum() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGA"),
            ("C", "TGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let site_lnl = site_log_likelihoods(&tree, &aln, &model).unwrap();
        let total = log_likelihood(&tree, &aln, &model).unwrap();

        let sum: f64 = site_lnl.iter().sum();
        assert!(
            (sum - total).abs() < 1e-10,
            "sum of site lnl ({}) should equal total ({})",
            sum,
            total
        );
    }

    #[test]
    fn test_four_taxa_jc69() {
        let aln = make_alignment(&[
            ("A", "AACCGG"),
            ("B", "AACCGG"),
            ("C", "CCAACC"),
            ("D", "CCAACC"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();
        let model = Jc69Model;
        let lnl = log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    // --- Topology comparison tests ---

    #[test]
    fn test_better_topology_higher_likelihood() {
        // Data has clear ((A,B),(C,D)) signal
        let aln = make_alignment(&[
            ("A", "AAAAACCCCC"),
            ("B", "AAAAACCCCC"),
            ("C", "CCCCCAAAAA"),
            ("D", "CCCCCAAAAA"),
        ]);

        let good_tree =
            parse_newick("((A:0.05,B:0.05):0.2,(C:0.05,D:0.05):0.2);").unwrap();
        let bad_tree =
            parse_newick("((A:0.05,C:0.05):0.2,(B:0.05,D:0.05):0.2);").unwrap();

        let model = Jc69Model;
        let lnl_good = log_likelihood(&good_tree, &aln, &model).unwrap();
        let lnl_bad = log_likelihood(&bad_tree, &aln, &model).unwrap();

        assert!(
            lnl_good > lnl_bad,
            "correct topology ({}) should have higher likelihood than wrong ({})",
            lnl_good,
            lnl_bad
        );
    }

    // --- Branch length effect tests ---

    #[test]
    fn test_branch_length_affects_likelihood() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGTACGT"),
            ("C", "TGCATGCA"),
        ]);
        let model = Jc69Model;

        let short = parse_newick("((A:0.01,B:0.01):0.01,C:0.01);").unwrap();
        let long = parse_newick("((A:1.0,B:1.0):1.0,C:1.0);").unwrap();

        let lnl_short = log_likelihood(&short, &aln, &model).unwrap();
        let lnl_long = log_likelihood(&long, &aln, &model).unwrap();

        // These should be different
        assert!(
            (lnl_short - lnl_long).abs() > 0.01,
            "different branch lengths should give different likelihoods"
        );
    }

    // --- F84 model tests ---

    #[test]
    fn test_f84_model_basic() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGA"),
            ("C", "TGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = F84Model::equal_freqs(2.0);
        let lnl = log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    #[test]
    fn test_f84_unequal_freqs() {
        let aln = make_alignment(&[
            ("A", "AAACCCTTTGGG"),
            ("B", "AAACCCTTTGGG"),
            ("C", "CCCAAAGGGCCC"),
        ]);
        let tree = parse_newick("((A:0.2,B:0.2):0.1,C:0.3);").unwrap();
        let freqs = BaseFrequencies::new(0.3, 0.2, 0.2, 0.3);
        let model = F84Model::new(freqs, 2.0);
        let lnl = log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    #[test]
    fn test_f84_reduces_to_jc69_for_likelihood() {
        // With equal freqs and ttratio=1, F84 should give same lnl as JC69
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGAACGA"),
            ("C", "TGCATGCA"),
        ]);
        let tree = parse_newick("((A:0.15,B:0.15):0.1,C:0.25);").unwrap();

        let jc69 = Jc69Model;
        let f84 = F84Model::equal_freqs(1.0);

        let lnl_jc69 = log_likelihood(&tree, &aln, &jc69).unwrap();
        let lnl_f84 = log_likelihood(&tree, &aln, &f84).unwrap();

        assert!(
            (lnl_jc69 - lnl_f84).abs() < 1e-8,
            "JC69 lnl ({}) should equal F84(ttratio=1) lnl ({})",
            lnl_jc69,
            lnl_f84
        );
    }

    // --- Ambiguity code tests ---

    #[test]
    fn test_ambiguity_codes() {
        // N is fully ambiguous: should act as if all states are possible
        let aln = make_alignment(&[
            ("A", "AN"),
            ("B", "AC"),
            ("C", "TC"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = Jc69Model;
        let lnl = log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    #[test]
    fn test_gap_handling() {
        // Gap treated as fully ambiguous (all states possible)
        let aln = make_alignment(&[
            ("A", "A-"),
            ("B", "AC"),
            ("C", "TC"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = Jc69Model;
        let lnl = log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    // --- Error handling tests ---

    #[test]
    fn test_taxa_mismatch() {
        let aln = make_alignment(&[
            ("X", "ACGT"),
            ("Y", "ACGT"),
            ("Z", "ACGT"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = Jc69Model;
        let result = log_likelihood(&tree, &aln, &model);
        assert!(result.is_err());
    }

    // --- Branch length optimization tests ---

    #[test]
    fn test_optimize_improves_likelihood() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGG"),
        ]);
        // Start with poor branch lengths
        let tree =
            parse_newick("((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);").unwrap();
        let model = Jc69Model;

        let lnl_before = log_likelihood(&tree, &aln, &model).unwrap();

        let mut opt_tree = tree.clone();
        let lnl_after =
            optimize_branch_lengths(&mut opt_tree, &aln, &model).unwrap();

        assert!(
            lnl_after >= lnl_before - 1e-6,
            "optimization should not decrease likelihood: before={}, after={}",
            lnl_before,
            lnl_after
        );
    }

    #[test]
    fn test_optimize_identical_sequences() {
        let aln = make_alignment(&[
            ("A", "ACGTACGTACGT"),
            ("B", "ACGTACGTACGT"),
            ("C", "ACGTACGTACGT"),
        ]);
        let tree = parse_newick("((A:0.5,B:0.5):0.5,C:0.5);").unwrap();
        let model = Jc69Model;

        let mut opt_tree = tree.clone();
        optimize_branch_lengths(&mut opt_tree, &aln, &model).unwrap();

        // With identical sequences, optimized branches should be very short
        for node in &opt_tree.nodes {
            if let Some(bl) = node.branch_length {
                assert!(
                    bl < 0.1,
                    "branch to {:?} should be short for identical seqs, got {}",
                    node.name,
                    bl
                );
            }
        }
    }

    #[test]
    fn test_evaluate_with_optimization() {
        let aln = make_alignment(&[
            ("A", "AACCGGTT"),
            ("B", "AACCGGTT"),
            ("C", "CCAATTGG"),
        ]);
        let tree = parse_newick("((A:0.5,B:0.5):0.5,C:0.5);").unwrap();
        let model = Jc69Model;

        let result_no_opt = evaluate(&tree, &aln, &model, false).unwrap();
        let result_opt = evaluate(&tree, &aln, &model, true).unwrap();

        assert!(
            result_opt.lnl >= result_no_opt.lnl - 1e-6,
            "optimized lnl ({}) should be >= unoptimized ({})",
            result_opt.lnl,
            result_no_opt.lnl
        );
        assert_eq!(result_opt.site_lnl.len(), 8);
    }

    // --- Default branch length tests ---

    #[test]
    fn test_tree_without_branch_lengths() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "TGCA"),
        ]);
        // Tree without branch lengths: should use defaults
        let tree = parse_newick("((A,B),C);").unwrap();
        let model = Jc69Model;
        let lnl = log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    // --- Multifurcating tree tests ---

    #[test]
    fn test_multifurcating_tree() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "TGCA"),
            ("D", "TGCA"),
        ]);
        let tree = parse_newick("(A:0.1,B:0.1,C:0.2,D:0.2);").unwrap();
        let model = Jc69Model;
        let lnl = log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    // --- Numerical property tests ---

    #[test]
    fn test_likelihood_monotone_with_divergence() {
        // More divergent data on the same tree should have lower likelihood
        // (with fixed branch lengths)
        let aln_close = make_alignment(&[
            ("A", "ACGTACGTACGTACGT"),
            ("B", "ACGTACGTACGTACGA"), // 1 difference
            ("C", "ACGTACGTACGTACGT"),
        ]);
        let aln_far = make_alignment(&[
            ("A", "ACGTACGTACGTACGT"),
            ("B", "TGCATGCATGCATGCA"), // all different
            ("C", "ACGTACGTACGTACGT"),
        ]);

        let tree = parse_newick("((A:0.01,B:0.01):0.01,C:0.01);").unwrap();
        let model = Jc69Model;

        let lnl_close = log_likelihood(&tree, &aln_close, &model).unwrap();
        let lnl_far = log_likelihood(&tree, &aln_far, &model).unwrap();

        assert!(
            lnl_close > lnl_far,
            "close sequences ({}) should have higher lnl than far ({}) on short tree",
            lnl_close,
            lnl_far
        );
    }

    #[test]
    fn test_five_taxa_f84() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAA"),
            ("B", "AACCGGTTAA"),
            ("C", "CCAATTGGCC"),
            ("D", "CCAATTGGCC"),
            ("E", "GGTTAACCGG"),
        ]);
        let tree = parse_newick(
            "((A:0.1,B:0.1):0.05,((C:0.1,D:0.1):0.05,E:0.15):0.05);",
        )
        .unwrap();
        let model = F84Model::equal_freqs(2.0);
        let lnl = log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    #[test]
    fn test_per_site_all_negative() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGAACGT"),
            ("C", "TGCATGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;
        let site_lnl = site_log_likelihoods(&tree, &aln, &model).unwrap();
        for (i, &sl) in site_lnl.iter().enumerate() {
            assert!(
                sl < 0.0,
                "site {} log-likelihood should be negative, got {}",
                i,
                sl
            );
        }
    }

    #[test]
    fn test_known_three_taxon_jc69() {
        // Hand-calculate a simple case:
        // Tree: ((A:t,B:t):0,C:t) with t=0.1
        // Site: A=A, B=A, C=A -> all same
        // P(same; 0.1) = 0.25 + 0.75*exp(-4/30) ≈ 0.9046
        // P(diff; 0.1) ≈ 0.0318
        //
        // At internal (A,B) node:
        //   L(A) = P(A->A;t)*1 * P(A->A;t)*1 = 0.9046^2
        //   L(C) = P(C->A;t)*1 * P(C->A;t)*1 = 0.0318^2
        //   L(G) = P(G->A;t)*1 * P(G->A;t)*1 = 0.0318^2
        //   L(T) = P(T->A;t)*1 * P(T->A;t)*1 = 0.0318^2
        //
        // At root (internal has branch_length 0, so P = identity):
        // Actually the internal node has branch 0 to root. With zero branch,
        // P = identity. So root cond = internal cond * P(root->C; t=0.1).
        // Let me use a simple tree where internal has nonzero branch.

        let aln = make_alignment(&[("A", "A"), ("B", "A"), ("C", "A")]);
        let tree = parse_newick("((A:0.1,B:0.1):0.1,C:0.1);").unwrap();
        let model = Jc69Model;

        let lnl = log_likelihood(&tree, &aln, &model).unwrap();
        // Should be close to ln(0.25) ≈ -1.386 since site is constant
        // but not exactly because of finite branch lengths
        assert!(lnl < 0.0);
        assert!(
            lnl > -3.0,
            "single constant site should have reasonable lnl, got {}",
            lnl
        );
    }

    #[test]
    fn test_optimize_f84() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
        ]);
        let tree = parse_newick("((A:0.5,B:0.5):0.5,C:0.5);").unwrap();
        let model = F84Model::equal_freqs(2.0);

        let lnl_before = log_likelihood(&tree, &aln, &model).unwrap();
        let mut opt_tree = tree.clone();
        let lnl_after =
            optimize_branch_lengths(&mut opt_tree, &aln, &model).unwrap();

        assert!(
            lnl_after >= lnl_before - 1e-6,
            "F84 optimization should not decrease lnl: {} -> {}",
            lnl_before,
            lnl_after
        );
    }
}

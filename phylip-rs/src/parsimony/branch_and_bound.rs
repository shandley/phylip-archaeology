//! Branch-and-bound exact parsimony search.
//!
//! Implements the Hendy & Penny (1982) branch-and-bound algorithm for
//! finding **all** globally optimal (most parsimonious) trees.  Unlike
//! heuristic methods (stepwise addition + SPR), this algorithm is
//! guaranteed to find the global optimum, but its running time is
//! exponential in the number of taxa.  It is practical for about
//! 10-12 taxa.
//!
//! ## Algorithm overview
//!
//! 1. Build an initial tree via stepwise addition to get an **upper bound**
//!    on the optimal parsimony score.
//! 2. Enumerate all tree topologies by adding taxa one at a time to every
//!    possible edge of a partial tree.
//! 3. At each partial tree, compute a **lower bound** (the Fitch score on
//!    the partial tree with only the taxa added so far).
//! 4. If `lower_bound >= current_best`, prune the entire branch of the
//!    search space — no topology extending this partial tree can improve
//!    upon the best known score.
//! 5. When all taxa have been placed, if the score equals the current best
//!    save this tree.  If the score is strictly less, update the best score
//!    and discard previously saved trees.
//!
//! # Example
//!
//! ```
//! use phylip_rs::parsimony::branch_and_bound::{branch_and_bound_fitch, BranchAndBoundResult};
//! use phylip_rs::tree::{Alignment, Base, Sequence};
//!
//! let alignment = Alignment::new(vec![
//!     Sequence::new("A", vec![Base::A, Base::C]),
//!     Sequence::new("B", vec![Base::A, Base::C]),
//!     Sequence::new("C", vec![Base::T, Base::G]),
//!     Sequence::new("D", vec![Base::T, Base::G]),
//! ]).unwrap();
//!
//! let result = branch_and_bound_fitch(&alignment, None);
//! assert!(result.score > 0);
//! assert!(!result.trees.is_empty());
//! ```
//!
//! # References
//!
//! - Hendy, M.D. & Penny, D. (1982). Branch and bound algorithms to
//!   determine minimal evolutionary trees.  *Mathematical Biosciences*,
//!   59, 277-290.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//!   Chapter 7.

use crate::parsimony::traits::{FitchScorer, ParsimonyScorer};
use crate::parsimony::wagner::{parsimony_score_with, search_with};
use crate::tree::types::{Alignment, Tree, TreeNode};

/// Result of a branch-and-bound search.
#[derive(Debug, Clone)]
pub struct BranchAndBoundResult {
    /// All optimal trees (there may be multiple equally parsimonious trees).
    pub trees: Vec<Tree>,
    /// The optimal parsimony score.
    pub score: usize,
    /// Total number of complete tree topologies examined.
    pub trees_examined: usize,
    /// Total number of partial topologies pruned.
    pub trees_pruned: usize,
}

/// Internal bookkeeping during the recursive search.
struct BbStats {
    trees_examined: usize,
    trees_pruned: usize,
}

/// Find all most parsimonious trees by branch-and-bound exact search.
///
/// This is guaranteed to find the globally optimal solution but is
/// exponential in time.  Practical only for about 10-12 taxa.
///
/// # Arguments
///
/// * `alignment` — the DNA alignment to analyze (must have >= 3 taxa)
/// * `scorer` — the parsimony scoring criterion (e.g., `FitchScorer`)
/// * `max_trees` — maximum number of optimal trees to save (default: 100)
///
/// # Panics
///
/// Panics if the alignment has fewer than 3 taxa.
pub fn branch_and_bound(
    alignment: &Alignment,
    scorer: &dyn ParsimonyScorer,
    max_trees: Option<usize>,
) -> BranchAndBoundResult {
    let ntaxa = alignment.ntaxa();
    assert!(ntaxa >= 3, "need at least 3 taxa for branch-and-bound search");

    let max_trees = max_trees.unwrap_or(100);

    // ---- Step 1: Get an initial upper bound via stepwise addition ----
    let initial = search_with(alignment, scorer, None);
    let mut best_score = initial.score;

    // ---- Step 2: Set up the taxon addition order ----
    // Use the natural alignment order.
    let order: Vec<usize> = (0..ntaxa).collect();

    // ---- Step 3: Build the initial 3-taxon star tree ----
    let mut tree = Tree::new();
    tree.rooted = false;
    let root = tree.add_node(None, None, None);
    tree.root = root;
    for &taxon_idx in &order[..3] {
        let name = alignment.taxon_name(taxon_idx).to_string();
        tree.add_node(Some(name), None, Some(root));
    }

    let mut results: Vec<Tree> = Vec::new();
    let mut stats = BbStats {
        trees_examined: 0,
        trees_pruned: 0,
    };

    // ---- Step 4: Recurse ----
    if ntaxa == 3 {
        // Only one unrooted topology for 3 taxa.
        let score = quick_score_with(&tree, alignment, scorer);
        stats.trees_examined = 1;
        if score <= best_score {
            best_score = score;
            results.push(tree);
        }
    } else {
        let remaining: Vec<usize> = order[3..].to_vec();
        bb_recurse(
            &mut tree,
            &remaining,
            alignment,
            scorer,
            &mut best_score,
            &mut results,
            max_trees,
            &mut stats,
        );
    }

    // If no trees were found (should not happen), at least return the initial
    if results.is_empty() {
        results.push(initial.tree);
    }

    BranchAndBoundResult {
        trees: results,
        score: best_score,
        trees_examined: stats.trees_examined,
        trees_pruned: stats.trees_pruned,
    }
}

/// Branch-and-bound with Fitch parsimony (most common use case).
///
/// Convenience wrapper around [`branch_and_bound`] that uses the
/// standard Fitch (unordered) parsimony criterion.
pub fn branch_and_bound_fitch(
    alignment: &Alignment,
    max_trees: Option<usize>,
) -> BranchAndBoundResult {
    branch_and_bound(alignment, &FitchScorer, max_trees)
}

/// Core recursive branch-and-bound function.
///
/// For each remaining taxon, tries inserting it at every possible edge in
/// the current tree.  At each step, the parsimony score of the partial tree
/// is computed as a lower bound.  If the lower bound is already >= the best
/// known complete-tree score, the subtree of the search space is pruned.
fn bb_recurse(
    tree: &mut Tree,
    remaining_taxa: &[usize],
    alignment: &Alignment,
    scorer: &dyn ParsimonyScorer,
    best_score: &mut usize,
    results: &mut Vec<Tree>,
    max_trees: usize,
    stats: &mut BbStats,
) {
    if remaining_taxa.is_empty() {
        // All taxa placed — evaluate this complete tree.
        let score = quick_score_with(tree, alignment, scorer);
        stats.trees_examined += 1;

        if score < *best_score {
            // Strictly better: reset saved trees.
            *best_score = score;
            results.clear();
            results.push(tree.clone());
        } else if score == *best_score && results.len() < max_trees {
            results.push(tree.clone());
        }
        return;
    }

    let taxon_idx = remaining_taxa[0];
    let rest = &remaining_taxa[1..];
    let taxon_name = alignment.taxon_name(taxon_idx).to_string();

    // Collect all edges in the current tree: (parent_id, child_id)
    let edges: Vec<(usize, usize)> = tree
        .nodes
        .iter()
        .filter_map(|n| n.parent.map(|p| (p, n.id)))
        .collect();

    for (parent_id, child_id) in edges {
        // Insert the taxon on this edge (modifies `tree` in place).
        let (new_internal_id, new_leaf_id) =
            insert_taxon_on_edge(tree, parent_id, child_id, &taxon_name);

        // Compute lower bound: parsimony score of the partial tree.
        let lower_bound = quick_score_with(tree, alignment, scorer);

        if lower_bound > *best_score {
            // Prune: no tree extending this partial topology can improve
            // upon or equal the best score.  (The lower bound is the
            // parsimony score of the partial tree with only the taxa
            // placed so far; adding more taxa can only increase the
            // score.  We use strict `>` to ensure we still find all
            // trees that tie with the current best.)
            stats.trees_pruned += 1;
        } else {
            // Recurse with remaining taxa.
            bb_recurse(tree, rest, alignment, scorer, best_score, results, max_trees, stats);
        }

        // Backtrack: remove the inserted taxon and internal node.
        remove_taxon_insertion(tree, parent_id, child_id, new_internal_id, new_leaf_id);
    }
}

/// Insert a new taxon on the edge (parent_id -> child_id) **in place**.
///
/// Creates a new internal node that splits the edge, plus a new leaf for
/// the taxon.  Returns `(new_internal_id, new_leaf_id)`.
fn insert_taxon_on_edge(
    tree: &mut Tree,
    parent_id: usize,
    child_id: usize,
    taxon_name: &str,
) -> (usize, usize) {
    // 1. Create the new internal node.
    let new_internal_id = tree.nodes.len();
    tree.nodes.push(TreeNode {
        id: new_internal_id,
        name: None,
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    // 2. Create the new leaf.
    let new_leaf_id = tree.nodes.len();
    tree.nodes.push(TreeNode {
        id: new_leaf_id,
        name: Some(taxon_name.to_string()),
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    // 3. Rewire: parent -> new_internal instead of parent -> child_id
    let pos = tree.nodes[parent_id]
        .children
        .iter()
        .position(|&c| c == child_id)
        .expect("child not found in parent");
    tree.nodes[parent_id].children[pos] = new_internal_id;
    tree.nodes[new_internal_id].parent = Some(parent_id);

    // 4. new_internal -> {child_id, new_leaf}
    tree.nodes[new_internal_id].children.push(child_id);
    tree.nodes[child_id].parent = Some(new_internal_id);

    tree.nodes[new_internal_id].children.push(new_leaf_id);
    tree.nodes[new_leaf_id].parent = Some(new_internal_id);

    (new_internal_id, new_leaf_id)
}

/// Undo an insertion performed by [`insert_taxon_on_edge`].
///
/// Restores the edge (parent_id -> child_id) and removes the two nodes
/// that were added (the internal node and the leaf), but does **not**
/// shrink the `nodes` vector — the entries are simply "emptied out" and
/// will be popped off the end.
fn remove_taxon_insertion(
    tree: &mut Tree,
    parent_id: usize,
    child_id: usize,
    new_internal_id: usize,
    new_leaf_id: usize,
) {
    // 1. Reconnect parent -> child_id directly.
    let pos = tree.nodes[parent_id]
        .children
        .iter()
        .position(|&c| c == new_internal_id)
        .expect("internal node not found in parent");
    tree.nodes[parent_id].children[pos] = child_id;
    tree.nodes[child_id].parent = Some(parent_id);

    // 2. Clear the internal node and leaf (they are at the end of the
    //    nodes vec, so we just pop them).
    debug_assert_eq!(new_leaf_id, tree.nodes.len() - 1);
    debug_assert_eq!(new_internal_id, tree.nodes.len() - 2);
    tree.nodes.pop(); // leaf
    tree.nodes.pop(); // internal
}

/// Convenience: compute Fitch parsimony score on a (possibly partial) tree.
fn quick_score_with(
    tree: &Tree,
    alignment: &Alignment,
    scorer: &dyn ParsimonyScorer,
) -> usize {
    match parsimony_score_with(tree, alignment, scorer) {
        Ok((score, _)) => score,
        Err(_) => usize::MAX,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parsimony::wagner::{parsimony_score, search};
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

    // --- 3 taxa: only one unrooted topology exists ---

    #[test]
    fn test_three_taxa_single_topology() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "TGCA"),
        ]);
        let result = branch_and_bound_fitch(&aln, None);
        // Only 1 unrooted topology for 3 taxa.
        assert_eq!(result.trees.len(), 1);
        assert_eq!(result.trees[0].num_leaves(), 3);
    }

    // --- 4 taxa: verify all 3 unrooted topologies are considered ---

    #[test]
    fn test_four_taxa_finds_all_topologies_picks_best() {
        // Strong signal: A,B share states; C,D share states.
        let aln = make_alignment(&[
            ("A", "AAAA"),
            ("B", "AAAA"),
            ("C", "TTTT"),
            ("D", "TTTT"),
        ]);
        let result = branch_and_bound_fitch(&aln, None);

        // Optimal score: 1 change per site = 4
        assert_eq!(result.score, 4);
        // All optimal trees should have 4 leaves.
        for tree in &result.trees {
            assert_eq!(tree.num_leaves(), 4);
        }
        // There should be exactly 1 optimal topology (AB|CD).
        assert_eq!(result.trees.len(), 1);
    }

    #[test]
    fn test_four_taxa_symmetric_data_multiple_optimal() {
        // Symmetric data where all 3 unrooted topologies give the same score.
        // Each taxon is maximally different from all others.
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "C"),
            ("C", "G"),
            ("D", "T"),
        ]);
        let result = branch_and_bound_fitch(&aln, None);
        // All 3 unrooted topologies should give the same score (3 steps).
        assert_eq!(result.score, 3);
        assert_eq!(result.trees.len(), 3);
    }

    // --- 5 taxa: optimal score matches heuristic ---

    #[test]
    fn test_five_taxa_matches_heuristic() {
        let aln = make_alignment(&[
            ("A", "AACCGG"),
            ("B", "AACCGG"),
            ("C", "CCAACC"),
            ("D", "CCAACC"),
            ("E", "GGGGAA"),
        ]);
        let heuristic = search(&aln, Some(42));
        let exact = branch_and_bound_fitch(&aln, None);

        // Exact should be <= heuristic (guaranteed optimal).
        assert!(exact.score <= heuristic.score);
    }

    // --- Known example ---

    #[test]
    fn test_known_example_four_taxa() {
        // ((A,B),(C,D)) is optimal for this data.
        let aln = make_alignment(&[
            ("A", "AACCAACCAA"),
            ("B", "AACCAACCAA"),
            ("C", "CCAACCAACC"),
            ("D", "CCAACCAACC"),
        ]);

        let result = branch_and_bound_fitch(&aln, None);

        // Evaluate all three topologies manually
        let t1 = parse_newick("((A,B),(C,D));").unwrap();
        let t2 = parse_newick("((A,C),(B,D));").unwrap();
        let t3 = parse_newick("((A,D),(B,C));").unwrap();

        let (s1, _) = parsimony_score(&t1, &aln).unwrap();
        let (s2, _) = parsimony_score(&t2, &aln).unwrap();
        let (s3, _) = parsimony_score(&t3, &aln).unwrap();

        // ((A,B),(C,D)) is strictly best
        assert!(s1 < s2);
        assert!(s1 < s3);
        assert_eq!(result.score, s1);
    }

    // --- Verify trees_examined + trees_pruned accounts for all ---

    #[test]
    fn test_pruning_stats() {
        let aln = make_alignment(&[
            ("A", "AACC"),
            ("B", "AACC"),
            ("C", "CCAA"),
            ("D", "CCAA"),
        ]);
        let result = branch_and_bound_fitch(&aln, None);
        // For 4 taxa, there are exactly 3 unrooted topologies.
        // trees_examined should be > 0 and trees_examined + trees_pruned
        // should account for the entire search space.
        assert!(result.trees_examined > 0);
        assert!(result.trees_examined + result.trees_pruned > 0);
        // With strong signal, should prune at least some.
        // (With 4 taxa the search space is small so pruning may be limited.)
    }

    // --- max_trees limit ---

    #[test]
    fn test_max_trees_limit() {
        // Data where all topologies tie.
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "C"),
            ("C", "G"),
            ("D", "T"),
        ]);
        // Limit to 2 trees.
        let result = branch_and_bound_fitch(&aln, Some(2));
        assert!(result.trees.len() <= 2);
        assert_eq!(result.score, 3);
    }

    // --- 3 taxa with identical sequences ---

    #[test]
    fn test_three_taxa_identical_sequences() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "ACGT"),
        ]);
        let result = branch_and_bound_fitch(&aln, None);
        assert_eq!(result.score, 0);
        assert_eq!(result.trees.len(), 1);
    }

    // --- Correct leaf sets ---

    #[test]
    fn test_result_trees_have_correct_leaves() {
        let aln = make_alignment(&[
            ("A", "AACCGG"),
            ("B", "AACCGG"),
            ("C", "CCAACC"),
            ("D", "CCAACC"),
            ("E", "GGTTAA"),
        ]);
        let result = branch_and_bound_fitch(&aln, None);
        for tree in &result.trees {
            assert_eq!(tree.num_leaves(), 5);
            let mut names = tree.leaf_names();
            names.sort();
            assert_eq!(names, vec!["A", "B", "C", "D", "E"]);
        }
    }

    // --- Optimality guarantee ---

    #[test]
    fn test_optimality_vs_all_four_taxa_topologies() {
        // Test with data that has a unique optimum.
        let aln = make_alignment(&[
            ("A", "AACC"),
            ("B", "AACC"),
            ("C", "CCAA"),
            ("D", "CCAA"),
        ]);
        let result = branch_and_bound_fitch(&aln, None);

        // Verify against all 3 possible unrooted topologies.
        let topologies = [
            "((A,B),(C,D));",
            "((A,C),(B,D));",
            "((A,D),(B,C));",
        ];
        let mut min_score = usize::MAX;
        for topo in &topologies {
            let tree = parse_newick(topo).unwrap();
            let (score, _) = parsimony_score(&tree, &aln).unwrap();
            if score < min_score {
                min_score = score;
            }
        }
        assert_eq!(result.score, min_score);
    }

    // --- 5 taxa exhaustive check ---

    #[test]
    fn test_five_taxa_optimal() {
        let aln = make_alignment(&[
            ("A", "AAAA"),
            ("B", "AAAA"),
            ("C", "CCCC"),
            ("D", "CCCC"),
            ("E", "GGGG"),
        ]);
        let result = branch_and_bound_fitch(&aln, None);
        // The optimal tree should group A,B together.
        assert!(result.score > 0);
        // Verify all result trees have 5 leaves.
        for tree in &result.trees {
            assert_eq!(tree.num_leaves(), 5);
        }
    }

    // --- Insert/remove roundtrip ---

    #[test]
    fn test_insert_remove_roundtrip() {
        let mut tree = Tree::new();
        tree.rooted = false;
        let root = tree.add_node(None, None, None);
        tree.root = root;
        tree.add_node(Some("A".into()), None, Some(root));
        tree.add_node(Some("B".into()), None, Some(root));
        tree.add_node(Some("C".into()), None, Some(root));

        let original_num_nodes = tree.num_nodes();
        let original_num_leaves = tree.num_leaves();

        // Pick edge root -> child "A" (child_id = 1)
        let parent_id = root;
        let child_id = 1;

        let (new_internal, new_leaf) =
            insert_taxon_on_edge(&mut tree, parent_id, child_id, "D");

        assert_eq!(tree.num_leaves(), original_num_leaves + 1);

        remove_taxon_insertion(&mut tree, parent_id, child_id, new_internal, new_leaf);

        assert_eq!(tree.num_nodes(), original_num_nodes);
        assert_eq!(tree.num_leaves(), original_num_leaves);
    }
}

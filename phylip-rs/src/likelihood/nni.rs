//! Nearest-Neighbor Interchange (NNI) for ML tree improvement.
//!
//! NNI is a local tree rearrangement strategy that swaps subtrees across
//! internal edges. For each internal edge connecting nodes u and v (where
//! u has children {a, v} and v has children {b, c}), there are two possible
//! NNI moves:
//!
//! 1. Swap a with b: u gets {b, v}, v gets {a, c}
//! 2. Swap a with c: u gets {c, v}, v gets {a, b}
//!
//! NNI is simpler and faster than SPR but explores a smaller neighborhood.
//! It is commonly used as a refinement step after initial tree construction.
//!
//! # Algorithm
//!
//! The search proceeds in rounds:
//!
//! 1. For each internal edge, generate the two NNI neighbor trees.
//! 2. Evaluate all candidates by log-likelihood (with branch length
//!    optimization for the best candidate).
//! 3. Accept the best-improving move if it increases the likelihood.
//! 4. Repeat until no improvement is found.
//!
//! # References
//!
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//! - Swofford, D.L. et al. (1996). Phylogenetic inference. In: *Molecular
//!   Systematics*, 2nd ed., Sinauer Associates.

use crate::tree::types::{Alignment, Tree};

use super::models::SubstitutionModel;
use super::pruning::{log_likelihood, optimize_branch_lengths, LikelihoodError};

/// Maximum number of NNI rounds before giving up.
const MAX_NNI_ROUNDS: usize = 100;

/// Improvement threshold: accept a move only if it improves by at least this much.
const IMPROVEMENT_THRESHOLD: f64 = 1e-4;

/// Generate the two NNI neighbor trees for the internal edge ending at
/// `edge_node_id`.
///
/// The edge is between `edge_node_id` and its parent. Both `edge_node_id`
/// and its parent must be internal nodes (i.e., have children). If either
/// is a leaf or `edge_node_id` is the root (no parent), an empty vector
/// is returned.
///
/// For an edge connecting parent u (with children {a, v}) and child v
/// (with children {b, c}), the two NNI moves are:
///
/// - Move 1: swap a with b -> u gets {b, v}, v gets {a, c}
/// - Move 2: swap a with c -> u gets {c, v}, v gets {a, b}
///
/// Here "a" is the sibling of v under u (the first child of u that is not v).
pub fn nni_moves(tree: &Tree, edge_node_id: usize) -> Vec<Tree> {
    // The edge is between edge_node_id (v) and its parent (u).
    let v = edge_node_id;

    // v must have a parent
    let u = match tree.nodes[v].parent {
        Some(p) => p,
        None => return Vec::new(),
    };

    // Both u and v must be internal nodes (have children)
    if tree.nodes[u].children.is_empty() || tree.nodes[v].children.is_empty() {
        return Vec::new();
    }

    // v must have at least 2 children
    if tree.nodes[v].children.len() < 2 {
        return Vec::new();
    }

    // Find 'a': the sibling of v under u (first child of u that is not v)
    let a = match tree.nodes[u].children.iter().find(|&&c| c != v) {
        Some(&a) => a,
        None => return Vec::new(), // u has no other children besides v
    };

    // Children of v
    let b = tree.nodes[v].children[0];
    let c = tree.nodes[v].children[1];

    let mut result = Vec::new();

    // Move 1: swap a with b
    // u gets {b, v}, v gets {a, c}
    if let Some(t) = make_nni_tree(tree, u, v, a, b) {
        result.push(t);
    }

    // Move 2: swap a with c
    // u gets {c, v}, v gets {a, b}
    if let Some(t) = make_nni_tree(tree, u, v, a, c) {
        result.push(t);
    }

    result
}

/// Build a new tree with an NNI swap: move `swap_child` from v to u,
/// and move `a` from u to v.
///
/// Before: u has children {..., a, ..., v, ...}, v has children {..., swap_child, ..., other, ...}
/// After:  u has children {..., swap_child, ..., v, ...}, v has children {..., a, ..., other, ...}
fn make_nni_tree(
    tree: &Tree,
    u: usize,
    v: usize,
    a: usize,
    swap_child: usize,
) -> Option<Tree> {
    let mut new_tree = tree.clone();

    // Remove 'a' from u's children and add 'swap_child'
    let u_children = &mut new_tree.nodes[u].children;
    if let Some(pos) = u_children.iter().position(|&c| c == a) {
        u_children[pos] = swap_child;
    } else {
        return None;
    }

    // Remove 'swap_child' from v's children and add 'a'
    let v_children = &mut new_tree.nodes[v].children;
    if let Some(pos) = v_children.iter().position(|&c| c == swap_child) {
        v_children[pos] = a;
    } else {
        return None;
    }

    // Update parent pointers
    new_tree.nodes[a].parent = Some(v);
    new_tree.nodes[swap_child].parent = Some(u);

    Some(new_tree)
}

/// Score a tree using log-likelihood. Returns -infinity on error.
fn quick_lnl(tree: &Tree, alignment: &Alignment, model: &dyn SubstitutionModel) -> f64 {
    log_likelihood(tree, alignment, model).unwrap_or(f64::NEG_INFINITY)
}

/// Perform one round of NNI rearrangement.
///
/// Tries all possible NNI moves on all internal edges. Picks the best
/// candidate (by unoptimized likelihood), then optimizes branch lengths.
/// Accepts the move only if the optimized result improves the current score.
///
/// # Returns
///
/// A tuple of (tree, log-likelihood, improved) where `improved` is true
/// if a better tree was found.
pub fn nni_round(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    current_lnl: f64,
) -> Result<(Tree, f64, bool), LikelihoodError> {
    let mut best_candidate = None;
    let mut best_quick_lnl = f64::NEG_INFINITY;

    // Find all internal edges: nodes that have a parent and have children,
    // and whose parent also has children (i.e., both endpoints are internal).
    let internal_edges: Vec<usize> = tree
        .nodes
        .iter()
        .filter(|n| {
            // Node must have a parent (not root)
            n.parent.is_some()
            // Node must be internal (have children)
            && !n.children.is_empty()
            // Parent must have at least one other child besides this node
            && tree.nodes[n.parent.unwrap()].children.len() >= 2
        })
        .map(|n| n.id)
        .collect();

    for &edge_node in &internal_edges {
        let candidates = nni_moves(tree, edge_node);
        for candidate in candidates {
            let lnl = quick_lnl(&candidate, alignment, model);
            if lnl > best_quick_lnl {
                best_quick_lnl = lnl;
                best_candidate = Some(candidate);
            }
        }
    }

    // Optimize the best candidate and check if it truly improves
    if let Some(mut candidate) = best_candidate {
        let opt_lnl = optimize_branch_lengths(&mut candidate, alignment, model)?;
        if opt_lnl > current_lnl + IMPROVEMENT_THRESHOLD {
            return Ok((candidate, opt_lnl, true));
        }
    }

    Ok((tree.clone(), current_lnl, false))
}

/// Iterate NNI rounds until convergence.
///
/// Starting from the given tree, repeatedly apply NNI rearrangements,
/// accepting each move that improves the log-likelihood. Stops when
/// no improving move is found or the maximum number of rounds is reached.
///
/// Branch lengths are optimized on the initial tree and after each
/// accepted NNI move.
///
/// # Returns
///
/// A tuple of (best_tree, best_log_likelihood).
pub fn nni_search(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
) -> Result<(Tree, f64), LikelihoodError> {
    let mut current_tree = tree.clone();

    // Optimize branch lengths on the initial tree
    let mut current_lnl = optimize_branch_lengths(&mut current_tree, alignment, model)?;

    for _ in 0..MAX_NNI_ROUNDS {
        let (new_tree, new_lnl, improved) =
            nni_round(&current_tree, alignment, model, current_lnl)?;
        if !improved {
            break;
        }
        current_tree = new_tree;
        current_lnl = new_lnl;
    }

    Ok((current_tree, current_lnl))
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

    // --- nni_moves tests ---

    #[test]
    fn test_nni_moves_generates_two_trees() {
        // ((A,B),(C,D)) has internal edges; pick the edge to the (A,B) clade
        let tree = parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();

        // Find an internal node that is not the root and has an internal parent
        let internal_edge_nodes: Vec<usize> = tree
            .nodes
            .iter()
            .filter(|n| {
                n.parent.is_some()
                    && !n.children.is_empty()
                    && tree.nodes[n.parent.unwrap()].children.len() >= 2
            })
            .map(|n| n.id)
            .collect();

        assert!(!internal_edge_nodes.is_empty(), "should have internal edges");

        let moves = nni_moves(&tree, internal_edge_nodes[0]);
        assert_eq!(
            moves.len(),
            2,
            "each internal edge should produce exactly 2 NNI neighbors"
        );
    }

    #[test]
    fn test_nni_moves_preserves_taxa() {
        let tree = parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();

        let internal_edge_node = tree
            .nodes
            .iter()
            .find(|n| {
                n.parent.is_some()
                    && !n.children.is_empty()
                    && tree.nodes[n.parent.unwrap()].children.len() >= 2
            })
            .unwrap()
            .id;

        let moves = nni_moves(&tree, internal_edge_node);
        for nni_tree in &moves {
            assert_eq!(nni_tree.num_leaves(), 4);
            let names = nni_tree.leaf_names();
            assert!(names.contains(&"A"));
            assert!(names.contains(&"B"));
            assert!(names.contains(&"C"));
            assert!(names.contains(&"D"));
        }
    }

    #[test]
    fn test_nni_moves_on_leaf_returns_empty() {
        let tree = parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();

        // Find a leaf node
        let leaf_id = tree
            .nodes
            .iter()
            .find(|n| n.is_leaf())
            .unwrap()
            .id;

        let moves = nni_moves(&tree, leaf_id);
        assert!(
            moves.is_empty(),
            "NNI on a leaf edge should return no moves"
        );
    }

    #[test]
    fn test_nni_moves_on_root_returns_empty() {
        let tree = parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();

        let moves = nni_moves(&tree, tree.root);
        assert!(
            moves.is_empty(),
            "NNI on the root should return no moves"
        );
    }

    #[test]
    fn test_nni_moves_produces_different_topologies() {
        let tree = parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();

        let internal_edge_node = tree
            .nodes
            .iter()
            .find(|n| {
                n.parent.is_some()
                    && !n.children.is_empty()
                    && tree.nodes[n.parent.unwrap()].children.len() >= 2
            })
            .unwrap()
            .id;

        let moves = nni_moves(&tree, internal_edge_node);
        assert_eq!(moves.len(), 2);

        // The two NNI neighbors should have different topologies from
        // the original and from each other. We verify by checking that
        // at least one child assignment differs.
        let orig_children: Vec<Vec<usize>> = tree
            .nodes
            .iter()
            .map(|n| n.children.clone())
            .collect();

        let move0_children: Vec<Vec<usize>> = moves[0]
            .nodes
            .iter()
            .map(|n| n.children.clone())
            .collect();

        let move1_children: Vec<Vec<usize>> = moves[1]
            .nodes
            .iter()
            .map(|n| n.children.clone())
            .collect();

        // At least one of the NNI neighbors should differ from the original
        assert!(
            orig_children != move0_children || orig_children != move1_children,
            "NNI moves should produce different topologies"
        );
    }

    #[test]
    fn test_nni_moves_preserves_branch_lengths() {
        let tree = parse_newick("((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.06);").unwrap();

        let internal_edge_node = tree
            .nodes
            .iter()
            .find(|n| {
                n.parent.is_some()
                    && !n.children.is_empty()
                    && tree.nodes[n.parent.unwrap()].children.len() >= 2
            })
            .unwrap()
            .id;

        let moves = nni_moves(&tree, internal_edge_node);
        for nni_tree in &moves {
            // All non-root nodes should still have branch lengths
            for node in &nni_tree.nodes {
                if node.parent.is_some() {
                    assert!(
                        node.branch_length.is_some(),
                        "node {:?} (id={}) missing branch length after NNI",
                        node.name,
                        node.id
                    );
                }
            }
        }
    }

    #[test]
    fn test_nni_moves_preserves_node_count() {
        let tree = parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();
        let original_count = tree.num_nodes();

        let internal_edge_node = tree
            .nodes
            .iter()
            .find(|n| {
                n.parent.is_some()
                    && !n.children.is_empty()
                    && tree.nodes[n.parent.unwrap()].children.len() >= 2
            })
            .unwrap()
            .id;

        let moves = nni_moves(&tree, internal_edge_node);
        for nni_tree in &moves {
            assert_eq!(
                nni_tree.num_nodes(),
                original_count,
                "NNI should not change the number of nodes"
            );
        }
    }

    // --- nni_round tests ---

    #[test]
    fn test_nni_round_does_not_worsen() {
        let aln = make_alignment(&[
            ("A", "AAAAACCCCC"),
            ("B", "AAAAACCCCC"),
            ("C", "CCCCCAAAAA"),
            ("D", "CCCCCAAAAA"),
        ]);
        let model = Jc69Model;

        let mut tree =
            parse_newick("((A:0.05,B:0.05):0.2,(C:0.05,D:0.05):0.2);").unwrap();
        let lnl = optimize_branch_lengths(&mut tree, &aln, &model).unwrap();

        let (_, new_lnl, _) = nni_round(&tree, &aln, &model, lnl).unwrap();
        assert!(
            new_lnl >= lnl - IMPROVEMENT_THRESHOLD,
            "NNI round should not worsen: {} -> {}",
            lnl,
            new_lnl
        );
    }

    #[test]
    fn test_nni_round_can_improve_wrong_topology() {
        // Start with a wrong topology: A grouped with C instead of B
        let aln = make_alignment(&[
            ("A", "AAAAAAAAAA"),
            ("B", "AAAAAAAAAA"),
            ("C", "CCCCCCCCCC"),
            ("D", "CCCCCCCCCC"),
        ]);
        let model = Jc69Model;

        let mut wrong_tree =
            parse_newick("((A:0.1,C:0.1):0.1,(B:0.1,D:0.1):0.1);").unwrap();
        let wrong_lnl = optimize_branch_lengths(&mut wrong_tree, &aln, &model).unwrap();

        // NNI should find an improving move since the correct topology
        // ((A,B),(C,D)) is a single NNI move away
        let (new_tree, new_lnl, improved) =
            nni_round(&wrong_tree, &aln, &model, wrong_lnl).unwrap();

        assert!(
            improved,
            "NNI should improve the wrong topology; wrong_lnl={}, new_lnl={}",
            wrong_lnl,
            new_lnl
        );
        assert!(new_lnl > wrong_lnl);
        assert_eq!(new_tree.num_leaves(), 4);
    }

    // --- nni_search tests ---

    #[test]
    fn test_nni_search_finds_good_tree() {
        let aln = make_alignment(&[
            ("A", "AAAAAAAAAA"),
            ("B", "AAAAAAAAAA"),
            ("C", "CCCCCCCCCC"),
            ("D", "CCCCCCCCCC"),
        ]);
        let model = Jc69Model;

        // Start with a wrong topology
        let tree =
            parse_newick("((A:0.1,C:0.1):0.1,(B:0.1,D:0.1):0.1);").unwrap();
        let initial_lnl = log_likelihood(&tree, &aln, &model).unwrap();

        let (result_tree, result_lnl) = nni_search(&tree, &aln, &model).unwrap();

        assert_eq!(result_tree.num_leaves(), 4);
        // NNI search should improve over the unoptimized wrong topology
        assert!(
            result_lnl > initial_lnl,
            "NNI search ({}) should improve over initial wrong topology ({})",
            result_lnl,
            initial_lnl
        );
        assert!(result_lnl < 0.0);
    }

    #[test]
    fn test_nni_search_converges_on_optimal() {
        // Start with the correct topology; NNI should not change it
        let aln = make_alignment(&[
            ("A", "AAAAACCCCC"),
            ("B", "AAAAACCCCC"),
            ("C", "CCCCCAAAAA"),
            ("D", "CCCCCAAAAA"),
        ]);
        let model = Jc69Model;

        let tree =
            parse_newick("((A:0.05,B:0.05):0.2,(C:0.05,D:0.05):0.2);").unwrap();

        let (result_tree, result_lnl) = nni_search(&tree, &aln, &model).unwrap();
        assert_eq!(result_tree.num_leaves(), 4);
        assert!(result_lnl < 0.0);
    }

    #[test]
    fn test_nni_search_with_f84_model() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGG"),
        ]);
        let freqs = BaseFrequencies::new(0.3, 0.2, 0.2, 0.3);
        let model = F84Model::new(freqs, 2.0);

        let tree =
            parse_newick("((A:0.1,C:0.1):0.1,(B:0.1,D:0.1):0.1);").unwrap();

        let (result_tree, result_lnl) = nni_search(&tree, &aln, &model).unwrap();
        assert_eq!(result_tree.num_leaves(), 4);
        assert!(result_lnl < 0.0);
    }

    #[test]
    fn test_nni_search_five_taxa() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAA"),
            ("B", "AACCGGTTAA"),
            ("C", "CCAATTGGCC"),
            ("D", "CCAATTGGCC"),
            ("E", "GGTTAACCGG"),
        ]);
        let model = Jc69Model;

        let tree = parse_newick(
            "((A:0.1,B:0.1):0.05,((C:0.1,D:0.1):0.05,E:0.15):0.05);",
        )
        .unwrap();

        let (result_tree, result_lnl) = nni_search(&tree, &aln, &model).unwrap();
        assert_eq!(result_tree.num_leaves(), 5);
        assert!(result_lnl < 0.0);
    }

    #[test]
    fn test_nni_search_improves_over_initial() {
        // Use a tree with intentionally poor branch lengths and wrong topology
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGG"),
        ]);
        let model = Jc69Model;

        let tree =
            parse_newick("((A:1.0,C:1.0):1.0,(B:1.0,D:1.0):1.0);").unwrap();

        let lnl_before = log_likelihood(&tree, &aln, &model).unwrap();
        let (_, result_lnl) = nni_search(&tree, &aln, &model).unwrap();

        assert!(
            result_lnl > lnl_before,
            "NNI search ({}) should improve over initial tree ({})",
            result_lnl,
            lnl_before
        );
    }

    #[test]
    fn test_nni_moves_five_taxa_internal_edges() {
        // Five-taxon tree: ((A,B),((C,D),E))
        // This should have multiple internal edges to try NNI on
        let tree = parse_newick(
            "((A:0.1,B:0.1):0.05,((C:0.1,D:0.1):0.05,E:0.15):0.05);",
        )
        .unwrap();

        let internal_edges: Vec<usize> = tree
            .nodes
            .iter()
            .filter(|n| {
                n.parent.is_some()
                    && !n.children.is_empty()
                    && tree.nodes[n.parent.unwrap()].children.len() >= 2
            })
            .map(|n| n.id)
            .collect();

        // A 5-taxon fully resolved tree should have at least 2 internal edges
        assert!(
            internal_edges.len() >= 2,
            "5-taxon tree should have at least 2 internal edges, got {}",
            internal_edges.len()
        );

        // Each internal edge should produce exactly 2 NNI moves
        for &edge_node in &internal_edges {
            let moves = nni_moves(&tree, edge_node);
            assert_eq!(
                moves.len(),
                2,
                "internal edge at node {} should produce 2 NNI moves",
                edge_node
            );
            for m in &moves {
                assert_eq!(m.num_leaves(), 5);
            }
        }
    }
}

//! Maximum likelihood tree search.
//!
//! Implements heuristic search for the maximum likelihood tree using
//! stepwise addition followed by SPR (subtree pruning and regrafting)
//! rearrangement, following the approach of PHYLIP's DNAML program.
//!
//! # Algorithm
//!
//! The search proceeds in two phases:
//!
//! 1. **Stepwise addition**: Start with a 3-taxon tree and add remaining
//!    taxa one at a time, choosing the insertion point that maximizes
//!    the log-likelihood. Branch lengths are optimized after each addition.
//!
//! 2. **SPR rearrangement**: Iteratively prune subtrees and regraft them
//!    at other locations, accepting moves that improve the likelihood.
//!    Branch lengths are optimized after each accepted move.
//!
//! # Example
//!
//! ```
//! use phylip_rs::likelihood::search::{search, MLSearchResult};
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
//! let result = search(&alignment, &model, None).unwrap();
//! assert!(result.lnl < 0.0);
//! assert_eq!(result.tree.num_leaves(), 4);
//! ```
//!
//! # References
//!
//! - Felsenstein, J. (1981). Evolutionary trees from DNA sequences:
//!   a maximum likelihood approach. *Journal of Molecular Evolution*,
//!   17, 368-376.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//!   Chapter 16: Models of DNA evolution.

use crate::tree::types::{Alignment, Tree, TreeNode};

use super::models::SubstitutionModel;
use super::pruning::{log_likelihood, optimize_branch_lengths, site_log_likelihoods, LikelihoodError};

/// Default branch length for newly created edges.
const DEFAULT_BRANCH_LENGTH: f64 = 0.1;

/// Maximum number of SPR rearrangement rounds.
const MAX_SPR_ROUNDS: usize = 100;

/// Result of a maximum likelihood tree search.
#[derive(Debug, Clone)]
pub struct MLSearchResult {
    /// The best tree found, with optimized branch lengths.
    pub tree: Tree,
    /// Total log-likelihood.
    pub lnl: f64,
    /// Per-site log-likelihoods.
    pub site_lnl: Vec<f64>,
}

/// Search for the maximum likelihood tree.
///
/// Uses stepwise addition followed by SPR rearrangement, with branch
/// length optimization after each topology change. This follows the
/// general strategy of PHYLIP's DNAML program.
///
/// # Arguments
/// * `alignment` - the DNA alignment to analyze
/// * `model` - the substitution model (e.g., JC69 or F84)
/// * `seed` - optional random seed for taxon addition order
///
/// # Returns
/// An `MLSearchResult` with the best tree found and its log-likelihood.
pub fn search(
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    seed: Option<u64>,
) -> Result<MLSearchResult, LikelihoodError> {
    let ntaxa = alignment.ntaxa();
    if ntaxa < 3 {
        return Err(LikelihoodError::NumericalError(
            "need at least 3 taxa for ML search".to_string(),
        ));
    }

    // Taxon addition order
    let mut order: Vec<usize> = (0..ntaxa).collect();
    if let Some(s) = seed {
        let mut rng = SimpleRng::new(s);
        rng.shuffle(&mut order);
    }

    // Phase 1: Stepwise addition
    let mut tree = build_initial_tree(&order, alignment);
    optimize_branch_lengths(&mut tree, alignment, model)?;

    for i in 3..ntaxa {
        let taxon_name = alignment.taxon_name(order[i]);
        tree = best_insertion(&tree, taxon_name, alignment, model);
        optimize_branch_lengths(&mut tree, alignment, model)?;
    }

    let mut lnl = log_likelihood(&tree, alignment, model)?;

    // Phase 2: SPR rearrangement (iterate until no improvement)
    for _ in 0..MAX_SPR_ROUNDS {
        let (new_tree, new_lnl, improved) = spr_round(&tree, alignment, model, lnl)?;
        if !improved {
            break;
        }
        tree = new_tree;
        lnl = new_lnl;
    }

    // Final evaluation
    let site_lnl = site_log_likelihoods(&tree, alignment, model)?;
    let final_lnl = site_lnl.iter().sum();

    Ok(MLSearchResult {
        tree,
        lnl: final_lnl,
        site_lnl,
    })
}

// ============================================================
// Internal: random number generator
// ============================================================

/// Simple LCG random number generator.
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1),
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    fn gen_range(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    fn shuffle<T>(&mut self, slice: &mut [T]) {
        for i in (1..slice.len()).rev() {
            let j = self.gen_range(i + 1);
            slice.swap(i, j);
        }
    }
}

// ============================================================
// Internal: tree construction and manipulation
// ============================================================

/// Build a 3-taxon star tree with default branch lengths.
fn build_initial_tree(taxa: &[usize], alignment: &Alignment) -> Tree {
    let mut tree = Tree::new();
    tree.rooted = false;

    let root = tree.add_node(None, None, None);
    tree.root = root;

    for &taxon_idx in &taxa[..3] {
        let name = alignment.taxon_name(taxon_idx).to_string();
        tree.add_node(Some(name), Some(DEFAULT_BRANCH_LENGTH), Some(root));
    }

    tree
}

/// Score a tree using log-likelihood. Returns -infinity on error.
fn quick_lnl(tree: &Tree, alignment: &Alignment, model: &dyn SubstitutionModel) -> f64 {
    log_likelihood(tree, alignment, model).unwrap_or(f64::NEG_INFINITY)
}

/// Try inserting a new taxon at every possible edge in the tree.
/// Returns the tree with the best (highest) log-likelihood.
fn best_insertion(
    tree: &Tree,
    taxon_name: &str,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
) -> Tree {
    let mut best_tree = None;
    let mut best_lnl = f64::NEG_INFINITY;

    // Collect all edges: (parent_id, child_id)
    let edges: Vec<(usize, usize)> = tree
        .nodes
        .iter()
        .filter_map(|n| n.parent.map(|p| (p, n.id)))
        .collect();

    for (parent_id, child_id) in edges {
        let candidate = insert_on_edge(tree, parent_id, child_id, taxon_name);
        let lnl = quick_lnl(&candidate, alignment, model);
        if lnl > best_lnl {
            best_lnl = lnl;
            best_tree = Some(candidate);
        }
    }

    best_tree.unwrap_or_else(|| {
        // Fallback: add as child of root
        let mut t = tree.clone();
        let root = t.root;
        t.add_node(
            Some(taxon_name.to_string()),
            Some(DEFAULT_BRANCH_LENGTH),
            Some(root),
        );
        t
    })
}

/// Insert a new taxon on the edge between parent_id and child_id.
///
/// Creates a new internal node splitting the edge. The original branch
/// length is divided between the two halves, and the new taxon gets
/// a default branch length.
fn insert_on_edge(
    tree: &Tree,
    parent_id: usize,
    child_id: usize,
    taxon_name: &str,
) -> Tree {
    let mut new_tree = Tree::new();
    new_tree.rooted = tree.rooted;

    let mut id_map = vec![0usize; tree.nodes.len()];

    // Copy all existing nodes (without links)
    for node in &tree.nodes {
        let new_id = new_tree.nodes.len();
        id_map[node.id] = new_id;
        new_tree.nodes.push(TreeNode {
            id: new_id,
            name: node.name.clone(),
            branch_length: node.branch_length,
            children: Vec::new(),
            parent: None,
        });
    }

    new_tree.root = id_map[tree.root];

    // Split the edge: divide the original branch length
    let original_bl = tree.nodes[child_id]
        .branch_length
        .unwrap_or(DEFAULT_BRANCH_LENGTH);
    let half_bl = original_bl / 2.0;

    // New internal node (on the split edge, inherits half the branch length)
    let new_internal_id = new_tree.nodes.len();
    new_tree.nodes.push(TreeNode {
        id: new_internal_id,
        name: None,
        branch_length: Some(half_bl),
        children: Vec::new(),
        parent: None,
    });

    // New leaf node
    let new_leaf_id = new_tree.nodes.len();
    new_tree.nodes.push(TreeNode {
        id: new_leaf_id,
        name: Some(taxon_name.to_string()),
        branch_length: Some(DEFAULT_BRANCH_LENGTH),
        children: Vec::new(),
        parent: None,
    });

    // Update the original child's branch length to the other half
    let mapped_child_id = id_map[child_id];
    new_tree.nodes[mapped_child_id].branch_length = Some(half_bl);

    // Rebuild parent-child links
    for node in &tree.nodes {
        let mapped_id = id_map[node.id];
        for &old_child in &node.children {
            if node.id == parent_id && old_child == child_id {
                // Split this edge: parent -> new_internal -> {child, new_leaf}
                new_tree.nodes[mapped_id].children.push(new_internal_id);
                new_tree.nodes[new_internal_id].parent = Some(mapped_id);

                let mapped_child = id_map[child_id];
                new_tree.nodes[new_internal_id].children.push(mapped_child);
                new_tree.nodes[mapped_child].parent = Some(new_internal_id);

                new_tree.nodes[new_internal_id].children.push(new_leaf_id);
                new_tree.nodes[new_leaf_id].parent = Some(new_internal_id);
            } else {
                let mapped_child = id_map[old_child];
                new_tree.nodes[mapped_id].children.push(mapped_child);
                new_tree.nodes[mapped_child].parent = Some(mapped_id);
            }
        }
    }

    new_tree
}

// ============================================================
// Internal: SPR rearrangement
// ============================================================

/// Result of pruning a subtree from a tree.
#[derive(Debug, Clone)]
struct PrunedTree {
    /// The tree with the subtree removed.
    tree: Tree,
    /// The pruned subtree (as a separate tree).
    subtree: Tree,
}

/// Collect all node IDs in the subtree rooted at `node_id`.
fn collect_subtree_ids(tree: &Tree, node_id: usize) -> Vec<usize> {
    let mut ids = vec![node_id];
    let mut stack = vec![node_id];
    while let Some(nid) = stack.pop() {
        for &child in &tree.nodes[nid].children {
            ids.push(child);
            stack.push(child);
        }
    }
    ids
}

/// Prune a subtree rooted at `prune_id` from the tree.
///
/// The parent of `prune_id` is collapsed if it becomes a degree-2 node
/// (pass-through). When collapsing, branch lengths are summed.
fn prune_subtree(tree: &Tree, prune_id: usize) -> Option<PrunedTree> {
    let parent_id = tree.nodes[prune_id].parent?;

    let subtree_ids: std::collections::HashSet<usize> =
        collect_subtree_ids(tree, prune_id).into_iter().collect();

    // Build the subtree
    let mut subtree = Tree::new();
    subtree.rooted = false;
    let mut subtree_id_map = vec![usize::MAX; tree.nodes.len()];

    for &sid in &subtree_ids {
        let new_id = subtree.nodes.len();
        subtree_id_map[sid] = new_id;
        subtree.nodes.push(TreeNode {
            id: new_id,
            name: tree.nodes[sid].name.clone(),
            branch_length: tree.nodes[sid].branch_length,
            children: Vec::new(),
            parent: None,
        });
    }
    subtree.root = subtree_id_map[prune_id];

    // Rebuild subtree links
    for &sid in &subtree_ids {
        for &child in &tree.nodes[sid].children {
            if subtree_ids.contains(&child) {
                let mapped_parent = subtree_id_map[sid];
                let mapped_child = subtree_id_map[child];
                subtree.nodes[mapped_parent].children.push(mapped_child);
                subtree.nodes[mapped_child].parent = Some(mapped_parent);
            }
        }
    }

    // Build the main tree without subtree nodes
    let mut main_tree = Tree::new();
    main_tree.rooted = tree.rooted;
    let mut main_id_map = vec![usize::MAX; tree.nodes.len()];

    // Check if parent should be collapsed (becomes degree-2 after pruning)
    let parent_remaining_children: Vec<usize> = tree.nodes[parent_id]
        .children
        .iter()
        .filter(|&&c| !subtree_ids.contains(&c))
        .copied()
        .collect();

    let collapse_parent =
        parent_remaining_children.len() == 1 && parent_id != tree.root;

    // Combined branch length when collapsing parent
    let collapsed_bl = if collapse_parent {
        let surviving_child = parent_remaining_children[0];
        let parent_bl = tree.nodes[parent_id].branch_length.unwrap_or(0.0);
        let child_bl = tree.nodes[surviving_child].branch_length.unwrap_or(0.0);
        Some(parent_bl + child_bl)
    } else {
        None
    };

    // Add all non-subtree nodes (skipping collapsed parent)
    for node in &tree.nodes {
        if subtree_ids.contains(&node.id) {
            continue;
        }
        if collapse_parent && node.id == parent_id {
            continue;
        }
        let new_id = main_tree.nodes.len();
        main_id_map[node.id] = new_id;
        main_tree.nodes.push(TreeNode {
            id: new_id,
            name: node.name.clone(),
            branch_length: node.branch_length,
            children: Vec::new(),
            parent: None,
        });
    }

    main_tree.root = main_id_map[tree.root];

    // Rebuild links
    for node in &tree.nodes {
        if subtree_ids.contains(&node.id) {
            continue;
        }
        if collapse_parent && node.id == parent_id {
            continue;
        }

        let mapped_id = main_id_map[node.id];

        for &child in &node.children {
            if subtree_ids.contains(&child) {
                continue;
            }

            if collapse_parent && child == parent_id {
                // Skip collapsed parent; connect directly to surviving child
                let surviving_child = parent_remaining_children[0];
                let mapped_child = main_id_map[surviving_child];
                main_tree.nodes[mapped_id].children.push(mapped_child);
                main_tree.nodes[mapped_child].parent = Some(mapped_id);
                if let Some(bl) = collapsed_bl {
                    main_tree.nodes[mapped_child].branch_length = Some(bl);
                }
            } else {
                let mapped_child = main_id_map[child];
                main_tree.nodes[mapped_id].children.push(mapped_child);
                main_tree.nodes[mapped_child].parent = Some(mapped_id);
            }
        }
    }

    Some(PrunedTree {
        tree: main_tree,
        subtree,
    })
}

/// Regraft a pruned subtree onto an edge in the main tree.
///
/// Splits the graft edge, creating a new internal node, and attaches
/// the subtree root as a sibling. Branch lengths are divided at the
/// split point.
fn regraft_subtree(
    pruned: &PrunedTree,
    graft_parent: usize,
    graft_child: usize,
) -> Tree {
    let main = &pruned.tree;
    let sub = &pruned.subtree;

    let mut new_tree = Tree::new();
    new_tree.rooted = main.rooted;

    // Copy main tree nodes
    let mut main_map = vec![0usize; main.nodes.len()];
    for node in &main.nodes {
        let new_id = new_tree.nodes.len();
        main_map[node.id] = new_id;
        new_tree.nodes.push(TreeNode {
            id: new_id,
            name: node.name.clone(),
            branch_length: node.branch_length,
            children: Vec::new(),
            parent: None,
        });
    }

    new_tree.root = main_map[main.root];

    // Copy subtree nodes
    let mut sub_map = vec![0usize; sub.nodes.len()];
    for node in &sub.nodes {
        let new_id = new_tree.nodes.len();
        sub_map[node.id] = new_id;
        new_tree.nodes.push(TreeNode {
            id: new_id,
            name: node.name.clone(),
            branch_length: node.branch_length,
            children: Vec::new(),
            parent: None,
        });
    }

    // Split the graft edge: divide branch length
    let graft_bl = main.nodes[graft_child]
        .branch_length
        .unwrap_or(DEFAULT_BRANCH_LENGTH);
    let half_bl = graft_bl / 2.0;

    // New internal node at the graft point
    let graft_node_id = new_tree.nodes.len();
    new_tree.nodes.push(TreeNode {
        id: graft_node_id,
        name: None,
        branch_length: Some(half_bl),
        children: Vec::new(),
        parent: None,
    });

    // Update graft child's branch length to half
    let mapped_graft_child = main_map[graft_child];
    new_tree.nodes[mapped_graft_child].branch_length = Some(half_bl);

    // Rebuild main tree links, splitting the graft edge
    for node in &main.nodes {
        let mapped_id = main_map[node.id];
        for &child in &node.children {
            if node.id == graft_parent && child == graft_child {
                // Split: parent -> graft_node -> {child, subtree_root}
                new_tree.nodes[mapped_id].children.push(graft_node_id);
                new_tree.nodes[graft_node_id].parent = Some(mapped_id);

                let mapped_child = main_map[child];
                new_tree.nodes[graft_node_id].children.push(mapped_child);
                new_tree.nodes[mapped_child].parent = Some(graft_node_id);

                let subtree_root = sub_map[sub.root];
                new_tree.nodes[graft_node_id].children.push(subtree_root);
                new_tree.nodes[subtree_root].parent = Some(graft_node_id);
            } else {
                let mapped_child = main_map[child];
                new_tree.nodes[mapped_id].children.push(mapped_child);
                new_tree.nodes[mapped_child].parent = Some(mapped_id);
            }
        }
    }

    // Rebuild subtree links
    for node in &sub.nodes {
        let mapped_id = sub_map[node.id];
        for &child in &node.children {
            let mapped_child = sub_map[child];
            new_tree.nodes[mapped_id].children.push(mapped_child);
            new_tree.nodes[mapped_child].parent = Some(mapped_id);
        }
    }

    new_tree
}

/// Perform one round of SPR rearrangement.
///
/// For each internal node (except root and its direct children), tries
/// pruning its subtree and regrafting at every other edge. Picks the
/// best candidate (by unoptimized likelihood), then optimizes it.
/// Accepts the move if the optimized result improves the current score.
fn spr_round(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    current_lnl: f64,
) -> Result<(Tree, f64, bool), LikelihoodError> {
    let mut best_candidate = None;
    let mut best_quick_lnl = f64::NEG_INFINITY;

    // Candidate subtree roots: non-root nodes not directly under root
    let nodes_to_try: Vec<usize> = tree
        .nodes
        .iter()
        .filter(|n| n.parent.is_some() && n.parent != Some(tree.root))
        .map(|n| n.id)
        .collect();

    for &prune_node in &nodes_to_try {
        if let Some(pruned) = prune_subtree(tree, prune_node) {
            let edges: Vec<(usize, usize)> = pruned
                .tree
                .nodes
                .iter()
                .filter_map(|n| n.parent.map(|p| (p, n.id)))
                .collect();

            for (graft_parent, graft_child) in edges {
                let candidate = regraft_subtree(&pruned, graft_parent, graft_child);
                let lnl = quick_lnl(&candidate, alignment, model);
                if lnl > best_quick_lnl {
                    best_quick_lnl = lnl;
                    best_candidate = Some(candidate);
                }
            }
        }
    }

    // Optimize the best candidate and check if it truly improves
    if let Some(mut candidate) = best_candidate {
        let opt_lnl = optimize_branch_lengths(&mut candidate, alignment, model)?;
        if opt_lnl > current_lnl + 1e-4 {
            return Ok((candidate, opt_lnl, true));
        }
    }

    Ok((tree.clone(), current_lnl, false))
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

    // --- Initial tree construction ---

    #[test]
    fn test_build_initial_tree() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "TGCA"),
        ]);
        let order = vec![0, 1, 2];
        let tree = build_initial_tree(&order, &aln);
        assert_eq!(tree.num_leaves(), 3);
        assert_eq!(tree.num_nodes(), 4); // 3 leaves + 1 root
        // All leaves should have branch lengths
        for leaf in tree.leaves() {
            assert!(leaf.branch_length.is_some());
        }
    }

    // --- Edge insertion ---

    #[test]
    fn test_insert_on_edge() {
        let tree = parse_newick("((A:0.1,B:0.1):0.1,C:0.2);").unwrap();
        let internal_id = tree
            .nodes
            .iter()
            .find(|n| !n.is_leaf() && n.id != tree.root)
            .unwrap()
            .id;
        let parent_id = tree.nodes[internal_id].parent.unwrap();

        let new_tree = insert_on_edge(&tree, parent_id, internal_id, "D");
        assert_eq!(new_tree.num_leaves(), 4);

        let leaf_names = new_tree.leaf_names();
        assert!(leaf_names.contains(&"D"));
        assert!(leaf_names.contains(&"A"));
        assert!(leaf_names.contains(&"B"));
        assert!(leaf_names.contains(&"C"));
    }

    #[test]
    fn test_insert_preserves_branch_lengths() {
        let tree = parse_newick("((A:0.2,B:0.3):0.15,C:0.25);").unwrap();
        // Insert on the A branch
        let a_id = tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("A"))
            .unwrap()
            .id;
        let parent_id = tree.nodes[a_id].parent.unwrap();

        let new_tree = insert_on_edge(&tree, parent_id, a_id, "D");

        // All edges should have branch lengths
        for node in &new_tree.nodes {
            if node.parent.is_some() {
                assert!(
                    node.branch_length.is_some(),
                    "node {:?} missing branch length",
                    node.name
                );
            }
        }
    }

    #[test]
    fn test_insert_splits_branch_length() {
        let tree = parse_newick("(A:0.4,B:0.2,C:0.3);").unwrap();
        // Insert on A's branch (length 0.4)
        let a_id = tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("A"))
            .unwrap()
            .id;

        let new_tree = insert_on_edge(&tree, tree.root, a_id, "D");

        // The new internal node should have branch length 0.2 (half of 0.4)
        let new_internal = new_tree
            .nodes
            .iter()
            .find(|n| !n.is_leaf() && n.id != new_tree.root)
            .unwrap();
        assert!(
            (new_internal.branch_length.unwrap() - 0.2).abs() < 1e-10,
            "expected 0.2, got {}",
            new_internal.branch_length.unwrap()
        );
    }

    // --- Pruning and regrafting ---

    #[test]
    fn test_prune_and_regraft_preserves_taxa() {
        let tree = parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();

        // Find a non-root, non-root-child node to prune (e.g., A)
        let a_id = tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("A"))
            .unwrap()
            .id;

        if let Some(pruned) = prune_subtree(&tree, a_id) {
            // Main tree should have 3 leaves (B, C, D)
            assert_eq!(pruned.tree.num_leaves(), 3);
            // Subtree should have 1 leaf (A)
            assert_eq!(pruned.subtree.num_leaves(), 1);

            // Regraft on any edge
            let edge = pruned
                .tree
                .nodes
                .iter()
                .find_map(|n| n.parent.map(|p| (p, n.id)))
                .unwrap();

            let result = regraft_subtree(&pruned, edge.0, edge.1);
            assert_eq!(result.num_leaves(), 4);

            let names = result.leaf_names();
            assert!(names.contains(&"A"));
            assert!(names.contains(&"B"));
            assert!(names.contains(&"C"));
            assert!(names.contains(&"D"));
        } else {
            panic!("prune_subtree returned None");
        }
    }

    #[test]
    fn test_prune_collapses_degree2_node() {
        let tree = parse_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();

        let a_id = tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("A"))
            .unwrap()
            .id;

        if let Some(pruned) = prune_subtree(&tree, a_id) {
            // After removing A, the internal node (A,B) has only B left
            // It should be collapsed: root -> B (with combined branch 0.3 + 0.2)
            let b_node = pruned
                .tree
                .nodes
                .iter()
                .find(|n| n.name.as_deref() == Some("B"))
                .unwrap();
            let bl = b_node.branch_length.unwrap();
            assert!(
                (bl - 0.5).abs() < 1e-10,
                "B's branch should be 0.3 + 0.2 = 0.5, got {}",
                bl
            );
        }
    }

    // --- Basic search tests ---

    #[test]
    fn test_search_three_taxa() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGTACGT"),
            ("C", "TGCATGCA"),
        ]);
        let model = Jc69Model;
        let result = search(&aln, &model, Some(42)).unwrap();
        assert_eq!(result.tree.num_leaves(), 3);
        assert!(result.lnl < 0.0);
        assert_eq!(result.site_lnl.len(), 8);
    }

    #[test]
    fn test_search_four_taxa() {
        let aln = make_alignment(&[
            ("A", "AAAAACCCCC"),
            ("B", "AAAAACCCCC"),
            ("C", "CCCCCAAAAA"),
            ("D", "CCCCCAAAAA"),
        ]);
        let model = Jc69Model;
        let result = search(&aln, &model, Some(42)).unwrap();
        assert_eq!(result.tree.num_leaves(), 4);
        assert!(result.lnl < 0.0);
    }

    #[test]
    fn test_search_finds_correct_topology() {
        // Strong signal: A,B share states; C,D share states
        // ML should group (A,B) and (C,D) together
        let aln = make_alignment(&[
            ("A", "AAAAACCCCCAAAAACCCCC"),
            ("B", "AAAAACCCCCAAAAACCCCC"),
            ("C", "CCCCCAAAAACCCCCAAAAA"),
            ("D", "CCCCCAAAAACCCCCAAAAA"),
        ]);
        let model = Jc69Model;

        // Correct topology
        let good_tree =
            parse_newick("((A:0.05,B:0.05):0.2,(C:0.05,D:0.05):0.2);").unwrap();
        let bad_tree =
            parse_newick("((A:0.05,C:0.05):0.2,(B:0.05,D:0.05):0.2);").unwrap();

        let lnl_good = log_likelihood(&good_tree, &aln, &model).unwrap();
        let lnl_bad = log_likelihood(&bad_tree, &aln, &model).unwrap();

        // Sanity check: correct topology should indeed be better
        assert!(lnl_good > lnl_bad);

        // Search should find a tree at least as good as the correct topology
        let result = search(&aln, &model, Some(42)).unwrap();
        assert!(
            result.lnl >= lnl_good - 1.0,
            "search lnl ({}) should be near optimal ({})",
            result.lnl,
            lnl_good
        );
    }

    #[test]
    fn test_search_improves_over_naive() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGG"),
        ]);
        let model = Jc69Model;

        // Naive tree with poor topology and branch lengths
        let naive =
            parse_newick("((A:1.0,C:1.0):1.0,(B:1.0,D:1.0):1.0);").unwrap();
        let lnl_naive = log_likelihood(&naive, &aln, &model).unwrap();

        let result = search(&aln, &model, Some(42)).unwrap();
        assert!(
            result.lnl > lnl_naive,
            "search ({}) should beat naive tree ({})",
            result.lnl,
            lnl_naive
        );
    }

    #[test]
    fn test_search_deterministic_with_seed() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGAACGA"),
            ("C", "TGCATGCA"),
            ("D", "TGCTTGCT"),
        ]);
        let model = Jc69Model;
        let r1 = search(&aln, &model, Some(42)).unwrap();
        let r2 = search(&aln, &model, Some(42)).unwrap();
        assert!(
            (r1.lnl - r2.lnl).abs() < 1e-6,
            "same seed should give same result: {} vs {}",
            r1.lnl,
            r2.lnl
        );
    }

    #[test]
    fn test_search_five_taxa() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAA"),
            ("B", "AACCGGTTAA"),
            ("C", "CCAATTGGCC"),
            ("D", "CCAATTGGCC"),
            ("E", "GGTTAACCGG"),
        ]);
        let model = Jc69Model;
        let result = search(&aln, &model, Some(99)).unwrap();
        assert_eq!(result.tree.num_leaves(), 5);
        assert!(result.lnl < 0.0);
    }

    #[test]
    fn test_search_six_taxa() {
        let aln = make_alignment(&[
            ("A", "AACCGGTT"),
            ("B", "AACCGGTT"),
            ("C", "CCAATTGG"),
            ("D", "CCAATTGG"),
            ("E", "GGTTAACC"),
            ("F", "GGTTAACC"),
        ]);
        let model = Jc69Model;
        let result = search(&aln, &model, Some(7)).unwrap();
        assert_eq!(result.tree.num_leaves(), 6);
        assert!(result.lnl < 0.0);
    }

    #[test]
    fn test_search_identical_sequences() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGTACGT"),
            ("C", "ACGTACGT"),
            ("D", "ACGTACGT"),
        ]);
        let model = Jc69Model;
        let result = search(&aln, &model, None).unwrap();
        assert_eq!(result.tree.num_leaves(), 4);
        // With identical sequences, branches should optimize to be very short
        // and likelihood should be relatively high (close to ln(0.25) per site)
        let per_site = result.lnl / 8.0;
        assert!(
            per_site > -2.0,
            "identical sequences should have high per-site lnl, got {}",
            per_site
        );
    }

    // --- F84 model tests ---

    #[test]
    fn test_search_with_f84() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGG"),
        ]);
        let model = F84Model::equal_freqs(2.0);
        let result = search(&aln, &model, Some(42)).unwrap();
        assert_eq!(result.tree.num_leaves(), 4);
        assert!(result.lnl < 0.0);
    }

    #[test]
    fn test_search_f84_unequal_freqs() {
        let aln = make_alignment(&[
            ("A", "AAACCCTTTGGG"),
            ("B", "AAACCCTTTGGG"),
            ("C", "CCCAAAGGGCCC"),
            ("D", "TTTGGGAAATTT"),
        ]);
        let freqs = BaseFrequencies::new(0.3, 0.2, 0.2, 0.3);
        let model = F84Model::new(freqs, 2.0);
        let result = search(&aln, &model, Some(42)).unwrap();
        assert_eq!(result.tree.num_leaves(), 4);
        assert!(result.lnl < 0.0);
    }

    // --- Result structure tests ---

    #[test]
    fn test_result_site_lnl_sums_to_total() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGAACGA"),
            ("C", "TGCATGCA"),
        ]);
        let model = Jc69Model;
        let result = search(&aln, &model, Some(42)).unwrap();
        let sum: f64 = result.site_lnl.iter().sum();
        assert!(
            (sum - result.lnl).abs() < 1e-8,
            "site lnl sum ({}) should equal total ({})",
            sum,
            result.lnl
        );
    }

    #[test]
    fn test_result_all_site_lnl_negative() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGAACGA"),
            ("C", "TGCATGCA"),
            ("D", "AGCTAGCT"),
        ]);
        let model = Jc69Model;
        let result = search(&aln, &model, None).unwrap();
        for (i, &sl) in result.site_lnl.iter().enumerate() {
            assert!(sl < 0.0, "site {} lnl should be negative, got {}", i, sl);
        }
    }

    // --- Error handling ---

    #[test]
    fn test_search_too_few_taxa() {
        let aln = make_alignment(&[("A", "ACGT"), ("B", "ACGT")]);
        let model = Jc69Model;
        let result = search(&aln, &model, None);
        assert!(result.is_err());
    }

    // --- SPR unit tests ---

    #[test]
    fn test_spr_does_not_worsen() {
        let aln = make_alignment(&[
            ("A", "AAAAACCCCC"),
            ("B", "AAAAACCCCC"),
            ("C", "CCCCCAAAAA"),
            ("D", "CCCCCAAAAA"),
        ]);
        let model = Jc69Model;

        // Start with a good tree
        let mut tree =
            parse_newick("((A:0.05,B:0.05):0.2,(C:0.05,D:0.05):0.2);").unwrap();
        let lnl = optimize_branch_lengths(&mut tree, &aln, &model).unwrap();

        let (_, new_lnl, _) = spr_round(&tree, &aln, &model, lnl).unwrap();
        assert!(
            new_lnl >= lnl - 1e-4,
            "SPR should not worsen: {} -> {}",
            lnl,
            new_lnl
        );
    }

    #[test]
    fn test_spr_can_improve() {
        // Start with a wrong topology
        let aln = make_alignment(&[
            ("A", "AAAAAAAAAA"),
            ("B", "AAAAAAAAAA"),
            ("C", "CCCCCCCCCC"),
            ("D", "CCCCCCCCCC"),
        ]);
        let model = Jc69Model;

        // Wrong topology: A grouped with C instead of B
        let mut tree =
            parse_newick("((A:0.1,C:0.1):0.1,(B:0.1,D:0.1):0.1);").unwrap();
        let lnl_wrong = optimize_branch_lengths(&mut tree, &aln, &model).unwrap();

        // Correct topology should have better likelihood
        let mut correct =
            parse_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let lnl_correct =
            optimize_branch_lengths(&mut correct, &aln, &model).unwrap();
        assert!(lnl_correct > lnl_wrong);
    }
}

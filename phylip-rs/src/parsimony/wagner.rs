//! Wagner parsimony for DNA sequences.
//!
//! Implements the Fitch (1971) algorithm for computing minimum-change
//! parsimony scores, and a heuristic tree search via stepwise addition
//! with SPR (subtree pruning and regrafting) rearrangement.
//!
//! The Fitch algorithm works site-by-site:
//! - At each leaf: the state set is the observed nucleotide
//! - At each internal node: take the intersection of children's sets.
//!   If the intersection is empty, take the union and add one step.
//!
//! This gives the minimum number of substitutions required to explain
//! the data on a given tree topology.
//!
//! # Example
//!
//! ```
//! use phylip_rs::parsimony::wagner::{parsimony_score, search, ParsimonyResult};
//! use phylip_rs::tree::{Alignment, Base, Sequence, Tree};
//! use phylip_rs::tree::newick::parse_newick;
//!
//! let alignment = Alignment::new(vec![
//!     Sequence::new("A", vec![Base::A, Base::C, Base::G, Base::T]),
//!     Sequence::new("B", vec![Base::A, Base::C, Base::G, Base::T]),
//!     Sequence::new("C", vec![Base::A, Base::C, Base::A, Base::T]),
//!     Sequence::new("D", vec![Base::A, Base::T, Base::A, Base::C]),
//! ]).unwrap();
//!
//! let result = search(&alignment, None);
//! assert!(result.score > 0);
//! assert_eq!(result.tree.num_leaves(), 4);
//! ```
//!
//! # References
//!
//! - Fitch, W.M. (1971). Toward defining the course of evolution: minimum
//!   change for a specific tree topology. *Systematic Zoology*, 20, 406-416.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//!   Chapter 7: Counting evolutionary changes on trees.

use std::fmt;

use crate::tree::types::{Alignment, Base, Tree};

/// A set of nucleotide states, represented as a 5-bit mask.
/// Bit 0 = A, 1 = C, 2 = G, 3 = T, 4 = gap (as 5th state).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StateSet(u8);

impl StateSet {
    pub const EMPTY: StateSet = StateSet(0);
    pub const A: StateSet = StateSet(1);
    pub const C: StateSet = StateSet(2);
    pub const G: StateSet = StateSet(4);
    pub const T: StateSet = StateSet(8);
    pub const GAP: StateSet = StateSet(16);
    pub const ALL: StateSet = StateSet(0b11111);

    /// Create a state set from a `Base`.
    pub fn from_base(base: Base) -> StateSet {
        match base {
            Base::A => StateSet::A,
            Base::C => StateSet::C,
            Base::G => StateSet::G,
            Base::T => StateSet::T,
            Base::Gap => StateSet::GAP,
            // Ambiguity codes map to unions of possible states
            Base::M => StateSet(0b00011), // A or C
            Base::R => StateSet(0b00101), // A or G
            Base::W => StateSet(0b01001), // A or T
            Base::S => StateSet(0b00110), // C or G
            Base::Y => StateSet(0b01010), // C or T
            Base::K => StateSet(0b01100), // G or T
            Base::V => StateSet(0b00111), // A, C, G
            Base::H => StateSet(0b01011), // A, C, T
            Base::D => StateSet(0b01101), // A, G, T
            Base::B => StateSet(0b01110), // C, G, T
            Base::N => StateSet(0b01111), // A, C, G, T
        }
    }

    /// Intersection of two state sets.
    pub fn intersection(self, other: StateSet) -> StateSet {
        StateSet(self.0 & other.0)
    }

    /// Union of two state sets.
    pub fn union(self, other: StateSet) -> StateSet {
        StateSet(self.0 | other.0)
    }

    /// True if this set is empty (no states).
    pub fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// Number of states in the set.
    pub fn count(self) -> u32 {
        self.0.count_ones()
    }

    /// Returns one representative base from the set (lowest-bit).
    pub fn representative(self) -> Option<Base> {
        if self.0 & 1 != 0 {
            Some(Base::A)
        } else if self.0 & 2 != 0 {
            Some(Base::C)
        } else if self.0 & 4 != 0 {
            Some(Base::G)
        } else if self.0 & 8 != 0 {
            Some(Base::T)
        } else if self.0 & 16 != 0 {
            Some(Base::Gap)
        } else {
            None
        }
    }
}

impl fmt::Display for StateSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{{")?;
        let mut first = true;
        for (bit, ch) in [(1, 'A'), (2, 'C'), (4, 'G'), (8, 'T'), (16, '-')] {
            if self.0 & bit != 0 {
                if !first {
                    write!(f, ",")?;
                }
                write!(f, "{}", ch)?;
                first = false;
            }
        }
        write!(f, "}}")
    }
}

/// Per-node Fitch algorithm data for all sites.
#[derive(Debug, Clone)]
struct NodeData {
    /// State set at this node for each site.
    states: Vec<StateSet>,
    /// Accumulated step count below this node for each site.
    steps: Vec<usize>,
}

impl NodeData {
    fn new(nsites: usize) -> Self {
        Self {
            states: vec![StateSet::EMPTY; nsites],
            steps: vec![0; nsites],
        }
    }
}

/// Result of a parsimony analysis.
#[derive(Debug, Clone)]
pub struct ParsimonyResult {
    /// The best tree found.
    pub tree: Tree,
    /// The parsimony score (total number of steps).
    pub score: usize,
    /// Per-site step counts.
    pub site_steps: Vec<usize>,
}

/// Errors that can occur during parsimony computation.
#[derive(Debug, Clone)]
pub enum ParsimonyError {
    /// Fewer than 3 taxa in alignment.
    TooFewTaxa(usize),
    /// Tree leaf names don't match alignment taxa.
    TaxaMismatch,
    /// Tree has no leaves.
    EmptyTree,
}

impl fmt::Display for ParsimonyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParsimonyError::TooFewTaxa(n) => {
                write!(f, "need at least 3 taxa for parsimony, got {}", n)
            }
            ParsimonyError::TaxaMismatch => {
                write!(f, "tree leaf names do not match alignment taxa")
            }
            ParsimonyError::EmptyTree => write!(f, "tree has no leaves"),
        }
    }
}

impl std::error::Error for ParsimonyError {}

/// Compute the parsimony score (total steps) for a given tree and alignment.
///
/// Uses the Fitch (1971) algorithm: postorder traversal computing
/// state-set intersections/unions at each internal node.
///
/// Returns the total score and per-site step counts.
pub fn parsimony_score(
    tree: &Tree,
    alignment: &Alignment,
) -> Result<(usize, Vec<usize>), ParsimonyError> {
    let nsites = alignment.nsites();
    let nleaves = tree.num_leaves();
    if nleaves == 0 {
        return Err(ParsimonyError::EmptyTree);
    }

    // Build name -> alignment index map
    let name_to_idx: std::collections::HashMap<&str, usize> = alignment
        .sequences
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.as_str(), i))
        .collect();

    // Initialize node data
    let nnodes = tree.num_nodes();
    let mut data = vec![NodeData::new(nsites); nnodes];

    // Set leaf states from alignment
    for node in &tree.nodes {
        if node.is_leaf() {
            if let Some(name) = &node.name {
                let taxon_idx = name_to_idx
                    .get(name.as_str())
                    .ok_or(ParsimonyError::TaxaMismatch)?;
                for site in 0..nsites {
                    data[node.id].states[site] =
                        StateSet::from_base(alignment.base(*taxon_idx, site));
                }
            }
        }
    }

    // Postorder traversal: Fitch algorithm
    let postorder = tree.postorder();
    for &node_id in &postorder {
        let children = &tree.nodes[node_id].children;
        if children.is_empty() {
            continue; // leaf: already initialized
        }

        // Combine all children using Fitch's rule
        let first_child = children[0];
        for site in 0..nsites {
            let mut current_states = data[first_child].states[site];
            let mut current_steps: usize =
                children.iter().map(|&c| data[c].steps[site]).sum();

            for &child in &children[1..] {
                let child_states = data[child].states[site];
                let isect = current_states.intersection(child_states);
                if isect.is_empty() {
                    // Union: add a step
                    current_states = current_states.union(child_states);
                    current_steps += 1;
                } else {
                    // Intersection: no extra step
                    current_states = isect;
                }
            }

            data[node_id].states[site] = current_states;
            data[node_id].steps[site] = current_steps;
        }
    }

    // Total score from root
    let root = tree.root;
    let site_steps = data[root].steps.clone();
    let total: usize = site_steps.iter().sum();

    Ok((total, site_steps))
}

/// Reconstruct ancestral states at internal nodes using the Fitch method.
///
/// First pass (postorder): compute state sets (already done by parsimony_score).
/// Second pass (preorder): resolve ambiguous state sets using parent's state.
///
/// Returns a vector indexed by node id, each containing the reconstructed
/// state for each site.
pub fn reconstruct_ancestors(
    tree: &Tree,
    alignment: &Alignment,
) -> Result<Vec<Vec<StateSet>>, ParsimonyError> {
    let nsites = alignment.nsites();
    let nnodes = tree.num_nodes();

    // Build name -> alignment index map
    let name_to_idx: std::collections::HashMap<&str, usize> = alignment
        .sequences
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.as_str(), i))
        .collect();

    // Initialize node data
    let mut data = vec![NodeData::new(nsites); nnodes];

    // Set leaf states
    for node in &tree.nodes {
        if node.is_leaf() {
            if let Some(name) = &node.name {
                let taxon_idx = name_to_idx
                    .get(name.as_str())
                    .ok_or(ParsimonyError::TaxaMismatch)?;
                for site in 0..nsites {
                    data[node.id].states[site] =
                        StateSet::from_base(alignment.base(*taxon_idx, site));
                }
            }
        }
    }

    // Postorder: compute state sets
    let postorder = tree.postorder();
    for &node_id in &postorder {
        let children = &tree.nodes[node_id].children;
        if children.is_empty() {
            continue;
        }

        let first_child = children[0];
        for site in 0..nsites {
            let mut current = data[first_child].states[site];
            for &child in &children[1..] {
                let child_states = data[child].states[site];
                let isect = current.intersection(child_states);
                if isect.is_empty() {
                    current = current.union(child_states);
                } else {
                    current = isect;
                }
            }
            data[node_id].states[site] = current;
        }
    }

    // Preorder: resolve ambiguities using parent's assigned state
    // The root keeps its full state set; we pick a representative.
    let mut resolved = vec![vec![StateSet::EMPTY; nsites]; nnodes];

    let preorder = tree.preorder();
    for &node_id in &preorder {
        for site in 0..nsites {
            let node_states = data[node_id].states[site];

            if tree.nodes[node_id].parent.is_none() {
                // Root: keep the postorder state set
                resolved[node_id][site] = node_states;
            } else {
                let parent_id = tree.nodes[node_id].parent.unwrap();
                let parent_state = resolved[parent_id][site];
                // If parent's state is in our set, use it (most parsimonious)
                let isect = node_states.intersection(parent_state);
                if !isect.is_empty() {
                    resolved[node_id][site] = isect;
                } else {
                    // Parent's state not in our set: keep our postorder set
                    resolved[node_id][site] = node_states;
                }
            }
        }
    }

    Ok(resolved)
}

// ============================================================
// Tree search: stepwise addition + SPR rearrangement
// ============================================================

/// Simple LCG random number generator (same as bootstrap module).
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
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.state
    }

    /// Random index in [0, n)
    fn gen_range(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    /// Shuffle a slice in place (Fisher-Yates).
    fn shuffle<T>(&mut self, slice: &mut [T]) {
        for i in (1..slice.len()).rev() {
            let j = self.gen_range(i + 1);
            slice.swap(i, j);
        }
    }
}

/// Build a 3-taxon star tree from the first 3 taxa.
fn build_initial_tree(taxa: &[usize], alignment: &Alignment) -> Tree {
    let mut tree = Tree::new();
    tree.rooted = false;

    // Root (internal node)
    let root = tree.add_node(None, None, None);
    tree.root = root;

    // Add first 3 taxa as children of root
    for &taxon_idx in &taxa[..3] {
        let name = alignment.taxon_name(taxon_idx).to_string();
        tree.add_node(Some(name), None, Some(root));
    }

    tree
}

/// Score a tree quickly using the Fitch algorithm.
/// Returns the total parsimony score.
fn quick_score(tree: &Tree, alignment: &Alignment) -> usize {
    match parsimony_score(tree, alignment) {
        Ok((score, _)) => score,
        Err(_) => usize::MAX,
    }
}

/// Try inserting a new taxon at every possible edge in the tree.
/// Returns the tree with the best (lowest) parsimony score.
fn best_insertion(
    tree: &Tree,
    taxon_name: &str,
    alignment: &Alignment,
) -> Tree {
    let mut best_tree = None;
    let mut best_score = usize::MAX;

    // Collect all edges: (parent_id, child_id)
    let edges: Vec<(usize, usize)> = tree
        .nodes
        .iter()
        .filter_map(|n| n.parent.map(|p| (p, n.id)))
        .collect();

    for (parent_id, child_id) in edges {
        let candidate = insert_on_edge(tree, parent_id, child_id, taxon_name);
        let score = quick_score(&candidate, alignment);
        if score < best_score {
            best_score = score;
            best_tree = Some(candidate);
        }
    }

    best_tree.unwrap_or_else(|| {
        // Fallback: add as child of root
        let mut t = tree.clone();
        let root = t.root;
        t.add_node(Some(taxon_name.to_string()), None, Some(root));
        t
    })
}

/// Insert a new taxon on the edge between parent_id and child_id.
/// Creates a new internal node that becomes parent of both child_id and the new taxon.
fn insert_on_edge(
    tree: &Tree,
    parent_id: usize,
    child_id: usize,
    taxon_name: &str,
) -> Tree {
    let mut new_tree = Tree::new();
    new_tree.rooted = tree.rooted;

    // Map old node IDs to new node IDs
    let mut id_map = vec![0usize; tree.nodes.len()];

    // First pass: copy all existing nodes (without children/parent links)
    for node in &tree.nodes {
        let new_id = new_tree.nodes.len();
        id_map[node.id] = new_id;
        new_tree.nodes.push(crate::tree::types::TreeNode {
            id: new_id,
            name: node.name.clone(),
            branch_length: node.branch_length,
            children: Vec::new(),
            parent: None,
        });
    }

    new_tree.root = id_map[tree.root];

    // Add the new internal node
    let new_internal_id = new_tree.nodes.len();
    new_tree.nodes.push(crate::tree::types::TreeNode {
        id: new_internal_id,
        name: None,
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    // Add the new leaf
    let new_leaf_id = new_tree.nodes.len();
    new_tree.nodes.push(crate::tree::types::TreeNode {
        id: new_leaf_id,
        name: Some(taxon_name.to_string()),
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    // Rebuild parent-child links
    for node in &tree.nodes {
        let mapped_id = id_map[node.id];
        for &old_child in &node.children {
            if node.id == parent_id && old_child == child_id {
                // This is the edge we're splitting:
                // parent -> new_internal -> {child, new_leaf}
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

/// Perform one round of SPR rearrangement on the tree.
///
/// For each internal node (except root), try pruning its subtree
/// and regrafting it at every other edge. If any move improves the score,
/// accept the best one and restart.
///
/// Returns the improved tree and whether any improvement was made.
fn spr_round(tree: &Tree, alignment: &Alignment, current_score: usize) -> (Tree, usize, bool) {
    let mut best_tree = tree.clone();
    let mut best_score = current_score;
    let mut improved = false;

    // Collect candidate subtree roots: all non-root nodes that are not
    // direct children of root (pruning a direct child of a trifurcating
    // root is handled specially).
    let nodes_to_try: Vec<usize> = tree
        .nodes
        .iter()
        .filter(|n| {
            n.parent.is_some() && n.parent != Some(tree.root)
        })
        .map(|n| n.id)
        .collect();

    for &prune_node in &nodes_to_try {
        // Try pruning this subtree and regrafting elsewhere
        if let Some(pruned) = prune_subtree(tree, prune_node) {
            // Collect regraft edges in the pruned tree
            let edges: Vec<(usize, usize)> = pruned
                .tree
                .nodes
                .iter()
                .filter_map(|n| n.parent.map(|p| (p, n.id)))
                .collect();

            for (graft_parent, graft_child) in edges {
                let candidate = regraft_subtree(
                    &pruned,
                    graft_parent,
                    graft_child,
                );
                let score = quick_score(&candidate, alignment);
                if score < best_score {
                    best_score = score;
                    best_tree = candidate;
                    improved = true;
                }
            }
        }
    }

    (best_tree, best_score, improved)
}

/// Result of pruning a subtree from a tree.
#[derive(Debug, Clone)]
struct PrunedTree {
    /// The tree with the subtree removed.
    tree: Tree,
    /// The pruned subtree (as a separate tree).
    subtree: Tree,
    /// Mapping from original node IDs to new IDs in the main tree.
    /// Nodes in the subtree map to usize::MAX.
    _id_map: Vec<usize>,
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
/// The parent of `prune_id` is collapsed if it becomes a degree-2 node.
fn prune_subtree(tree: &Tree, prune_id: usize) -> Option<PrunedTree> {
    let parent_id = tree.nodes[prune_id].parent?;

    // Collect subtree node IDs
    let subtree_ids: std::collections::HashSet<usize> =
        collect_subtree_ids(tree, prune_id).into_iter().collect();

    // Build the subtree
    let mut subtree = Tree::new();
    subtree.rooted = false;
    let mut subtree_id_map = vec![usize::MAX; tree.nodes.len()];

    // Add subtree nodes
    for &sid in &subtree_ids {
        let new_id = subtree.nodes.len();
        subtree_id_map[sid] = new_id;
        subtree.nodes.push(crate::tree::types::TreeNode {
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

    // Determine which nodes to skip: subtree nodes, and also the parent
    // if it would become a pass-through (degree-2) node.
    let parent_remaining_children: Vec<usize> = tree.nodes[parent_id]
        .children
        .iter()
        .filter(|&&c| !subtree_ids.contains(&c))
        .copied()
        .collect();

    let collapse_parent = parent_remaining_children.len() == 1
        && parent_id != tree.root;

    // Add all non-subtree nodes (possibly skipping collapsed parent)
    for node in &tree.nodes {
        if subtree_ids.contains(&node.id) {
            continue;
        }
        if collapse_parent && node.id == parent_id {
            continue;
        }
        let new_id = main_tree.nodes.len();
        main_id_map[node.id] = new_id;
        main_tree.nodes.push(crate::tree::types::TreeNode {
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
                continue; // pruned away
            }

            if collapse_parent && child == parent_id {
                // Skip the collapsed parent; connect directly to its remaining child
                let surviving_child = parent_remaining_children[0];
                let mapped_child = main_id_map[surviving_child];
                main_tree.nodes[mapped_id].children.push(mapped_child);
                main_tree.nodes[mapped_child].parent = Some(mapped_id);
            } else if collapse_parent && node.id == parent_id {
                // This shouldn't happen since we skip parent_id
                continue;
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
        _id_map: main_id_map,
    })
}

/// Regraft a pruned subtree onto an edge in the main tree.
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
        new_tree.nodes.push(crate::tree::types::TreeNode {
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
        new_tree.nodes.push(crate::tree::types::TreeNode {
            id: new_id,
            name: node.name.clone(),
            branch_length: node.branch_length,
            children: Vec::new(),
            parent: None,
        });
    }

    // Add new internal node for the graft point
    let graft_node_id = new_tree.nodes.len();
    new_tree.nodes.push(crate::tree::types::TreeNode {
        id: graft_node_id,
        name: None,
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    // Rebuild main tree links, splitting the graft edge
    for node in &main.nodes {
        let mapped_id = main_map[node.id];
        for &child in &node.children {
            if node.id == graft_parent && child == graft_child {
                // Split this edge: parent -> graft_node -> {child, subtree_root}
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

/// Search for the most parsimonious tree using stepwise addition
/// followed by SPR rearrangement.
///
/// # Arguments
/// * `alignment` — the DNA alignment to analyze
/// * `seed` — optional random seed for taxon addition order. If None,
///   taxa are added in alignment order.
///
/// # Returns
/// A `ParsimonyResult` with the best tree found and its score.
pub fn search(alignment: &Alignment, seed: Option<u64>) -> ParsimonyResult {
    let ntaxa = alignment.ntaxa();
    assert!(ntaxa >= 3, "need at least 3 taxa for tree search");

    // Taxon addition order
    let mut order: Vec<usize> = (0..ntaxa).collect();
    if let Some(s) = seed {
        let mut rng = SimpleRng::new(s);
        rng.shuffle(&mut order);
    }

    // Phase 1: Stepwise addition
    let mut tree = build_initial_tree(&order, alignment);

    for i in 3..ntaxa {
        let taxon_name = alignment.taxon_name(order[i]);
        tree = best_insertion(&tree, taxon_name, alignment);
    }

    let mut score = quick_score(&tree, alignment);

    // Phase 2: SPR rearrangement (iterate until no improvement)
    loop {
        let (new_tree, new_score, improved) = spr_round(&tree, alignment, score);
        if !improved {
            break;
        }
        tree = new_tree;
        score = new_score;
    }

    // Get final per-site scores
    let (final_score, site_steps) = parsimony_score(&tree, alignment)
        .expect("final scoring should not fail");

    ParsimonyResult {
        tree,
        score: final_score,
        site_steps,
    }
}

/// Evaluate a user-provided tree on the given alignment.
///
/// Returns a `ParsimonyResult` with the tree and its parsimony score.
pub fn evaluate(
    tree: &Tree,
    alignment: &Alignment,
) -> Result<ParsimonyResult, ParsimonyError> {
    let (score, site_steps) = parsimony_score(tree, alignment)?;
    Ok(ParsimonyResult {
        tree: tree.clone(),
        score,
        site_steps,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
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

    // --- StateSet tests ---

    #[test]
    fn test_state_set_from_base() {
        assert_eq!(StateSet::from_base(Base::A), StateSet::A);
        assert_eq!(StateSet::from_base(Base::C), StateSet::C);
        assert_eq!(StateSet::from_base(Base::G), StateSet::G);
        assert_eq!(StateSet::from_base(Base::T), StateSet::T);
        assert_eq!(StateSet::from_base(Base::Gap), StateSet::GAP);
    }

    #[test]
    fn test_state_set_ambiguity() {
        // M = A or C
        let m = StateSet::from_base(Base::M);
        assert!(!m.intersection(StateSet::A).is_empty());
        assert!(!m.intersection(StateSet::C).is_empty());
        assert!(m.intersection(StateSet::G).is_empty());

        // N = A, C, G, T
        let n = StateSet::from_base(Base::N);
        assert_eq!(n.count(), 4);
        assert!(!n.intersection(StateSet::A).is_empty());
        assert!(!n.intersection(StateSet::T).is_empty());
    }

    #[test]
    fn test_state_set_intersection_union() {
        let a = StateSet::A;
        let c = StateSet::C;
        let ac = a.union(c);

        assert!(a.intersection(c).is_empty());
        assert_eq!(ac.count(), 2);
        assert_eq!(ac.intersection(a), a);
        assert_eq!(ac.intersection(c), c);
    }

    #[test]
    fn test_state_set_representative() {
        assert_eq!(StateSet::A.representative(), Some(Base::A));
        assert_eq!(StateSet::C.representative(), Some(Base::C));
        assert_eq!(StateSet::G.representative(), Some(Base::G));
        assert_eq!(StateSet::T.representative(), Some(Base::T));
        assert_eq!(StateSet::EMPTY.representative(), None);
        // Union should return lowest-bit
        assert_eq!(StateSet::A.union(StateSet::C).representative(), Some(Base::A));
    }

    #[test]
    fn test_state_set_display() {
        assert_eq!(format!("{}", StateSet::A), "{A}");
        assert_eq!(format!("{}", StateSet::A.union(StateSet::T)), "{A,T}");
        assert_eq!(format!("{}", StateSet::EMPTY), "{}");
    }

    // --- Parsimony score tests ---

    #[test]
    fn test_identical_sequences_zero_steps() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "ACGT"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let (score, site_steps) = parsimony_score(&tree, &aln).unwrap();
        assert_eq!(score, 0);
        assert_eq!(site_steps, vec![0, 0, 0, 0]);
    }

    #[test]
    fn test_one_change_one_site() {
        // A and B share 'A' at site 0, C has 'T' -> 1 step
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "A"),
            ("C", "T"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let (score, _) = parsimony_score(&tree, &aln).unwrap();
        assert_eq!(score, 1);
    }

    #[test]
    fn test_all_different_three_taxa() {
        // A=A, B=C, C=G at one site: on ((A,B),C), need 2 steps
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "C"),
            ("C", "G"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let (score, _) = parsimony_score(&tree, &aln).unwrap();
        assert_eq!(score, 2);
    }

    #[test]
    fn test_four_taxa_parsimony() {
        // ((A,B),(C,D)) with data:
        // Site 0: A=A, B=A, C=C, D=C -> {A}&{A}={A}, {C}&{C}={C}, {A}&{C}=empty -> 1 step
        // Site 1: all same -> 0 steps
        let aln = make_alignment(&[
            ("A", "AG"),
            ("B", "AG"),
            ("C", "CG"),
            ("D", "CG"),
        ]);
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let (score, site_steps) = parsimony_score(&tree, &aln).unwrap();
        assert_eq!(score, 1);
        assert_eq!(site_steps, vec![1, 0]);
    }

    #[test]
    fn test_four_taxa_two_changes() {
        // ((A,B),(C,D)) with A=A, B=C, C=G, D=T at one site
        // {A}&{C}=empty -> union {A,C}, step=1
        // {G}&{T}=empty -> union {G,T}, step=1
        // {A,C}&{G,T}=empty -> union {A,C,G,T}, step=1
        // Total = 3 steps
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "C"),
            ("C", "G"),
            ("D", "T"),
        ]);
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let (score, _) = parsimony_score(&tree, &aln).unwrap();
        assert_eq!(score, 3);
    }

    #[test]
    fn test_topology_matters() {
        // Data: A=A, B=A, C=T, D=T
        // On ((A,B),(C,D)): 1 step
        // On ((A,C),(B,D)): 2 steps
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "A"),
            ("C", "T"),
            ("D", "T"),
        ]);

        let tree1 = parse_newick("((A,B),(C,D));").unwrap();
        let tree2 = parse_newick("((A,C),(B,D));").unwrap();

        let (s1, _) = parsimony_score(&tree1, &aln).unwrap();
        let (s2, _) = parsimony_score(&tree2, &aln).unwrap();

        assert_eq!(s1, 1);
        assert_eq!(s2, 2);
    }

    #[test]
    fn test_ambiguity_codes_in_parsimony() {
        // M = A or C. If A has M and B has A, intersection = {A} -> 0 steps
        let aln = make_alignment(&[
            ("A", "M"),
            ("B", "A"),
            ("C", "A"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let (score, _) = parsimony_score(&tree, &aln).unwrap();
        assert_eq!(score, 0);
    }

    #[test]
    fn test_gap_as_fifth_state() {
        // Gap is treated as a 5th state, different from ACGT
        let aln = make_alignment(&[
            ("A", "-"),
            ("B", "-"),
            ("C", "A"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let (score, _) = parsimony_score(&tree, &aln).unwrap();
        assert_eq!(score, 1);
    }

    #[test]
    fn test_taxa_mismatch_error() {
        let aln = make_alignment(&[
            ("X", "A"),
            ("Y", "A"),
            ("Z", "A"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let result = parsimony_score(&tree, &aln);
        assert!(result.is_err());
    }

    #[test]
    fn test_multifurcating_tree() {
        // Star tree: (A,B,C,D) — all children of root
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "A"),
            ("C", "T"),
            ("D", "T"),
        ]);
        let tree = parse_newick("(A,B,C,D);").unwrap();
        let (score, _) = parsimony_score(&tree, &aln).unwrap();
        // On a star: {A}&{A}={A}, then {A}&{T}=empty->union, +1,
        //   then {A,T}&{T}={T}, no extra step -> total = 1
        assert_eq!(score, 1);
    }

    #[test]
    fn test_multisite_scoring() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGTACGT"),
            ("C", "TGCATGCA"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let (score, site_steps) = parsimony_score(&tree, &aln).unwrap();
        // Each site where A/B differ from C costs 1 step
        // A=A,B=A,C=T -> 1; A=C,B=C,C=G -> 1; etc.
        assert_eq!(site_steps.len(), 8);
        assert_eq!(score, site_steps.iter().sum::<usize>());
    }

    // --- Ancestral reconstruction tests ---

    #[test]
    fn test_reconstruct_identical() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "ACGT"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let ancestors = reconstruct_ancestors(&tree, &aln).unwrap();

        // All internal nodes should have definite states matching the input
        let root = tree.root;
        for site in 0..4 {
            let root_state = ancestors[root][site];
            assert_eq!(root_state.count(), 1);
        }
    }

    #[test]
    fn test_reconstruct_resolves_with_parent() {
        // A=A, B=C, C=A: on ((A,B),C)
        // Internal (A,B): intersection of {A} and {C} is empty -> union {A,C}
        // Root: intersection of {A,C} and {A} = {A}
        // Preorder: root gets {A}, internal gets {A,C} & parent {A} = {A}
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "C"),
            ("C", "A"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let ancestors = reconstruct_ancestors(&tree, &aln).unwrap();

        let root = tree.root;
        assert_eq!(ancestors[root][0], StateSet::A);
    }

    // --- Tree search tests ---

    #[test]
    fn test_search_three_taxa() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "TGCA"),
        ]);
        let result = search(&aln, Some(42));
        assert_eq!(result.tree.num_leaves(), 3);
        assert!(result.score > 0);
    }

    #[test]
    fn test_search_four_taxa() {
        let aln = make_alignment(&[
            ("A", "AAAA"),
            ("B", "AAAA"),
            ("C", "TTTT"),
            ("D", "TTTT"),
        ]);
        let result = search(&aln, Some(42));
        assert_eq!(result.tree.num_leaves(), 4);
        // Best tree groups A,B together and C,D together -> 1 step per site
        assert_eq!(result.score, 4);
    }

    #[test]
    fn test_search_finds_optimal() {
        // Clear phylogenetic signal: (A,B) and (C,D) share states
        let aln = make_alignment(&[
            ("A", "AAAACCCC"),
            ("B", "AAAACCCC"),
            ("C", "CCCCAAAA"),
            ("D", "CCCCAAAA"),
        ]);
        let result = search(&aln, Some(123));
        // Optimal score: 1 change per site = 8
        assert_eq!(result.score, 8);
    }

    #[test]
    fn test_search_five_taxa() {
        let aln = make_alignment(&[
            ("A", "AACCGG"),
            ("B", "AACCGG"),
            ("C", "CCAACC"),
            ("D", "CCAACC"),
            ("E", "GGGGAA"),
        ]);
        let result = search(&aln, Some(99));
        assert_eq!(result.tree.num_leaves(), 5);
        assert!(result.score > 0);
    }

    #[test]
    fn test_search_deterministic_with_seed() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGAACGA"),
            ("C", "TGCATGCA"),
            ("D", "TGCTTGCT"),
        ]);
        let r1 = search(&aln, Some(42));
        let r2 = search(&aln, Some(42));
        assert_eq!(r1.score, r2.score);
    }

    #[test]
    fn test_search_identical_sequences() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "ACGT"),
            ("D", "ACGT"),
        ]);
        let result = search(&aln, None);
        assert_eq!(result.score, 0);
    }

    #[test]
    fn test_evaluate_user_tree() {
        let aln = make_alignment(&[
            ("A", "AACC"),
            ("B", "AACC"),
            ("C", "CCAA"),
            ("D", "CCAA"),
        ]);
        let good_tree = parse_newick("((A,B),(C,D));").unwrap();
        let bad_tree = parse_newick("((A,C),(B,D));").unwrap();

        let good_result = evaluate(&good_tree, &aln).unwrap();
        let bad_result = evaluate(&bad_tree, &aln).unwrap();

        assert!(good_result.score < bad_result.score);
    }

    #[test]
    fn test_search_matches_evaluate() {
        let aln = make_alignment(&[
            ("A", "AACCGG"),
            ("B", "AACCGG"),
            ("C", "CCAACC"),
            ("D", "CCAACC"),
        ]);
        let result = search(&aln, Some(42));
        let eval_result = evaluate(&result.tree, &aln).unwrap();
        assert_eq!(result.score, eval_result.score);
    }

    #[test]
    fn test_parsimony_known_optimal_four_taxa() {
        // Classic example: when AB share state and CD share state,
        // grouping AB|CD requires fewer changes than AC|BD or AD|BC.
        let aln = make_alignment(&[
            ("A", "AACCAACCAA"),
            ("B", "AACCAACCAA"),
            ("C", "CCAACCAACC"),
            ("D", "CCAACCAACC"),
        ]);

        // ((A,B),(C,D)) should be optimal
        let optimal = parse_newick("((A,B),(C,D));").unwrap();
        let suboptimal1 = parse_newick("((A,C),(B,D));").unwrap();
        let suboptimal2 = parse_newick("((A,D),(B,C));").unwrap();

        let (s_opt, _) = parsimony_score(&optimal, &aln).unwrap();
        let (s_sub1, _) = parsimony_score(&suboptimal1, &aln).unwrap();
        let (s_sub2, _) = parsimony_score(&suboptimal2, &aln).unwrap();

        assert!(s_opt < s_sub1);
        assert!(s_opt < s_sub2);
        assert_eq!(s_sub1, s_sub2); // symmetric data

        // Search should find the optimal
        let result = search(&aln, Some(42));
        assert_eq!(result.score, s_opt);
    }

    #[test]
    fn test_insert_on_edge() {
        let tree = parse_newick("((A,B),C);").unwrap();
        // Find edge from root to the (A,B) internal node
        let internal_id = tree
            .nodes
            .iter()
            .find(|n| !n.is_leaf() && n.id != tree.root)
            .unwrap()
            .id;
        let parent_id = tree.nodes[internal_id].parent.unwrap();

        let new_tree = insert_on_edge(&tree, parent_id, internal_id, "D");
        assert_eq!(new_tree.num_leaves(), 4);

        // Verify D is in the tree
        let leaf_names: Vec<&str> = new_tree.leaf_names().into_iter().collect();
        assert!(leaf_names.contains(&"D"));
    }

    #[test]
    fn test_six_taxa_search() {
        let aln = make_alignment(&[
            ("A", "AACCGGTT"),
            ("B", "AACCGGTT"),
            ("C", "CCAATTGG"),
            ("D", "CCAATTGG"),
            ("E", "GGTTAACC"),
            ("F", "GGTTAACC"),
        ]);
        let result = search(&aln, Some(7));
        assert_eq!(result.tree.num_leaves(), 6);
        // With 3 pairs of identical sequences, optimal groups them
        // Score should be 2 changes per site = 16
        assert_eq!(result.score, 16);
    }

    #[test]
    fn test_per_site_steps() {
        let aln = make_alignment(&[
            ("A", "AC"),
            ("B", "AC"),
            ("C", "CA"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let (_, site_steps) = parsimony_score(&tree, &aln).unwrap();
        // Site 0: A,A,C -> 1 step; Site 1: C,C,A -> 1 step
        assert_eq!(site_steps, vec![1, 1]);
    }
}

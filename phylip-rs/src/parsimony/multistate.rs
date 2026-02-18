//! Multistate parsimony for discrete characters.
//!
//! Generalizes parsimony from the 5-state DNA alphabet to arbitrary
//! discrete state alphabets of up to 32 states.  This corresponds to
//! PHYLIP's `pars` program for discrete multistate data.
//!
//! ## Key types
//!
//! - [`MultiStateSet`] — a 32-bit bitmask representing a set of discrete
//!   character states, analogous to [`StateSet`](super::wagner::StateSet)
//!   for DNA.
//! - [`MultiStateAlignment`] — a taxon-by-character matrix where each cell
//!   holds a state index (0..n_states-1) or `None` for missing data.
//! - [`StepMatrix`] — a cost matrix specifying the penalty for each pair
//!   of state changes.
//!
//! ## Parsimony criteria
//!
//! Two standard criteria are provided:
//!
//! - **Unordered (Fitch)**: all state changes cost 1.  Built via
//!   [`StepMatrix::unordered`].
//! - **Ordered (Wagner)**: states form a linear series 0–1–2–...–(n-1) and
//!   the cost of a change is `|i - j|`.  Built via [`StepMatrix::ordered`].
//!
//! # Example
//!
//! ```
//! use phylip_rs::parsimony::multistate::{
//!     MultiStateAlignment, MultiStateSet, StepMatrix,
//!     multistate_parsimony_score, multistate_search,
//! };
//! use phylip_rs::tree::newick::parse_newick;
//!
//! let data = MultiStateAlignment {
//!     taxa: vec!["A".into(), "B".into(), "C".into()],
//!     characters: vec![
//!         vec![Some(0), Some(1)],
//!         vec![Some(0), Some(1)],
//!         vec![Some(1), Some(0)],
//!     ],
//!     n_states: 2,
//!     n_chars: 2,
//!     state_labels: vec!['0', '1'],
//! };
//!
//! let step_matrix = StepMatrix::unordered(2);
//! let tree = parse_newick("((A,B),C);").unwrap();
//! let score = multistate_parsimony_score(&tree, &data, &step_matrix);
//! assert_eq!(score, 2); // 1 step per character
//! ```
//!
//! # References
//!
//! - Fitch, W.M. (1971). Toward defining the course of evolution: minimum
//!   change for a specific tree topology. *Systematic Zoology*, 20, 406-416.
//! - Kluge, A.G. & Farris, J.S. (1969). Quantitative phyletics and the
//!   evolution of anurans. *Systematic Zoology*, 18, 1-32.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.

use std::fmt;

use crate::tree::types::{Tree, TreeNode};

/// A set of discrete character states, represented as a 32-bit mask.
///
/// Supports up to 32 states (0..31).  Each bit position corresponds to
/// a state index: bit `i` is set if and only if state `i` is in the set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MultiStateSet(u32);

impl MultiStateSet {
    /// The empty set (no states).
    pub const EMPTY: MultiStateSet = MultiStateSet(0);

    /// A set containing exactly one state.
    pub fn singleton(state: u8) -> MultiStateSet {
        debug_assert!(state < 32, "state index must be < 32");
        MultiStateSet(1u32 << state)
    }

    /// Construct a set directly from a bitmask.
    pub fn from_bits(bits: u32) -> MultiStateSet {
        MultiStateSet(bits)
    }

    /// The raw bitmask.
    pub fn bits(self) -> u32 {
        self.0
    }

    /// Set containing all states from 0 to `n_states - 1`.
    pub fn all(n_states: usize) -> MultiStateSet {
        debug_assert!(n_states <= 32);
        if n_states == 32 {
            MultiStateSet(u32::MAX)
        } else {
            MultiStateSet((1u32 << n_states) - 1)
        }
    }

    /// Intersection of two sets.
    pub fn intersection(self, other: MultiStateSet) -> MultiStateSet {
        MultiStateSet(self.0 & other.0)
    }

    /// Union of two sets.
    pub fn union(self, other: MultiStateSet) -> MultiStateSet {
        MultiStateSet(self.0 | other.0)
    }

    /// True if this set contains no states.
    pub fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// Number of states in the set.
    pub fn count(self) -> u32 {
        self.0.count_ones()
    }

    /// True if the set contains the given state.
    pub fn contains(self, state: u8) -> bool {
        debug_assert!(state < 32);
        self.0 & (1u32 << state) != 0
    }

    /// Return the lowest state index in the set, or `None` if empty.
    pub fn lowest(self) -> Option<u8> {
        if self.is_empty() {
            None
        } else {
            Some(self.0.trailing_zeros() as u8)
        }
    }
}

impl fmt::Display for MultiStateSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{{")?;
        let mut first = true;
        for i in 0..32u8 {
            if self.contains(i) {
                if !first {
                    write!(f, ",")?;
                }
                write!(f, "{}", i)?;
                first = false;
            }
        }
        write!(f, "}}")
    }
}

/// A discrete character matrix with arbitrary state labels.
#[derive(Debug, Clone)]
pub struct MultiStateAlignment {
    /// Taxon names.
    pub taxa: Vec<String>,
    /// Character matrix: `characters[taxon][site]` = state index
    /// (0..n_states-1), or `None` for missing data.
    pub characters: Vec<Vec<Option<u8>>>,
    /// Number of possible states.
    pub n_states: usize,
    /// Number of characters (sites).
    pub n_chars: usize,
    /// Human-readable label for each state (e.g., `['0','1','2','3']`).
    pub state_labels: Vec<char>,
}

impl MultiStateAlignment {
    /// Number of taxa.
    pub fn ntaxa(&self) -> usize {
        self.taxa.len()
    }

    /// Taxon name by index.
    pub fn taxon_name(&self, i: usize) -> &str {
        &self.taxa[i]
    }
}

/// Cost matrix for state changes.
///
/// `costs[i][j]` = cost of changing from state `i` to state `j`.
#[derive(Debug, Clone)]
pub struct StepMatrix {
    /// Number of states.
    pub n_states: usize,
    /// Cost matrix: `costs[i][j]` is the cost of i -> j.
    pub costs: Vec<Vec<usize>>,
}

impl StepMatrix {
    /// Unordered (Fitch-like): all changes cost 1, self-transitions cost 0.
    pub fn unordered(n_states: usize) -> StepMatrix {
        let mut costs = vec![vec![1usize; n_states]; n_states];
        for i in 0..n_states {
            costs[i][i] = 0;
        }
        StepMatrix { n_states, costs }
    }

    /// Ordered (Wagner-like): cost = |i - j|.
    ///
    /// States are assumed to form a linear series 0-1-2-...(n-1).
    pub fn ordered(n_states: usize) -> StepMatrix {
        let mut costs = vec![vec![0usize; n_states]; n_states];
        for i in 0..n_states {
            for j in 0..n_states {
                costs[i][j] = if i >= j { i - j } else { j - i };
            }
        }
        StepMatrix { n_states, costs }
    }
}

/// Per-node data for the multistate Fitch-like postorder pass.
#[derive(Debug, Clone)]
struct MsNodeData {
    /// State set at this node for each character.
    states: Vec<MultiStateSet>,
    /// Accumulated step count below this node for each character.
    steps: Vec<usize>,
}

impl MsNodeData {
    fn new(n_chars: usize) -> Self {
        Self {
            states: vec![MultiStateSet::EMPTY; n_chars],
            steps: vec![0; n_chars],
        }
    }
}

/// Compute the multistate parsimony score for a tree and alignment.
///
/// Uses a Fitch-like postorder algorithm generalized for multiple states.
/// For unordered characters (all changes cost 1), this is exactly the
/// Fitch algorithm.  For ordered or weighted step matrices, this uses a
/// heuristic that may slightly overestimate the cost for complex step
/// matrices, but is exact for unordered and ordered (linear) cases.
///
/// # Arguments
///
/// * `tree`        — the phylogenetic tree to score
/// * `data`        — the multistate character matrix
/// * `step_matrix` — the step-cost matrix
///
/// # Returns
///
/// The total parsimony score (sum over all characters).
pub fn multistate_parsimony_score(
    tree: &Tree,
    data: &MultiStateAlignment,
    step_matrix: &StepMatrix,
) -> usize {
    let n_chars = data.n_chars;
    let n_states = data.n_states;

    // Build name -> taxon index map.
    let name_to_idx: std::collections::HashMap<&str, usize> = data
        .taxa
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let nnodes = tree.num_nodes();
    let mut node_data = vec![MsNodeData::new(n_chars); nnodes];

    // Initialize leaf states.
    for node in &tree.nodes {
        if node.is_leaf() {
            if let Some(name) = &node.name {
                if let Some(&taxon_idx) = name_to_idx.get(name.as_str()) {
                    for ch in 0..n_chars {
                        node_data[node.id].states[ch] = match data.characters[taxon_idx][ch] {
                            Some(s) => MultiStateSet::singleton(s),
                            None => MultiStateSet::all(n_states), // missing = all states
                        };
                    }
                }
            }
        }
    }

    // Check if the step matrix is "unordered" (all off-diagonal costs = 1).
    let is_unordered = is_unordered_step_matrix(step_matrix);

    // Postorder traversal.
    let postorder = tree.postorder();
    for &node_id in &postorder {
        let children = &tree.nodes[node_id].children;
        if children.is_empty() {
            continue;
        }

        let first_child = children[0];
        for ch in 0..n_chars {
            let mut current_states = node_data[first_child].states[ch];
            let mut current_steps: usize =
                children.iter().map(|&c| node_data[c].steps[ch]).sum();

            for &child in &children[1..] {
                let child_states = node_data[child].states[ch];

                if is_unordered {
                    // Standard Fitch: intersection/union.
                    let isect = current_states.intersection(child_states);
                    if isect.is_empty() {
                        current_states = current_states.union(child_states);
                        current_steps += 1;
                    } else {
                        current_states = isect;
                    }
                } else {
                    // Weighted / ordered: use minimum-cost assignment.
                    let (combined, cost) =
                        weighted_combine(current_states, child_states, step_matrix);
                    current_states = combined;
                    current_steps += cost;
                }
            }

            node_data[node_id].states[ch] = current_states;
            node_data[node_id].steps[ch] = current_steps;
        }
    }

    // Sum over all characters at the root.
    let root = tree.root;
    node_data[root].steps.iter().sum()
}

/// Combine two state sets under a weighted step matrix.
///
/// For each possible parent state, find the minimum cost of connecting
/// to both child state sets.  The optimal parent set is the set of states
/// achieving this minimum cost.
fn weighted_combine(
    left: MultiStateSet,
    right: MultiStateSet,
    step_matrix: &StepMatrix,
) -> (MultiStateSet, usize) {
    let n = step_matrix.n_states;
    let mut min_cost = usize::MAX;
    let mut best_states = 0u32;

    for parent_state in 0..n {
        // Minimum cost from parent_state to any state in left set.
        let mut left_min = usize::MAX;
        for s in 0..n {
            if left.contains(s as u8) {
                let c = step_matrix.costs[parent_state][s];
                if c < left_min {
                    left_min = c;
                }
            }
        }

        // Minimum cost from parent_state to any state in right set.
        let mut right_min = usize::MAX;
        for s in 0..n {
            if right.contains(s as u8) {
                let c = step_matrix.costs[parent_state][s];
                if c < right_min {
                    right_min = c;
                }
            }
        }

        let total = left_min.saturating_add(right_min);
        if total < min_cost {
            min_cost = total;
            best_states = 1u32 << parent_state;
        } else if total == min_cost {
            best_states |= 1u32 << parent_state;
        }
    }

    // The "extra cost" at this node is the minimum cost minus the cost
    // of having parent in one of the child sets (which is zero when there
    // is an intersection).
    (MultiStateSet::from_bits(best_states), min_cost)
}

/// Check whether a step matrix is unordered (Fitch-like).
fn is_unordered_step_matrix(sm: &StepMatrix) -> bool {
    for i in 0..sm.n_states {
        for j in 0..sm.n_states {
            if i == j {
                if sm.costs[i][j] != 0 {
                    return false;
                }
            } else if sm.costs[i][j] != 1 {
                return false;
            }
        }
    }
    true
}

/// Result of a multistate parsimony search.
#[derive(Debug, Clone)]
pub struct MultiStateParsimonyResult {
    /// The best tree found.
    pub tree: Tree,
    /// The parsimony score (total steps across all characters).
    pub score: usize,
}

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

/// Search for the most parsimonious tree for multistate data.
///
/// Uses stepwise addition followed by SPR rearrangement, scoring
/// trees with the provided step matrix.
///
/// # Arguments
///
/// * `data`        — the multistate character matrix
/// * `step_matrix` — the step-cost matrix
/// * `seed`        — optional random seed for taxon addition order
pub fn multistate_search(
    data: &MultiStateAlignment,
    step_matrix: &StepMatrix,
    seed: Option<u64>,
) -> MultiStateParsimonyResult {
    let ntaxa = data.ntaxa();
    assert!(ntaxa >= 3, "need at least 3 taxa for tree search");

    // Taxon addition order.
    let mut order: Vec<usize> = (0..ntaxa).collect();
    if let Some(s) = seed {
        let mut rng = SimpleRng::new(s);
        rng.shuffle(&mut order);
    }

    // Phase 1: Stepwise addition — start with a 3-taxon star.
    let mut tree = Tree::new();
    tree.rooted = false;
    let root = tree.add_node(None, None, None);
    tree.root = root;
    for &taxon_idx in &order[..3] {
        let name = data.taxon_name(taxon_idx).to_string();
        tree.add_node(Some(name), None, Some(root));
    }

    // Add remaining taxa one by one.
    for i in 3..ntaxa {
        let taxon_name = data.taxon_name(order[i]);
        tree = ms_best_insertion(&tree, taxon_name, data, step_matrix);
    }

    let mut score = multistate_parsimony_score(&tree, data, step_matrix);

    // Phase 2: SPR rearrangement.
    loop {
        let (new_tree, new_score, improved) =
            ms_spr_round(&tree, data, step_matrix, score);
        if !improved {
            break;
        }
        tree = new_tree;
        score = new_score;
    }

    MultiStateParsimonyResult { tree, score }
}

/// Try inserting a new taxon at every possible edge; return the tree
/// with the best parsimony score.
fn ms_best_insertion(
    tree: &Tree,
    taxon_name: &str,
    data: &MultiStateAlignment,
    step_matrix: &StepMatrix,
) -> Tree {
    let mut best_tree = None;
    let mut best_score = usize::MAX;

    let edges: Vec<(usize, usize)> = tree
        .nodes
        .iter()
        .filter_map(|n| n.parent.map(|p| (p, n.id)))
        .collect();

    for (parent_id, child_id) in edges {
        let candidate = ms_insert_on_edge(tree, parent_id, child_id, taxon_name);
        let score = multistate_parsimony_score(&candidate, data, step_matrix);
        if score < best_score {
            best_score = score;
            best_tree = Some(candidate);
        }
    }

    best_tree.unwrap_or_else(|| {
        let mut t = tree.clone();
        let root = t.root;
        t.add_node(Some(taxon_name.to_string()), None, Some(root));
        t
    })
}

/// Insert a new taxon on the edge between parent_id and child_id.
fn ms_insert_on_edge(
    tree: &Tree,
    parent_id: usize,
    child_id: usize,
    taxon_name: &str,
) -> Tree {
    let mut new_tree = Tree::new();
    new_tree.rooted = tree.rooted;

    let mut id_map = vec![0usize; tree.nodes.len()];

    // Copy existing nodes.
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

    // New internal node.
    let new_internal_id = new_tree.nodes.len();
    new_tree.nodes.push(TreeNode {
        id: new_internal_id,
        name: None,
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    // New leaf.
    let new_leaf_id = new_tree.nodes.len();
    new_tree.nodes.push(TreeNode {
        id: new_leaf_id,
        name: Some(taxon_name.to_string()),
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    // Rebuild parent-child links.
    for node in &tree.nodes {
        let mapped_id = id_map[node.id];
        for &old_child in &node.children {
            if node.id == parent_id && old_child == child_id {
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

/// One round of SPR rearrangement for multistate data.
fn ms_spr_round(
    tree: &Tree,
    data: &MultiStateAlignment,
    step_matrix: &StepMatrix,
    current_score: usize,
) -> (Tree, usize, bool) {
    let mut best_tree = tree.clone();
    let mut best_score = current_score;
    let mut improved = false;

    let nodes_to_try: Vec<usize> = tree
        .nodes
        .iter()
        .filter(|n| n.parent.is_some() && n.parent != Some(tree.root))
        .map(|n| n.id)
        .collect();

    for &prune_node in &nodes_to_try {
        if let Some(pruned) = ms_prune_subtree(tree, prune_node) {
            let edges: Vec<(usize, usize)> = pruned
                .main_tree
                .nodes
                .iter()
                .filter_map(|n| n.parent.map(|p| (p, n.id)))
                .collect();

            for (graft_parent, graft_child) in edges {
                let candidate = ms_regraft_subtree(&pruned, graft_parent, graft_child);
                let score = multistate_parsimony_score(&candidate, data, step_matrix);
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

/// Result of pruning a subtree.
struct MsPrunedTree {
    main_tree: Tree,
    subtree: Tree,
}

/// Collect all node IDs in the subtree rooted at `node_id`.
fn ms_collect_subtree_ids(tree: &Tree, node_id: usize) -> Vec<usize> {
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

/// Prune a subtree rooted at `prune_id`.
fn ms_prune_subtree(tree: &Tree, prune_id: usize) -> Option<MsPrunedTree> {
    let parent_id = tree.nodes[prune_id].parent?;

    let subtree_ids: std::collections::HashSet<usize> =
        ms_collect_subtree_ids(tree, prune_id).into_iter().collect();

    // Build the subtree.
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

    for &sid in &subtree_ids {
        for &child in &tree.nodes[sid].children {
            if subtree_ids.contains(&child) {
                let mp = subtree_id_map[sid];
                let mc = subtree_id_map[child];
                subtree.nodes[mp].children.push(mc);
                subtree.nodes[mc].parent = Some(mp);
            }
        }
    }

    // Build the main tree (without subtree nodes).
    let mut main_tree = Tree::new();
    main_tree.rooted = tree.rooted;
    let mut main_id_map = vec![usize::MAX; tree.nodes.len()];

    let parent_remaining_children: Vec<usize> = tree.nodes[parent_id]
        .children
        .iter()
        .filter(|&&c| !subtree_ids.contains(&c))
        .copied()
        .collect();

    let collapse_parent =
        parent_remaining_children.len() == 1 && parent_id != tree.root;

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
                let surviving_child = parent_remaining_children[0];
                let mapped_child = main_id_map[surviving_child];
                main_tree.nodes[mapped_id].children.push(mapped_child);
                main_tree.nodes[mapped_child].parent = Some(mapped_id);
            } else {
                let mapped_child = main_id_map[child];
                main_tree.nodes[mapped_id].children.push(mapped_child);
                main_tree.nodes[mapped_child].parent = Some(mapped_id);
            }
        }
    }

    Some(MsPrunedTree {
        main_tree,
        subtree,
    })
}

/// Regraft a pruned subtree onto an edge in the main tree.
fn ms_regraft_subtree(
    pruned: &MsPrunedTree,
    graft_parent: usize,
    graft_child: usize,
) -> Tree {
    let main = &pruned.main_tree;
    let sub = &pruned.subtree;

    let mut new_tree = Tree::new();
    new_tree.rooted = main.rooted;

    // Copy main tree nodes.
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

    // Copy subtree nodes.
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

    // Graft node.
    let graft_node_id = new_tree.nodes.len();
    new_tree.nodes.push(TreeNode {
        id: graft_node_id,
        name: None,
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    // Rebuild main links, splitting the graft edge.
    for node in &main.nodes {
        let mapped_id = main_map[node.id];
        for &child in &node.children {
            if node.id == graft_parent && child == graft_child {
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

    // Rebuild subtree links.
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

/// Unordered multistate parsimony scorer that implements `ParsimonyScorer`.
///
/// Uses the standard Fitch algorithm: intersection if non-empty, else
/// union with +1 step.  This scorer works with the DNA-based
/// `parsimony_score_with` function from `wagner.rs` and maps each of the
/// n_states to the first n bits of a `StateSet`.
///
/// NOTE: This scorer is provided for API completeness but the primary
/// multistate scoring path is [`multistate_parsimony_score`], which
/// supports weighted step matrices.
pub struct MultiStateFitchScorer {
    /// Number of discrete states.
    pub n_states: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::newick::parse_newick;

    // --- MultiStateSet tests ---

    #[test]
    fn test_singleton() {
        let s = MultiStateSet::singleton(0);
        assert!(s.contains(0));
        assert!(!s.contains(1));
        assert_eq!(s.count(), 1);
    }

    #[test]
    fn test_singleton_higher() {
        let s = MultiStateSet::singleton(5);
        assert!(s.contains(5));
        assert!(!s.contains(0));
        assert!(!s.contains(4));
        assert_eq!(s.count(), 1);
    }

    #[test]
    fn test_intersection() {
        let a = MultiStateSet::singleton(0).union(MultiStateSet::singleton(1));
        let b = MultiStateSet::singleton(1).union(MultiStateSet::singleton(2));
        let isect = a.intersection(b);
        assert!(isect.contains(1));
        assert!(!isect.contains(0));
        assert!(!isect.contains(2));
        assert_eq!(isect.count(), 1);
    }

    #[test]
    fn test_union() {
        let a = MultiStateSet::singleton(0);
        let b = MultiStateSet::singleton(2);
        let u = a.union(b);
        assert!(u.contains(0));
        assert!(u.contains(2));
        assert!(!u.contains(1));
        assert_eq!(u.count(), 2);
    }

    #[test]
    fn test_empty() {
        assert!(MultiStateSet::EMPTY.is_empty());
        assert_eq!(MultiStateSet::EMPTY.count(), 0);
    }

    #[test]
    fn test_from_bits() {
        let s = MultiStateSet::from_bits(0b1010);
        assert!(s.contains(1));
        assert!(s.contains(3));
        assert!(!s.contains(0));
        assert!(!s.contains(2));
    }

    #[test]
    fn test_all() {
        let s = MultiStateSet::all(4);
        for i in 0..4 {
            assert!(s.contains(i));
        }
        assert!(!s.contains(4));
        assert_eq!(s.count(), 4);
    }

    #[test]
    fn test_lowest() {
        assert_eq!(MultiStateSet::EMPTY.lowest(), None);
        assert_eq!(MultiStateSet::singleton(3).lowest(), Some(3));
        let s = MultiStateSet::singleton(2).union(MultiStateSet::singleton(5));
        assert_eq!(s.lowest(), Some(2));
    }

    #[test]
    fn test_display() {
        let s = MultiStateSet::singleton(0).union(MultiStateSet::singleton(3));
        assert_eq!(format!("{}", s), "{0,3}");
        assert_eq!(format!("{}", MultiStateSet::EMPTY), "{}");
    }

    // --- StepMatrix tests ---

    #[test]
    fn test_unordered_step_matrix() {
        let sm = StepMatrix::unordered(4);
        assert_eq!(sm.n_states, 4);
        for i in 0..4 {
            for j in 0..4 {
                if i == j {
                    assert_eq!(sm.costs[i][j], 0);
                } else {
                    assert_eq!(sm.costs[i][j], 1);
                }
            }
        }
    }

    #[test]
    fn test_ordered_step_matrix() {
        let sm = StepMatrix::ordered(4);
        assert_eq!(sm.n_states, 4);
        assert_eq!(sm.costs[0][0], 0);
        assert_eq!(sm.costs[0][1], 1);
        assert_eq!(sm.costs[0][2], 2);
        assert_eq!(sm.costs[0][3], 3);
        assert_eq!(sm.costs[1][3], 2);
        assert_eq!(sm.costs[3][0], 3);
        assert_eq!(sm.costs[2][1], 1);
    }

    // --- Scoring tests ---

    #[test]
    fn test_identical_data_zero_score() {
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into()],
            characters: vec![
                vec![Some(0), Some(1)],
                vec![Some(0), Some(1)],
                vec![Some(0), Some(1)],
            ],
            n_states: 2,
            n_chars: 2,
            state_labels: vec!['0', '1'],
        };
        let sm = StepMatrix::unordered(2);
        let tree = parse_newick("((A,B),C);").unwrap();
        let score = multistate_parsimony_score(&tree, &data, &sm);
        assert_eq!(score, 0);
    }

    #[test]
    fn test_two_state_one_change() {
        // A=0, B=0, C=1: one change needed.
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into()],
            characters: vec![
                vec![Some(0)],
                vec![Some(0)],
                vec![Some(1)],
            ],
            n_states: 2,
            n_chars: 1,
            state_labels: vec!['0', '1'],
        };
        let sm = StepMatrix::unordered(2);
        let tree = parse_newick("((A,B),C);").unwrap();
        let score = multistate_parsimony_score(&tree, &data, &sm);
        assert_eq!(score, 1);
    }

    #[test]
    fn test_two_state_matches_dna_parsimony() {
        // Encode binary data: 0 -> state 0, 1 -> state 1.
        // This should match Fitch parsimony on the equivalent DNA data.
        // Data: A=0,0; B=0,0; C=1,1; D=1,1
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into(), "D".into()],
            characters: vec![
                vec![Some(0), Some(0)],
                vec![Some(0), Some(0)],
                vec![Some(1), Some(1)],
                vec![Some(1), Some(1)],
            ],
            n_states: 2,
            n_chars: 2,
            state_labels: vec!['0', '1'],
        };
        let sm = StepMatrix::unordered(2);
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let score = multistate_parsimony_score(&tree, &data, &sm);
        // Same as Fitch: 1 step per character = 2.
        assert_eq!(score, 2);
    }

    #[test]
    fn test_four_state_matches_dna() {
        // 4-state characters matching A,C,G,T as states 0,1,2,3.
        // A=0, B=0, C=1, D=1 at one character -> 1 step.
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into(), "D".into()],
            characters: vec![
                vec![Some(0)],
                vec![Some(0)],
                vec![Some(1)],
                vec![Some(1)],
            ],
            n_states: 4,
            n_chars: 1,
            state_labels: vec!['0', '1', '2', '3'],
        };
        let sm = StepMatrix::unordered(4);
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let score = multistate_parsimony_score(&tree, &data, &sm);
        assert_eq!(score, 1);
    }

    #[test]
    fn test_eight_state_characters() {
        // 8 states. A has 0, B has 0, C has 7 -> 1 step.
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into()],
            characters: vec![
                vec![Some(0)],
                vec![Some(0)],
                vec![Some(7)],
            ],
            n_states: 8,
            n_chars: 1,
            state_labels: vec!['0', '1', '2', '3', '4', '5', '6', '7'],
        };
        let sm = StepMatrix::unordered(8);
        let tree = parse_newick("((A,B),C);").unwrap();
        let score = multistate_parsimony_score(&tree, &data, &sm);
        assert_eq!(score, 1);
    }

    #[test]
    fn test_eight_state_all_different() {
        // 3 taxa, all different states -> 2 steps (unordered).
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into()],
            characters: vec![
                vec![Some(0)],
                vec![Some(3)],
                vec![Some(7)],
            ],
            n_states: 8,
            n_chars: 1,
            state_labels: vec!['0', '1', '2', '3', '4', '5', '6', '7'],
        };
        let sm = StepMatrix::unordered(8);
        let tree = parse_newick("((A,B),C);").unwrap();
        let score = multistate_parsimony_score(&tree, &data, &sm);
        assert_eq!(score, 2);
    }

    // --- Ordered vs unordered ---

    #[test]
    fn test_ordered_vs_unordered_different_scores() {
        // A=0, B=0, C=3 on an unordered matrix: 1 step.
        // On an ordered matrix (4 states): cost = |0-3| = 3 steps.
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into()],
            characters: vec![
                vec![Some(0)],
                vec![Some(0)],
                vec![Some(3)],
            ],
            n_states: 4,
            n_chars: 1,
            state_labels: vec!['0', '1', '2', '3'],
        };
        let sm_unordered = StepMatrix::unordered(4);
        let sm_ordered = StepMatrix::ordered(4);
        let tree = parse_newick("((A,B),C);").unwrap();

        let score_unord = multistate_parsimony_score(&tree, &data, &sm_unordered);
        let score_ord = multistate_parsimony_score(&tree, &data, &sm_ordered);

        assert_eq!(score_unord, 1);
        assert_eq!(score_ord, 3);
        assert!(score_ord > score_unord);
    }

    #[test]
    fn test_ordered_adjacent_states() {
        // A=0, B=1, C=2: ordered cost = min path.
        // On ((A,B),C): parent of (A,B) can be 0 or 1.
        // If parent = 0: cost to A = 0, cost to B = 1 -> 1
        //   root with state set for parent={0} and C={2}:
        //   best root = 1: cost 1 to {0} + cost 1 to {2} = 2.
        // If parent = 1: cost to A = 1, cost to B = 0 -> 1
        //   best root = 1: cost 0 + cost 1 to {2} = 1. Total = 1+1 = 2.
        // Actually, let's compute correctly:
        // (A,B): A=0, B=1. For each potential parent state p:
        //   p=0: min(|0-0|,|0-1|) = 0 (left) + ... wait, this is Fitch generalized.
        //
        // The weighted_combine for (A={0}, B={1}):
        //   p=0: cost to A(0)=0, cost to B(1)=1 -> total 1
        //   p=1: cost to A(0)=1, cost to B(1)=0 -> total 1
        //   p=2: cost to A(0)=2, cost to B(1)=1 -> total 3
        //   p=3: cost to A(0)=3, cost to B(1)=2 -> total 5
        //   min = 1, states = {0,1}
        //
        // Root combine({0,1}, C={2}):
        //   p=0: min(0,1)=0 + 2 = 2
        //   p=1: min(1,0)=0 + 1 = 1
        //   p=2: min(2,1)=1 + 0 = 1
        //   p=3: min(3,2)=2 + 1 = 3
        //   min = 1, states = {1,2}
        //
        // Total = 1 (from AB) + 1 (from root) = 2.
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into()],
            characters: vec![
                vec![Some(0)],
                vec![Some(1)],
                vec![Some(2)],
            ],
            n_states: 4,
            n_chars: 1,
            state_labels: vec!['0', '1', '2', '3'],
        };
        let sm = StepMatrix::ordered(4);
        let tree = parse_newick("((A,B),C);").unwrap();
        let score = multistate_parsimony_score(&tree, &data, &sm);
        assert_eq!(score, 2);
    }

    #[test]
    fn test_ordered_same_as_unordered_for_adjacent() {
        // A=0, B=0, C=1: cost 1 in both ordered and unordered.
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into()],
            characters: vec![
                vec![Some(0)],
                vec![Some(0)],
                vec![Some(1)],
            ],
            n_states: 4,
            n_chars: 1,
            state_labels: vec!['0', '1', '2', '3'],
        };
        let sm_u = StepMatrix::unordered(4);
        let sm_o = StepMatrix::ordered(4);
        let tree = parse_newick("((A,B),C);").unwrap();

        let su = multistate_parsimony_score(&tree, &data, &sm_u);
        let so = multistate_parsimony_score(&tree, &data, &sm_o);
        assert_eq!(su, 1);
        assert_eq!(so, 1);
    }

    // --- Missing data ---

    #[test]
    fn test_missing_data() {
        // Missing data (None) should be compatible with any state.
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into()],
            characters: vec![
                vec![Some(0)],
                vec![Some(0)],
                vec![None], // missing
            ],
            n_states: 2,
            n_chars: 1,
            state_labels: vec!['0', '1'],
        };
        let sm = StepMatrix::unordered(2);
        let tree = parse_newick("((A,B),C);").unwrap();
        let score = multistate_parsimony_score(&tree, &data, &sm);
        assert_eq!(score, 0); // missing data doesn't force a change
    }

    // --- Search tests ---

    #[test]
    fn test_search_basic() {
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into(), "D".into()],
            characters: vec![
                vec![Some(0), Some(0)],
                vec![Some(0), Some(0)],
                vec![Some(1), Some(1)],
                vec![Some(1), Some(1)],
            ],
            n_states: 2,
            n_chars: 2,
            state_labels: vec!['0', '1'],
        };
        let sm = StepMatrix::unordered(2);
        let result = multistate_search(&data, &sm, Some(42));
        assert_eq!(result.tree.num_leaves(), 4);
        assert_eq!(result.score, 2); // 1 step per character
    }

    #[test]
    fn test_search_five_taxa() {
        let data = MultiStateAlignment {
            taxa: vec![
                "A".into(),
                "B".into(),
                "C".into(),
                "D".into(),
                "E".into(),
            ],
            characters: vec![
                vec![Some(0), Some(0), Some(0)],
                vec![Some(0), Some(0), Some(0)],
                vec![Some(1), Some(1), Some(1)],
                vec![Some(1), Some(1), Some(1)],
                vec![Some(2), Some(2), Some(2)],
            ],
            n_states: 3,
            n_chars: 3,
            state_labels: vec!['0', '1', '2'],
        };
        let sm = StepMatrix::unordered(3);
        let result = multistate_search(&data, &sm, Some(7));
        assert_eq!(result.tree.num_leaves(), 5);
        assert!(result.score > 0);
    }

    #[test]
    fn test_search_ordered() {
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into()],
            characters: vec![
                vec![Some(0)],
                vec![Some(0)],
                vec![Some(2)],
            ],
            n_states: 3,
            n_chars: 1,
            state_labels: vec!['0', '1', '2'],
        };
        let sm = StepMatrix::ordered(3);
        let result = multistate_search(&data, &sm, None);
        assert_eq!(result.tree.num_leaves(), 3);
        // Ordered: |0 - 2| = 2.
        assert_eq!(result.score, 2);
    }

    // --- Topology comparison ---

    #[test]
    fn test_topology_matters_multistate() {
        // A=0, B=0, C=1, D=1.
        // ((A,B),(C,D)) = 1 step; ((A,C),(B,D)) = 2 steps.
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into(), "D".into()],
            characters: vec![
                vec![Some(0)],
                vec![Some(0)],
                vec![Some(1)],
                vec![Some(1)],
            ],
            n_states: 2,
            n_chars: 1,
            state_labels: vec!['0', '1'],
        };
        let sm = StepMatrix::unordered(2);

        let t1 = parse_newick("((A,B),(C,D));").unwrap();
        let t2 = parse_newick("((A,C),(B,D));").unwrap();

        let s1 = multistate_parsimony_score(&t1, &data, &sm);
        let s2 = multistate_parsimony_score(&t2, &data, &sm);

        assert_eq!(s1, 1);
        assert_eq!(s2, 2);
    }

    // --- Weighted combine correctness ---

    #[test]
    fn test_weighted_combine_unordered() {
        let sm = StepMatrix::unordered(3);
        let left = MultiStateSet::singleton(0);
        let right = MultiStateSet::singleton(0);
        let (combined, cost) = weighted_combine(left, right, &sm);
        assert!(combined.contains(0));
        assert_eq!(cost, 0);

        let right2 = MultiStateSet::singleton(1);
        let (combined2, cost2) = weighted_combine(left, right2, &sm);
        // All states have cost 1 (except the one matching).
        // p=0: 0 + 1 = 1; p=1: 1 + 0 = 1; p=2: 1 + 1 = 2
        assert_eq!(cost2, 1);
        assert!(combined2.contains(0));
        assert!(combined2.contains(1));
    }

    #[test]
    fn test_weighted_combine_ordered() {
        let sm = StepMatrix::ordered(4);
        let left = MultiStateSet::singleton(0);
        let right = MultiStateSet::singleton(3);
        let (combined, cost) = weighted_combine(left, right, &sm);
        // p=0: 0 + 3 = 3; p=1: 1 + 2 = 3; p=2: 2 + 1 = 3; p=3: 3 + 0 = 3.
        // All tie at cost 3.
        assert_eq!(cost, 3);
        assert_eq!(combined.count(), 4); // all states are optimal
    }

    // --- Multiple characters ---

    #[test]
    fn test_multiple_characters() {
        let data = MultiStateAlignment {
            taxa: vec!["A".into(), "B".into(), "C".into()],
            characters: vec![
                vec![Some(0), Some(1), Some(2)],
                vec![Some(0), Some(1), Some(2)],
                vec![Some(1), Some(0), Some(0)],
            ],
            n_states: 3,
            n_chars: 3,
            state_labels: vec!['0', '1', '2'],
        };
        let sm = StepMatrix::unordered(3);
        let tree = parse_newick("((A,B),C);").unwrap();
        let score = multistate_parsimony_score(&tree, &data, &sm);
        // Each character: A,B share, C differs -> 1 step each = 3.
        assert_eq!(score, 3);
    }

    // --- is_unordered_step_matrix ---

    #[test]
    fn test_is_unordered() {
        assert!(is_unordered_step_matrix(&StepMatrix::unordered(4)));
        assert!(!is_unordered_step_matrix(&StepMatrix::ordered(4)));
        assert!(is_unordered_step_matrix(&StepMatrix::unordered(2)));
        // Ordered with 2 states is also unordered (|0-1| = 1).
        assert!(is_unordered_step_matrix(&StepMatrix::ordered(2)));
    }
}

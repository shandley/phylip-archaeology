//! DNA compatibility analysis (PHYLIP's dnacomp).
//!
//! For each tree topology, counts how many sites are "compatible" with the
//! tree — i.e., their Fitch parsimony score equals the minimum possible
//! number of changes: (number of distinct states at that site - 1).
//!
//! A site is compatible with a tree if it can be explained without any
//! homoplasy: each state transition occurs exactly once.
//!
//! The search algorithm uses stepwise addition with SPR rearrangement,
//! optimizing for the maximum number of compatible sites rather than
//! minimum parsimony score.
//!
//! # References
//!
//! - Felsenstein, J. (1981). A likelihood approach to character weighting
//!   and what it tells us about parsimony and compatibility. *Biological
//!   Journal of the Linnean Society*, 16, 183-196.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//!   Chapter 10: Compatibility.

use crate::parsimony::wagner::{parsimony_score, StateSet};
use crate::tree::types::{Alignment, Tree};

/// Result of a DNA compatibility analysis.
#[derive(Debug, Clone)]
pub struct DnaCompatResult {
    /// The tree that maximizes the number of compatible sites.
    pub tree: Tree,
    /// Number of sites compatible with the tree.
    pub compatible_sites: usize,
    /// Total number of sites in the alignment.
    pub total_sites: usize,
    /// Fraction of sites that are compatible.
    pub compatibility_fraction: f64,
}

/// Count the number of distinct nucleotide states observed at a given site.
///
/// Only counts unambiguous bases (A, C, G, T). Gaps are treated as a
/// separate state. Ambiguity codes are ignored.
fn count_distinct_states(alignment: &Alignment, site: usize) -> usize {
    let mut seen = 0u8; // bits: 0=A, 1=C, 2=G, 3=T, 4=gap

    for taxon in 0..alignment.ntaxa() {
        let base = alignment.base(taxon, site);
        let ss = StateSet::from_base(base);
        // Only count single-state (unambiguous) observations
        if ss.count() == 1 {
            seen |= ss_to_bits(ss);
        }
    }

    seen.count_ones() as usize
}

/// Convert a single-state StateSet to a bit mask.
fn ss_to_bits(ss: StateSet) -> u8 {
    let mut bits = 0u8;
    if !ss.intersection(StateSet::A).is_empty() {
        bits |= 1;
    }
    if !ss.intersection(StateSet::C).is_empty() {
        bits |= 2;
    }
    if !ss.intersection(StateSet::G).is_empty() {
        bits |= 4;
    }
    if !ss.intersection(StateSet::T).is_empty() {
        bits |= 8;
    }
    if !ss.intersection(StateSet::GAP).is_empty() {
        bits |= 16;
    }
    bits
}

/// Count sites compatible with a given tree.
///
/// A site is compatible with the tree if its Fitch parsimony score equals
/// (number of distinct states at that site - 1). This means the character
/// can be explained on this tree without any homoplasy.
///
/// Invariant sites (only one state observed) are always compatible.
pub fn count_compatible_sites(tree: &Tree, alignment: &Alignment) -> usize {
    let (_, site_steps) = match parsimony_score(tree, alignment) {
        Ok(result) => result,
        Err(_) => return 0,
    };

    let nsites = alignment.nsites();
    let mut compatible = 0;

    for site in 0..nsites {
        let distinct = count_distinct_states(alignment, site);
        let min_changes = if distinct > 0 { distinct - 1 } else { 0 };
        if site_steps[site] == min_changes {
            compatible += 1;
        }
    }

    compatible
}

/// Search for the tree that maximizes the number of compatible sites.
///
/// Uses stepwise addition followed by SPR rearrangement, analogous to
/// the Wagner parsimony search but optimizing for compatible site count
/// instead of minimum total parsimony score.
///
/// # Arguments
/// * `alignment` — the DNA alignment to analyze
/// * `seed` — optional random seed for taxon addition order
///
/// # Returns
/// A `DnaCompatResult` with the best tree found and its compatibility score.
pub fn dna_compat_search(
    alignment: &Alignment,
    seed: Option<u64>,
) -> DnaCompatResult {
    let ntaxa = alignment.ntaxa();
    assert!(ntaxa >= 3, "need at least 3 taxa for tree search");

    let nsites = alignment.nsites();

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
        tree = best_insertion_compat(&tree, taxon_name, alignment);
    }

    let mut compat_count = count_compatible_sites(&tree, alignment);

    // Phase 2: SPR rearrangement (iterate until no improvement)
    loop {
        let (new_tree, new_count, improved) =
            spr_round_compat(&tree, alignment, compat_count);
        if !improved {
            break;
        }
        tree = new_tree;
        compat_count = new_count;
    }

    // Final count
    let final_compat = count_compatible_sites(&tree, alignment);
    let fraction = if nsites > 0 {
        final_compat as f64 / nsites as f64
    } else {
        0.0
    };

    DnaCompatResult {
        tree,
        compatible_sites: final_compat,
        total_sites: nsites,
        compatibility_fraction: fraction,
    }
}

// ============================================================
// Tree search internals (adapted from wagner.rs for compat)
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

/// Build a 3-taxon star tree from the first 3 taxa.
fn build_initial_tree(taxa: &[usize], alignment: &Alignment) -> Tree {
    let mut tree = Tree::new();
    tree.rooted = false;

    let root = tree.add_node(None, None, None);
    tree.root = root;

    for &taxon_idx in &taxa[..3] {
        let name = alignment.taxon_name(taxon_idx).to_string();
        tree.add_node(Some(name), None, Some(root));
    }

    tree
}

/// Insert a new taxon on the edge between parent_id and child_id.
fn insert_on_edge(
    tree: &Tree,
    parent_id: usize,
    child_id: usize,
    taxon_name: &str,
) -> Tree {
    let mut new_tree = Tree::new();
    new_tree.rooted = tree.rooted;

    let mut id_map = vec![0usize; tree.nodes.len()];

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

    let new_internal_id = new_tree.nodes.len();
    new_tree.nodes.push(crate::tree::types::TreeNode {
        id: new_internal_id,
        name: None,
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    let new_leaf_id = new_tree.nodes.len();
    new_tree.nodes.push(crate::tree::types::TreeNode {
        id: new_leaf_id,
        name: Some(taxon_name.to_string()),
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    for node in &tree.nodes {
        let mapped_id = id_map[node.id];
        for &old_child in &node.children {
            if node.id == parent_id && old_child == child_id {
                new_tree.nodes[mapped_id]
                    .children
                    .push(new_internal_id);
                new_tree.nodes[new_internal_id].parent = Some(mapped_id);

                let mapped_child = id_map[child_id];
                new_tree.nodes[new_internal_id]
                    .children
                    .push(mapped_child);
                new_tree.nodes[mapped_child].parent = Some(new_internal_id);

                new_tree.nodes[new_internal_id]
                    .children
                    .push(new_leaf_id);
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

/// Try inserting a new taxon at every possible edge, maximizing compatible sites.
fn best_insertion_compat(
    tree: &Tree,
    taxon_name: &str,
    alignment: &Alignment,
) -> Tree {
    let mut best_tree = None;
    let mut best_compat = 0;
    let mut first = true;

    let edges: Vec<(usize, usize)> = tree
        .nodes
        .iter()
        .filter_map(|n| n.parent.map(|p| (p, n.id)))
        .collect();

    for (parent_id, child_id) in edges {
        let candidate = insert_on_edge(tree, parent_id, child_id, taxon_name);
        let compat = count_compatible_sites(&candidate, alignment);
        if first || compat > best_compat {
            best_compat = compat;
            best_tree = Some(candidate);
            first = false;
        }
    }

    best_tree.unwrap_or_else(|| {
        let mut t = tree.clone();
        let root = t.root;
        t.add_node(Some(taxon_name.to_string()), None, Some(root));
        t
    })
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

/// Result of pruning a subtree.
#[derive(Debug, Clone)]
struct PrunedTree {
    tree: Tree,
    subtree: Tree,
}

/// Prune a subtree rooted at `prune_id` from the tree.
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
        subtree.nodes.push(crate::tree::types::TreeNode {
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
        main_tree.nodes.push(crate::tree::types::TreeNode {
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

    Some(PrunedTree {
        tree: main_tree,
        subtree,
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

    let graft_node_id = new_tree.nodes.len();
    new_tree.nodes.push(crate::tree::types::TreeNode {
        id: graft_node_id,
        name: None,
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    for node in &main.nodes {
        let mapped_id = main_map[node.id];
        for &child in &node.children {
            if node.id == graft_parent && child == graft_child {
                new_tree.nodes[mapped_id].children.push(graft_node_id);
                new_tree.nodes[graft_node_id].parent = Some(mapped_id);

                let mapped_child = main_map[child];
                new_tree.nodes[graft_node_id]
                    .children
                    .push(mapped_child);
                new_tree.nodes[mapped_child].parent = Some(graft_node_id);

                let subtree_root = sub_map[sub.root];
                new_tree.nodes[graft_node_id]
                    .children
                    .push(subtree_root);
                new_tree.nodes[subtree_root].parent = Some(graft_node_id);
            } else {
                let mapped_child = main_map[child];
                new_tree.nodes[mapped_id].children.push(mapped_child);
                new_tree.nodes[mapped_child].parent = Some(mapped_id);
            }
        }
    }

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

/// Perform one round of SPR rearrangement optimizing for compatible sites.
fn spr_round_compat(
    tree: &Tree,
    alignment: &Alignment,
    current_compat: usize,
) -> (Tree, usize, bool) {
    let mut best_tree = tree.clone();
    let mut best_compat = current_compat;
    let mut improved = false;

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
                let candidate =
                    regraft_subtree(&pruned, graft_parent, graft_child);
                let compat = count_compatible_sites(&candidate, alignment);
                if compat > best_compat {
                    best_compat = compat;
                    best_tree = candidate;
                    improved = true;
                }
            }
        }
    }

    (best_tree, best_compat, improved)
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

    // --- count_distinct_states tests ---

    #[test]
    fn test_count_distinct_states_single() {
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "A"),
            ("C", "A"),
        ]);
        assert_eq!(count_distinct_states(&aln, 0), 1);
    }

    #[test]
    fn test_count_distinct_states_two() {
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "T"),
            ("C", "A"),
        ]);
        assert_eq!(count_distinct_states(&aln, 0), 2);
    }

    #[test]
    fn test_count_distinct_states_four() {
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "C"),
            ("C", "G"),
            ("D", "T"),
        ]);
        assert_eq!(count_distinct_states(&aln, 0), 4);
    }

    #[test]
    fn test_count_distinct_states_with_ambiguity() {
        // N is ambiguous, should not be counted
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "N"),
            ("C", "T"),
        ]);
        assert_eq!(count_distinct_states(&aln, 0), 2); // only A and T
    }

    // --- count_compatible_sites tests ---

    #[test]
    fn test_invariant_sites_all_compatible() {
        let aln = make_alignment(&[
            ("A", "AAAA"),
            ("B", "AAAA"),
            ("C", "AAAA"),
            ("D", "AAAA"),
        ]);
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let compat = count_compatible_sites(&tree, &aln);
        assert_eq!(compat, 4); // all sites invariant -> all compatible
    }

    #[test]
    fn test_two_state_compatible_on_correct_tree() {
        // A,B share one state; C,D share another
        // On ((A,B),(C,D)): 1 change per site = n_distinct_states - 1 = 1
        let aln = make_alignment(&[
            ("A", "AAAA"),
            ("B", "AAAA"),
            ("C", "TTTT"),
            ("D", "TTTT"),
        ]);
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let compat = count_compatible_sites(&tree, &aln);
        assert_eq!(compat, 4); // all compatible on this tree
    }

    #[test]
    fn test_two_state_incompatible_on_wrong_tree() {
        // A,B share one state; C,D share another
        // On ((A,C),(B,D)): 2 changes needed, but minimum is 1 -> incompatible
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "A"),
            ("C", "T"),
            ("D", "T"),
        ]);
        let tree = parse_newick("((A,C),(B,D));").unwrap();
        let compat = count_compatible_sites(&tree, &aln);
        assert_eq!(compat, 0); // not compatible on this topology
    }

    #[test]
    fn test_homoplastic_site_incompatible() {
        // All 4 bases at a site on ((A,B),(C,D)):
        // Fitch score = 3, but minimum = 4-1 = 3 -> actually compatible!
        // (4 states, each appears once: min = 3 steps, Fitch gives 3)
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "C"),
            ("C", "G"),
            ("D", "T"),
        ]);
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let compat = count_compatible_sites(&tree, &aln);
        assert_eq!(compat, 1); // 4 states, min changes = 3, Fitch = 3 -> compatible

        // Now test a truly homoplastic case: 2 states, but needs 2 changes
        let aln2 = make_alignment(&[
            ("A", "A"),
            ("B", "T"),
            ("C", "A"),
            ("D", "T"),
        ]);
        // On ((A,B),(C,D)): {A}&{T}=empty->union,+1; {A}&{T}=empty->union,+1;
        // {A,T}&{A,T}={A,T},+0 -> total=2 steps
        // Distinct states = 2, min = 1. Fitch = 2 > 1 -> incompatible
        let compat2 = count_compatible_sites(&tree, &aln2);
        assert_eq!(compat2, 0);
    }

    #[test]
    fn test_mixed_sites_compatible_and_incompatible() {
        let aln = make_alignment(&[
            // site 0: A,A,T,T -> compatible on ((A,B),(C,D))
            // site 1: A,T,A,T -> incompatible on ((A,B),(C,D))
            // site 2: A,A,A,A -> invariant, compatible
            ("A", "AAA"),
            ("B", "ATA"),
            ("C", "TAA"),
            ("D", "TTA"),
        ]);
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let compat = count_compatible_sites(&tree, &aln);
        assert_eq!(compat, 2); // sites 0 and 2 are compatible
    }

    // --- dna_compat_search tests ---

    #[test]
    fn test_search_finds_tree_with_max_compat() {
        // Clear signal: A,B share one state; C,D share another
        let aln = make_alignment(&[
            ("A", "AAAACCCC"),
            ("B", "AAAACCCC"),
            ("C", "CCCCAAAA"),
            ("D", "CCCCAAAA"),
        ]);
        let result = dna_compat_search(&aln, Some(42));
        assert_eq!(result.tree.num_leaves(), 4);
        // All 8 sites should be compatible on the correct tree
        assert_eq!(result.compatible_sites, 8);
        assert_eq!(result.total_sites, 8);
        assert!((result.compatibility_fraction - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_search_three_taxa() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "TGCA"),
        ]);
        let result = dna_compat_search(&aln, Some(42));
        assert_eq!(result.tree.num_leaves(), 3);
        // With 3 taxa, any binary tree gives the same parsimony score
        // All 2-state sites should be compatible
        assert!(result.compatible_sites > 0);
    }

    #[test]
    fn test_search_all_identical() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "ACGT"),
            ("D", "ACGT"),
        ]);
        let result = dna_compat_search(&aln, None);
        // All sites invariant -> all compatible
        assert_eq!(result.compatible_sites, 4);
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
        let result = dna_compat_search(&aln, Some(99));
        assert_eq!(result.tree.num_leaves(), 5);
        assert!(result.compatible_sites > 0);
        assert!(result.compatibility_fraction >= 0.0);
        assert!(result.compatibility_fraction <= 1.0);
    }

    #[test]
    fn test_search_deterministic_with_seed() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGAACGA"),
            ("C", "TGCATGCA"),
            ("D", "TGCTTGCT"),
        ]);
        let r1 = dna_compat_search(&aln, Some(42));
        let r2 = dna_compat_search(&aln, Some(42));
        assert_eq!(r1.compatible_sites, r2.compatible_sites);
    }

    #[test]
    fn test_search_prefers_more_compatible_tree() {
        // A,B share states at most sites; C,D share states at most sites
        // Search should find ((A,B),(C,D)) which makes more sites compatible
        let aln = make_alignment(&[
            ("A", "AAAAAA"),
            ("B", "AAAAAA"),
            ("C", "TTTTTT"),
            ("D", "TTTTTT"),
        ]);
        let result = dna_compat_search(&aln, Some(42));
        assert_eq!(result.compatible_sites, 6);

        // Compare with the "wrong" tree
        let wrong_tree = parse_newick("((A,C),(B,D));").unwrap();
        let wrong_compat = count_compatible_sites(&wrong_tree, &aln);
        assert!(result.compatible_sites >= wrong_compat);
    }
}

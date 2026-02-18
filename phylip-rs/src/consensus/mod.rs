//! Consensus tree construction methods.
//!
//! Implements strict consensus, majority-rule consensus, extended majority-rule
//! consensus, and threshold consensus for summarizing sets of phylogenetic trees.
//!
//! These methods take a collection of trees (typically bootstrap replicates) and
//! produce a single summary tree that represents the areas of agreement.
//!
//! # Consensus methods
//!
//! - **Strict**: Only includes splits present in *every* input tree
//! - **Majority rule**: Includes splits present in >50% of input trees
//! - **Extended majority rule** (default): Majority rule plus compatible splits
//!   added greedily in decreasing frequency order
//! - **Threshold**: Includes splits present in >= threshold fraction of trees
//!
//! # Example
//!
//! ```
//! use phylip_rs::consensus::{consensus_tree, ConsensusMethod};
//! use phylip_rs::tree::newick;
//!
//! let trees: Vec<_> = vec![
//!     "((A,B),(C,D));",
//!     "((A,B),(C,D));",
//!     "((A,C),(B,D));",
//! ].into_iter().map(|s| newick::parse_newick(s).unwrap()).collect();
//!
//! let result = consensus_tree(&trees, &ConsensusMethod::MajorityRule).unwrap();
//! ```
//!
//! Reference: Felsenstein, J. (2004). Inferring Phylogenies, Chapter 30.

use crate::tree::splits::Split;
use crate::tree::Tree;
use std::collections::{HashMap, HashSet};
use std::fmt;

// ============================================================
// Consensus method configuration
// ============================================================

/// Method for computing the consensus tree.
#[derive(Debug, Clone)]
pub enum ConsensusMethod {
    /// Include only splits present in ALL input trees.
    Strict,
    /// Include splits present in >50% of input trees.
    MajorityRule,
    /// Majority rule plus compatible splits in decreasing frequency.
    /// This is PHYLIP's default and produces a fully resolved tree when possible.
    ExtendedMajorityRule,
    /// Include splits present in >= threshold fraction of trees.
    /// Threshold must be in (0.0, 1.0].
    Threshold(f64),
}

/// Errors during consensus tree computation.
#[derive(Debug, Clone)]
pub enum ConsensusError {
    /// No input trees provided.
    NoTrees,
    /// Trees have different leaf sets.
    IncompatibleTaxa,
    /// Invalid threshold value.
    InvalidThreshold(f64),
}

impl fmt::Display for ConsensusError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ConsensusError::NoTrees => write!(f, "no input trees provided"),
            ConsensusError::IncompatibleTaxa => write!(f, "trees have different leaf sets"),
            ConsensusError::InvalidThreshold(t) => {
                write!(f, "invalid threshold {}: must be in (0.0, 1.0]", t)
            }
        }
    }
}

impl std::error::Error for ConsensusError {}

// Split type is imported from crate::tree::splits::Split

// ============================================================
// Split frequency counting
// ============================================================

/// A split with its frequency count and optional average branch length.
#[derive(Debug, Clone)]
struct SplitInfo {
    split: Split,
    count: usize,
    branch_length_sum: f64,
}

/// Extract all non-trivial splits from a tree, given a taxon-to-index mapping.
///
/// Delegates to the shared `crate::tree::splits::extract_splits` function,
/// converting the taxon index map into the sorted taxon order it expects.
fn consensus_extract_splits(
    tree: &Tree,
    taxon_index: &HashMap<String, usize>,
    ntaxa: usize,
) -> Vec<(Split, Option<f64>)> {
    // Build sorted taxon order from the index map
    let mut taxon_order = vec![String::new(); ntaxa];
    for (name, &idx) in taxon_index {
        taxon_order[idx] = name.clone();
    }
    crate::tree::splits::extract_splits(tree, &taxon_order)
}

/// Count split frequencies across multiple trees.
fn count_splits(
    trees: &[Tree],
    taxon_index: &HashMap<String, usize>,
    ntaxa: usize,
) -> Vec<SplitInfo> {
    let mut split_map: HashMap<Split, (usize, f64)> = HashMap::new();

    for tree in trees {
        let tree_splits = consensus_extract_splits(tree, taxon_index, ntaxa);
        // Deduplicate within a single tree (same split from both children of root)
        let mut seen = HashSet::new();
        for (split, bl) in tree_splits {
            if seen.insert(split.clone()) {
                let entry = split_map.entry(split).or_insert((0, 0.0));
                entry.0 += 1;
                if let Some(len) = bl {
                    entry.1 += len;
                }
            }
        }
    }

    split_map
        .into_iter()
        .map(|(split, (count, bl_sum))| SplitInfo {
            split,
            count,
            branch_length_sum: bl_sum,
        })
        .collect()
}

// ============================================================
// Consensus tree building
// ============================================================

/// Result of consensus tree computation.
#[derive(Debug)]
pub struct ConsensusResult {
    /// The consensus tree.
    pub tree: Tree,
    /// Number of input trees.
    pub ntrees: usize,
    /// Split frequencies: (taxon names on one side, count, proportion).
    pub split_frequencies: Vec<(Vec<String>, usize, f64)>,
}

/// Compute a consensus tree from a set of input trees.
///
/// All input trees must have the same leaf set (taxon names).
pub fn consensus_tree(
    trees: &[Tree],
    method: &ConsensusMethod,
) -> Result<ConsensusResult, ConsensusError> {
    if trees.is_empty() {
        return Err(ConsensusError::NoTrees);
    }

    if let ConsensusMethod::Threshold(t) = method {
        if *t <= 0.0 || *t > 1.0 {
            return Err(ConsensusError::InvalidThreshold(*t));
        }
    }

    // Get taxon names from the first tree
    let mut taxa: Vec<String> = trees[0]
        .leaves()
        .iter()
        .filter_map(|n| n.name.clone())
        .collect();
    taxa.sort();

    let ntaxa = taxa.len();

    // Build taxon-to-index mapping
    let taxon_index: HashMap<String, usize> = taxa
        .iter()
        .enumerate()
        .map(|(i, name)| (name.clone(), i))
        .collect();

    // Verify all trees have the same taxa
    for tree in trees.iter().skip(1) {
        let mut tree_taxa: Vec<String> = tree
            .leaves()
            .iter()
            .filter_map(|n| n.name.clone())
            .collect();
        tree_taxa.sort();
        if tree_taxa != taxa {
            return Err(ConsensusError::IncompatibleTaxa);
        }
    }

    let ntrees = trees.len();

    // Count split frequencies
    let mut split_infos = count_splits(trees, &taxon_index, ntaxa);

    // Sort by frequency (descending), then by split size for stability
    split_infos.sort_by(|a, b| {
        b.count
            .cmp(&a.count)
            .then_with(|| a.split.count().cmp(&b.split.count()))
    });

    // Filter splits based on consensus method
    let accepted = select_splits(&split_infos, ntrees, method);

    // Build split_frequencies output
    let split_frequencies: Vec<(Vec<String>, usize, f64)> = accepted
        .iter()
        .map(|info| {
            let names: Vec<String> = info
                .split
                .taxa_indices()
                .iter()
                .map(|&i| taxa[i].clone())
                .collect();
            (names, info.count, info.count as f64 / ntrees as f64)
        })
        .collect();

    // Build the consensus tree from accepted splits
    let tree = build_tree_from_splits(&taxa, &accepted, ntrees);

    Ok(ConsensusResult {
        tree,
        ntrees,
        split_frequencies,
    })
}

/// Select splits based on the consensus method.
fn select_splits(
    split_infos: &[SplitInfo],
    ntrees: usize,
    method: &ConsensusMethod,
) -> Vec<SplitInfo> {
    match method {
        ConsensusMethod::Strict => split_infos
            .iter()
            .filter(|s| s.count == ntrees)
            .cloned()
            .collect(),

        ConsensusMethod::MajorityRule => split_infos
            .iter()
            .filter(|s| s.count as f64 > ntrees as f64 / 2.0)
            .cloned()
            .collect(),

        ConsensusMethod::Threshold(threshold) => split_infos
            .iter()
            .filter(|s| s.count as f64 >= *threshold * ntrees as f64)
            .cloned()
            .collect(),

        ConsensusMethod::ExtendedMajorityRule => {
            // Start with majority rule splits
            let mut accepted: Vec<SplitInfo> = split_infos
                .iter()
                .filter(|s| s.count as f64 > ntrees as f64 / 2.0)
                .cloned()
                .collect();

            // Greedily add compatible splits in decreasing frequency
            for info in split_infos {
                if info.count as f64 > ntrees as f64 / 2.0 {
                    continue; // Already included
                }

                // Check compatibility with all already-accepted splits
                let compatible = accepted
                    .iter()
                    .all(|a| a.split.is_compatible(&info.split));

                if compatible {
                    accepted.push(info.clone());
                }
            }

            accepted
        }
    }
}

/// Build a tree from a set of accepted splits.
///
/// Uses a top-down approach: start with a star tree (all taxa connected
/// to root), then refine by adding splits one at a time. Each split
/// groups a subset of taxa into a new internal node.
fn build_tree_from_splits(taxa: &[String], splits: &[SplitInfo], ntrees: usize) -> Tree {
    let ntaxa = taxa.len();
    let mut tree = Tree::new();
    tree.rooted = false;

    // Create root
    let root = tree.add_node(None, None, None);
    tree.root = root;

    // Create leaf nodes, all initially children of root
    let mut leaf_ids = Vec::with_capacity(ntaxa);
    for name in taxa {
        let leaf = tree.add_node(Some(name.clone()), None, Some(root));
        leaf_ids.push(leaf);
    }

    // Sort splits by size (smallest first) for proper nesting
    let mut sorted_splits: Vec<&SplitInfo> = splits.iter().collect();
    sorted_splits.sort_by_key(|s| s.split.count());

    // For each split, create an internal node that groups the relevant taxa
    for info in &sorted_splits {
        let taxa_in_split: Vec<usize> = info.split.taxa_indices();

        // Find the lowest common ancestor node of these taxa in the current tree.
        // Because we process smallest splits first, the taxa in this split should
        // currently all be children of the same parent.
        // Find which nodes (leaves or internal) need to be grouped.
        let nodes_to_group: Vec<usize> = find_nodes_to_group(&tree, &leaf_ids, &taxa_in_split);

        if nodes_to_group.len() < 2 {
            continue; // Nothing to group
        }

        // Find their common parent
        let parent = tree.nodes[nodes_to_group[0]].parent;
        if parent.is_none() {
            continue;
        }
        let parent_id = parent.unwrap();

        // Verify all nodes share the same parent
        if !nodes_to_group
            .iter()
            .all(|&n| tree.nodes[n].parent == Some(parent_id))
        {
            continue; // Split not compatible with current tree topology
        }

        // Create a new internal node
        let freq_label = format!("{:.1}", info.count as f64 / ntrees as f64 * 100.0);
        let avg_bl = if info.branch_length_sum > 0.0 && info.count > 0 {
            Some(info.branch_length_sum / info.count as f64)
        } else {
            None
        };
        let new_node = tree.add_node(Some(freq_label), avg_bl, Some(parent_id));

        // Move the grouped nodes under the new internal node
        for &node_id in &nodes_to_group {
            // Remove from parent's children
            tree.nodes[parent_id]
                .children
                .retain(|&c| c != node_id);
            // Add to new node
            tree.nodes[new_node].children.push(node_id);
            tree.nodes[node_id].parent = Some(new_node);
        }
    }

    tree
}

/// Find which current tree nodes correspond to the taxa in a split.
///
/// For each taxon index in the split, we need the "effective node" — either
/// the leaf itself, or the internal node that already groups a subset containing it.
fn find_nodes_to_group(tree: &Tree, leaf_ids: &[usize], taxa_indices: &[usize]) -> Vec<usize> {
    // For each taxon in the split, walk up to find the child of its current ancestor
    // that contains it. We need to find the set of distinct immediate children of some
    // common ancestor that collectively contain exactly these taxa.

    // Simple approach: collect the set of nodes that are ancestors of these taxa
    // at the lowest level where they all share a parent.
    let target_leaves: HashSet<usize> = taxa_indices.iter().map(|&i| leaf_ids[i]).collect();

    // For each target leaf, find the highest ancestor that is entirely within the target set
    let mut effective_nodes: HashSet<usize> = HashSet::new();

    for &leaf in &target_leaves {
        let mut node = leaf;
        loop {
            let parent = match tree.nodes[node].parent {
                Some(p) => p,
                None => break,
            };
            // Check if all leaves under `parent` are in target_leaves
            let parent_leaves = collect_leaves(tree, parent);
            if parent_leaves.iter().all(|l| target_leaves.contains(l)) {
                // Parent is entirely within target, go up
                node = parent;
            } else {
                break;
            }
        }
        effective_nodes.insert(node);
    }

    effective_nodes.into_iter().collect()
}

/// Collect all leaf node IDs under a given node.
fn collect_leaves(tree: &Tree, node_id: usize) -> Vec<usize> {
    let mut leaves = Vec::new();
    let mut stack = vec![node_id];
    while let Some(nid) = stack.pop() {
        let node = &tree.nodes[nid];
        if node.is_leaf() {
            leaves.push(nid);
        } else {
            for &child in &node.children {
                stack.push(child);
            }
        }
    }
    leaves
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::newick;

    fn parse(s: &str) -> Tree {
        newick::parse_newick(s).unwrap()
    }

    fn leaf_set(tree: &Tree) -> Vec<String> {
        let mut names: Vec<String> = tree
            .leaves()
            .iter()
            .filter_map(|n| n.name.clone())
            .collect();
        names.sort();
        names
    }

    // --- Split tests ---

    #[test]
    fn test_split_set_get() {
        let mut s = Split::new(8);
        s.set(0);
        s.set(3);
        s.set(7);
        assert!(s.get(0));
        assert!(!s.get(1));
        assert!(s.get(3));
        assert!(s.get(7));
        assert_eq!(s.count(), 3);
    }

    #[test]
    fn test_split_canonicalize() {
        let mut s1 = Split::new(4);
        s1.set(0);
        s1.set(1);
        s1.canonicalize();
        // Taxon 0 was set, so it gets complemented: now {2, 3}
        assert!(!s1.get(0));
        assert!(!s1.get(1));
        assert!(s1.get(2));
        assert!(s1.get(3));
    }

    #[test]
    fn test_split_canonical_equality() {
        // {0,1} | {2,3} and {2,3} | {0,1} should be equal after canonicalization
        let mut s1 = Split::new(4);
        s1.set(0);
        s1.set(1);
        s1.canonicalize();

        let mut s2 = Split::new(4);
        s2.set(2);
        s2.set(3);
        s2.canonicalize();

        assert_eq!(s1, s2);
    }

    #[test]
    fn test_split_complement() {
        let mut s = Split::new(4);
        s.set(1);
        s.complement();
        assert!(s.get(0));
        assert!(!s.get(1));
        assert!(s.get(2));
        assert!(s.get(3));
        assert_eq!(s.count(), 3);
    }

    #[test]
    fn test_split_trivial() {
        let mut s = Split::new(6);
        s.set(2);
        assert!(s.is_trivial()); // Only 1 taxon

        let mut s2 = Split::new(6);
        for i in 0..5 {
            s2.set(i);
        }
        assert!(s2.is_trivial()); // 5 out of 6
    }

    #[test]
    fn test_split_compatible_disjoint() {
        // {1,2} | {0,3,4} and {3,4} | {0,1,2} — compatible (A∩B = {}, A'∩B' = {0})
        let mut s1 = Split::new(5);
        s1.set(1);
        s1.set(2);
        s1.canonicalize();

        let mut s2 = Split::new(5);
        s2.set(3);
        s2.set(4);
        s2.canonicalize();

        assert!(s1.is_compatible(&s2));
    }

    #[test]
    fn test_split_compatible_nested() {
        // {1} | {0,2,3,4} and {1,2} | {0,3,4} — compatible ({1} ⊂ {1,2})
        let mut s1 = Split::new(5);
        s1.set(1);
        s1.canonicalize();

        let mut s2 = Split::new(5);
        s2.set(1);
        s2.set(2);
        s2.canonicalize();

        assert!(s1.is_compatible(&s2));
    }

    #[test]
    fn test_split_incompatible() {
        // {1,2} | {0,3,4} and {2,3} | {0,1,4} — not compatible
        // A∩B = {2}, A∩B' = {1}, A'∩B = {3}, A'∩B' = {0,4}
        let mut s1 = Split::new(5);
        s1.set(1);
        s1.set(2);
        s1.canonicalize();

        let mut s2 = Split::new(5);
        s2.set(2);
        s2.set(3);
        s2.canonicalize();

        assert!(!s1.is_compatible(&s2));
    }

    // --- Extract splits tests ---

    #[test]
    fn test_extract_splits_4taxa() {
        let tree = parse("((A,B),(C,D));");
        let taxa: Vec<String> = vec!["A".into(), "B".into(), "C".into(), "D".into()];
        let taxon_index: HashMap<String, usize> =
            taxa.iter().enumerate().map(|(i, n)| (n.clone(), i)).collect();

        let splits = consensus_extract_splits(&tree, &taxon_index, 4);
        // Should have 1 unique non-trivial split (after dedup, may see 2 raw)
        let unique: HashSet<Split> = splits.into_iter().map(|(s, _)| s).collect();
        assert_eq!(unique.len(), 1, "4-taxon symmetric tree has 1 non-trivial split");
    }

    #[test]
    fn test_extract_splits_5taxa() {
        let tree = parse("(((A,B),C),(D,E));");
        let mut taxa: Vec<String> = tree
            .leaves()
            .iter()
            .filter_map(|n| n.name.clone())
            .collect();
        taxa.sort();
        let taxon_index: HashMap<String, usize> =
            taxa.iter().enumerate().map(|(i, n)| (n.clone(), i)).collect();

        let splits = consensus_extract_splits(&tree, &taxon_index, 5);
        let unique: HashSet<Split> = splits.into_iter().map(|(s, _)| s).collect();
        assert!(unique.len() >= 2, "5-taxon tree should have at least 2 non-trivial splits");
    }

    // --- Consensus tree tests ---

    #[test]
    fn test_strict_consensus_identical() {
        let trees: Vec<Tree> = (0..10).map(|_| parse("((A,B),(C,D));")).collect();
        let result = consensus_tree(&trees, &ConsensusMethod::Strict).unwrap();

        assert_eq!(result.ntrees, 10);
        assert_eq!(leaf_set(&result.tree), vec!["A", "B", "C", "D"]);
        assert_eq!(result.split_frequencies.len(), 1);
        assert_eq!(result.split_frequencies[0].1, 10); // count = 10
    }

    #[test]
    fn test_strict_consensus_conflicting() {
        let mut trees = Vec::new();
        for _ in 0..5 {
            trees.push(parse("((A,B),(C,D));"));
        }
        for _ in 0..5 {
            trees.push(parse("((A,C),(B,D));"));
        }

        let result = consensus_tree(&trees, &ConsensusMethod::Strict).unwrap();
        // No split is in all 10 trees
        assert_eq!(result.split_frequencies.len(), 0);
    }

    #[test]
    fn test_majority_rule_basic() {
        let mut trees = Vec::new();
        for _ in 0..7 {
            trees.push(parse("((A,B),(C,D));"));
        }
        for _ in 0..3 {
            trees.push(parse("((A,C),(B,D));"));
        }

        let result = consensus_tree(&trees, &ConsensusMethod::MajorityRule).unwrap();
        assert_eq!(result.split_frequencies.len(), 1);
        assert_eq!(result.split_frequencies[0].1, 7);
    }

    #[test]
    fn test_majority_rule_both_above_50() {
        // 5 taxa, two independent splits both above 50%
        let mut trees = Vec::new();
        // ((A,B),(C,(D,E))) - has splits {A,B} and {D,E}
        for _ in 0..8 {
            trees.push(parse("(((A,B),C),(D,E));"));
        }
        // ((A,C),(B,(D,E))) - has split {D,E} but not {A,B}
        for _ in 0..2 {
            trees.push(parse("(((A,C),B),(D,E));"));
        }

        let result = consensus_tree(&trees, &ConsensusMethod::MajorityRule).unwrap();
        // {A,B} appears 8/10 = 80%, {D,E} appears 10/10 = 100%
        assert_eq!(result.split_frequencies.len(), 2);
    }

    #[test]
    fn test_threshold_consensus() {
        let mut trees = Vec::new();
        for _ in 0..9 {
            trees.push(parse("((A,B),(C,D));"));
        }
        trees.push(parse("((A,C),(B,D));"));

        // 90% threshold: {A,B}|{C,D} appears 90% -> included
        let result = consensus_tree(&trees, &ConsensusMethod::Threshold(0.9)).unwrap();
        assert_eq!(result.split_frequencies.len(), 1);

        // 95% threshold: {A,B}|{C,D} appears 90% -> excluded
        let result2 = consensus_tree(&trees, &ConsensusMethod::Threshold(0.95)).unwrap();
        assert_eq!(result2.split_frequencies.len(), 0);
    }

    #[test]
    fn test_extended_majority_rule() {
        // 5 taxa: majority rule split + a weaker compatible split
        let mut trees = Vec::new();
        // {A,B} at 80%, {D,E} at 40% (compatible with {A,B})
        for _ in 0..8 {
            trees.push(parse("(((A,B),C),(D,E));"));
        }
        for _ in 0..2 {
            trees.push(parse("(((A,C),D),(B,E));"));
        }

        let mr_result = consensus_tree(&trees, &ConsensusMethod::MajorityRule).unwrap();
        let mre_result = consensus_tree(&trees, &ConsensusMethod::ExtendedMajorityRule).unwrap();

        // MRe should have at least as many splits as MR
        assert!(
            mre_result.split_frequencies.len() >= mr_result.split_frequencies.len(),
            "MRe should have >= MR splits"
        );
    }

    #[test]
    fn test_extended_majority_rule_rejects_incompatible() {
        // Ensure MRe doesn't add incompatible splits
        let mut trees = Vec::new();
        // 6 trees: ((A,B),(C,D)); -> {A,B}|{C,D} at 60%
        for _ in 0..6 {
            trees.push(parse("((A,B),(C,D));"));
        }
        // 4 trees: ((A,C),(B,D)); -> {A,C}|{B,D} at 40%
        for _ in 0..4 {
            trees.push(parse("((A,C),(B,D));"));
        }

        let result = consensus_tree(&trees, &ConsensusMethod::ExtendedMajorityRule).unwrap();
        // {A,B}|{C,D} and {A,C}|{B,D} are incompatible, only the majority one should be kept
        assert_eq!(
            result.split_frequencies.len(),
            1,
            "incompatible minority split should be rejected by MRe"
        );
    }

    #[test]
    fn test_consensus_tree_has_correct_leaves() {
        let trees: Vec<Tree> = (0..5).map(|_| parse("((A,B),(C,D));")).collect();
        let result = consensus_tree(&trees, &ConsensusMethod::Strict).unwrap();
        assert_eq!(leaf_set(&result.tree), vec!["A", "B", "C", "D"]);
    }

    #[test]
    fn test_consensus_empty_input() {
        let trees: Vec<Tree> = Vec::new();
        let result = consensus_tree(&trees, &ConsensusMethod::Strict);
        assert!(result.is_err());
    }

    #[test]
    fn test_consensus_incompatible_taxa() {
        let trees = vec![
            parse("((A,B),(C,D));"),
            parse("((A,B),(C,E));"), // E instead of D
        ];
        let result = consensus_tree(&trees, &ConsensusMethod::Strict);
        assert!(result.is_err());
    }

    #[test]
    fn test_consensus_star_tree() {
        // All conflicting -> strict consensus is a star tree
        let trees = vec![
            parse("((A,B),(C,D));"),
            parse("((A,C),(B,D));"),
            parse("((A,D),(B,C));"),
        ];
        let result = consensus_tree(&trees, &ConsensusMethod::Strict).unwrap();
        assert_eq!(result.split_frequencies.len(), 0, "fully conflicting trees -> star tree");
        assert_eq!(leaf_set(&result.tree), vec!["A", "B", "C", "D"]);
    }

    #[test]
    fn test_split_frequency_proportions() {
        let mut trees = Vec::new();
        for _ in 0..75 {
            trees.push(parse("((A,B),(C,D));"));
        }
        for _ in 0..25 {
            trees.push(parse("((A,C),(B,D));"));
        }

        let result = consensus_tree(&trees, &ConsensusMethod::MajorityRule).unwrap();
        assert_eq!(result.split_frequencies.len(), 1);
        let (_, count, proportion) = &result.split_frequencies[0];
        assert_eq!(*count, 75);
        assert!((proportion - 0.75).abs() < 1e-10);
    }

    #[test]
    fn test_invalid_threshold() {
        let trees = vec![parse("((A,B),(C,D));")];
        assert!(consensus_tree(&trees, &ConsensusMethod::Threshold(0.0)).is_err());
        assert!(consensus_tree(&trees, &ConsensusMethod::Threshold(1.5)).is_err());
        assert!(consensus_tree(&trees, &ConsensusMethod::Threshold(-0.1)).is_err());
    }

    #[test]
    fn test_consensus_5taxa_multiple_splits() {
        let trees: Vec<Tree> = (0..100).map(|_| parse("(((A,B),C),(D,E));")).collect();
        let result = consensus_tree(&trees, &ConsensusMethod::Strict).unwrap();

        // Should have 2 splits: {A,B} and {D,E}
        assert_eq!(result.split_frequencies.len(), 2);
        for (_, count, _) in &result.split_frequencies {
            assert_eq!(*count, 100);
        }
    }

    #[test]
    fn test_consensus_6taxa() {
        let trees: Vec<Tree> = (0..50)
            .map(|_| parse("(((A,B),(C,D)),(E,F));"))
            .collect();
        let result = consensus_tree(&trees, &ConsensusMethod::Strict).unwrap();

        assert_eq!(leaf_set(&result.tree), vec!["A", "B", "C", "D", "E", "F"]);
        // Should have 3 splits: {A,B}, {C,D}, {E,F}
        assert_eq!(result.split_frequencies.len(), 3);
    }
}

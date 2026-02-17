//! Bootstrap support values for phylogenetic tree branches.
//!
//! Computes bootstrap proportions by comparing bipartitions (splits)
//! between the original tree and replicate trees. Each branch in the
//! original tree defines a bipartition of the leaf set; the bootstrap
//! support for that branch is the fraction of replicate trees that
//! contain the same bipartition.

use crate::tree::Tree;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// A bipartition (split) of the leaf set induced by removing a branch.
///
/// Represented as the sorted set of leaf names on the smaller side
/// of the split. This canonical form ensures that equivalent splits
/// are represented identically regardless of rooting.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Bipartition {
    /// Leaf names on the smaller side of the split (sorted).
    taxa: Vec<String>,
}

impl Bipartition {
    /// Create a bipartition from a set of leaf names.
    ///
    /// Canonicalizes by always storing the smaller side of the split.
    /// If `taxa` contains more than half the total leaves, the
    /// complement is stored instead.
    fn new(taxa: HashSet<String>, all_taxa: &[String]) -> Self {
        let complement: HashSet<String> = all_taxa
            .iter()
            .filter(|t| !taxa.contains(*t))
            .cloned()
            .collect();

        // Always store the smaller side for canonical representation.
        // When sides are equal size, pick the one whose lexicographic
        // minimum is smaller, ensuring {A,B}|{C,D} and {C,D}|{A,B}
        // produce the same canonical form.
        let chosen = if taxa.len() < complement.len() {
            taxa
        } else if taxa.len() > complement.len() {
            complement
        } else {
            // Equal size: pick side with lexicographically smaller minimum element
            let min_taxa = taxa.iter().min().cloned();
            let min_comp = complement.iter().min().cloned();
            if min_taxa <= min_comp {
                taxa
            } else {
                complement
            }
        };

        let mut sorted: Vec<String> = chosen.into_iter().collect();
        sorted.sort();
        Self { taxa: sorted }
    }

    /// Returns true if this is a trivial split (single leaf or all-but-one).
    pub fn is_trivial(&self) -> bool {
        self.taxa.len() <= 1
    }
}

impl fmt::Display for Bipartition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{{{}}}", self.taxa.join(", "))
    }
}

/// Bootstrap support values for each non-trivial branch in a tree.
#[derive(Debug, Clone)]
pub struct BootstrapSupport {
    /// Map from bipartition to (support count, total replicates).
    support: HashMap<Bipartition, (usize, usize)>,
}

impl BootstrapSupport {
    /// Compute the bootstrap support proportion for a bipartition.
    pub fn get(&self, bipartition: &Bipartition) -> Option<f64> {
        self.support
            .get(bipartition)
            .map(|&(count, total)| count as f64 / total as f64)
    }

    /// Returns all bipartitions with their support values.
    pub fn all_supports(&self) -> Vec<(&Bipartition, f64)> {
        self.support
            .iter()
            .map(|(bp, &(count, total))| (bp, count as f64 / total as f64))
            .collect()
    }

    /// Returns the number of non-trivial bipartitions tracked.
    pub fn num_bipartitions(&self) -> usize {
        self.support.len()
    }
}

/// Extract all non-trivial bipartitions from a tree.
///
/// Each internal branch defines a bipartition of the leaf set.
/// We skip trivial bipartitions (pendant branches leading to single leaves)
/// and the root split (which is an artifact of rooting).
pub fn tree_bipartitions(tree: &Tree) -> Vec<Bipartition> {
    let all_taxa: Vec<String> = tree
        .leaves()
        .iter()
        .filter_map(|n| n.name.clone())
        .collect();

    if all_taxa.len() < 4 {
        // Trees with <4 taxa have no non-trivial bipartitions
        return Vec::new();
    }

    // For each node, compute the set of leaf descendants
    let mut descendant_leaves: Vec<HashSet<String>> = vec![HashSet::new(); tree.num_nodes()];

    for &node_id in &tree.postorder() {
        let node = &tree.nodes[node_id];
        if node.is_leaf() {
            if let Some(ref name) = node.name {
                descendant_leaves[node_id].insert(name.clone());
            }
        } else {
            let child_leaves: HashSet<String> = node
                .children
                .iter()
                .flat_map(|&c| descendant_leaves[c].iter().cloned())
                .collect();
            descendant_leaves[node_id] = child_leaves;
        }
    }

    // Extract bipartitions from internal branches, deduplicating.
    // In a rooted tree, both children of the root define the same bipartition,
    // so we use a HashSet to collect unique splits.
    let mut seen = HashSet::new();
    let mut bipartitions = Vec::new();
    let total_leaves = all_taxa.len();

    for node_id in 0..tree.num_nodes() {
        let node = &tree.nodes[node_id];

        // Skip leaves (trivial splits) and root
        if node.is_leaf() || node_id == tree.root {
            continue;
        }

        let desc = &descendant_leaves[node_id];
        let desc_count = desc.len();

        // Skip trivial (single-leaf) and total (all-leaf) splits
        if desc_count <= 1 || desc_count >= total_leaves - 1 {
            continue;
        }

        let bp = Bipartition::new(desc.clone(), &all_taxa);
        if seen.insert(bp.clone()) {
            bipartitions.push(bp);
        }
    }

    bipartitions
}

/// Compute bootstrap support for a reference tree given replicate trees.
///
/// For each non-trivial bipartition in the reference tree, counts how
/// many replicate trees contain the same bipartition. The support value
/// is the fraction of replicates containing the split.
///
/// # Arguments
///
/// * `reference` - The tree to annotate with support values
/// * `replicates` - Trees inferred from bootstrap replicate datasets
///
/// # Returns
///
/// A `BootstrapSupport` mapping bipartitions to their support proportions.
pub fn compute_support(reference: &Tree, replicates: &[Tree]) -> BootstrapSupport {
    let ref_bipartitions = tree_bipartitions(reference);
    let nreps = replicates.len();

    // Count how often each reference bipartition appears in replicates
    let mut counts: HashMap<Bipartition, usize> = HashMap::new();
    for bp in &ref_bipartitions {
        counts.insert(bp.clone(), 0);
    }

    let ref_set: HashSet<&Bipartition> = ref_bipartitions.iter().collect();

    for rep_tree in replicates {
        let rep_bps = tree_bipartitions(rep_tree);
        let rep_set: HashSet<Bipartition> = rep_bps.into_iter().collect();

        for bp in &ref_set {
            if rep_set.contains(*bp) {
                *counts.get_mut(*bp).unwrap() += 1;
            }
        }
    }

    let support = counts
        .into_iter()
        .map(|(bp, count)| (bp, (count, nreps)))
        .collect();

    BootstrapSupport { support }
}

/// Map bootstrap support values onto tree branch labels.
///
/// Returns a new tree with bootstrap support percentages (0-100)
/// stored as branch labels on internal nodes.
pub fn label_tree_with_support(tree: &Tree, support: &BootstrapSupport) -> Tree {
    let all_taxa: Vec<String> = tree
        .leaves()
        .iter()
        .filter_map(|n| n.name.clone())
        .collect();

    // Compute descendant leaves for each node
    let mut descendant_leaves: Vec<HashSet<String>> = vec![HashSet::new(); tree.num_nodes()];
    for &node_id in &tree.postorder() {
        let node = &tree.nodes[node_id];
        if node.is_leaf() {
            if let Some(ref name) = node.name {
                descendant_leaves[node_id].insert(name.clone());
            }
        } else {
            let child_leaves: HashSet<String> = node
                .children
                .iter()
                .flat_map(|&c| descendant_leaves[c].iter().cloned())
                .collect();
            descendant_leaves[node_id] = child_leaves;
        }
    }

    let mut new_tree = tree.clone();

    for node_id in 0..new_tree.num_nodes() {
        let node = &new_tree.nodes[node_id];
        if node.is_leaf() || node_id == new_tree.root {
            continue;
        }

        let desc = &descendant_leaves[node_id];
        if desc.len() <= 1 || desc.len() >= all_taxa.len() - 1 {
            continue;
        }

        let bp = Bipartition::new(desc.clone(), &all_taxa);
        if let Some(pct) = support.get(&bp) {
            let label = format!("{:.0}", pct * 100.0);
            new_tree.nodes[node_id].name = Some(label);
        }
    }

    new_tree
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

    #[test]
    fn test_bipartition_canonical() {
        let all = vec!["A".to_string(), "B".to_string(), "C".to_string(), "D".to_string()];

        let set1: HashSet<String> = ["A", "B"].iter().map(|s| s.to_string()).collect();
        let set2: HashSet<String> = ["C", "D"].iter().map(|s| s.to_string()).collect();

        let bp1 = Bipartition::new(set1, &all);
        let bp2 = Bipartition::new(set2, &all);

        assert_eq!(bp1, bp2, "complementary splits should be canonical equals");
    }

    #[test]
    fn test_bipartition_trivial() {
        let all = vec!["A".to_string(), "B".to_string(), "C".to_string(), "D".to_string()];
        let single: HashSet<String> = ["A"].iter().map(|s| s.to_string()).collect();
        let bp = Bipartition::new(single, &all);
        assert!(bp.is_trivial());
    }

    #[test]
    fn test_bipartition_nontrivial() {
        let all = vec!["A".to_string(), "B".to_string(), "C".to_string(), "D".to_string()];
        let pair: HashSet<String> = ["A", "B"].iter().map(|s| s.to_string()).collect();
        let bp = Bipartition::new(pair, &all);
        assert!(!bp.is_trivial());
    }

    #[test]
    fn test_tree_bipartitions_4taxa() {
        // ((A,B),(C,D)); has one non-trivial bipartition: {A,B} | {C,D}
        let tree = parse("((A,B),(C,D));");
        let bps = tree_bipartitions(&tree);
        assert_eq!(bps.len(), 1, "4-taxon symmetric tree has 1 non-trivial split");
        assert_eq!(bps[0].taxa, vec!["A".to_string(), "B".to_string()]);
    }

    #[test]
    fn test_tree_bipartitions_5taxa() {
        // (((A,B),C),(D,E)); has 2 non-trivial bipartitions:
        // {A,B} | {C,D,E} and {A,B,C} | {D,E}
        let tree = parse("(((A,B),C),(D,E));");
        let bps = tree_bipartitions(&tree);
        assert_eq!(bps.len(), 2, "5-taxon tree should have 2 non-trivial splits");
    }

    #[test]
    fn test_tree_bipartitions_3taxa() {
        // (A,B,C); — only 3 taxa, no non-trivial bipartitions
        let tree = parse("(A,B,C);");
        let bps = tree_bipartitions(&tree);
        assert_eq!(bps.len(), 0);
    }

    #[test]
    fn test_compute_support_identical() {
        // All replicates have the same topology as reference -> 100% support
        let ref_tree = parse("((A,B),(C,D));");
        let replicates: Vec<Tree> = (0..100).map(|_| parse("((A,B),(C,D));")).collect();

        let support = compute_support(&ref_tree, &replicates);
        let values: Vec<f64> = support.all_supports().iter().map(|(_, v)| *v).collect();

        assert_eq!(values.len(), 1);
        assert!((values[0] - 1.0).abs() < 1e-10, "identical replicates should give 100% support");
    }

    #[test]
    fn test_compute_support_zero() {
        // All replicates have different topology -> 0% support for ref split
        let ref_tree = parse("((A,B),(C,D));");
        // {A,C} | {B,D} is a different split
        let replicates: Vec<Tree> = (0..100).map(|_| parse("((A,C),(B,D));")).collect();

        let support = compute_support(&ref_tree, &replicates);
        let values: Vec<f64> = support.all_supports().iter().map(|(_, v)| *v).collect();

        assert_eq!(values.len(), 1);
        assert!(values[0].abs() < 1e-10, "conflicting replicates should give 0% support");
    }

    #[test]
    fn test_compute_support_mixed() {
        let ref_tree = parse("((A,B),(C,D));");

        let mut replicates = Vec::new();
        // 70 replicates with same topology
        for _ in 0..70 {
            replicates.push(parse("((A,B),(C,D));"));
        }
        // 30 replicates with different topology
        for _ in 0..30 {
            replicates.push(parse("((A,C),(B,D));"));
        }

        let support = compute_support(&ref_tree, &replicates);
        let values: Vec<f64> = support.all_supports().iter().map(|(_, v)| *v).collect();

        assert_eq!(values.len(), 1);
        assert!((values[0] - 0.70).abs() < 1e-10, "expected 70% support, got {}", values[0]);
    }

    #[test]
    fn test_compute_support_5taxa() {
        let ref_tree = parse("(((A,B),C),(D,E));");
        // Bipartitions: {A,B}|{C,D,E} and {D,E}|{A,B,C}

        let mut replicates = Vec::new();
        // 80 replicates keep {A,B} but scramble {D,E}
        for _ in 0..80 {
            replicates.push(parse("(((A,B),D),(C,E));"));
        }
        // 20 replicates keep both
        for _ in 0..20 {
            replicates.push(parse("(((A,B),C),(D,E));"));
        }

        let support = compute_support(&ref_tree, &replicates);

        // {A,B} should appear in all 100 replicates
        let all = vec!["A".to_string(), "B".to_string(), "C".to_string(), "D".to_string(), "E".to_string()];
        let ab_set: HashSet<String> = ["A", "B"].iter().map(|s| s.to_string()).collect();
        let ab_bp = Bipartition::new(ab_set, &all);
        let ab_support = support.get(&ab_bp).unwrap();
        assert!((ab_support - 1.0).abs() < 1e-10, "{{A,B}} should have 100% support");

        // {D,E} should appear in only 20 replicates
        let de_set: HashSet<String> = ["D", "E"].iter().map(|s| s.to_string()).collect();
        let de_bp = Bipartition::new(de_set, &all);
        let de_support = support.get(&de_bp).unwrap();
        assert!((de_support - 0.20).abs() < 1e-10, "{{D,E}} should have 20% support, got {}", de_support);
    }

    #[test]
    fn test_label_tree_with_support() {
        let ref_tree = parse("((A:1,B:1):0.5,(C:1,D:1):0.5);");
        let replicates: Vec<Tree> = (0..100).map(|_| parse("((A,B),(C,D));")).collect();

        let support = compute_support(&ref_tree, &replicates);
        let labeled = label_tree_with_support(&ref_tree, &support);

        // The internal node grouping {A,B} should be labeled "100"
        let mut found_label = false;
        for node in &labeled.nodes {
            if !node.is_leaf() && node.id != labeled.root {
                if let Some(ref name) = node.name {
                    if name == "100" {
                        found_label = true;
                    }
                }
            }
        }
        assert!(found_label, "should find a '100' bootstrap label on an internal node");
    }

    #[test]
    fn test_support_with_branch_lengths() {
        // Branch lengths shouldn't affect bipartition matching
        let ref_tree = parse("((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.05);");
        let replicates: Vec<Tree> = (0..50).map(|_| parse("((A:0.5,B:0.5):0.1,(C:0.5,D:0.5):0.1);")).collect();

        let support = compute_support(&ref_tree, &replicates);
        let values: Vec<f64> = support.all_supports().iter().map(|(_, v)| *v).collect();
        assert_eq!(values.len(), 1);
        assert!((values[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_empty_replicates() {
        let ref_tree = parse("((A,B),(C,D));");
        let replicates: Vec<Tree> = Vec::new();

        let support = compute_support(&ref_tree, &replicates);
        // With 0 replicates, support is 0/0 = NaN, but we should handle gracefully
        assert_eq!(support.num_bipartitions(), 1);
    }
}

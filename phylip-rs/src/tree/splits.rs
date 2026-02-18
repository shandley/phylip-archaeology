//! Bipartition (split) representation and extraction from phylogenetic trees.
//!
//! A split divides the taxon set into two disjoint subsets. Each internal
//! branch of an unrooted tree induces a bipartition of the leaf set.
//! Splits are represented as bit vectors for efficient comparison.
//!
//! # Canonical form
//!
//! Each split has two equivalent representations (the two sides of the
//! bipartition). We canonicalize by ensuring taxon 0 is always on the
//! "0" side of the bit vector.
//!
//! # References
//!
//! - Robinson, D.F. & Foulds, L.R. (1981). Comparison of phylogenetic trees.
//!   *Mathematical Biosciences*, 53:131-147.
//! - Penny, D. & Hendy, M.D. (1985). The use of tree comparison metrics.
//!   *Systematic Zoology*, 34:75-82.

use std::fmt;

use super::Tree;

/// A bipartition (split) of the taxon set, stored as a bit vector.
///
/// Each bit position corresponds to a taxon (by index). A set bit means
/// the taxon is on one side of the split. The canonical form always stores
/// the side that does NOT contain taxon 0 (index 0 is always 0 in the
/// canonical form), unless the split is trivial.
#[derive(Clone)]
pub struct Split {
    bits: Vec<u64>,
    ntaxa: usize,
}

impl Split {
    /// Create a new empty split for `ntaxa` taxa.
    pub fn new(ntaxa: usize) -> Self {
        let nwords = (ntaxa + 63) / 64;
        Self {
            bits: vec![0u64; nwords],
            ntaxa,
        }
    }

    /// Set taxon `i` as being on the "1" side of the split.
    pub fn set(&mut self, i: usize) {
        self.bits[i / 64] |= 1u64 << (i % 64);
    }

    /// Check if taxon `i` is on the "1" side.
    pub fn get(&self, i: usize) -> bool {
        (self.bits[i / 64] >> (i % 64)) & 1 == 1
    }

    /// Count the number of taxa on the "1" side.
    pub fn count(&self) -> usize {
        self.bits.iter().map(|w| w.count_ones() as usize).sum()
    }

    /// Canonicalize: ensure taxon 0 is always on the "0" side.
    /// This means complementary representations map to the same canonical form.
    pub fn canonicalize(&mut self) {
        if self.get(0) {
            self.complement();
        }
    }

    /// Flip all bits (complement the split).
    pub fn complement(&mut self) {
        for w in &mut self.bits {
            *w = !*w;
        }
        // Clear bits beyond ntaxa
        let remainder = self.ntaxa % 64;
        if remainder > 0 {
            if let Some(last) = self.bits.last_mut() {
                *last &= (1u64 << remainder) - 1;
            }
        }
    }

    /// Check if this split is trivial (0 or 1 taxa on either side).
    pub fn is_trivial(&self) -> bool {
        let c = self.count();
        c <= 1 || c >= self.ntaxa - 1
    }

    /// Check if two splits are compatible.
    ///
    /// Two splits A|A' and B|B' are compatible if one of the four
    /// intersections A n B, A n B', A' n B, A' n B' is empty.
    pub fn is_compatible(&self, other: &Split) -> bool {
        // A n B
        let ab: usize = self
            .bits
            .iter()
            .zip(&other.bits)
            .map(|(a, b)| (a & b).count_ones() as usize)
            .sum();

        // A n B' (bits in self but not other)
        let ab_prime: usize = self
            .bits
            .iter()
            .zip(&other.bits)
            .map(|(a, b)| (a & !b).count_ones() as usize)
            .sum();

        // A' n B (bits in other but not self)
        let a_prime_b: usize = self
            .bits
            .iter()
            .zip(&other.bits)
            .map(|(a, b)| (!a & b).count_ones() as usize)
            .sum();

        // A' n B' = ntaxa - ab - ab_prime - a_prime_b
        let a_prime_b_prime = self.ntaxa - ab - ab_prime - a_prime_b;

        ab == 0 || ab_prime == 0 || a_prime_b == 0 || a_prime_b_prime == 0
    }

    /// Compute a hash value for this split.
    ///
    /// Uses a simple FNV-1a-style hash over the bit vector words.
    pub fn hash_value(&self) -> u64 {
        let mut h: u64 = 0xcbf29ce484222325;
        for &w in &self.bits {
            h ^= w;
            h = h.wrapping_mul(0x100000001b3);
        }
        h
    }

    /// Get the set of taxon indices on the "1" side.
    pub fn taxa_indices(&self) -> Vec<usize> {
        (0..self.ntaxa).filter(|&i| self.get(i)).collect()
    }

    /// Number of taxa this split is defined over.
    pub fn ntaxa(&self) -> usize {
        self.ntaxa
    }
}

impl PartialEq for Split {
    fn eq(&self, other: &Self) -> bool {
        self.ntaxa == other.ntaxa && self.bits == other.bits
    }
}

impl Eq for Split {}

impl std::hash::Hash for Split {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.ntaxa.hash(state);
        self.bits.hash(state);
    }
}

impl fmt::Debug for Split {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let indices: Vec<usize> = self.taxa_indices();
        write!(f, "Split({:?})", indices)
    }
}

/// Build a sorted list of taxon names from a tree.
///
/// Returns all leaf node names in sorted order, providing a canonical
/// ordering for split construction.
pub fn build_taxon_index(tree: &Tree) -> Vec<String> {
    let mut names: Vec<String> = tree
        .leaves()
        .iter()
        .filter_map(|n| n.name.clone())
        .collect();
    names.sort();
    names
}

/// Extract all non-trivial splits from a tree, given a taxon ordering.
///
/// Performs a postorder traversal to compute the descendant leaf set for
/// each node, then extracts the induced bipartition for each internal
/// branch. Splits are returned in canonical form with their associated
/// branch lengths.
///
/// # Arguments
///
/// * `tree` - The phylogenetic tree to extract splits from
/// * `taxon_order` - Sorted taxon names defining the bit positions
///
/// # Returns
///
/// A vector of (split, optional_branch_length) pairs for each non-trivial
/// internal branch.
pub fn extract_splits(tree: &Tree, taxon_order: &[String]) -> Vec<(Split, Option<f64>)> {
    let ntaxa = taxon_order.len();

    // Build name -> index mapping
    let taxon_index: std::collections::HashMap<&str, usize> = taxon_order
        .iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i))
        .collect();

    // Compute descendant taxon sets via postorder traversal
    let mut desc_bits: Vec<Split> = vec![Split::new(ntaxa); tree.num_nodes()];

    for &node_id in &tree.postorder() {
        let node = &tree.nodes[node_id];
        if node.is_leaf() {
            if let Some(ref name) = node.name {
                if let Some(&idx) = taxon_index.get(name.as_str()) {
                    desc_bits[node_id].set(idx);
                }
            }
        } else {
            for &child in &node.children {
                for (w_idx, &w) in desc_bits[child].bits.clone().iter().enumerate() {
                    desc_bits[node_id].bits[w_idx] |= w;
                }
            }
        }
    }

    let mut splits = Vec::new();
    for node_id in 0..tree.num_nodes() {
        let node = &tree.nodes[node_id];
        if node.is_leaf() || node_id == tree.root {
            continue;
        }

        let mut split = desc_bits[node_id].clone();
        split.canonicalize();

        if !split.is_trivial() {
            splits.push((split, node.branch_length));
        }
    }

    splits
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::newick;

    fn parse(s: &str) -> Tree {
        newick::parse_newick(s).unwrap()
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

    #[test]
    fn test_split_hash_value() {
        let mut s1 = Split::new(4);
        s1.set(1);
        s1.set(2);

        let mut s2 = Split::new(4);
        s2.set(1);
        s2.set(2);

        assert_eq!(s1.hash_value(), s2.hash_value());

        let mut s3 = Split::new(4);
        s3.set(0);
        s3.set(3);
        // Different splits should (very likely) have different hashes
        assert_ne!(s1.hash_value(), s3.hash_value());
    }

    #[test]
    fn test_split_eq() {
        let mut s1 = Split::new(4);
        s1.set(1);
        s1.set(2);

        let mut s2 = Split::new(4);
        s2.set(1);
        s2.set(2);

        assert_eq!(s1, s2);

        let mut s3 = Split::new(4);
        s3.set(0);
        s3.set(3);

        assert_ne!(s1, s3);
    }

    #[test]
    fn test_build_taxon_index() {
        let tree = parse("((C,A),(D,B));");
        let index = build_taxon_index(&tree);
        assert_eq!(index, vec!["A", "B", "C", "D"]);
    }

    // --- Extract splits tests ---

    #[test]
    fn test_extract_splits_4taxa() {
        let tree = parse("((A,B),(C,D));");
        let taxa = build_taxon_index(&tree);
        let splits = extract_splits(&tree, &taxa);
        // Deduplicate: may get the same split from both children of root
        let unique: std::collections::HashSet<Split> =
            splits.into_iter().map(|(s, _)| s).collect();
        assert_eq!(unique.len(), 1, "4-taxon symmetric tree has 1 non-trivial split");
    }

    #[test]
    fn test_extract_splits_5taxa() {
        let tree = parse("(((A,B),C),(D,E));");
        let taxa = build_taxon_index(&tree);
        let splits = extract_splits(&tree, &taxa);
        let unique: std::collections::HashSet<Split> =
            splits.into_iter().map(|(s, _)| s).collect();
        assert!(
            unique.len() >= 2,
            "5-taxon tree should have at least 2 non-trivial splits"
        );
    }

    #[test]
    fn test_extract_splits_with_branch_lengths() {
        let tree = parse("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);");
        let taxa = build_taxon_index(&tree);
        let splits = extract_splits(&tree, &taxa);
        // At least one split should have a branch length
        let has_bl = splits.iter().any(|(_, bl)| bl.is_some());
        assert!(has_bl, "splits from tree with branch lengths should carry them");
    }

    #[test]
    fn test_extract_splits_same_topology_canonical() {
        // Two different Newick representations of the same unrooted topology
        let tree1 = parse("((A,B),(C,D));");
        let tree2 = parse("((C,D),(A,B));");
        let taxa1 = build_taxon_index(&tree1);
        let taxa2 = build_taxon_index(&tree2);
        assert_eq!(taxa1, taxa2);

        let splits1: std::collections::HashSet<Split> =
            extract_splits(&tree1, &taxa1).into_iter().map(|(s, _)| s).collect();
        let splits2: std::collections::HashSet<Split> =
            extract_splits(&tree2, &taxa2).into_iter().map(|(s, _)| s).collect();

        assert_eq!(splits1, splits2);
    }
}

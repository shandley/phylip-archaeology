//! Tree distance metrics for comparing phylogenetic tree topologies.
//!
//! Provides three distance measures between phylogenetic trees:
//!
//! - **Robinson-Foulds (RF) distance**: counts the number of splits
//!   (bipartitions) present in one tree but not the other. This is the
//!   most widely used topological distance.
//!
//! - **Normalized Robinson-Foulds distance**: RF divided by the maximum
//!   possible RF distance, giving a value in [0, 1].
//!
//! - **Branch Score Distance (BSD)**: the Euclidean distance in split-space,
//!   where each split is a coordinate axis and branch lengths define the
//!   position. This captures both topological and branch-length differences.
//!
//! # References
//!
//! - Robinson, D.F. & Foulds, L.R. (1981). Comparison of phylogenetic trees.
//!   *Mathematical Biosciences*, 53:131-147.
//! - Kuhner, M.K. & Felsenstein, J. (1994). A simulation comparison of
//!   phylogeny algorithms under equal and unequal evolutionary rates.
//!   *Molecular Biology and Evolution*, 11:459-468.

use std::collections::HashSet;

use super::splits::{build_taxon_index, extract_splits, Split};
use super::Tree;

/// Compute the Robinson-Foulds distance between two trees.
///
/// The RF distance is the size of the symmetric difference of the two
/// non-trivial split sets: RF = |S1 \ S2| + |S2 \ S1|.
///
/// Both trees must have the same set of leaf taxa.
///
/// # Errors
///
/// Returns an error if the trees have different leaf sets.
pub fn robinson_foulds(tree1: &Tree, tree2: &Tree) -> Result<usize, String> {
    let taxa1 = build_taxon_index(tree1);
    let taxa2 = build_taxon_index(tree2);
    if taxa1 != taxa2 {
        return Err("trees have different leaf sets".to_string());
    }

    let splits1: HashSet<Split> = extract_splits(tree1, &taxa1)
        .into_iter()
        .map(|(s, _)| s)
        .collect();
    let splits2: HashSet<Split> = extract_splits(tree2, &taxa2)
        .into_iter()
        .map(|(s, _)| s)
        .collect();

    // Symmetric difference: splits in 1 but not 2, plus splits in 2 but not 1
    let only_in_1 = splits1.difference(&splits2).count();
    let only_in_2 = splits2.difference(&splits1).count();

    Ok(only_in_1 + only_in_2)
}

/// Compute the normalized Robinson-Foulds distance between two trees.
///
/// The normalized RF distance is RF / (2 * (n - 3)) for n taxa, which
/// gives a value in [0, 1]. For n <= 3 there are no non-trivial splits
/// and the distance is always 0.
///
/// # Errors
///
/// Returns an error if the trees have different leaf sets.
pub fn normalized_robinson_foulds(tree1: &Tree, tree2: &Tree) -> Result<f64, String> {
    let rf = robinson_foulds(tree1, tree2)?;
    let n = tree1.num_leaves();
    if n <= 3 {
        return Ok(0.0);
    }
    let max_rf = 2 * (n - 3);
    Ok(rf as f64 / max_rf as f64)
}

/// Compute the Branch Score Distance between two trees.
///
/// The BSD (Kuhner & Felsenstein, 1994) is the Euclidean distance in
/// split-space:
///
/// ```text
/// BSD = sqrt( sum over all splits s in S1 union S2 of (b1(s) - b2(s))^2 )
/// ```
///
/// where b_i(s) is the branch length associated with split s in tree i
/// (0 if the split is absent from tree i).
///
/// Both trees must have the same set of leaf taxa.
///
/// # Errors
///
/// Returns an error if the trees have different leaf sets.
pub fn branch_score_distance(tree1: &Tree, tree2: &Tree) -> Result<f64, String> {
    let taxa1 = build_taxon_index(tree1);
    let taxa2 = build_taxon_index(tree2);
    if taxa1 != taxa2 {
        return Err("trees have different leaf sets".to_string());
    }

    let raw_splits1 = extract_splits(tree1, &taxa1);
    let raw_splits2 = extract_splits(tree2, &taxa2);

    // Build maps: split -> branch_length for each tree
    // If a split appears multiple times (e.g. both children of root),
    // take the first branch length encountered.
    let mut map1: std::collections::HashMap<Split, f64> = std::collections::HashMap::new();
    for (split, bl) in &raw_splits1 {
        map1.entry(split.clone())
            .or_insert_with(|| bl.unwrap_or(0.0));
    }

    let mut map2: std::collections::HashMap<Split, f64> = std::collections::HashMap::new();
    for (split, bl) in &raw_splits2 {
        map2.entry(split.clone())
            .or_insert_with(|| bl.unwrap_or(0.0));
    }

    // Collect all unique splits across both trees
    let all_splits: HashSet<&Split> = map1.keys().chain(map2.keys()).collect();

    let mut sum_sq = 0.0;
    for split in all_splits {
        let b1 = map1.get(split).copied().unwrap_or(0.0);
        let b2 = map2.get(split).copied().unwrap_or(0.0);
        let diff = b1 - b2;
        sum_sq += diff * diff;
    }

    Ok(sum_sq.sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::newick;

    fn parse(s: &str) -> Tree {
        newick::parse_newick(s).unwrap()
    }

    // --- Robinson-Foulds tests ---

    #[test]
    fn test_rf_identical_trees() {
        let t1 = parse("((A,B),(C,D));");
        let t2 = parse("((A,B),(C,D));");
        assert_eq!(robinson_foulds(&t1, &t2).unwrap(), 0);
    }

    #[test]
    fn test_rf_different_topologies_4taxa() {
        // ((A,B),(C,D)) has split {A,B}|{C,D} (canonicalized to {C,D})
        // ((A,C),(B,D)) has split {A,C}|{B,D} (canonicalized to {B,D})
        // These are different splits, so RF = 2 (one in each but not the other)
        let t1 = parse("((A,B),(C,D));");
        let t2 = parse("((A,C),(B,D));");
        assert_eq!(robinson_foulds(&t1, &t2).unwrap(), 2);
    }

    #[test]
    fn test_rf_symmetry() {
        let t1 = parse("((A,B),(C,D));");
        let t2 = parse("((A,C),(B,D));");
        let d12 = robinson_foulds(&t1, &t2).unwrap();
        let d21 = robinson_foulds(&t2, &t1).unwrap();
        assert_eq!(d12, d21);
    }

    #[test]
    fn test_rf_different_leaf_sets_error() {
        let t1 = parse("((A,B),(C,D));");
        let t2 = parse("((A,B),(C,E));");
        assert!(robinson_foulds(&t1, &t2).is_err());
    }

    #[test]
    fn test_rf_identical_5taxa() {
        let t1 = parse("(((A,B),C),(D,E));");
        let t2 = parse("(((A,B),C),(D,E));");
        assert_eq!(robinson_foulds(&t1, &t2).unwrap(), 0);
    }

    #[test]
    fn test_rf_5taxa_one_split_different() {
        // t1: (((A,B),C),(D,E)) has splits {A,B} and {D,E}
        // t2: (((A,C),B),(D,E)) has splits {A,C} and {D,E}
        // {D,E} is shared, {A,B} vs {A,C} differ -> RF = 2
        let t1 = parse("(((A,B),C),(D,E));");
        let t2 = parse("(((A,C),B),(D,E));");
        assert_eq!(robinson_foulds(&t1, &t2).unwrap(), 2);
    }

    #[test]
    fn test_rf_5taxa_all_splits_different() {
        // t1: (((A,B),C),(D,E)) has splits {A,B} and {D,E}
        // t2: (((A,D),E),(B,C)) has splits {B,C} and {A,D}
        // No shared splits -> RF = 4
        let t1 = parse("(((A,B),C),(D,E));");
        let t2 = parse("(((A,D),E),(B,C));");
        assert_eq!(robinson_foulds(&t1, &t2).unwrap(), 4);
    }

    #[test]
    fn test_rf_3taxa_always_zero() {
        // 3-taxon trees have no non-trivial splits
        let t1 = parse("((A,B),C);");
        let t2 = parse("((A,C),B);");
        assert_eq!(robinson_foulds(&t1, &t2).unwrap(), 0);
    }

    #[test]
    fn test_rf_with_branch_lengths_topology_only() {
        // RF should not depend on branch lengths, only topology
        let t1 = parse("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);");
        let t2 = parse("((A:0.5,B:0.6):0.7,(C:0.8,D:0.9):1.0);");
        assert_eq!(robinson_foulds(&t1, &t2).unwrap(), 0);
    }

    #[test]
    fn test_rf_same_topology_doubled_lengths() {
        let t1 = parse("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);");
        let t2 = parse("((A:0.2,B:0.4):0.6,(C:0.8,D:1.0):1.2);");
        assert_eq!(robinson_foulds(&t1, &t2).unwrap(), 0);
    }

    // --- Normalized RF tests ---

    #[test]
    fn test_nrf_identical() {
        let t1 = parse("((A,B),(C,D));");
        let t2 = parse("((A,B),(C,D));");
        let nrf = normalized_robinson_foulds(&t1, &t2).unwrap();
        assert!(nrf.abs() < 1e-10);
    }

    #[test]
    fn test_nrf_maximally_different_4taxa() {
        // 4 taxa: max RF = 2*(4-3) = 2, RF = 2 -> nRF = 1.0
        let t1 = parse("((A,B),(C,D));");
        let t2 = parse("((A,C),(B,D));");
        let nrf = normalized_robinson_foulds(&t1, &t2).unwrap();
        assert!((nrf - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_nrf_3taxa_zero() {
        let t1 = parse("((A,B),C);");
        let t2 = parse("((A,C),B);");
        let nrf = normalized_robinson_foulds(&t1, &t2).unwrap();
        assert!(nrf.abs() < 1e-10);
    }

    #[test]
    fn test_nrf_range() {
        let t1 = parse("(((A,B),C),(D,E));");
        let t2 = parse("(((A,C),B),(D,E));");
        let nrf = normalized_robinson_foulds(&t1, &t2).unwrap();
        assert!(nrf >= 0.0 && nrf <= 1.0);
    }

    // --- Branch Score Distance tests ---

    #[test]
    fn test_bsd_identical_trees() {
        let t1 = parse("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);");
        let t2 = parse("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);");
        let bsd = branch_score_distance(&t1, &t2).unwrap();
        assert!(bsd.abs() < 1e-10);
    }

    #[test]
    fn test_bsd_same_topology_different_lengths() {
        // Same topology, different branch lengths -> RF=0 but BSD>0
        let t1 = parse("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);");
        let t2 = parse("((A:0.2,B:0.4):0.6,(C:0.8,D:1.0):1.2);");
        let rf = robinson_foulds(&t1, &t2).unwrap();
        let bsd = branch_score_distance(&t1, &t2).unwrap();
        assert_eq!(rf, 0);
        assert!(bsd > 0.0, "BSD should be positive for different branch lengths");
    }

    #[test]
    fn test_bsd_symmetry() {
        let t1 = parse("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);");
        let t2 = parse("((A:0.2,B:0.4):0.6,(C:0.8,D:1.0):1.2);");
        let d12 = branch_score_distance(&t1, &t2).unwrap();
        let d21 = branch_score_distance(&t2, &t1).unwrap();
        assert!((d12 - d21).abs() < 1e-10);
    }

    #[test]
    fn test_bsd_different_topologies() {
        // Different topologies with branch lengths
        let t1 = parse("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);");
        let t2 = parse("((A:0.1,C:0.2):0.3,(B:0.4,D:0.5):0.6);");
        let bsd = branch_score_distance(&t1, &t2).unwrap();
        assert!(bsd > 0.0);
    }

    #[test]
    fn test_bsd_no_branch_lengths() {
        // Trees without branch lengths should use 0.0 for all splits
        let t1 = parse("((A,B),(C,D));");
        let t2 = parse("((A,B),(C,D));");
        let bsd = branch_score_distance(&t1, &t2).unwrap();
        assert!(bsd.abs() < 1e-10);
    }

    #[test]
    fn test_bsd_5taxa_known_value() {
        // Both trees have the same 2 non-trivial splits: {A,B} and {D,E}
        // We verify by extracting splits, then check BSD is as expected.
        let t1 = parse("(((A:0.1,B:0.1):0.5,C:0.2),(D:0.1,E:0.1):0.3);");
        let t2 = parse("(((A:0.1,B:0.1):0.7,C:0.2),(D:0.1,E:0.1):0.1);");

        // Verify same topology
        assert_eq!(robinson_foulds(&t1, &t2).unwrap(), 0);

        // BSD should be positive due to different internal branch lengths
        let bsd = branch_score_distance(&t1, &t2).unwrap();
        assert!(bsd > 0.0, "BSD should be positive for different branch lengths");

        // Symmetry check
        let bsd_rev = branch_score_distance(&t2, &t1).unwrap();
        assert!(
            (bsd - bsd_rev).abs() < 1e-10,
            "BSD should be symmetric: {} vs {}",
            bsd,
            bsd_rev
        );
    }

    #[test]
    fn test_bsd_4taxa_known_value() {
        // 4 taxa: one non-trivial split {C,D} with branch length
        // t1: ((A,B),(C,D)):  split {C,D} from (C,D) internal node with bl=0.3
        //     and split {A,B} from (A,B) internal node with bl=0.3
        //     But {A,B} = complement of {C,D}, so both are the same split.
        //     Only one unique split, with length 0.3 (first one found).
        // t2: same split with length 0.5
        // BSD = |0.3 - 0.5| = 0.2
        let t1 = parse("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);");
        let t2 = parse("((A:0.1,B:0.1):0.5,(C:0.1,D:0.1):0.5);");
        let bsd = branch_score_distance(&t1, &t2).unwrap();
        assert_eq!(robinson_foulds(&t1, &t2).unwrap(), 0);
        assert!(
            (bsd - 0.2).abs() < 1e-10,
            "expected BSD=0.2, got {}",
            bsd
        );
    }

    #[test]
    fn test_bsd_different_leaf_sets_error() {
        let t1 = parse("((A,B),(C,D));");
        let t2 = parse("((A,B),(C,E));");
        assert!(branch_score_distance(&t1, &t2).is_err());
    }
}

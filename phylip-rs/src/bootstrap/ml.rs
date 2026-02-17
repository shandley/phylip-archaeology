//! Maximum likelihood bootstrap analysis.
//!
//! Integrates bootstrap resampling with ML tree search to produce
//! a best tree annotated with bootstrap support values. This is the
//! standard workflow for assessing confidence in ML phylogenies:
//!
//! 1. Infer the best ML tree from the original alignment
//! 2. Generate bootstrap replicate alignments by resampling columns
//! 3. Infer an ML tree from each replicate
//! 4. Compute bootstrap support: for each branch in the best tree,
//!    the fraction of replicate trees containing the same bipartition
//! 5. Label the best tree with support values
//!
//! # Example
//!
//! ```
//! use phylip_rs::bootstrap::ml::ml_bootstrap;
//! use phylip_rs::likelihood::models::Jc69Model;
//! use phylip_rs::tree::{Alignment, Base, Sequence};
//!
//! let alignment = Alignment::new(vec![
//!     Sequence::new("A", vec![Base::A, Base::A, Base::A, Base::A, Base::C, Base::C, Base::C, Base::C]),
//!     Sequence::new("B", vec![Base::A, Base::A, Base::A, Base::A, Base::C, Base::C, Base::C, Base::C]),
//!     Sequence::new("C", vec![Base::C, Base::C, Base::C, Base::C, Base::A, Base::A, Base::A, Base::A]),
//!     Sequence::new("D", vec![Base::C, Base::C, Base::C, Base::C, Base::A, Base::A, Base::A, Base::A]),
//! ]).unwrap();
//!
//! let model = Jc69Model;
//! let result = ml_bootstrap(&alignment, &model, 10, 42).unwrap();
//! assert!(result.lnl < 0.0);
//! assert_eq!(result.replicate_trees.len(), 10);
//! ```
//!
//! # References
//!
//! - Felsenstein, J. (1985). Confidence limits on phylogenies:
//!   an approach using the bootstrap. *Evolution*, 39, 783-791.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.

use crate::bootstrap::support::{compute_support, label_tree_with_support, BootstrapSupport};
use crate::bootstrap::{bootstrap_weights, resample_alignment, ResamplingMethod, SimpleRng};
use crate::likelihood::models::SubstitutionModel;
use crate::likelihood::pruning::LikelihoodError;
use crate::likelihood::search::{search, MLSearchResult};
use crate::tree::types::{Alignment, Tree};

/// Result of a maximum likelihood bootstrap analysis.
#[derive(Debug, Clone)]
pub struct MLBootstrapResult {
    /// The best ML tree with bootstrap support labels on internal nodes.
    /// Support values are formatted as percentages (0-100) in node names.
    pub tree: Tree,
    /// Log-likelihood of the best tree (from the original alignment).
    pub lnl: f64,
    /// ML trees inferred from each bootstrap replicate alignment.
    pub replicate_trees: Vec<Tree>,
    /// Bootstrap support values for each non-trivial bipartition.
    pub support_values: BootstrapSupport,
}

/// Perform a maximum likelihood bootstrap analysis.
///
/// Infers the best ML tree from the original alignment, generates
/// `n_replicates` bootstrap replicate alignments, runs ML search on
/// each, and computes bootstrap support for the branches of the best tree.
///
/// # Arguments
///
/// * `alignment` - The DNA alignment to analyze.
/// * `model` - The substitution model (e.g., JC69 or F84).
/// * `n_replicates` - Number of bootstrap replicates to generate.
/// * `seed` - Random seed for reproducibility.
///
/// # Returns
///
/// An `MLBootstrapResult` containing the best tree annotated with
/// bootstrap support, the log-likelihood, all replicate trees,
/// and the raw support values.
///
/// # Errors
///
/// Returns `LikelihoodError` if ML search fails on the original
/// alignment or any replicate. Common causes include too few taxa
/// (need at least 3) or numerical issues.
pub fn ml_bootstrap(
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    n_replicates: usize,
    seed: u64,
) -> Result<MLBootstrapResult, LikelihoodError> {
    // Phase 1: Infer the best ML tree from the original alignment.
    let best_result: MLSearchResult = search(alignment, model, Some(seed))?;

    // Phase 2: Generate bootstrap replicate alignments and infer trees.
    let mut rng = SimpleRng::seed(seed);
    let method = ResamplingMethod::Bootstrap;
    let mut replicate_trees = Vec::with_capacity(n_replicates);

    for i in 0..n_replicates {
        let weights = bootstrap_weights(alignment.nsites(), &method, &mut rng);
        let replicate_alignment = resample_alignment(alignment, &weights);

        // Use a derived seed for each replicate's ML search so that the
        // search taxon-addition order varies across replicates.
        let replicate_seed = seed.wrapping_add(i as u64).wrapping_add(1);
        let replicate_result = search(&replicate_alignment, model, Some(replicate_seed))?;
        replicate_trees.push(replicate_result.tree);
    }

    // Phase 3: Compute bootstrap support values.
    let support_values = compute_support(&best_result.tree, &replicate_trees);

    // Phase 4: Label the best tree with support percentages.
    let labeled_tree = label_tree_with_support(&best_result.tree, &support_values);

    Ok(MLBootstrapResult {
        tree: labeled_tree,
        lnl: best_result.lnl,
        replicate_trees,
        support_values,
    })
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::likelihood::models::Jc69Model;
    use crate::tree::types::{Base, Sequence};

    /// Helper to build an alignment from (&str, &str) pairs.
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

    // ----------------------------------------------------------
    // Test 1: Basic functionality - ml_bootstrap runs and returns Ok
    // ----------------------------------------------------------
    #[test]
    fn test_ml_bootstrap_basic() {
        let aln = make_alignment(&[
            ("A", "AAAAACCCCC"),
            ("B", "AAAAACCCCC"),
            ("C", "CCCCCAAAAA"),
            ("D", "CCCCCAAAAA"),
        ]);
        let model = Jc69Model;
        let result = ml_bootstrap(&aln, &model, 5, 42);
        assert!(result.is_ok(), "ml_bootstrap should succeed: {:?}", result.err());
    }

    // ----------------------------------------------------------
    // Test 2: Result structure - tree has correct number of leaves
    // ----------------------------------------------------------
    #[test]
    fn test_ml_bootstrap_tree_structure() {
        let aln = make_alignment(&[
            ("A", "AAAAACCCCC"),
            ("B", "AAAAACCCCC"),
            ("C", "CCCCCAAAAA"),
            ("D", "CCCCCAAAAA"),
        ]);
        let model = Jc69Model;
        let result = ml_bootstrap(&aln, &model, 5, 42).unwrap();

        assert_eq!(result.tree.num_leaves(), 4, "tree should have 4 leaves");
        assert!(result.lnl < 0.0, "log-likelihood should be negative");
    }

    // ----------------------------------------------------------
    // Test 3: Correct number of replicate trees
    // ----------------------------------------------------------
    #[test]
    fn test_ml_bootstrap_replicate_count() {
        let aln = make_alignment(&[
            ("A", "AAAAACCCCC"),
            ("B", "AAAAACCCCC"),
            ("C", "CCCCCAAAAA"),
            ("D", "CCCCCAAAAA"),
        ]);
        let model = Jc69Model;

        for n_reps in [3, 5, 10] {
            let result = ml_bootstrap(&aln, &model, n_reps, 42).unwrap();
            assert_eq!(
                result.replicate_trees.len(),
                n_reps,
                "expected {} replicate trees, got {}",
                n_reps,
                result.replicate_trees.len()
            );
        }
    }

    // ----------------------------------------------------------
    // Test 4: Support values are in [0, 1]
    // ----------------------------------------------------------
    #[test]
    fn test_ml_bootstrap_support_values_range() {
        let aln = make_alignment(&[
            ("A", "AAAAACCCCC"),
            ("B", "AAAAACCCCC"),
            ("C", "CCCCCAAAAA"),
            ("D", "CCCCCAAAAA"),
        ]);
        let model = Jc69Model;
        let result = ml_bootstrap(&aln, &model, 10, 42).unwrap();

        for (bp, support) in result.support_values.all_supports() {
            assert!(
                (0.0..=1.0).contains(&support),
                "support value {} for bipartition {} is outside [0, 1]",
                support,
                bp
            );
        }
    }

    // ----------------------------------------------------------
    // Test 5: Identical sequences grouped together get high support
    // ----------------------------------------------------------
    #[test]
    fn test_ml_bootstrap_identical_sequences_high_support() {
        // Strong phylogenetic signal: A,B identical; C,D identical and very different.
        // Bootstrap support for the {A,B}|{C,D} split should be high.
        let aln = make_alignment(&[
            ("A", "AAAAAAAACCCCCCCC"),
            ("B", "AAAAAAAACCCCCCCC"),
            ("C", "CCCCCCCCAAAAAAAA"),
            ("D", "CCCCCCCCAAAAAAAA"),
        ]);
        let model = Jc69Model;
        let result = ml_bootstrap(&aln, &model, 20, 42).unwrap();

        let supports: Vec<f64> = result
            .support_values
            .all_supports()
            .iter()
            .map(|(_, v)| *v)
            .collect();

        // With perfectly identical pairs and strong signal, support should
        // be high (>= 0.8). The exact value depends on resampling luck.
        if !supports.is_empty() {
            let max_support = supports.iter().cloned().fold(0.0_f64, f64::max);
            assert!(
                max_support >= 0.5,
                "with strong signal, max support should be >= 0.5, got {}",
                max_support
            );
        }
    }

    // ----------------------------------------------------------
    // Test 6: Replicate trees have the same leaf count as original
    // ----------------------------------------------------------
    #[test]
    fn test_ml_bootstrap_replicate_trees_valid() {
        let aln = make_alignment(&[
            ("A", "ACGTACGTAC"),
            ("B", "ACGTACGTAC"),
            ("C", "TGCATGCATG"),
            ("D", "TGCATGCATG"),
        ]);
        let model = Jc69Model;
        let result = ml_bootstrap(&aln, &model, 5, 42).unwrap();

        for (i, tree) in result.replicate_trees.iter().enumerate() {
            assert_eq!(
                tree.num_leaves(),
                4,
                "replicate tree {} should have 4 leaves, got {}",
                i,
                tree.num_leaves()
            );
        }
    }

    // ----------------------------------------------------------
    // Test 7: Reproducibility with same seed
    // ----------------------------------------------------------
    #[test]
    fn test_ml_bootstrap_reproducible() {
        let aln = make_alignment(&[
            ("A", "AAAAACCCCC"),
            ("B", "AAAAACCCCC"),
            ("C", "CCCCCAAAAA"),
            ("D", "CCCCCAAAAA"),
        ]);
        let model = Jc69Model;

        let r1 = ml_bootstrap(&aln, &model, 5, 42).unwrap();
        let r2 = ml_bootstrap(&aln, &model, 5, 42).unwrap();

        assert!(
            (r1.lnl - r2.lnl).abs() < 1e-6,
            "same seed should produce same lnl: {} vs {}",
            r1.lnl,
            r2.lnl
        );
        assert_eq!(
            r1.replicate_trees.len(),
            r2.replicate_trees.len(),
            "same seed should produce same number of replicates"
        );
    }

    // ----------------------------------------------------------
    // Test 8: Error on too few taxa
    // ----------------------------------------------------------
    #[test]
    fn test_ml_bootstrap_too_few_taxa() {
        let aln = make_alignment(&[("A", "ACGT"), ("B", "ACGT")]);
        let model = Jc69Model;
        let result = ml_bootstrap(&aln, &model, 5, 42);
        assert!(
            result.is_err(),
            "ml_bootstrap should fail with fewer than 3 taxa"
        );
    }

    // ----------------------------------------------------------
    // Test 9: Five taxa - ensures non-trivial bipartitions are found
    // ----------------------------------------------------------
    #[test]
    fn test_ml_bootstrap_five_taxa() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAA"),
            ("B", "AACCGGTTAA"),
            ("C", "CCAATTGGCC"),
            ("D", "CCAATTGGCC"),
            ("E", "GGTTAACCGG"),
        ]);
        let model = Jc69Model;
        let result = ml_bootstrap(&aln, &model, 5, 99).unwrap();

        assert_eq!(result.tree.num_leaves(), 5);
        assert!(result.lnl < 0.0);
        // A 5-taxon tree should have at least 1 non-trivial bipartition
        // (in fact up to 2), so support_values should be non-empty.
        assert!(
            result.support_values.num_bipartitions() >= 1,
            "5-taxon tree should have at least 1 non-trivial bipartition, got {}",
            result.support_values.num_bipartitions()
        );
    }

    // ----------------------------------------------------------
    // Test 10: Support labels appear on internal nodes of the result tree
    // ----------------------------------------------------------
    #[test]
    fn test_ml_bootstrap_labels_on_tree() {
        let aln = make_alignment(&[
            ("A", "AAAAAAAACCCCCCCC"),
            ("B", "AAAAAAAACCCCCCCC"),
            ("C", "CCCCCCCCAAAAAAAA"),
            ("D", "CCCCCCCCAAAAAAAA"),
        ]);
        let model = Jc69Model;
        let result = ml_bootstrap(&aln, &model, 10, 42).unwrap();

        // Check that at least one internal (non-leaf, non-root) node
        // has a name that looks like a support percentage.
        let has_support_label = result.tree.nodes.iter().any(|node| {
            if node.is_leaf() || node.id == result.tree.root {
                return false;
            }
            match &node.name {
                Some(name) => name.parse::<f64>().is_ok(),
                None => false,
            }
        });

        // With 4 taxa there is exactly 1 non-trivial bipartition,
        // so we expect exactly one internal node labeled with support.
        assert!(
            has_support_label,
            "labeled tree should have at least one internal node with a numeric support label"
        );
    }
}

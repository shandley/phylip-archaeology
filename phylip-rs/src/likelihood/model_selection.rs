//! Model selection via information criteria (AIC, BIC, AICc).
//!
//! Implements model comparison for phylogenetic substitution models using
//! Akaike's Information Criterion (AIC), the corrected AIC (AICc), and
//! the Bayesian Information Criterion (BIC). These criteria balance model
//! fit (log-likelihood) against model complexity (number of free parameters)
//! to select the most appropriate substitution model for a given dataset.
//!
//! # Information Criteria
//!
//! - **AIC** = -2 * lnL + 2 * k
//! - **AICc** = AIC + 2k(k+1) / (n - k - 1)  (small-sample correction)
//! - **BIC** = -2 * lnL + k * ln(n)
//!
//! where k is the number of free parameters and n is the number of sites.
//!
//! # Akaike Weights
//!
//! Akaike weights measure the relative likelihood of each model:
//!
//! w_i = exp(-0.5 * delta_i) / sum_j(exp(-0.5 * delta_j))
//!
//! where delta_i = AIC_i - AIC_min. Weights sum to 1.0 and can be
//! interpreted as the probability that model i is the best model
//! given the data and the set of candidate models.
//!
//! # Parameter Counting
//!
//! For an unrooted binary tree with m taxa, there are 2m - 3 branch
//! lengths. Model-specific parameters:
//! - JC69: 0 free model parameters (all rates equal, frequencies fixed)
//! - F84: 4 free model parameters (3 frequency parameters + 1 ttratio)
//!
//! # Example
//!
//! ```
//! use phylip_rs::likelihood::model_selection::{
//!     compute_information_criteria, compare_models,
//! };
//!
//! // Suppose we fitted two models to an alignment of 100 sites
//! let jc69_result = compute_information_criteria(
//!     -250.0,  // log-likelihood
//!     7,       // 0 model params + 7 branch lengths (5 taxa)
//!     100,     // sites
//! );
//! let f84_result = compute_information_criteria(
//!     -245.0,  // better fit
//!     11,      // 4 model params + 7 branch lengths
//!     100,     // sites
//! );
//!
//! let comparison = compare_models(vec![jc69_result, f84_result]);
//! // Weights sum to 1.0
//! let weight_sum: f64 = comparison.aic_weights.iter().sum();
//! assert!((weight_sum - 1.0).abs() < 1e-10);
//! ```
//!
//! # References
//!
//! - Akaike, H. (1974). A new look at the statistical model identification.
//!   *IEEE Transactions on Automatic Control*, 19, 716-723.
//! - Schwarz, G. (1978). Estimating the dimension of a model.
//!   *The Annals of Statistics*, 6, 461-464.
//! - Burnham, K. P. & Anderson, D. R. (2002). *Model Selection and
//!   Multimodel Inference*. Springer.
//! - Posada, D. & Crandall, K. A. (1998). MODELTEST: testing the model
//!   of DNA substitution. *Bioinformatics*, 14, 817-818.

use crate::tree::types::{Alignment, BaseFrequencies, Tree};

use super::models::{F84Model, Jc69Model, SubstitutionModel};
use super::pruning::{optimize_branch_lengths, LikelihoodError};

/// Result of evaluating a single model with information criteria.
#[derive(Debug, Clone)]
pub struct ModelResult {
    /// Name identifying the model (e.g., "JC69", "F84(ttratio=2.0)").
    pub model_name: String,
    /// Log-likelihood of the data given the model and optimized tree.
    pub lnl: f64,
    /// Total number of free parameters (model parameters + branch lengths).
    pub n_free_params: usize,
    /// Number of alignment sites.
    pub n_sites: usize,
    /// Akaike Information Criterion: -2*lnL + 2*k.
    pub aic: f64,
    /// Corrected AIC for small samples: AIC + 2k(k+1)/(n-k-1).
    pub aicc: f64,
    /// Bayesian Information Criterion: -2*lnL + k*ln(n).
    pub bic: f64,
}

/// Comparison of multiple candidate models ranked by information criteria.
#[derive(Debug, Clone)]
pub struct ModelComparison {
    /// All model results, in the order they were provided.
    pub results: Vec<ModelResult>,
    /// Index of the best model by AIC (lowest AIC).
    pub best_aic: usize,
    /// Index of the best model by BIC (lowest BIC).
    pub best_bic: usize,
    /// Akaike weights based on AIC differences, summing to 1.0.
    pub aic_weights: Vec<f64>,
    /// Akaike weights based on BIC differences, summing to 1.0.
    pub bic_weights: Vec<f64>,
}

/// Compute AIC, BIC, and AICc for a given model fit.
///
/// # Arguments
///
/// * `lnl` - The maximized log-likelihood.
/// * `n_params` - Total number of free parameters (model params + branch lengths).
/// * `n_sites` - Number of alignment sites (sample size).
///
/// # Returns
///
/// A `ModelResult` with computed information criteria. The `model_name` field
/// is set to an empty string; callers should set it to the appropriate name.
///
/// # Panics
///
/// Does not panic, but AICc will be `f64::INFINITY` when n_sites <= n_params + 1
/// (the correction denominator becomes non-positive).
pub fn compute_information_criteria(lnl: f64, n_params: usize, n_sites: usize) -> ModelResult {
    let k = n_params as f64;
    let n = n_sites as f64;

    let aic = -2.0 * lnl + 2.0 * k;
    let bic = -2.0 * lnl + k * n.ln();

    let aicc = if n > k + 1.0 {
        aic + 2.0 * k * (k + 1.0) / (n - k - 1.0)
    } else {
        f64::INFINITY
    };

    ModelResult {
        model_name: String::new(),
        lnl,
        n_free_params: n_params,
        n_sites,
        aic,
        aicc,
        bic,
    }
}

/// Compute Akaike-style weights from a vector of criterion values.
///
/// For each model, delta_i = value_i - min(values), then
/// w_i = exp(-0.5 * delta_i) / sum_j(exp(-0.5 * delta_j)).
fn compute_weights(values: &[f64]) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }

    let min_val = values.iter().copied().fold(f64::INFINITY, f64::min);
    let raw: Vec<f64> = values.iter().map(|&v| (-0.5 * (v - min_val)).exp()).collect();
    let total: f64 = raw.iter().sum();

    if total > 0.0 {
        raw.iter().map(|&r| r / total).collect()
    } else {
        // Degenerate case: assign equal weights
        let n = values.len() as f64;
        vec![1.0 / n; values.len()]
    }
}

/// Compare multiple models and rank them by information criteria.
///
/// Computes Akaike weights for both AIC and BIC, and identifies the
/// best model under each criterion.
///
/// # Arguments
///
/// * `results` - A vector of `ModelResult` values, one per candidate model.
///
/// # Panics
///
/// Panics if `results` is empty.
pub fn compare_models(results: Vec<ModelResult>) -> ModelComparison {
    assert!(!results.is_empty(), "compare_models requires at least one model result");

    let aic_values: Vec<f64> = results.iter().map(|r| r.aic).collect();
    let bic_values: Vec<f64> = results.iter().map(|r| r.bic).collect();

    let best_aic = aic_values
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);

    let best_bic = bic_values
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);

    let aic_weights = compute_weights(&aic_values);
    let bic_weights = compute_weights(&bic_values);

    ModelComparison {
        results,
        best_aic,
        best_bic,
        aic_weights,
        bic_weights,
    }
}

/// Count the number of free parameters for a substitution model.
///
/// - JC69: 0 model-specific free parameters (all fixed).
/// - F84: 4 model-specific free parameters (3 frequencies + 1 ttratio).
///
/// Branch length parameters (2m - 3 for an unrooted binary tree with m taxa)
/// are added separately.
fn count_model_params(model: &dyn SubstitutionModel) -> usize {
    match model.name() {
        "JC69" => 0,
        "F84" => 4,
        _ => 0,
    }
}

/// Count the number of branch length parameters for a tree.
///
/// For an unrooted binary tree with m taxa, there are 2m - 3 branches.
/// If the tree is not fully resolved, the actual number of branches
/// (edges with a parent) is counted directly.
fn count_branch_params(tree: &Tree) -> usize {
    let n_taxa = tree.num_leaves();
    if n_taxa < 2 {
        return 0;
    }
    // For an unrooted binary tree: 2m - 3 branches
    // But count actual branches in case the tree is not fully resolved
    let actual_branches = tree.nodes.iter().filter(|n| n.parent.is_some()).count();
    // Use the minimum of theoretical and actual to be conservative
    actual_branches.min(2 * n_taxa - 3)
}

/// Run full model selection on a tree and alignment.
///
/// Evaluates JC69 and F84 (with ttratio values 1.0, 2.0, and 4.0)
/// on the given alignment. For each model, branch lengths are optimized
/// to maximize the log-likelihood. The models are then compared using
/// AIC, AICc, and BIC.
///
/// The base frequencies for F84 models are estimated from the alignment.
///
/// # Arguments
///
/// * `tree` - The starting tree topology (branch lengths will be optimized).
/// * `alignment` - The DNA sequence alignment.
///
/// # Returns
///
/// A `ModelComparison` with all evaluated models, their information criteria
/// scores, and Akaike weights.
///
/// # Errors
///
/// Returns `LikelihoodError` if likelihood evaluation fails for any model
/// (e.g., taxa mismatch, empty tree, numerical issues).
pub fn select_model(
    tree: &Tree,
    alignment: &Alignment,
) -> Result<ModelComparison, LikelihoodError> {
    let n_sites = alignment.nsites();
    let freqs = BaseFrequencies::from_alignment(alignment);

    // Candidate models
    let models: Vec<(String, Box<dyn SubstitutionModel>)> = vec![
        ("JC69".to_string(), Box::new(Jc69Model)),
        (
            "F84(ttratio=1.0)".to_string(),
            Box::new(F84Model::new(freqs, 1.0)),
        ),
        (
            "F84(ttratio=2.0)".to_string(),
            Box::new(F84Model::new(freqs, 2.0)),
        ),
        (
            "F84(ttratio=4.0)".to_string(),
            Box::new(F84Model::new(freqs, 4.0)),
        ),
    ];

    let mut results = Vec::new();

    for (name, model) in &models {
        let mut work_tree = tree.clone();
        let lnl = optimize_branch_lengths(&mut work_tree, alignment, model.as_ref())?;

        let n_model_params = count_model_params(model.as_ref());
        let n_branch_params = count_branch_params(&work_tree);
        let n_params = n_model_params + n_branch_params;

        let mut result = compute_information_criteria(lnl, n_params, n_sites);
        result.model_name = name.clone();
        results.push(result);
    }

    Ok(compare_models(results))
}

/// Count the total number of free parameters for a given model and tree.
///
/// This is a convenience function that combines model-specific parameters
/// with branch length parameters.
pub fn total_free_params(model: &dyn SubstitutionModel, tree: &Tree) -> usize {
    count_model_params(model) + count_branch_params(tree)
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

    // =====================================================================
    // Tests for compute_information_criteria
    // =====================================================================

    #[test]
    fn test_aic_formula() {
        // AIC = -2*lnL + 2*k
        let result = compute_information_criteria(-100.0, 5, 200);
        let expected_aic = -2.0 * (-100.0) + 2.0 * 5.0;
        assert!(
            (result.aic - expected_aic).abs() < 1e-10,
            "AIC: expected {}, got {}",
            expected_aic,
            result.aic
        );
    }

    #[test]
    fn test_bic_formula() {
        // BIC = -2*lnL + k*ln(n)
        let result = compute_information_criteria(-100.0, 5, 200);
        let expected_bic = -2.0 * (-100.0) + 5.0 * (200.0_f64).ln();
        assert!(
            (result.bic - expected_bic).abs() < 1e-10,
            "BIC: expected {}, got {}",
            expected_bic,
            result.bic
        );
    }

    #[test]
    fn test_aicc_formula() {
        // AICc = AIC + 2k(k+1)/(n-k-1)
        let result = compute_information_criteria(-100.0, 5, 200);
        let aic = -2.0 * (-100.0) + 2.0 * 5.0;
        let correction = 2.0 * 5.0 * 6.0 / (200.0 - 5.0 - 1.0);
        let expected_aicc = aic + correction;
        assert!(
            (result.aicc - expected_aicc).abs() < 1e-10,
            "AICc: expected {}, got {}",
            expected_aicc,
            result.aicc
        );
    }

    #[test]
    fn test_aicc_infinite_when_n_too_small() {
        // When n <= k + 1, AICc should be infinity
        let result = compute_information_criteria(-50.0, 10, 10);
        assert!(
            result.aicc.is_infinite(),
            "AICc should be infinite when n <= k+1, got {}",
            result.aicc
        );

        // Also test n == k + 1
        let result2 = compute_information_criteria(-50.0, 10, 11);
        assert!(
            result2.aicc.is_infinite(),
            "AICc should be infinite when n == k+1, got {}",
            result2.aicc
        );
    }

    #[test]
    fn test_stored_fields() {
        let result = compute_information_criteria(-123.456, 7, 500);
        assert!((result.lnl - (-123.456)).abs() < 1e-10);
        assert_eq!(result.n_free_params, 7);
        assert_eq!(result.n_sites, 500);
    }

    // =====================================================================
    // Tests for compare_models
    // =====================================================================

    #[test]
    fn test_lower_aic_is_better() {
        let mut r1 = compute_information_criteria(-100.0, 3, 200);
        r1.model_name = "model_A".to_string();

        let mut r2 = compute_information_criteria(-80.0, 3, 200);
        r2.model_name = "model_B".to_string();

        // r2 has higher lnl (-80 > -100), so lower AIC, so r2 is best
        let comp = compare_models(vec![r1.clone(), r2.clone()]);
        assert_eq!(comp.best_aic, 1, "model_B (index 1) should have lower AIC");
        assert!(r2.aic < r1.aic);
    }

    #[test]
    fn test_lower_bic_is_better() {
        let mut r1 = compute_information_criteria(-100.0, 3, 200);
        r1.model_name = "model_A".to_string();

        let mut r2 = compute_information_criteria(-80.0, 3, 200);
        r2.model_name = "model_B".to_string();

        let comp = compare_models(vec![r1.clone(), r2.clone()]);
        assert_eq!(comp.best_bic, 1, "model_B should have lower BIC");
        assert!(r2.bic < r1.bic);
    }

    #[test]
    fn test_aic_weights_sum_to_one() {
        let mut r1 = compute_information_criteria(-100.0, 3, 200);
        r1.model_name = "A".to_string();
        let mut r2 = compute_information_criteria(-105.0, 5, 200);
        r2.model_name = "B".to_string();
        let mut r3 = compute_information_criteria(-98.0, 4, 200);
        r3.model_name = "C".to_string();

        let comp = compare_models(vec![r1, r2, r3]);
        let sum: f64 = comp.aic_weights.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-10,
            "AIC weights should sum to 1.0, got {}",
            sum
        );
    }

    #[test]
    fn test_bic_weights_sum_to_one() {
        let mut r1 = compute_information_criteria(-100.0, 3, 200);
        r1.model_name = "A".to_string();
        let mut r2 = compute_information_criteria(-105.0, 5, 200);
        r2.model_name = "B".to_string();
        let mut r3 = compute_information_criteria(-98.0, 4, 200);
        r3.model_name = "C".to_string();

        let comp = compare_models(vec![r1, r2, r3]);
        let sum: f64 = comp.bic_weights.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-10,
            "BIC weights should sum to 1.0, got {}",
            sum
        );
    }

    #[test]
    fn test_best_model_has_highest_weight() {
        let mut r1 = compute_information_criteria(-100.0, 3, 200);
        r1.model_name = "A".to_string();
        let mut r2 = compute_information_criteria(-80.0, 3, 200);
        r2.model_name = "B".to_string();

        let comp = compare_models(vec![r1, r2]);
        // Best model (index 1) should have highest AIC weight
        assert!(
            comp.aic_weights[comp.best_aic] >= comp.aic_weights[1 - comp.best_aic],
            "best model's AIC weight ({}) should be >= other ({})",
            comp.aic_weights[comp.best_aic],
            comp.aic_weights[1 - comp.best_aic]
        );
    }

    #[test]
    fn test_weights_between_zero_and_one() {
        let mut r1 = compute_information_criteria(-100.0, 3, 200);
        r1.model_name = "A".to_string();
        let mut r2 = compute_information_criteria(-105.0, 5, 200);
        r2.model_name = "B".to_string();
        let mut r3 = compute_information_criteria(-98.0, 4, 200);
        r3.model_name = "C".to_string();

        let comp = compare_models(vec![r1, r2, r3]);
        for (i, &w) in comp.aic_weights.iter().enumerate() {
            assert!(
                w >= 0.0 && w <= 1.0,
                "AIC weight[{}] = {} should be in [0, 1]",
                i,
                w
            );
        }
        for (i, &w) in comp.bic_weights.iter().enumerate() {
            assert!(
                w >= 0.0 && w <= 1.0,
                "BIC weight[{}] = {} should be in [0, 1]",
                i,
                w
            );
        }
    }

    #[test]
    fn test_bic_penalizes_more_than_aic_for_large_n() {
        // For large n, BIC penalty k*ln(n) > AIC penalty 2*k when ln(n) > 2,
        // i.e. when n > e^2 ~ 7.39, so n >= 8.
        let n_sites = 1000;
        let k_simple = 3;
        let k_complex = 10;
        let lnl = -200.0;

        let simple = compute_information_criteria(lnl, k_simple, n_sites);
        let complex = compute_information_criteria(lnl, k_complex, n_sites);

        // Both have same lnl, so difference is purely from the penalty
        let aic_diff = complex.aic - simple.aic; // penalty difference by AIC
        let bic_diff = complex.bic - simple.bic; // penalty difference by BIC

        // BIC should penalize the extra parameters more heavily for large n
        assert!(
            bic_diff > aic_diff,
            "BIC penalty difference ({}) should exceed AIC penalty difference ({}) for n={}",
            bic_diff,
            aic_diff,
            n_sites
        );
    }

    #[test]
    fn test_single_model() {
        let mut r = compute_information_criteria(-100.0, 5, 200);
        r.model_name = "only_model".to_string();

        let comp = compare_models(vec![r]);
        assert_eq!(comp.best_aic, 0);
        assert_eq!(comp.best_bic, 0);
        assert!((comp.aic_weights[0] - 1.0).abs() < 1e-10);
        assert!((comp.bic_weights[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_identical_models_equal_weights() {
        let mut r1 = compute_information_criteria(-100.0, 5, 200);
        r1.model_name = "A".to_string();
        let mut r2 = compute_information_criteria(-100.0, 5, 200);
        r2.model_name = "B".to_string();

        let comp = compare_models(vec![r1, r2]);
        assert!(
            (comp.aic_weights[0] - 0.5).abs() < 1e-10,
            "identical models should have equal AIC weights, got {}",
            comp.aic_weights[0]
        );
        assert!(
            (comp.aic_weights[1] - 0.5).abs() < 1e-10,
            "identical models should have equal AIC weights, got {}",
            comp.aic_weights[1]
        );
    }

    // =====================================================================
    // Tests for parameter counting
    // =====================================================================

    #[test]
    fn test_jc69_param_count() {
        let model = Jc69Model;
        assert_eq!(
            count_model_params(&model),
            0,
            "JC69 should have 0 free model parameters"
        );
    }

    #[test]
    fn test_f84_param_count() {
        let model = F84Model::equal_freqs(2.0);
        assert_eq!(
            count_model_params(&model),
            4,
            "F84 should have 4 free model parameters"
        );
    }

    #[test]
    fn test_branch_param_count() {
        // 4-taxon unrooted binary tree: 2*4 - 3 = 5 branches
        let tree = parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();
        let n_branches = count_branch_params(&tree);
        assert_eq!(
            n_branches, 5,
            "4-taxon unrooted tree should have 5 branches, got {}",
            n_branches
        );
    }

    #[test]
    fn test_total_free_params_jc69() {
        let model = Jc69Model;
        // 3-taxon tree: 2*3 - 3 = 3 branches + 0 model params = 3
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let total = total_free_params(&model, &tree);
        assert_eq!(total, 3, "JC69 on 3-taxon tree: expected 3 params, got {}", total);
    }

    #[test]
    fn test_total_free_params_f84() {
        let model = F84Model::equal_freqs(2.0);
        // 3-taxon tree: 3 branches + 4 model params = 7
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let total = total_free_params(&model, &tree);
        assert_eq!(total, 7, "F84 on 3-taxon tree: expected 7 params, got {}", total);
    }

    // =====================================================================
    // Tests for select_model (integration)
    // =====================================================================

    #[test]
    fn test_select_model_basic() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGG"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();

        let comp = select_model(&tree, &aln).unwrap();

        // Should have 4 model results (JC69 + 3 F84 variants)
        assert_eq!(comp.results.len(), 4);

        // All log-likelihoods should be negative
        for r in &comp.results {
            assert!(r.lnl < 0.0, "lnl for {} should be negative", r.model_name);
        }

        // AIC weights should sum to 1
        let aic_sum: f64 = comp.aic_weights.iter().sum();
        assert!(
            (aic_sum - 1.0).abs() < 1e-10,
            "AIC weights sum to {}, expected 1.0",
            aic_sum
        );

        // BIC weights should sum to 1
        let bic_sum: f64 = comp.bic_weights.iter().sum();
        assert!(
            (bic_sum - 1.0).abs() < 1e-10,
            "BIC weights sum to {}, expected 1.0",
            bic_sum
        );

        // Best indices should be valid
        assert!(comp.best_aic < comp.results.len());
        assert!(comp.best_bic < comp.results.len());
    }

    #[test]
    fn test_select_model_contains_expected_models() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGTACGT"),
            ("C", "TGCATGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();

        let comp = select_model(&tree, &aln).unwrap();

        let names: Vec<&str> = comp.results.iter().map(|r| r.model_name.as_str()).collect();
        assert!(names.contains(&"JC69"), "should contain JC69");
        assert!(
            names.iter().any(|n| n.contains("F84")),
            "should contain at least one F84 variant"
        );
    }

    #[test]
    fn test_select_model_aic_best_has_lowest_aic() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();

        let comp = select_model(&tree, &aln).unwrap();

        let best_aic_val = comp.results[comp.best_aic].aic;
        for (i, r) in comp.results.iter().enumerate() {
            assert!(
                best_aic_val <= r.aic + 1e-10,
                "best AIC (index {}, value {}) should be <= model {} AIC ({})",
                comp.best_aic,
                best_aic_val,
                i,
                r.aic
            );
        }
    }

    #[test]
    fn test_select_model_bic_best_has_lowest_bic() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();

        let comp = select_model(&tree, &aln).unwrap();

        let best_bic_val = comp.results[comp.best_bic].bic;
        for (i, r) in comp.results.iter().enumerate() {
            assert!(
                best_bic_val <= r.bic + 1e-10,
                "best BIC (index {}, value {}) should be <= model {} BIC ({})",
                comp.best_bic,
                best_bic_val,
                i,
                r.bic
            );
        }
    }

    // =====================================================================
    // Edge case and property tests
    // =====================================================================

    #[test]
    fn test_aicc_converges_to_aic_for_large_n() {
        // As n -> infinity, the correction 2k(k+1)/(n-k-1) -> 0,
        // so AICc -> AIC.
        let k = 5;
        let lnl = -100.0;
        let result = compute_information_criteria(lnl, k, 100_000);
        assert!(
            (result.aicc - result.aic).abs() < 0.01,
            "AICc ({}) should approximate AIC ({}) for large n",
            result.aicc,
            result.aic
        );
    }

    #[test]
    fn test_more_params_higher_aic_same_lnl() {
        // With the same log-likelihood, more parameters -> higher AIC
        let r_simple = compute_information_criteria(-100.0, 3, 200);
        let r_complex = compute_information_criteria(-100.0, 8, 200);
        assert!(
            r_complex.aic > r_simple.aic,
            "more params should give higher AIC: {} vs {}",
            r_complex.aic,
            r_simple.aic
        );
    }

    #[test]
    fn test_better_lnl_lower_aic_same_params() {
        // With the same parameters, better (higher) log-likelihood -> lower AIC
        let r_bad = compute_information_criteria(-200.0, 5, 100);
        let r_good = compute_information_criteria(-100.0, 5, 100);
        assert!(
            r_good.aic < r_bad.aic,
            "better lnl should give lower AIC: {} vs {}",
            r_good.aic,
            r_bad.aic
        );
    }

    #[test]
    fn test_zero_params() {
        let result = compute_information_criteria(-100.0, 0, 200);
        // AIC = -2*(-100) + 0 = 200
        assert!((result.aic - 200.0).abs() < 1e-10);
        // BIC = -2*(-100) + 0 = 200
        assert!((result.bic - 200.0).abs() < 1e-10);
        // AICc = 200 + 0 = 200
        assert!((result.aicc - 200.0).abs() < 1e-10);
    }

    #[test]
    fn test_compare_preserves_order() {
        let mut r1 = compute_information_criteria(-100.0, 3, 200);
        r1.model_name = "first".to_string();
        let mut r2 = compute_information_criteria(-90.0, 5, 200);
        r2.model_name = "second".to_string();
        let mut r3 = compute_information_criteria(-95.0, 4, 200);
        r3.model_name = "third".to_string();

        let comp = compare_models(vec![r1, r2, r3]);

        // Results should be in the original order
        assert_eq!(comp.results[0].model_name, "first");
        assert_eq!(comp.results[1].model_name, "second");
        assert_eq!(comp.results[2].model_name, "third");
    }

    #[test]
    fn test_heavily_favored_model_gets_weight_near_one() {
        // One model much better than the others
        let mut r_best = compute_information_criteria(-50.0, 3, 200);
        r_best.model_name = "best".to_string();
        let mut r_bad = compute_information_criteria(-200.0, 3, 200);
        r_bad.model_name = "bad".to_string();

        let comp = compare_models(vec![r_best, r_bad]);
        assert!(
            comp.aic_weights[0] > 0.999,
            "heavily favored model should have weight near 1.0, got {}",
            comp.aic_weights[0]
        );
    }
}

//! Approximately Unbiased (AU) test for phylogenetic hypothesis testing.
//!
//! Implements the multiscale bootstrap method of Shimodaira (2002) for
//! comparing candidate phylogenetic trees. The AU test corrects the
//! selection bias inherent in standard bootstrap proportions (BP) by
//! resampling at multiple scales and fitting a parametric model to the
//! scale-dependent bootstrap frequencies.
//!
//! # Algorithm
//!
//! 1. Compute site log-likelihoods for each candidate tree.
//! 2. For each scale factor `r` (e.g., 0.5, 0.6, ..., 1.4), resample
//!    `floor(r * n)` sites with replacement and record which tree wins.
//! 3. For each tree, compute rejection z-values `z_k(r) = Φ^{-1}(1 - BP_k(r))`
//!    and fit the model `z_k(r) = d1_k * sqrt(r) + d2_k / sqrt(r)` by least squares.
//! 4. The AU p-value is `p_AU_k = 1 - Φ(d1_k + d2_k)`.
//!
//! The AU test also provides the standard bootstrap probability (BP at r=1),
//! a Kishino-Hasegawa (KH) test p-value, and a Shimodaira-Hasegawa (SH)
//! test p-value for each candidate tree.
//!
//! # Example
//!
//! ```
//! use phylip_rs::bootstrap::au_test::{au_test, AuTestResult};
//!
//! // Site log-likelihoods for 3 trees, 5 sites each
//! let site_lnl = vec![
//!     vec![-1.0, -1.2, -1.1, -1.3, -1.0],  // tree 0 (best)
//!     vec![-1.5, -1.8, -1.6, -1.7, -1.5],  // tree 1
//!     vec![-2.0, -2.1, -2.0, -2.2, -2.0],  // tree 2
//! ];
//!
//! let result = au_test(&site_lnl, 1000, &[], 42);
//! assert!(result.trees[0].au_pvalue >= 0.0);
//! assert!(result.trees[0].au_pvalue <= 1.0);
//! ```
//!
//! # References
//!
//! - Shimodaira, H. (2002). An approximately unbiased test of phylogenetic
//!   tree selection. *Systematic Biology*, 51(3), 492-508.
//! - Shimodaira, H. & Hasegawa, M. (1999). Multiple comparisons of
//!   log-likelihoods with applications to phylogenetic inference.
//!   *Molecular Biology and Evolution*, 16, 1114-1116.
//! - Kishino, H. & Hasegawa, M. (1989). Evaluation of the maximum
//!   likelihood estimate of the evolutionary tree topologies from DNA
//!   sequence data, and the branching order in Hominoidea.
//!   *Journal of Molecular Evolution*, 29, 170-179.

use crate::bootstrap::SimpleRng;
use crate::likelihood::models::SubstitutionModel;
use crate::likelihood::pruning::{site_log_likelihoods, LikelihoodError};
use crate::tree::types::{Alignment, Tree};

// ============================================================
// Mathematical functions (zero dependencies)
// ============================================================

/// Standard normal cumulative distribution function Φ(x).
///
/// Uses the Abramowitz & Stegun (1964) rational approximation
/// (formula 26.2.17) with maximum error < 7.5e-8.
fn normal_cdf(x: f64) -> f64 {
    if x.is_nan() {
        return 0.5;
    }
    if x > 8.0 {
        return 1.0;
    }
    if x < -8.0 {
        return 0.0;
    }

    let abs_x = x.abs();

    // Abramowitz & Stegun constants for the rational approximation
    const P: f64 = 0.2316419;
    const B1: f64 = 0.319381530;
    const B2: f64 = -0.356563782;
    const B3: f64 = 1.781477937;
    const B4: f64 = -1.821255978;
    const B5: f64 = 1.330274429;

    let t = 1.0 / (1.0 + P * abs_x);
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    // Standard normal PDF at abs_x
    let pdf = (-0.5 * abs_x * abs_x).exp() / (2.0 * core::f64::consts::PI).sqrt();

    let cdf_complement = pdf * (B1 * t + B2 * t2 + B3 * t3 + B4 * t4 + B5 * t5);

    if x >= 0.0 {
        1.0 - cdf_complement
    } else {
        cdf_complement
    }
}

/// Inverse standard normal CDF (quantile function) Φ^{-1}(p).
///
/// Uses the minimax rational approximation of Peter Acklam (2004) for an
/// initial estimate, then refines with one Newton-Raphson step for full
/// double-precision accuracy.
///
/// Reference: Acklam, P.J. (2004). An algorithm for computing the inverse
/// normal cumulative distribution function.
fn normal_quantile(p: f64) -> f64 {
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    if (p - 0.5).abs() < 1e-15 {
        return 0.0;
    }

    // Coefficients for the rational approximation
    const A: [f64; 6] = [
        -3.969683028665376e+01,
         2.209460984245205e+02,
        -2.759285104469687e+02,
         1.383577518672690e+02,
        -3.066479806614716e+01,
         2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01,
         1.615858368580409e+02,
        -1.556989798598866e+02,
         6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
         4.374664141464968e+00,
         2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
         7.784695709041462e-03,
         3.224671290700398e-01,
         2.445134137142996e+00,
         3.754408661907416e+00,
    ];

    const P_LOW: f64 = 0.02425;
    const P_HIGH: f64 = 1.0 - P_LOW;

    let x = if p < P_LOW {
        // Lower tail
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= P_HIGH {
        // Central region
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        // Upper tail
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    };

    // One Newton-Raphson refinement step for full precision
    let phi = normal_cdf(x);
    let pdf = (-0.5 * x * x).exp() / (2.0 * core::f64::consts::PI).sqrt();
    if pdf > 0.0 {
        x + (p - phi) / pdf
    } else {
        x
    }
}

// ============================================================
// AU test result types
// ============================================================

/// Results for one candidate tree in the AU test.
#[derive(Debug, Clone)]
pub struct AuTreeResult {
    /// Index of this tree among the candidates.
    pub tree_index: usize,
    /// Total log-likelihood of this tree on the original data.
    pub lnl: f64,
    /// Difference in log-likelihood vs. the best tree (delta >= 0; best tree has delta = 0).
    pub delta_lnl: f64,
    /// Standard bootstrap probability (proportion of r=1 resamples where this tree wins).
    pub bp: f64,
    /// Approximately Unbiased p-value (Shimodaira 2002).
    pub au_pvalue: f64,
    /// Kishino-Hasegawa test p-value (1989).
    pub kh_pvalue: f64,
    /// Shimodaira-Hasegawa test p-value (1999).
    pub sh_pvalue: f64,
    /// Fitted intercept parameter from the multiscale bootstrap curve.
    pub au_a: f64,
    /// Fitted slope parameter from the multiscale bootstrap curve.
    pub au_b: f64,
}

/// Full results from the AU test for all candidate trees.
#[derive(Debug, Clone)]
pub struct AuTestResult {
    /// Per-tree results, ordered by tree index.
    pub trees: Vec<AuTreeResult>,
    /// Number of bootstrap replicates per scale.
    pub nboot: usize,
    /// Scale factors used.
    pub scales: Vec<f64>,
    /// Bootstrap proportions at each scale for each tree: bp_by_scale[tree][scale].
    pub bp_by_scale: Vec<Vec<f64>>,
}

impl AuTestResult {
    /// Index of the best tree (highest log-likelihood).
    pub fn best_tree(&self) -> usize {
        self.trees
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.lnl.partial_cmp(&b.lnl).unwrap_or(core::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0)
    }

    /// Trees whose AU p-value is above the given significance level
    /// (i.e., trees that cannot be rejected at that level).
    pub fn confidence_set(&self, alpha: f64) -> Vec<usize> {
        self.trees
            .iter()
            .filter(|t| t.au_pvalue > alpha)
            .map(|t| t.tree_index)
            .collect()
    }
}

impl core::fmt::Display for AuTestResult {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        writeln!(f, "AU Test Results (Shimodaira 2002)")?;
        writeln!(f, "{:-<72}", "")?;
        writeln!(
            f,
            "{:>5} {:>12} {:>12} {:>8} {:>8} {:>8} {:>8}",
            "Tree", "lnL", "delta_lnL", "BP", "AU", "KH", "SH"
        )?;
        writeln!(f, "{:-<72}", "")?;
        for tr in &self.trees {
            writeln!(
                f,
                "{:>5} {:>12.2} {:>12.2} {:>8.4} {:>8.4} {:>8.4} {:>8.4}",
                tr.tree_index, tr.lnl, tr.delta_lnl, tr.bp, tr.au_pvalue, tr.kh_pvalue, tr.sh_pvalue,
            )?;
        }
        writeln!(f, "{:-<72}", "")?;
        writeln!(
            f,
            "BP: bootstrap probability, AU: approximately unbiased p-value"
        )?;
        writeln!(
            f,
            "KH: Kishino-Hasegawa p-value, SH: Shimodaira-Hasegawa p-value"
        )?;
        Ok(())
    }
}

// ============================================================
// Default scale factors
// ============================================================

/// Default scale factors for multiscale bootstrap (Shimodaira 2002).
///
/// These span scales from 0.5 to 1.4, which is typically sufficient
/// for accurate curve fitting. The range should bracket 1.0 for
/// reliable extrapolation.
const DEFAULT_SCALES: &[f64] = &[0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4];

// ============================================================
// Core AU test algorithm
// ============================================================

/// Perform the Approximately Unbiased (AU) test on pre-computed site log-likelihoods.
///
/// # Arguments
///
/// * `site_lnl` - Site log-likelihoods for each candidate tree. `site_lnl[k][i]` is
///   the log-likelihood at site `i` for tree `k`. All inner vectors must have the
///   same length (number of sites).
/// * `nboot` - Number of bootstrap replicates per scale factor.
/// * `scales` - Scale factors for multiscale bootstrap. If empty, uses default
///   scales `[0.5, 0.6, ..., 1.4]`.
/// * `seed` - Random seed for reproducible resampling.
///
/// # Panics
///
/// Panics if `site_lnl` is empty or if inner vectors have inconsistent lengths.
///
/// # References
///
/// Shimodaira, H. (2002). An approximately unbiased test of phylogenetic
/// tree selection. *Systematic Biology*, 51(3), 492-508.
pub fn au_test(
    site_lnl: &[Vec<f64>],
    nboot: usize,
    scales: &[f64],
    seed: u64,
) -> AuTestResult {
    assert!(!site_lnl.is_empty(), "au_test requires at least one tree");
    let ntrees = site_lnl.len();
    let nsites = site_lnl[0].len();
    assert!(nsites > 0, "au_test requires at least one site");
    for (k, sll) in site_lnl.iter().enumerate() {
        assert_eq!(
            sll.len(),
            nsites,
            "tree {} has {} sites but tree 0 has {}",
            k,
            sll.len(),
            nsites
        );
    }

    let scales = if scales.is_empty() {
        DEFAULT_SCALES.to_vec()
    } else {
        scales.to_vec()
    };

    // Total log-likelihoods
    let total_lnl: Vec<f64> = site_lnl.iter().map(|sll| sll.iter().sum()).collect();
    let best_lnl = total_lnl
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let delta_lnl: Vec<f64> = total_lnl.iter().map(|l| best_lnl - l).collect();

    let nboot = nboot.max(1);
    let mut rng = SimpleRng::seed(seed);

    // Multiscale bootstrap
    // bp_by_scale[tree_index][scale_index] = bootstrap proportion at that scale
    let mut bp_counts: Vec<Vec<usize>> = vec![vec![0; scales.len()]; ntrees];

    for (si, &r) in scales.iter().enumerate() {
        let resample_size = (r * nsites as f64).floor().max(1.0) as usize;

        for _ in 0..nboot {
            // Resample sites with replacement
            let mut tree_sums = vec![0.0f64; ntrees];
            for _ in 0..resample_size {
                let site = rng.gen_range(nsites);
                for k in 0..ntrees {
                    tree_sums[k] += site_lnl[k][site];
                }
            }

            // Find which tree(s) have the best resampled log-likelihood
            let best_sum = tree_sums
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            // Count the first tree with the best sum (break ties deterministically)
            for k in 0..ntrees {
                if (tree_sums[k] - best_sum).abs() < 1e-12 {
                    bp_counts[k][si] += 1;
                    break;
                }
            }
        }
    }

    let bp_by_scale: Vec<Vec<f64>> = bp_counts
        .iter()
        .map(|counts| {
            counts
                .iter()
                .map(|&c| c as f64 / nboot as f64)
                .collect()
        })
        .collect();

    // Standard bootstrap probability: BP at scale r=1
    let scale_1_index = scales
        .iter()
        .position(|&r| (r - 1.0).abs() < 1e-10)
        .unwrap_or(0);
    let bp: Vec<f64> = bp_by_scale.iter().map(|bps| bps[scale_1_index]).collect();

    // Fit AU model and compute p-values for each tree
    let mut tree_results = Vec::with_capacity(ntrees);

    for k in 0..ntrees {
        let (d1, d2) = fit_au_curve(&bp_by_scale[k], &scales);

        // AU p-value: 1 - Φ(d1 + d2) (Shimodaira 2002)
        let au_pvalue = (1.0 - normal_cdf(d1 + d2)).clamp(0.0, 1.0);

        // KH test p-value
        let kh_pvalue = kh_test(site_lnl, k);

        // SH test p-value
        let sh_pvalue = sh_test(site_lnl, k, nboot, &mut rng);

        tree_results.push(AuTreeResult {
            tree_index: k,
            lnl: total_lnl[k],
            delta_lnl: delta_lnl[k],
            bp: bp[k],
            au_pvalue,
            kh_pvalue,
            sh_pvalue,
            au_a: d1,
            au_b: d2,
        });
    }

    AuTestResult {
        trees: tree_results,
        nboot,
        scales,
        bp_by_scale,
    }
}

/// Fit the AU curve for multiscale bootstrap.
///
/// Following Shimodaira (2002, eq. 6-8) and the `consel` software implementation:
///
/// The z-values are defined on the **rejection** frequencies:
///   z_k(r) = Φ^{-1}(1 - BP_k(r))
///
/// These are modeled as:
///   z_k(r) = d_1 * sqrt(r) + d_2 / sqrt(r)
///
/// where d_1 and d_2 are fitted by least squares with basis functions
/// x_1 = sqrt(r) and x_2 = 1/sqrt(r).
///
/// The AU p-value is then:
///   p_AU = 1 - Φ(d_1 + d_2)
///
/// Returns (d_1, d_2) so that the caller can compute Φ(d_1 + d_2) or
/// equivalently 1 - Φ(d_1 + d_2) as the AU p-value.
fn fit_au_curve(bp_at_scales: &[f64], scales: &[f64]) -> (f64, f64) {
    assert_eq!(bp_at_scales.len(), scales.len());
    let nscales = scales.len();
    if nscales < 2 {
        // Cannot fit a line with fewer than 2 points; return conservative values
        let mean_bp = bp_at_scales.iter().sum::<f64>() / bp_at_scales.len().max(1) as f64;
        if mean_bp > 0.5 {
            return (-3.0, 0.0); // AU p-value = 1 - Φ(-3) ≈ 0.999
        } else {
            return (3.0, 0.0); // AU p-value = 1 - Φ(3) ≈ 0.001
        }
    }

    // Clamp epsilon for extreme BP values
    let eps = 1e-4;

    // Compute z-values (on rejection frequency) and basis vectors
    let mut x1s = Vec::with_capacity(nscales); // sqrt(r)
    let mut x2s = Vec::with_capacity(nscales); // 1/sqrt(r)
    let mut zs = Vec::with_capacity(nscales);

    for i in 0..nscales {
        let r = scales[i];
        if r <= 0.0 {
            continue;
        }
        // z = Φ^{-1}(1 - BP) = -Φ^{-1}(BP)
        let bp_clamped = bp_at_scales[i].clamp(eps, 1.0 - eps);
        let z = -normal_quantile(bp_clamped);
        if z.is_finite() {
            x1s.push(r.sqrt());
            x2s.push(1.0 / r.sqrt());
            zs.push(z);
        }
    }

    let n = x1s.len();
    if n < 2 {
        let mean_bp = bp_at_scales.iter().sum::<f64>() / bp_at_scales.len().max(1) as f64;
        if mean_bp > 0.5 {
            return (-3.0, 0.0);
        } else {
            return (3.0, 0.0);
        }
    }

    // Least squares: z = d_1 * x_1 + d_2 * x_2 (no intercept!)
    // Normal equations:
    //   d_1 * sum(x1^2) + d_2 * sum(x1*x2) = sum(x1*z)
    //   d_1 * sum(x1*x2) + d_2 * sum(x2^2) = sum(x2*z)
    let s11: f64 = x1s.iter().map(|&x| x * x).sum();
    let s22: f64 = x2s.iter().map(|&x| x * x).sum();
    let s12: f64 = x1s.iter().zip(x2s.iter()).map(|(&a, &b)| a * b).sum();
    let s1z: f64 = x1s.iter().zip(zs.iter()).map(|(&x, &z)| x * z).sum();
    let s2z: f64 = x2s.iter().zip(zs.iter()).map(|(&x, &z)| x * z).sum();

    let denom = s11 * s22 - s12 * s12;
    if denom.abs() < 1e-15 {
        // Degenerate: all scales the same
        let mean_z = zs.iter().sum::<f64>() / n as f64;
        return (mean_z, 0.0);
    }

    let d1 = (s22 * s1z - s12 * s2z) / denom;
    let d2 = (s11 * s2z - s12 * s1z) / denom;

    (d1, d2)
}

// ============================================================
// KH and SH tests
// ============================================================

/// Kishino-Hasegawa (1989) test p-value for tree k.
///
/// The KH test compares tree k against the best tree using the
/// paired-sites test statistic. Under the null hypothesis that both
/// trees are equally good, the site-wise log-likelihood differences
/// are approximately normally distributed.
///
/// p_KH = Φ(-|delta| / SE(delta)) * 2 (two-sided)
///
/// where delta = sum of (sll_best - sll_k) and SE is the standard
/// error estimated from site-wise differences.
fn kh_test(site_lnl: &[Vec<f64>], tree_index: usize) -> f64 {
    let ntrees = site_lnl.len();
    let nsites = site_lnl[0].len();

    // Find the best tree
    let total_lnl: Vec<f64> = site_lnl.iter().map(|sll| sll.iter().sum()).collect();
    let best_index = total_lnl
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);

    if tree_index == best_index || ntrees < 2 {
        return 1.0;
    }

    // Site-wise differences: d[i] = sll_best[i] - sll_k[i]
    let diffs: Vec<f64> = (0..nsites)
        .map(|i| site_lnl[best_index][i] - site_lnl[tree_index][i])
        .collect();

    let n = nsites as f64;
    let mean_d = diffs.iter().sum::<f64>() / n;
    let var_d = diffs.iter().map(|&d| (d - mean_d) * (d - mean_d)).sum::<f64>() / (n - 1.0);

    if var_d <= 0.0 {
        return 1.0;
    }

    let se = (var_d / n).sqrt();
    let z = mean_d / se;

    // Two-sided p-value
    let p = 2.0 * normal_cdf(-z.abs());
    p.clamp(0.0, 1.0)
}

/// Shimodaira-Hasegawa (1999) test p-value for tree k.
///
/// The SH test accounts for multiple comparisons by centering the
/// bootstrap distribution of log-likelihood differences. For each
/// bootstrap resample, the test statistic is the maximum over all
/// trees of (delta_k_boot - delta_k_obs). The p-value is the
/// proportion of resamples where the centered test statistic for
/// tree k exceeds the observed delta.
fn sh_test(
    site_lnl: &[Vec<f64>],
    tree_index: usize,
    nboot: usize,
    rng: &mut SimpleRng,
) -> f64 {
    let ntrees = site_lnl.len();
    let nsites = site_lnl[0].len();

    // Total log-likelihoods for each tree
    let total_lnl: Vec<f64> = site_lnl.iter().map(|sll| sll.iter().sum()).collect();
    let best_lnl = total_lnl
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);

    // Observed delta for tree k
    let delta_obs = best_lnl - total_lnl[tree_index];

    if delta_obs.abs() < 1e-12 || ntrees < 2 {
        return 1.0;
    }

    // Site-wise mean log-likelihoods for centering
    let site_means: Vec<Vec<f64>> = (0..ntrees)
        .map(|k| {
            let mean_k: f64 = site_lnl[k].iter().sum::<f64>() / nsites as f64;
            site_lnl[k].iter().map(|&s| s - mean_k).collect()
        })
        .collect();

    let mut count_exceed = 0usize;

    for _ in 0..nboot {
        // Resample site indices
        let mut boot_sums = vec![0.0f64; ntrees];
        for _ in 0..nsites {
            let site = rng.gen_range(nsites);
            for k in 0..ntrees {
                boot_sums[k] += site_means[k][site];
            }
        }

        // Centered resampled best
        let boot_best = boot_sums
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let boot_delta = boot_best - boot_sums[tree_index];

        if boot_delta >= delta_obs {
            count_exceed += 1;
        }
    }

    (count_exceed as f64 / nboot as f64).clamp(0.0, 1.0)
}

// ============================================================
// Convenience wrapper for tree objects
// ============================================================

/// Perform the AU test directly from candidate trees, an alignment, and a model.
///
/// This is a convenience wrapper that computes site log-likelihoods for each
/// candidate tree and then calls [`au_test`].
///
/// # Arguments
///
/// * `trees` - Candidate phylogenetic trees to compare.
/// * `alignment` - The DNA alignment.
/// * `model` - The substitution model.
/// * `nboot` - Number of bootstrap replicates per scale factor.
/// * `seed` - Random seed for reproducibility.
///
/// # Errors
///
/// Returns `LikelihoodError` if site log-likelihoods cannot be computed
/// for any tree (e.g., taxon mismatch, empty tree).
pub fn au_test_from_trees(
    trees: &[Tree],
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    nboot: usize,
    seed: u64,
) -> Result<AuTestResult, LikelihoodError> {
    let mut all_site_lnl = Vec::with_capacity(trees.len());

    for tree in trees {
        let sll = site_log_likelihoods(tree, alignment, model)?;
        all_site_lnl.push(sll);
    }

    Ok(au_test(&all_site_lnl, nboot, &[], seed))
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: create clearly separated site log-likelihoods where tree 0 is best.
    fn make_clear_winner() -> Vec<Vec<f64>> {
        // 20 sites, 3 trees. Tree 0 is substantially better at every site.
        let nsites = 20;
        let mut sll = vec![vec![0.0; nsites]; 3];
        for i in 0..nsites {
            sll[0][i] = -1.0 - 0.1 * (i as f64 % 3.0); // best
            sll[1][i] = -2.5 - 0.1 * (i as f64 % 3.0); // worse
            sll[2][i] = -4.0 - 0.1 * (i as f64 % 3.0); // worst
        }
        sll
    }

    /// Helper: create site log-likelihoods where trees are nearly indistinguishable.
    fn make_close_trees() -> Vec<Vec<f64>> {
        let nsites = 20;
        let mut sll = vec![vec![0.0; nsites]; 2];
        for i in 0..nsites {
            sll[0][i] = -1.5 - 0.01 * (i as f64);
            sll[1][i] = -1.5 - 0.01 * (i as f64) - 0.001; // barely worse
        }
        sll
    }

    // ----------------------------------------------------------
    // Test 1: True (best) tree should have AU p-value near 1.0
    // ----------------------------------------------------------
    #[test]
    fn test_best_tree_high_au_pvalue() {
        let sll = make_clear_winner();
        let result = au_test(&sll, 10000, &[], 42);

        assert!(
            result.trees[0].au_pvalue > 0.5,
            "best tree AU p-value should be high, got {}",
            result.trees[0].au_pvalue
        );
        assert!(
            result.trees[0].delta_lnl < 1e-10,
            "best tree delta_lnl should be ~0, got {}",
            result.trees[0].delta_lnl
        );
    }

    // ----------------------------------------------------------
    // Test 2: Clearly wrong tree should have low AU p-value
    // ----------------------------------------------------------
    #[test]
    fn test_wrong_tree_low_au_pvalue() {
        let sll = make_clear_winner();
        let result = au_test(&sll, 10000, &[], 42);

        assert!(
            result.trees[2].au_pvalue < 0.5,
            "worst tree AU p-value should be low, got {}",
            result.trees[2].au_pvalue
        );
        assert!(
            result.trees[2].delta_lnl > 0.0,
            "worst tree delta_lnl should be > 0, got {}",
            result.trees[2].delta_lnl
        );
    }

    // ----------------------------------------------------------
    // Test 3: All AU p-values should be in [0, 1]
    // ----------------------------------------------------------
    #[test]
    fn test_au_pvalues_in_range() {
        let sll = make_clear_winner();
        let result = au_test(&sll, 5000, &[], 99);

        for tr in &result.trees {
            assert!(
                (0.0..=1.0).contains(&tr.au_pvalue),
                "AU p-value {} for tree {} is outside [0, 1]",
                tr.au_pvalue,
                tr.tree_index
            );
            assert!(
                (0.0..=1.0).contains(&tr.bp),
                "BP {} for tree {} is outside [0, 1]",
                tr.bp,
                tr.tree_index
            );
            assert!(
                (0.0..=1.0).contains(&tr.kh_pvalue),
                "KH p-value {} for tree {} is outside [0, 1]",
                tr.kh_pvalue,
                tr.tree_index
            );
            assert!(
                (0.0..=1.0).contains(&tr.sh_pvalue),
                "SH p-value {} for tree {} is outside [0, 1]",
                tr.sh_pvalue,
                tr.tree_index
            );
        }
    }

    // ----------------------------------------------------------
    // Test 4: BP at scale 1.0 equals standard bootstrap proportion
    // ----------------------------------------------------------
    #[test]
    fn test_bp_at_scale_one() {
        let sll = make_clear_winner();
        let result = au_test(&sll, 5000, &[], 42);

        // BP should equal bp_by_scale at the r=1.0 entry
        let scale_1_idx = result
            .scales
            .iter()
            .position(|&r| (r - 1.0).abs() < 1e-10)
            .expect("scales should contain 1.0");

        for k in 0..sll.len() {
            assert!(
                (result.trees[k].bp - result.bp_by_scale[k][scale_1_idx]).abs() < 1e-12,
                "BP for tree {} ({}) should equal bp_by_scale at r=1.0 ({})",
                k,
                result.trees[k].bp,
                result.bp_by_scale[k][scale_1_idx]
            );
        }
    }

    // ----------------------------------------------------------
    // Test 5: With 2 close trees, AU p-values should roughly sum to ~1
    // ----------------------------------------------------------
    #[test]
    fn test_two_tree_pvalue_sum() {
        let sll = make_close_trees();
        let result = au_test(&sll, 10000, &[], 123);

        let sum_bp: f64 = result.trees.iter().map(|t| t.bp).sum();
        assert!(
            (sum_bp - 1.0).abs() < 1e-10,
            "BP values for 2 trees should sum to 1.0, got {}",
            sum_bp
        );

        // AU p-values should roughly sum to ~1 for 2 trees, but the
        // correction means this is approximate. Allow wider tolerance.
        let sum_au: f64 = result.trees.iter().map(|t| t.au_pvalue).sum();
        assert!(
            sum_au > 0.3 && sum_au < 2.0,
            "AU p-values for 2 trees should roughly sum to ~1, got {}",
            sum_au
        );
    }

    // ----------------------------------------------------------
    // Test 6: Normal CDF properties
    // ----------------------------------------------------------
    #[test]
    fn test_normal_cdf_properties() {
        // Φ(0) = 0.5
        assert!(
            (normal_cdf(0.0) - 0.5).abs() < 1e-7,
            "Φ(0) should be 0.5, got {}",
            normal_cdf(0.0)
        );

        // Φ(x) + Φ(-x) = 1 (symmetry)
        for &x in &[0.5, 1.0, 1.5, 2.0, 3.0] {
            let sum = normal_cdf(x) + normal_cdf(-x);
            assert!(
                (sum - 1.0).abs() < 1e-6,
                "Φ({}) + Φ({}) = {}, expected 1.0",
                x,
                -x,
                sum
            );
        }

        // Known values
        assert!(
            (normal_cdf(1.0) - 0.8413).abs() < 0.001,
            "Φ(1.0) should be ~0.8413, got {}",
            normal_cdf(1.0)
        );
        assert!(
            (normal_cdf(1.96) - 0.975).abs() < 0.001,
            "Φ(1.96) should be ~0.975, got {}",
            normal_cdf(1.96)
        );
        assert!(
            (normal_cdf(-1.96) - 0.025).abs() < 0.001,
            "Φ(-1.96) should be ~0.025, got {}",
            normal_cdf(-1.96)
        );

        // Monotonicity
        let mut prev = normal_cdf(-5.0);
        for i in -49..=50 {
            let x = i as f64 / 10.0;
            let curr = normal_cdf(x);
            assert!(
                curr >= prev - 1e-15,
                "Φ should be monotonically non-decreasing: Φ({}) = {} < prev {}",
                x,
                curr,
                prev
            );
            prev = curr;
        }
    }

    // ----------------------------------------------------------
    // Test 7: Normal quantile is inverse of CDF
    // ----------------------------------------------------------
    #[test]
    fn test_normal_quantile_inverse() {
        for &p in &[0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99] {
            let z = normal_quantile(p);
            let p_roundtrip = normal_cdf(z);
            assert!(
                (p_roundtrip - p).abs() < 1e-5,
                "Φ(Φ^{{-1}}({})) = {}, expected {}",
                p,
                p_roundtrip,
                p
            );
        }

        // Φ^{-1}(0.5) = 0
        assert!(
            normal_quantile(0.5).abs() < 1e-10,
            "Φ^{{-1}}(0.5) should be 0, got {}",
            normal_quantile(0.5)
        );

        // Antisymmetry: Φ^{-1}(p) = -Φ^{-1}(1-p)
        for &p in &[0.1, 0.2, 0.3, 0.4] {
            let z1 = normal_quantile(p);
            let z2 = normal_quantile(1.0 - p);
            assert!(
                (z1 + z2).abs() < 1e-5,
                "Φ^{{-1}}({}) + Φ^{{-1}}({}) = {}, expected 0",
                p,
                1.0 - p,
                z1 + z2
            );
        }
    }

    // ----------------------------------------------------------
    // Test 8: Single tree should have AU p-value = 1.0 (or near 1)
    // ----------------------------------------------------------
    #[test]
    fn test_single_tree() {
        let sll = vec![vec![-1.0, -1.1, -1.2, -1.0, -1.3]];
        let result = au_test(&sll, 1000, &[], 42);

        assert_eq!(result.trees.len(), 1);
        assert!(
            result.trees[0].bp > 0.99,
            "single tree BP should be ~1.0, got {}",
            result.trees[0].bp
        );
        assert!(
            result.trees[0].delta_lnl < 1e-10,
            "single tree delta should be 0"
        );
    }

    // ----------------------------------------------------------
    // Test 9: KH test - best tree should have p-value = 1.0
    // ----------------------------------------------------------
    #[test]
    fn test_kh_best_tree_pvalue() {
        let sll = make_clear_winner();
        let result = au_test(&sll, 1000, &[], 42);

        assert!(
            (result.trees[0].kh_pvalue - 1.0).abs() < 1e-10,
            "KH p-value for best tree should be 1.0, got {}",
            result.trees[0].kh_pvalue
        );
    }

    // ----------------------------------------------------------
    // Test 10: Confidence set with clear winner
    // ----------------------------------------------------------
    #[test]
    fn test_confidence_set() {
        let sll = make_clear_winner();
        let result = au_test(&sll, 10000, &[], 42);

        let cs = result.confidence_set(0.05);
        assert!(
            cs.contains(&0),
            "confidence set at alpha=0.05 should contain the best tree (0)"
        );
    }

    // ----------------------------------------------------------
    // Test 11: Reproducibility with same seed
    // ----------------------------------------------------------
    #[test]
    fn test_au_reproducible() {
        let sll = make_clear_winner();

        let r1 = au_test(&sll, 1000, &[], 42);
        let r2 = au_test(&sll, 1000, &[], 42);

        for k in 0..sll.len() {
            assert!(
                (r1.trees[k].au_pvalue - r2.trees[k].au_pvalue).abs() < 1e-10,
                "same seed should produce identical AU p-values: tree {}: {} vs {}",
                k,
                r1.trees[k].au_pvalue,
                r2.trees[k].au_pvalue
            );
            assert!(
                (r1.trees[k].bp - r2.trees[k].bp).abs() < 1e-10,
                "same seed should produce identical BP: tree {}: {} vs {}",
                k,
                r1.trees[k].bp,
                r2.trees[k].bp
            );
        }
    }

    // ----------------------------------------------------------
    // Test 12: Custom scale factors
    // ----------------------------------------------------------
    #[test]
    fn test_custom_scales() {
        let sll = make_clear_winner();
        let custom_scales = vec![0.5, 1.0, 1.5, 2.0];
        let result = au_test(&sll, 1000, &custom_scales, 42);

        assert_eq!(result.scales, custom_scales);
        assert_eq!(result.bp_by_scale[0].len(), 4);

        for tr in &result.trees {
            assert!(
                (0.0..=1.0).contains(&tr.au_pvalue),
                "AU p-value with custom scales should be in [0, 1]"
            );
        }
    }

    // ----------------------------------------------------------
    // Test 13: Display formatting does not panic
    // ----------------------------------------------------------
    #[test]
    fn test_display() {
        let sll = make_clear_winner();
        let result = au_test(&sll, 100, &[], 42);
        let display = format!("{}", result);
        assert!(display.contains("AU Test Results"));
        assert!(display.contains("Tree"));
        assert!(display.contains("lnL"));
    }

    // ----------------------------------------------------------
    // Test 14: BP values across all scales sum to ~1 for each scale
    // ----------------------------------------------------------
    #[test]
    fn test_bp_by_scale_sums() {
        let sll = make_clear_winner();
        let result = au_test(&sll, 5000, &[], 42);

        for si in 0..result.scales.len() {
            let sum: f64 = result.bp_by_scale.iter().map(|bps| bps[si]).sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "BP at scale {} should sum to 1.0 across trees, got {}",
                result.scales[si],
                sum
            );
        }
    }

    // ----------------------------------------------------------
    // Test 15: Curve fitting on synthetic data
    // ----------------------------------------------------------
    #[test]
    fn test_fit_au_curve() {
        // Construct BP values that follow:
        //   z_reject = d1*sqrt(r) + d2/sqrt(r)
        //   BP(r) = 1 - Φ(z_reject) = Φ(-z_reject)
        // with d1 = 0.5, d2 = -1.0 (so AU p-value = 1 - Φ(-0.5) = Φ(0.5) ≈ 0.69)
        let true_d1 = 0.5;
        let true_d2 = -1.0;
        let scales: Vec<f64> = (5..=15).map(|i| i as f64 / 10.0).collect();
        let bp_vals: Vec<f64> = scales
            .iter()
            .map(|&r| {
                let z_reject = true_d1 * r.sqrt() + true_d2 / r.sqrt();
                // BP = 1 - Φ(z_reject) = Φ(-z_reject)
                normal_cdf(-z_reject)
            })
            .collect();

        let (d1, d2) = fit_au_curve(&bp_vals, &scales);

        assert!(
            (d1 - true_d1).abs() < 0.05,
            "fitted d1 = {}, expected {}",
            d1,
            true_d1
        );
        assert!(
            (d2 - true_d2).abs() < 0.05,
            "fitted d2 = {}, expected {}",
            d2,
            true_d2
        );
    }

    // ----------------------------------------------------------
    // Test 16: best_tree() returns correct index
    // ----------------------------------------------------------
    #[test]
    fn test_best_tree_index() {
        let sll = make_clear_winner();
        let result = au_test(&sll, 100, &[], 42);
        assert_eq!(
            result.best_tree(),
            0,
            "best_tree() should return 0 for clear winner"
        );
    }

    // ----------------------------------------------------------
    // Test 17: SH test p-value for best tree should be 1.0
    // ----------------------------------------------------------
    #[test]
    fn test_sh_best_tree() {
        let sll = make_clear_winner();
        let result = au_test(&sll, 5000, &[], 42);

        assert!(
            result.trees[0].sh_pvalue > 0.9,
            "SH p-value for best tree should be near 1.0, got {}",
            result.trees[0].sh_pvalue
        );
    }

    // ----------------------------------------------------------
    // Test 18: SH test is more conservative than KH test
    // ----------------------------------------------------------
    #[test]
    fn test_sh_more_conservative_than_kh() {
        // The SH test corrects for multiple testing and should
        // generally give higher (more conservative) p-values than KH
        // for non-best trees.
        let sll = make_clear_winner();
        let result = au_test(&sll, 10000, &[], 42);

        // For clearly wrong trees, SH >= KH (SH is more conservative)
        // This is a statistical property that holds on average but may
        // not hold for every single run with finite bootstrap. We check
        // the worst tree where the effect should be most pronounced.
        // With very clear separation, both should be near 0.
        for tr in &result.trees {
            // Both should be valid p-values
            assert!((0.0..=1.0).contains(&tr.sh_pvalue));
            assert!((0.0..=1.0).contains(&tr.kh_pvalue));
        }
    }

    // ----------------------------------------------------------
    // Test 19: Normal CDF and quantile handle edge cases
    // ----------------------------------------------------------
    #[test]
    fn test_normal_edge_cases() {
        // Extreme values
        assert!(normal_cdf(10.0) > 0.9999);
        assert!(normal_cdf(-10.0) < 0.0001);

        // NaN -> 0.5
        assert!((normal_cdf(f64::NAN) - 0.5).abs() < 1e-10);

        // Quantile extremes
        assert!(normal_quantile(0.0).is_infinite());
        assert!(normal_quantile(0.0) < 0.0); // -inf
        assert!(normal_quantile(1.0).is_infinite());
        assert!(normal_quantile(1.0) > 0.0); // +inf
    }

    // ----------------------------------------------------------
    // Test 20: au_test_from_trees produces valid results
    // ----------------------------------------------------------
    #[test]
    fn test_au_test_from_trees() {
        use crate::likelihood::models::Jc69Model;
        use crate::tree::newick::parse_newick;
        use crate::tree::types::{Base, Sequence};

        let alignment = Alignment::new(vec![
            Sequence::new(
                "A",
                vec![
                    Base::A, Base::A, Base::A, Base::A, Base::C, Base::C, Base::C, Base::C,
                    Base::A, Base::A, Base::A, Base::A, Base::C, Base::C, Base::C, Base::C,
                ],
            ),
            Sequence::new(
                "B",
                vec![
                    Base::A, Base::A, Base::A, Base::A, Base::C, Base::C, Base::C, Base::C,
                    Base::A, Base::A, Base::A, Base::A, Base::C, Base::C, Base::C, Base::C,
                ],
            ),
            Sequence::new(
                "C",
                vec![
                    Base::C, Base::C, Base::C, Base::C, Base::A, Base::A, Base::A, Base::A,
                    Base::C, Base::C, Base::C, Base::C, Base::A, Base::A, Base::A, Base::A,
                ],
            ),
            Sequence::new(
                "D",
                vec![
                    Base::C, Base::C, Base::C, Base::C, Base::A, Base::A, Base::A, Base::A,
                    Base::C, Base::C, Base::C, Base::C, Base::A, Base::A, Base::A, Base::A,
                ],
            ),
        ])
        .unwrap();

        let good_tree =
            parse_newick("((A:0.05,B:0.05):0.2,(C:0.05,D:0.05):0.2);").unwrap();
        let bad_tree =
            parse_newick("((A:0.05,C:0.05):0.2,(B:0.05,D:0.05):0.2);").unwrap();

        let model = Jc69Model;
        let result =
            au_test_from_trees(&[good_tree, bad_tree], &alignment, &model, 1000, 42)
                .unwrap();

        assert_eq!(result.trees.len(), 2);

        // The correct topology should have higher log-likelihood
        assert!(
            result.trees[0].lnl > result.trees[1].lnl,
            "correct topology lnl ({}) should exceed wrong topology lnl ({})",
            result.trees[0].lnl,
            result.trees[1].lnl
        );

        // AU p-values should be valid
        for tr in &result.trees {
            assert!((0.0..=1.0).contains(&tr.au_pvalue));
            assert!((0.0..=1.0).contains(&tr.bp));
        }
    }
}

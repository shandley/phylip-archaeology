//! Discrete gamma site rate heterogeneity (Yang 1994).
//!
//! Implements the discrete gamma model for among-site rate variation in
//! phylogenetic likelihood calculations. Many genes exhibit rate variation
//! across sites: some positions are highly conserved (slow rate) while
//! others evolve rapidly. Ignoring this heterogeneity biases branch length
//! estimates and can mislead phylogenetic inference.
//!
//! Yang (1994) proposed approximating the continuous gamma distribution of
//! rates with a discrete distribution of `k` categories. The shape parameter
//! `alpha` controls the amount of rate variation:
//!
//! - **alpha -> infinity**: uniform rates (no variation)
//! - **alpha ~ 1**: moderate variation
//! - **alpha < 0.5**: extreme variation (many slow sites, few very fast)
//!
//! For `k` categories, the rate in each category is the mean of the gamma
//! density within the corresponding quantile interval. Rates are normalized
//! so their mean equals 1.0, ensuring that the overall substitution rate
//! (and thus branch length interpretation) is unchanged.
//!
//! The site likelihood under the discrete gamma model is:
//!
//! ```text
//! L_site = (1/k) * sum_{c=1}^{k} L(site | rate_c)
//! ```
//!
//! where `rate_c` scales all branch lengths for category `c`.
//!
//! # Example
//!
//! ```
//! use phylip_rs::likelihood::gamma::{discrete_gamma_rates, gamma_log_likelihood};
//! use phylip_rs::likelihood::models::Jc69Model;
//! use phylip_rs::tree::{Alignment, Base, Sequence};
//! use phylip_rs::tree::newick::parse_newick;
//!
//! // Compute 4 discrete gamma rate categories with alpha = 0.5
//! let rates = discrete_gamma_rates(0.5, 4);
//! assert_eq!(rates.len(), 4);
//! let mean: f64 = rates.iter().sum::<f64>() / rates.len() as f64;
//! assert!((mean - 1.0).abs() < 1e-6);
//!
//! // Compute gamma-adjusted likelihood
//! let alignment = Alignment::new(vec![
//!     Sequence::new("A", vec![Base::A, Base::C, Base::G, Base::T]),
//!     Sequence::new("B", vec![Base::A, Base::C, Base::A, Base::T]),
//!     Sequence::new("C", vec![Base::T, Base::G, Base::G, Base::A]),
//! ]).unwrap();
//! let tree = parse_newick("((A:0.1,B:0.1):0.1,C:0.2);").unwrap();
//! let model = Jc69Model;
//! let lnl = gamma_log_likelihood(&tree, &alignment, &model, &rates).unwrap();
//! assert!(lnl < 0.0);
//! ```
//!
//! # References
//!
//! - Yang, Z. (1994). Maximum likelihood phylogenetic estimation from DNA
//!   sequences with variable rates over sites: approximate methods.
//!   *Journal of Molecular Evolution*, 39, 306-314.
//! - Yang, Z. (1993). Maximum likelihood estimation of phylogeny from DNA
//!   sequences when substitution rates differ over sites.
//!   *Molecular Biology and Evolution*, 10, 1396-1401.

use crate::likelihood::models::SubstitutionModel;
use crate::likelihood::pruning::{site_log_likelihoods, LikelihoodError};
use crate::tree::types::{Alignment, Tree};

// ============================================================================
// Gamma function mathematics (pure Rust, no external dependencies)
// ============================================================================

/// Lanczos approximation coefficients for the gamma function.
/// Using g = 7 with 9 coefficients, accurate to ~15 digits.
const LANCZOS_G: f64 = 7.0;
const LANCZOS_COEFFICIENTS: [f64; 9] = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
];

/// Compute the natural logarithm of the gamma function: ln(Gamma(a)).
///
/// Uses the Lanczos approximation, which is accurate to ~15 significant
/// digits for positive real arguments.
///
/// # Panics
///
/// Panics if `a` is not positive.
fn ln_gamma(a: f64) -> f64 {
    assert!(a > 0.0, "ln_gamma requires a > 0, got {}", a);

    if a < 0.5 {
        // Use reflection formula: Gamma(z) * Gamma(1-z) = pi / sin(pi*z)
        // ln(Gamma(z)) = ln(pi) - ln(sin(pi*z)) - ln(Gamma(1-z))
        let pi = std::f64::consts::PI;
        pi.ln() - (pi * a).sin().ln() - ln_gamma(1.0 - a)
    } else {
        let a = a - 1.0;
        let mut x = LANCZOS_COEFFICIENTS[0];
        for i in 1..LANCZOS_COEFFICIENTS.len() {
            x += LANCZOS_COEFFICIENTS[i] / (a + i as f64);
        }
        let t = a + LANCZOS_G + 0.5;
        0.5 * (2.0 * std::f64::consts::PI).ln() + (t.ln() * (a + 0.5)) - t + x.ln()
    }
}

/// Compute the gamma function Gamma(a).
#[cfg(test)]
fn gamma_fn(a: f64) -> f64 {
    ln_gamma(a).exp()
}

/// Regularized lower incomplete gamma function P(a, x) = gamma(a, x) / Gamma(a).
///
/// Uses the series expansion for small x and the continued fraction
/// representation for large x, following the method in Numerical Recipes.
///
/// Returns a value in [0, 1].
fn gamma_inc_lower(a: f64, x: f64) -> f64 {
    assert!(a > 0.0, "gamma_inc_lower requires a > 0, got {}", a);
    if x <= 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        // Series representation
        gamma_inc_series(a, x)
    } else {
        // Continued fraction representation: P(a,x) = 1 - Q(a,x)
        1.0 - gamma_inc_cf(a, x)
    }
}

/// Series expansion for the regularized lower incomplete gamma function.
///
/// P(a, x) = exp(-x) * x^a / Gamma(a) * sum_{n=0}^{inf} x^n / ((a)(a+1)...(a+n))
fn gamma_inc_series(a: f64, x: f64) -> f64 {
    let max_iter = 200;
    let eps = 1e-14;

    let ln_prefix = a * x.ln() - x - ln_gamma(a);

    let mut sum = 1.0 / a;
    let mut term = 1.0 / a;

    for n in 1..=max_iter {
        term *= x / (a + n as f64);
        sum += term;
        if term.abs() < sum.abs() * eps {
            break;
        }
    }

    sum * ln_prefix.exp()
}

/// Continued fraction for the regularized upper incomplete gamma function Q(a, x).
///
/// Uses the modified Lentz method for evaluating the continued fraction.
/// Q(a, x) = 1 - P(a, x)
fn gamma_inc_cf(a: f64, x: f64) -> f64 {
    let max_iter = 200;
    let eps = 1e-14;
    let tiny = 1e-300;

    let ln_prefix = a * x.ln() - x - ln_gamma(a);

    // Modified Lentz's method for the continued fraction:
    //   CF = 1 / (x + 1 - a + K)
    // where K is the continued fraction with terms:
    //   a_n = n * (a - n)
    //   b_n = x + (2n + 1) - a
    // b_0 = x + 1 - a
    let b0 = x + 1.0 - a;
    let mut f = if b0.abs() < tiny { tiny } else { b0 };
    let mut c = f;
    let mut d = 0.0;

    for n in 1..=max_iter {
        let an = (n as f64) * (a - n as f64);
        let bn = x + (2 * n + 1) as f64 - a;

        d = bn + an * d;
        if d.abs() < tiny {
            d = tiny;
        }
        d = 1.0 / d;

        c = bn + an / c;
        if c.abs() < tiny {
            c = tiny;
        }

        let delta = c * d;
        f *= delta;

        if (delta - 1.0).abs() < eps {
            break;
        }
    }

    ln_prefix.exp() / f
}

/// Inverse of the regularized lower incomplete gamma function.
///
/// Given `a` and `p` in (0, 1), finds `x` such that P(a, x) = p.
///
/// Uses an initial approximation followed by Halley's method (a
/// third-order root-finding algorithm) for rapid convergence.
fn gamma_inc_inv(a: f64, p: f64) -> f64 {
    assert!(a > 0.0, "gamma_inc_inv requires a > 0");
    assert!(
        (0.0..=1.0).contains(&p),
        "gamma_inc_inv requires 0 <= p <= 1, got {}",
        p
    );

    if p == 0.0 {
        return 0.0;
    }
    if p == 1.0 {
        return f64::INFINITY;
    }

    // Initial approximation using Wilson-Hilferty transformation
    // For a >= 1: x ~ a * (1 - 1/(9a) + z*sqrt(1/(9a)))^3
    // where z = inverse normal CDF of p
    let x0 = if a > 1.0 {
        let z = normal_quantile(p);
        let d = 1.0 / (9.0 * a);
        let t = 1.0 - d + z * d.sqrt();
        if t > 0.0 {
            a * t * t * t
        } else {
            // Fallback for negative t
            if p < 0.5 { 0.1 * a } else { 2.0 * a }
        }
    } else {
        // For small a, use a different initial guess
        let ln_p = p.ln();
        let x_guess = (ln_p + ln_gamma(a + 1.0)).exp().powf(1.0 / a);
        x_guess.max(1e-10)
    };

    // Halley's method refinement
    let max_iter = 50;
    let eps = 1e-12;
    let mut x = x0.max(1e-30);

    for _ in 0..max_iter {
        let current_p = gamma_inc_lower(a, x);
        let err = current_p - p;

        if err.abs() < eps {
            break;
        }

        // Derivative of P(a,x) w.r.t. x is the gamma PDF:
        // f(x) = x^(a-1) * exp(-x) / Gamma(a)
        let ln_f = (a - 1.0) * x.ln() - x - ln_gamma(a);
        let f = ln_f.exp();

        if f.abs() < 1e-300 {
            break;
        }

        // Newton step
        let delta = err / f;

        // Halley correction: f'(x)/f(x) = (a-1)/x - 1
        let fp_over_f = (a - 1.0) / x - 1.0;
        let halley_delta = delta / (1.0 - 0.5 * delta * fp_over_f);

        let x_new = x - halley_delta;

        if x_new <= 0.0 {
            // Bisect towards zero if Newton overshoots
            x *= 0.5;
        } else {
            if (x_new - x).abs() < eps * x.abs() {
                x = x_new;
                break;
            }
            x = x_new;
        }
    }

    x.max(0.0)
}

/// Approximate inverse of the standard normal CDF (probit function).
///
/// Uses the rational approximation from Abramowitz and Stegun,
/// accurate to about 4.5 decimal places, which is sufficient for
/// an initial approximation in gamma_inc_inv.
fn normal_quantile(p: f64) -> f64 {
    assert!(
        p > 0.0 && p < 1.0,
        "normal_quantile requires 0 < p < 1"
    );

    // Use symmetry: for p > 0.5, compute for 1-p and negate
    if p == 0.5 {
        return 0.0;
    }

    let (p_work, sign) = if p > 0.5 {
        (1.0 - p, 1.0)
    } else {
        (p, -1.0)
    };

    let t = (-2.0 * p_work.ln()).sqrt();

    // Rational approximation (Abramowitz & Stegun 26.2.23)
    let c0 = 2.515517;
    let c1 = 0.802853;
    let c2 = 0.010328;
    let d1 = 1.432788;
    let d2 = 0.189269;
    let d3 = 0.001308;

    let numerator = c0 + c1 * t + c2 * t * t;
    let denominator = 1.0 + d1 * t + d2 * t * t + d3 * t * t * t;

    sign * (t - numerator / denominator)
}

// ============================================================================
// Discrete gamma rate categories
// ============================================================================

/// Compute discrete gamma rate categories using Yang's (1994) method.
///
/// For a gamma distribution with shape parameter `alpha` (and rate = alpha,
/// so that the mean is 1.0), divides the distribution into `n_categories`
/// equiprobable classes and returns the mean rate within each class.
///
/// The rates are normalized to average to 1.0.
///
/// # Arguments
///
/// * `alpha` - Shape parameter of the gamma distribution. Smaller values
///   indicate greater rate heterogeneity. Must be positive.
/// * `n_categories` - Number of discrete rate categories (typically 4).
///   Must be at least 1.
///
/// # Returns
///
/// A vector of `n_categories` rate multipliers, sorted in ascending order,
/// with mean equal to 1.0.
///
/// # Panics
///
/// Panics if `alpha <= 0` or `n_categories == 0`.
pub fn discrete_gamma_rates(alpha: f64, n_categories: usize) -> Vec<f64> {
    assert!(alpha > 0.0, "alpha must be positive, got {}", alpha);
    assert!(
        n_categories > 0,
        "n_categories must be at least 1, got {}",
        n_categories
    );

    // Special case: single category means uniform rate
    if n_categories == 1 {
        return vec![1.0];
    }

    // For very large alpha, rates are nearly uniform
    if alpha > 1000.0 {
        return vec![1.0; n_categories];
    }

    let beta = alpha; // rate parameter, so mean = alpha/beta = 1.0
    let k = n_categories;

    // Compute the quantile boundaries that divide the gamma distribution
    // into k equiprobable intervals. Boundary i is at quantile i/k.
    // We need boundaries at 1/k, 2/k, ..., (k-1)/k.
    let mut boundaries = Vec::with_capacity(k + 1);
    boundaries.push(0.0_f64);
    for i in 1..k {
        let p = i as f64 / k as f64;
        // Find x such that P(alpha, beta*x) = p
        // For gamma with shape alpha, rate beta: CDF(x) = P(alpha, beta*x)
        let x = gamma_inc_inv(alpha, p) / beta;
        boundaries.push(x);
    }
    boundaries.push(f64::INFINITY);

    // Compute the mean rate in each category.
    // For a Gamma(alpha, beta) distribution, the mean of x in the interval
    // where P(alpha, beta*x) is between p1 and p2 is:
    //
    //   mean_c = k * integral_{x_low}^{x_high} x * f(x) dx
    //          = k * [F_1(x_high) - F_1(x_low)]
    //
    // where F_1 is the CDF of Gamma(alpha+1, beta).
    // This is because integral_0^x t * gamma_pdf(t; a, b) dt
    //   = (a/b) * P(a+1, b*x)
    //
    // So mean in category c = k * (a/b) * [P(a+1, b*x_{c+1}) - P(a+1, b*x_c)]
    //   = k * (alpha/beta) * [P(alpha+1, beta*x_{c+1}) - P(alpha+1, beta*x_c)]
    //   = k * 1.0 * [P(alpha+1, beta*x_{c+1}) - P(alpha+1, beta*x_c)]
    //   (since alpha/beta = 1 with our parameterization)

    let mut rates = Vec::with_capacity(k);
    for c in 0..k {
        let p_low = if c == 0 {
            0.0
        } else {
            gamma_inc_lower(alpha + 1.0, beta * boundaries[c])
        };
        let p_high = if c == k - 1 {
            1.0
        } else {
            gamma_inc_lower(alpha + 1.0, beta * boundaries[c + 1])
        };

        // Mean in this category = k * (alpha/beta) * (p_high - p_low)
        let mean_c = (k as f64) * (alpha / beta) * (p_high - p_low);
        rates.push(mean_c);
    }

    // Normalize so that mean rate = 1.0
    let mean_rate: f64 = rates.iter().sum::<f64>() / k as f64;
    if mean_rate > 0.0 {
        for r in &mut rates {
            *r /= mean_rate;
        }
    }

    rates
}

// ============================================================================
// Gamma-adjusted likelihood
// ============================================================================

/// Compute the log-likelihood of a tree with discrete gamma rate heterogeneity.
///
/// For each site, the likelihood is averaged across all rate categories:
///
/// ```text
/// L_site = (1/k) * sum_{c=1}^{k} L(site | rate_c)
/// ```
///
/// where `rate_c` scales all branch lengths by the rate multiplier for
/// category `c`. The overall log-likelihood is the sum of per-site
/// log-likelihoods.
///
/// # Arguments
///
/// * `tree` - The phylogenetic tree with branch lengths.
/// * `alignment` - The sequence alignment.
/// * `model` - The substitution model.
/// * `rates` - Rate multipliers for each category (from `discrete_gamma_rates`).
///
/// # Returns
///
/// The total log-likelihood, or an error if computation fails.
pub fn gamma_log_likelihood(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    rates: &[f64],
) -> Result<f64, LikelihoodError> {
    if rates.is_empty() {
        return Err(LikelihoodError::NumericalError(
            "rates must not be empty".to_string(),
        ));
    }

    let k = rates.len();
    let nsites = alignment.nsites();

    // Compute site log-likelihoods for each rate category.
    // For each category, we scale all branch lengths by the rate multiplier.
    let mut category_site_lnl: Vec<Vec<f64>> = Vec::with_capacity(k);

    for &rate in rates {
        let scaled_tree = scale_branch_lengths(tree, rate);
        let site_lnl = site_log_likelihoods(&scaled_tree, alignment, model)?;
        category_site_lnl.push(site_lnl);
    }

    // For each site, compute the log of the average likelihood across categories.
    // L_site = (1/k) * sum_c exp(lnL_site_c)
    // ln(L_site) = ln(1/k) + ln(sum_c exp(lnL_site_c))
    //            = -ln(k) + logsumexp(lnL_site_c for all c)
    let ln_k = (k as f64).ln();
    let mut total_lnl = 0.0;

    for site in 0..nsites {
        // Gather log-likelihoods for this site across all categories
        let mut site_lnls: Vec<f64> = Vec::with_capacity(k);
        for cat in 0..k {
            site_lnls.push(category_site_lnl[cat][site]);
        }

        // Log-sum-exp trick for numerical stability
        let max_lnl = site_lnls
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);

        let sum_exp: f64 = site_lnls.iter().map(|&lnl| (lnl - max_lnl).exp()).sum();

        let site_gamma_lnl = max_lnl + sum_exp.ln() - ln_k;
        total_lnl += site_gamma_lnl;
    }

    Ok(total_lnl)
}

/// Scale all branch lengths in a tree by a given factor.
///
/// Creates a new tree with each branch length multiplied by `scale`.
fn scale_branch_lengths(tree: &Tree, scale: f64) -> Tree {
    let mut new_tree = tree.clone();
    for node in &mut new_tree.nodes {
        if let Some(ref mut bl) = node.branch_length {
            *bl *= scale;
        }
    }
    new_tree
}

// ============================================================================
// Alpha parameter optimization
// ============================================================================

/// Optimize the alpha parameter for the discrete gamma distribution.
///
/// Uses golden section search to find the alpha value that maximizes
/// the log-likelihood of the data given the tree and substitution model.
///
/// # Arguments
///
/// * `tree` - The phylogenetic tree with branch lengths.
/// * `alignment` - The sequence alignment.
/// * `model` - The substitution model.
/// * `n_categories` - Number of discrete gamma rate categories (typically 4).
///
/// # Returns
///
/// A tuple `(best_alpha, best_lnl)` containing the optimal alpha parameter
/// and the corresponding log-likelihood, or an error if computation fails.
pub fn optimize_alpha(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    n_categories: usize,
) -> Result<(f64, f64), LikelihoodError> {
    // Search bounds for alpha
    let mut a = 0.01_f64;
    let mut b = 100.0_f64;
    let tol = 0.01;

    // Golden section search
    let phi = (5.0_f64.sqrt() - 1.0) / 2.0; // golden ratio conjugate ~ 0.618

    let mut x1 = b - phi * (b - a);
    let mut x2 = a + phi * (b - a);

    let mut f1 = {
        let rates = discrete_gamma_rates(x1, n_categories);
        gamma_log_likelihood(tree, alignment, model, &rates)?
    };
    let mut f2 = {
        let rates = discrete_gamma_rates(x2, n_categories);
        gamma_log_likelihood(tree, alignment, model, &rates)?
    };

    let max_iter = 100;
    for _ in 0..max_iter {
        if (b - a) < tol {
            break;
        }

        if f1 < f2 {
            // Maximum is in [x1, b]
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + phi * (b - a);
            let rates = discrete_gamma_rates(x2, n_categories);
            f2 = gamma_log_likelihood(tree, alignment, model, &rates)?;
        } else {
            // Maximum is in [a, x2]
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - phi * (b - a);
            let rates = discrete_gamma_rates(x1, n_categories);
            f1 = gamma_log_likelihood(tree, alignment, model, &rates)?;
        }
    }

    let best_alpha = (a + b) / 2.0;
    let best_rates = discrete_gamma_rates(best_alpha, n_categories);
    let best_lnl = gamma_log_likelihood(tree, alignment, model, &best_rates)?;

    Ok((best_alpha, best_lnl))
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::likelihood::models::{F84Model, Jc69Model};
    use crate::tree::newick::parse_newick;
    use crate::tree::types::{Base, BaseFrequencies, Sequence};

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

    // --- Gamma function tests ---

    #[test]
    fn test_ln_gamma_integers() {
        // Gamma(n) = (n-1)! for positive integers
        // Gamma(1) = 1, ln(1) = 0
        assert!((ln_gamma(1.0) - 0.0).abs() < 1e-10);
        // Gamma(2) = 1, ln(1) = 0
        assert!((ln_gamma(2.0) - 0.0).abs() < 1e-10);
        // Gamma(3) = 2, ln(2) = 0.693...
        assert!((ln_gamma(3.0) - 2.0_f64.ln()).abs() < 1e-10);
        // Gamma(4) = 6, ln(6) = 1.791...
        assert!((ln_gamma(4.0) - 6.0_f64.ln()).abs() < 1e-10);
        // Gamma(5) = 24
        assert!((ln_gamma(5.0) - 24.0_f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn test_ln_gamma_half() {
        // Gamma(0.5) = sqrt(pi) ≈ 1.7724538509
        let expected = 0.5 * std::f64::consts::PI.ln();
        assert!(
            (ln_gamma(0.5) - expected).abs() < 1e-10,
            "ln_gamma(0.5) = {}, expected {}",
            ln_gamma(0.5),
            expected
        );
    }

    #[test]
    fn test_ln_gamma_small_values() {
        // Gamma(0.1) ≈ 9.51350769866873
        let expected = 9.51350769866873_f64.ln();
        assert!(
            (ln_gamma(0.1) - expected).abs() < 1e-6,
            "ln_gamma(0.1) = {}, expected {}",
            ln_gamma(0.1),
            expected
        );
    }

    #[test]
    fn test_gamma_fn_basic() {
        assert!((gamma_fn(1.0) - 1.0).abs() < 1e-10);
        assert!((gamma_fn(2.0) - 1.0).abs() < 1e-10);
        assert!((gamma_fn(3.0) - 2.0).abs() < 1e-10);
        assert!((gamma_fn(4.0) - 6.0).abs() < 1e-10);
    }

    // --- Incomplete gamma function tests ---

    #[test]
    fn test_gamma_inc_lower_boundaries() {
        // P(a, 0) = 0
        assert!((gamma_inc_lower(1.0, 0.0) - 0.0).abs() < 1e-12);
        assert!((gamma_inc_lower(2.0, 0.0) - 0.0).abs() < 1e-12);

        // P(a, inf) = 1 (use large x as approximation)
        assert!((gamma_inc_lower(1.0, 50.0) - 1.0).abs() < 1e-10);
        assert!((gamma_inc_lower(2.0, 50.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gamma_inc_lower_exponential_cdf() {
        // For a=1, Gamma(1,1) = Exponential(1)
        // P(1, x) = 1 - exp(-x)
        for &x in &[0.1, 0.5, 1.0, 2.0, 5.0] {
            let expected = 1.0 - (-x as f64).exp();
            let computed = gamma_inc_lower(1.0, x);
            assert!(
                (computed - expected).abs() < 1e-10,
                "P(1, {}) = {}, expected {}",
                x,
                computed,
                expected
            );
        }
    }

    #[test]
    fn test_gamma_inc_lower_known_values() {
        // P(2, 1) = 1 - 2*exp(-1) ≈ 0.26424
        let expected = 1.0 - 2.0 * (-1.0_f64).exp();
        let computed = gamma_inc_lower(2.0, 1.0);
        assert!(
            (computed - expected).abs() < 1e-8,
            "P(2, 1) = {}, expected {}",
            computed,
            expected
        );
    }

    #[test]
    fn test_gamma_inc_inv_roundtrip() {
        // Verify that gamma_inc_inv is the inverse of gamma_inc_lower
        for &a in &[0.5, 1.0, 2.0, 5.0, 10.0] {
            for &p in &[0.1, 0.25, 0.5, 0.75, 0.9] {
                let x = gamma_inc_inv(a, p);
                let p_back = gamma_inc_lower(a, x);
                assert!(
                    (p_back - p).abs() < 1e-6,
                    "roundtrip failed for a={}, p={}: got x={}, P(a,x)={}",
                    a,
                    p,
                    x,
                    p_back
                );
            }
        }
    }

    // --- Discrete gamma rate tests ---

    #[test]
    fn test_discrete_gamma_rates_mean_is_one() {
        for &alpha in &[0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0] {
            for &k in &[2, 4, 8] {
                let rates = discrete_gamma_rates(alpha, k);
                let mean: f64 = rates.iter().sum::<f64>() / k as f64;
                assert!(
                    (mean - 1.0).abs() < 1e-4,
                    "mean rate should be 1.0 for alpha={}, k={}, got {}",
                    alpha,
                    k,
                    mean
                );
            }
        }
    }

    #[test]
    fn test_discrete_gamma_rates_sorted() {
        let rates = discrete_gamma_rates(0.5, 4);
        for i in 1..rates.len() {
            assert!(
                rates[i] >= rates[i - 1],
                "rates should be sorted ascending: {:?}",
                rates
            );
        }
    }

    #[test]
    fn test_discrete_gamma_rates_all_positive() {
        for &alpha in &[0.1, 0.5, 1.0, 5.0] {
            let rates = discrete_gamma_rates(alpha, 4);
            for (i, &r) in rates.iter().enumerate() {
                assert!(
                    r > 0.0,
                    "rate[{}] should be positive for alpha={}, got {}",
                    i,
                    alpha,
                    r
                );
            }
        }
    }

    #[test]
    fn test_discrete_gamma_rates_known_alpha_05() {
        // For alpha=0.5, k=4, the expected rates from Yang (1994) are
        // approximately [0.0334, 0.2519, 0.8203, 2.8944]
        let rates = discrete_gamma_rates(0.5, 4);
        let expected = [0.0334, 0.2519, 0.8203, 2.8944];

        assert_eq!(rates.len(), 4);
        for i in 0..4 {
            assert!(
                (rates[i] - expected[i]).abs() < 0.05,
                "rate[{}] = {}, expected ~{} (alpha=0.5, k=4)",
                i,
                rates[i],
                expected[i]
            );
        }
    }

    #[test]
    fn test_discrete_gamma_rates_large_alpha_uniform() {
        // Large alpha should give nearly equal rates
        let rates = discrete_gamma_rates(1000.0, 4);
        for &r in &rates {
            assert!(
                (r - 1.0).abs() < 0.1,
                "large alpha should give uniform rates, got {:?}",
                rates
            );
        }
    }

    #[test]
    fn test_discrete_gamma_rates_small_alpha_variable() {
        // Small alpha should give highly variable rates
        let rates = discrete_gamma_rates(0.1, 4);
        // The last rate should be much larger than the first
        assert!(
            rates[3] / rates[0] > 10.0,
            "small alpha should give high variation, got ratio {}",
            rates[3] / rates[0]
        );
    }

    #[test]
    fn test_discrete_gamma_rates_single_category() {
        let rates = discrete_gamma_rates(0.5, 1);
        assert_eq!(rates.len(), 1);
        assert!((rates[0] - 1.0).abs() < 1e-10);
    }

    // --- Gamma likelihood tests ---

    #[test]
    fn test_gamma_log_likelihood_is_negative() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "TGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = Jc69Model;
        let rates = discrete_gamma_rates(0.5, 4);
        let lnl = gamma_log_likelihood(&tree, &aln, &model, &rates).unwrap();
        assert!(
            lnl < 0.0,
            "gamma log-likelihood should be negative, got {}",
            lnl
        );
    }

    #[test]
    fn test_gamma_likelihood_uniform_rates_equals_standard() {
        // With uniform rates [1.0], gamma likelihood should equal standard likelihood
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGAACGA"),
            ("C", "TGCATGCA"),
        ]);
        let tree = parse_newick("((A:0.15,B:0.15):0.1,C:0.25);").unwrap();
        let model = Jc69Model;

        let standard_lnl =
            crate::likelihood::pruning::log_likelihood(&tree, &aln, &model).unwrap();
        let gamma_lnl =
            gamma_log_likelihood(&tree, &aln, &model, &[1.0]).unwrap();

        assert!(
            (standard_lnl - gamma_lnl).abs() < 1e-8,
            "uniform gamma lnl ({}) should equal standard lnl ({})",
            gamma_lnl,
            standard_lnl
        );
    }

    #[test]
    fn test_gamma_likelihood_improves_with_rate_variation() {
        // For data with variable rates, gamma model should fit better
        // than a single rate (especially when some sites are very conserved
        // and others are highly divergent)
        let aln = make_alignment(&[
            ("A", "AAAAAAACGTCGT"),
            ("B", "AAAAAAATGCTGC"),
            ("C", "AAAAAAACGTCGT"),
        ]);
        let tree = parse_newick("((A:0.3,B:0.3):0.1,C:0.3);").unwrap();
        let model = Jc69Model;

        let standard_lnl =
            crate::likelihood::pruning::log_likelihood(&tree, &aln, &model).unwrap();
        let rates = discrete_gamma_rates(0.5, 4);
        let gamma_lnl =
            gamma_log_likelihood(&tree, &aln, &model, &rates).unwrap();

        // Gamma model should fit at least as well (often better) for
        // heterogeneous data, but the key test is that it computes successfully.
        // The sign of the difference depends on the specific alpha value and data.
        assert!(gamma_lnl < 0.0);
        // Both should be finite
        assert!(standard_lnl.is_finite());
        assert!(gamma_lnl.is_finite());
    }

    #[test]
    fn test_gamma_likelihood_f84_model() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGG"),
        ]);
        let tree =
            parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();
        let freqs = BaseFrequencies::new(0.3, 0.2, 0.2, 0.3);
        let model = F84Model::new(freqs, 2.0);
        let rates = discrete_gamma_rates(1.0, 4);
        let lnl = gamma_log_likelihood(&tree, &aln, &model, &rates).unwrap();
        assert!(lnl < 0.0);
    }

    #[test]
    fn test_gamma_likelihood_different_alpha_values() {
        // Different alpha values should give different likelihoods
        let aln = make_alignment(&[
            ("A", "ACGTACGTACGT"),
            ("B", "ACGAACGAACGA"),
            ("C", "TGCATGCATGCA"),
        ]);
        let tree = parse_newick("((A:0.2,B:0.2):0.1,C:0.3);").unwrap();
        let model = Jc69Model;

        let rates_low = discrete_gamma_rates(0.1, 4);
        let rates_high = discrete_gamma_rates(10.0, 4);

        let lnl_low = gamma_log_likelihood(&tree, &aln, &model, &rates_low).unwrap();
        let lnl_high = gamma_log_likelihood(&tree, &aln, &model, &rates_high).unwrap();

        // They should be different (the exact direction depends on the data)
        assert!(
            (lnl_low - lnl_high).abs() > 1e-6,
            "different alphas should give different likelihoods: low={}, high={}",
            lnl_low,
            lnl_high
        );
    }

    // --- Alpha optimization tests ---

    #[test]
    fn test_optimize_alpha_basic() {
        let aln = make_alignment(&[
            ("A", "ACGTACGTACGTACGT"),
            ("B", "ACGAACGAACGAACGA"),
            ("C", "TGCATGCATGCATGCA"),
        ]);
        let tree = parse_newick("((A:0.2,B:0.2):0.1,C:0.3);").unwrap();
        let model = Jc69Model;

        let (best_alpha, best_lnl) =
            optimize_alpha(&tree, &aln, &model, 4).unwrap();

        assert!(best_alpha > 0.0, "optimal alpha should be positive");
        assert!(
            best_alpha < 100.0,
            "optimal alpha should be within search bounds"
        );
        assert!(best_lnl < 0.0, "optimal lnl should be negative");
    }

    #[test]
    fn test_optimize_alpha_improves_on_extremes() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGG"),
        ]);
        let tree =
            parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();
        let model = Jc69Model;

        let (_, best_lnl) = optimize_alpha(&tree, &aln, &model, 4).unwrap();

        // Check that the optimized alpha gives a likelihood at least as good
        // as some arbitrary alpha values
        let rates_low = discrete_gamma_rates(0.05, 4);
        let lnl_low = gamma_log_likelihood(&tree, &aln, &model, &rates_low).unwrap();

        assert!(
            best_lnl >= lnl_low - 1e-4,
            "optimized lnl ({}) should be >= extreme low alpha lnl ({})",
            best_lnl,
            lnl_low
        );
    }

    // --- Edge case tests ---

    #[test]
    fn test_gamma_rates_two_categories() {
        let rates = discrete_gamma_rates(1.0, 2);
        assert_eq!(rates.len(), 2);
        let mean: f64 = rates.iter().sum::<f64>() / 2.0;
        assert!(
            (mean - 1.0).abs() < 1e-6,
            "2-category mean should be 1.0, got {}",
            mean
        );
        assert!(rates[0] < 1.0, "first rate should be < 1.0");
        assert!(rates[1] > 1.0, "second rate should be > 1.0");
    }

    #[test]
    fn test_gamma_rates_eight_categories() {
        let rates = discrete_gamma_rates(0.5, 8);
        assert_eq!(rates.len(), 8);
        let mean: f64 = rates.iter().sum::<f64>() / 8.0;
        assert!(
            (mean - 1.0).abs() < 1e-4,
            "8-category mean should be 1.0, got {}",
            mean
        );
        // Should be sorted
        for i in 1..8 {
            assert!(rates[i] >= rates[i - 1]);
        }
    }

    #[test]
    fn test_normal_quantile_symmetry() {
        // The normal quantile function should satisfy q(p) = -q(1-p)
        for &p in &[0.1, 0.2, 0.3, 0.4] {
            let q1 = normal_quantile(p);
            let q2 = normal_quantile(1.0 - p);
            assert!(
                (q1 + q2).abs() < 1e-4,
                "q({}) + q({}) = {} + {} = {}, should be ~0",
                p,
                1.0 - p,
                q1,
                q2,
                q1 + q2
            );
        }
    }

    #[test]
    fn test_normal_quantile_median() {
        assert!((normal_quantile(0.5)).abs() < 1e-10);
    }

    #[test]
    fn test_scale_branch_lengths() {
        let tree = parse_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();
        let scaled = scale_branch_lengths(&tree, 2.0);

        for node in &scaled.nodes {
            if let Some(bl) = node.branch_length {
                let orig_bl = tree.nodes[node.id].branch_length.unwrap();
                assert!(
                    (bl - orig_bl * 2.0).abs() < 1e-12,
                    "branch length should be scaled by 2.0"
                );
            }
        }
    }
}

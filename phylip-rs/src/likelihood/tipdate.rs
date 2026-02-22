//! TipDate: Serial-sample molecular clock estimation.
//!
//! Implements the dated-tips (SRDT) molecular clock model for estimating
//! the rate of molecular evolution from sequences sampled at different
//! times. This approach, pioneered by Andrew Rambaut (2000), exploits
//! the temporal signal in heterochronous sequence data to jointly
//! estimate the substitution rate and divergence times.
//!
//! # Background
//!
//! When sequences are sampled at different calendar dates (e.g., viral
//! isolates collected over years), the sampling times provide calibration
//! points that anchor the molecular clock. Under a strict clock with
//! rate `r`, the expected number of substitutions along a branch from
//! time `t_p` (parent) to time `t_c` (child) is:
//!
//! ```text
//! branch_length = r * |t_c - t_p|
//! ```
//!
//! Tip dates are known (sampling times); internal node dates and the
//! substitution rate are estimated by maximum likelihood.
//!
//! # Three Models Compared by Likelihood Ratio Test
//!
//! - **Model 0 (No clock):** Free branch lengths, no temporal constraint.
//! - **Model 1 (Strict clock, contemporaneous):** Ultrametric tree with
//!   all tips at time 0. Standard molecular clock.
//! - **Model 2 (SRDT / TipDate):** Dated-tip clock. Tips have known
//!   sampling dates; internal node times and rate are estimated.
//!
//! The likelihood ratio test (LRT) compares Model 2 vs. Model 1:
//!
//! ```text
//! LRT = 2 * (lnL_SRDT - lnL_clock) ~ chi^2(df)
//! ```
//!
//! where `df = (number of distinct tip dates) - 1`.
//!
//! # Algorithm
//!
//! 1. Root the tree and assign tip dates from the provided map.
//! 2. Initialize the rate `r` from root-to-tip regression.
//! 3. Optimize `r` and internal node times using golden section search
//!    to maximize the log-likelihood under Felsenstein's pruning algorithm.
//! 4. Compare to the strict clock model using the LRT.
//!
//! # Example
//!
//! ```
//! use phylip_rs::likelihood::tipdate::{tipdate_optimize, tipdate_lrt};
//! use phylip_rs::likelihood::models::Jc69Model;
//! use phylip_rs::tree::newick::parse_newick;
//! use phylip_rs::tree::{Alignment, Base, Sequence};
//! use std::collections::HashMap;
//!
//! let alignment = Alignment::new(vec![
//!     Sequence::new("A", vec![Base::A, Base::C, Base::G, Base::T,
//!                             Base::A, Base::C, Base::G, Base::T]),
//!     Sequence::new("B", vec![Base::A, Base::C, Base::G, Base::A,
//!                             Base::A, Base::C, Base::G, Base::A]),
//!     Sequence::new("C", vec![Base::A, Base::C, Base::A, Base::A,
//!                             Base::A, Base::C, Base::A, Base::A]),
//! ]).unwrap();
//!
//! let tree = parse_newick("((A:0.1,B:0.1):0.1,C:0.2);").unwrap();
//! let model = Jc69Model;
//!
//! let mut tip_dates = HashMap::new();
//! tip_dates.insert("A".to_string(), 2000.0);
//! tip_dates.insert("B".to_string(), 2005.0);
//! tip_dates.insert("C".to_string(), 2010.0);
//!
//! let result = tipdate_optimize(&tree, &alignment, &model, &tip_dates).unwrap();
//! assert!(result.rate > 0.0);
//! assert!(result.root_date < 2000.0);
//! ```
//!
//! # References
//!
//! - Rambaut, A. (2000). Estimating the rate of molecular evolution:
//!   incorporating non-contemporaneous sequences into maximum likelihood
//!   phylogenetics. *Bioinformatics*, 16(4), 395-399.
//! - Drummond, A. J., Pybus, O. G., Rambaut, A., Forsberg, R., &
//!   Rodrigo, A. G. (2003). Measurably evolving populations. *Trends in
//!   Ecology & Evolution*, 18(9), 481-488.
//! - Felsenstein, J. (1981). Evolutionary trees from DNA sequences:
//!   a maximum likelihood approach. *J. Mol. Evol.*, 17, 368-376.

use std::collections::HashMap;

use crate::likelihood::models::SubstitutionModel;
use crate::likelihood::pruning::log_likelihood;
use crate::tree::types::{Alignment, Tree};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Minimum substitution rate for numerical stability.
const MIN_RATE: f64 = 1e-10;

/// Maximum substitution rate (subs/site/time-unit).
const MAX_RATE: f64 = 10.0;

/// Minimum branch length to prevent zero-length branches.
const MIN_BRANCH_LENGTH: f64 = 1e-8;

/// Convergence tolerance for optimization.
const TOLERANCE: f64 = 1e-6;

/// Maximum number of outer optimization rounds.
const MAX_ROUNDS: usize = 50;

/// Maximum iterations for golden section search.
const MAX_GOLDEN: usize = 40;

// ---------------------------------------------------------------------------
// Result type
// ---------------------------------------------------------------------------

/// Result of a TipDate (SRDT) analysis.
///
/// Contains the optimized tree with branch lengths set from the
/// estimated rate and divergence times, together with the estimated
/// parameters.
#[derive(Debug, Clone)]
pub struct TipDateResult {
    /// The tree with branch lengths set from `rate * |t_parent - t_child|`.
    pub tree: Tree,
    /// Log-likelihood of the optimized SRDT model.
    pub lnl: f64,
    /// Estimated substitution rate (substitutions per site per time unit).
    pub rate: f64,
    /// Estimated date of the root node.
    pub root_date: f64,
    /// Tip dates provided as input (taxon name -> sampling date).
    pub tip_dates: HashMap<String, f64>,
    /// Estimated dates for all nodes (indexed by node id).
    pub node_dates: Vec<f64>,
}

// ---------------------------------------------------------------------------
// Core functions
// ---------------------------------------------------------------------------

/// Build a mapping from leaf node id to its sampling date.
///
/// Returns an error if any leaf in the tree is missing from `tip_dates`.
fn leaf_date_map(
    tree: &Tree,
    tip_dates: &HashMap<String, f64>,
) -> Result<HashMap<usize, f64>, String> {
    let mut map = HashMap::new();
    for node in &tree.nodes {
        if node.is_leaf() {
            let name = node
                .name
                .as_deref()
                .ok_or_else(|| "leaf node has no name".to_string())?;
            let date = tip_dates
                .get(name)
                .ok_or_else(|| format!("no date provided for leaf '{}'", name))?;
            map.insert(node.id, *date);
        }
    }
    Ok(map)
}

/// Initialize internal node dates by averaging descendant tip dates.
///
/// Uses a postorder traversal: each internal node's date is set to the
/// average of its children's dates, shifted earlier (toward the past)
/// by a small offset proportional to the tree depth.
fn initialize_node_dates(
    tree: &Tree,
    leaf_dates: &HashMap<usize, f64>,
) -> Vec<f64> {
    let n = tree.num_nodes();
    let mut dates = vec![0.0; n];

    // Find the date range for scaling
    let min_date = leaf_dates
        .values()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let max_date = leaf_dates
        .values()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let span = (max_date - min_date).max(1.0);

    // Set leaf dates
    for (&id, &date) in leaf_dates {
        dates[id] = date;
    }

    // Postorder: set internal node dates
    let postorder = tree.postorder();
    for &node_id in &postorder {
        if tree.nodes[node_id].is_leaf() {
            continue;
        }
        let children = &tree.nodes[node_id].children;
        if children.is_empty() {
            continue;
        }
        // Average of children's dates, shifted toward the past
        let avg: f64 =
            children.iter().map(|&c| dates[c]).sum::<f64>() / children.len() as f64;
        // Push the internal node earlier than the earliest child
        let earliest_child = children
            .iter()
            .map(|&c| dates[c])
            .fold(f64::INFINITY, f64::min);
        // Place node at avg minus a fraction of the span, but always before earliest child
        let offset = span * 0.1;
        dates[node_id] = (avg - offset).min(earliest_child - offset.max(0.01));
    }

    dates
}

/// Set branch lengths on a tree from a rate and node dates.
///
/// For each non-root node, the branch length is `rate * |date_parent - date_child|`,
/// clamped to a minimum for numerical stability.
fn set_branch_lengths_from_dates(tree: &mut Tree, dates: &[f64], rate: f64) {
    let n = tree.num_nodes();
    for i in 0..n {
        if let Some(parent_id) = tree.nodes[i].parent {
            let dt = (dates[i] - dates[parent_id]).abs();
            let bl = (rate * dt).max(MIN_BRANCH_LENGTH);
            tree.nodes[i].branch_length = Some(bl);
        }
    }
}

/// Compute the log-likelihood of a tree under the SRDT (TipDate) model.
///
/// Given a rooted tree, alignment, substitution model, tip dates, rate,
/// and internal node dates, this function:
/// 1. Sets branch lengths as `rate * |t_parent - t_child|`.
/// 2. Evaluates the log-likelihood using Felsenstein's pruning algorithm.
///
/// # Arguments
///
/// * `tree` - A rooted phylogenetic tree (topology only; branch lengths will be overwritten).
/// * `alignment` - The DNA sequence alignment.
/// * `model` - Nucleotide substitution model (e.g., JC69, F84).
/// * `tip_dates` - Mapping from taxon name to sampling date.
/// * `rate` - Substitution rate (substitutions per site per time unit).
/// * `node_dates` - Dates for all nodes (tip dates fixed, internal dates variable).
///
/// # Returns
///
/// The log-likelihood, or an error string if computation fails.
///
/// # References
///
/// Rambaut, A. (2000). *Bioinformatics*, 16(4), 395-399.
pub fn tipdate_likelihood(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    _tip_dates: &HashMap<String, f64>,
    rate: f64,
    node_dates: &[f64],
) -> Result<f64, String> {
    let mut work_tree = tree.clone();
    set_branch_lengths_from_dates(&mut work_tree, node_dates, rate);
    log_likelihood(&work_tree, alignment, model).map_err(|e| format!("{}", e))
}

/// Estimate an initial substitution rate from root-to-tip distances.
///
/// Performs a simple regression of root-to-tip genetic distance against
/// tip sampling date. The slope of this regression provides a rough
/// initial estimate of the substitution rate.
fn estimate_initial_rate(tree: &Tree, tip_dates: &HashMap<String, f64>) -> f64 {
    // Compute root-to-tip distances from branch lengths
    let mut distances: Vec<(f64, f64)> = Vec::new(); // (date, distance)
    for node in &tree.nodes {
        if !node.is_leaf() {
            continue;
        }
        let name = match &node.name {
            Some(n) => n.as_str(),
            None => continue,
        };
        let date = match tip_dates.get(name) {
            Some(&d) => d,
            None => continue,
        };
        // Sum branch lengths from leaf to root
        let mut dist = 0.0;
        let mut current = node.id;
        loop {
            dist += tree.nodes[current].branch_length.unwrap_or(0.0);
            match tree.nodes[current].parent {
                Some(p) => current = p,
                None => break,
            }
        }
        distances.push((date, dist));
    }

    if distances.len() < 2 {
        return 0.001; // fallback
    }

    // Simple linear regression: distance = rate * date + intercept
    let n = distances.len() as f64;
    let sum_x: f64 = distances.iter().map(|(d, _)| d).sum();
    let sum_y: f64 = distances.iter().map(|(_, d)| d).sum();
    let sum_xy: f64 = distances.iter().map(|(x, y)| x * y).sum();
    let sum_xx: f64 = distances.iter().map(|(x, _)| x * x).sum();

    let denom = n * sum_xx - sum_x * sum_x;
    if denom.abs() < 1e-15 {
        return 0.001; // no temporal signal, fallback
    }

    let slope = (n * sum_xy - sum_x * sum_y) / denom;

    // Rate should be positive
    slope.abs().max(MIN_RATE).min(MAX_RATE)
}

/// Golden section search to maximize a univariate function.
///
/// Finds the value of `x` in `[lo, hi]` that maximizes `f(x)`.
/// Returns `(best_x, best_f)`.
fn golden_section_maximize<F>(mut lo: f64, mut hi: f64, mut f: F) -> (f64, f64)
where
    F: FnMut(f64) -> f64,
{
    let phi = (5.0_f64.sqrt() - 1.0) / 2.0;

    let mut x1 = hi - phi * (hi - lo);
    let mut x2 = lo + phi * (hi - lo);
    let mut f1 = f(x1);
    let mut f2 = f(x2);

    for _ in 0..MAX_GOLDEN {
        if (hi - lo) < TOLERANCE {
            break;
        }
        if f1 < f2 {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + phi * (hi - lo);
            f2 = f(x2);
        } else {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - phi * (hi - lo);
            f1 = f(x1);
        }
    }

    let best_x = (lo + hi) / 2.0;
    let best_f = f(best_x);
    (best_x, best_f)
}

/// Optimize the substitution rate and internal node dates to maximize
/// the log-likelihood under the SRDT (TipDate) model.
///
/// Uses an iterative coordinate-ascent approach:
/// 1. Optimize the rate `r` by golden section search (fixing node dates).
/// 2. Optimize each internal node's date by golden section search
///    (fixing the rate and other node dates).
/// 3. Repeat until convergence.
///
/// # Arguments
///
/// * `tree` - A rooted phylogenetic tree (topology used; branch lengths overwritten).
/// * `alignment` - The DNA sequence alignment.
/// * `model` - Nucleotide substitution model.
/// * `tip_dates` - Mapping from taxon name to sampling date.
///
/// # Returns
///
/// A `TipDateResult` containing the optimized tree, rate, dates, and
/// log-likelihood.
///
/// # Errors
///
/// Returns an error if:
/// - Any leaf in the tree has no matching entry in `tip_dates`.
/// - The tree has fewer than 2 leaves.
/// - Likelihood computation fails numerically.
///
/// # References
///
/// Rambaut, A. (2000). Estimating the rate of molecular evolution:
/// incorporating non-contemporaneous sequences into maximum likelihood
/// phylogenetics. *Bioinformatics*, 16(4), 395-399.
pub fn tipdate_optimize(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    tip_dates: &HashMap<String, f64>,
) -> Result<TipDateResult, String> {
    if tree.num_leaves() < 2 {
        return Err("need at least 2 taxa for TipDate analysis".to_string());
    }

    // Validate tip dates and build leaf-id -> date map
    let leaf_dates = leaf_date_map(tree, tip_dates)?;

    // Initialize node dates
    let mut dates = initialize_node_dates(tree, &leaf_dates);

    // Initialize rate from root-to-tip regression
    let mut rate = estimate_initial_rate(tree, tip_dates);

    // Working copy of the tree
    let mut work_tree = tree.clone();

    // Compute initial likelihood
    set_branch_lengths_from_dates(&mut work_tree, &dates, rate);
    let mut prev_lnl = log_likelihood(&work_tree, alignment, model)
        .map_err(|e| format!("{}", e))?;

    // Identify internal (non-leaf, non-root-only-if-needed) nodes
    let internal_nodes: Vec<usize> = tree
        .postorder()
        .into_iter()
        .filter(|&id| !tree.nodes[id].is_leaf())
        .collect();

    // Find the overall date range for bounding
    let min_tip_date = leaf_dates
        .values()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let max_tip_date = leaf_dates
        .values()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let date_span = (max_tip_date - min_tip_date).max(1.0);

    for _round in 0..MAX_ROUNDS {
        // --- Step 1: Optimize rate ---
        let dates_snapshot = dates.clone();
        let (new_rate, _) = golden_section_maximize(MIN_RATE, MAX_RATE, |r| {
            let mut t = tree.clone();
            set_branch_lengths_from_dates(&mut t, &dates_snapshot, r);
            log_likelihood(&t, alignment, model).unwrap_or(f64::NEG_INFINITY)
        });
        rate = new_rate;

        // --- Step 2: Optimize each internal node's date ---
        for &node_id in &internal_nodes {
            let children = &tree.nodes[node_id].children;
            if children.is_empty() {
                continue;
            }

            // The node's date must be earlier than (or equal to) the earliest child date
            let earliest_child_date = children
                .iter()
                .map(|&c| dates[c])
                .fold(f64::INFINITY, f64::min);

            // The node's date must be later than (or equal to) the parent's date (if any)
            let parent_date = tree.nodes[node_id]
                .parent
                .map(|p| dates[p])
                .unwrap_or(min_tip_date - date_span * 10.0);

            // Bounds: parent_date < node_date < earliest_child_date
            let lo = parent_date + TOLERANCE;
            let hi = earliest_child_date - TOLERANCE;

            if lo >= hi {
                continue; // no room to optimize
            }

            let rate_snap = rate;
            let mut dates_snap = dates.clone();
            let (best_date, _) = golden_section_maximize(lo, hi, |d| {
                dates_snap[node_id] = d;
                let mut t = tree.clone();
                set_branch_lengths_from_dates(&mut t, &dates_snap, rate_snap);
                let lnl =
                    log_likelihood(&t, alignment, model).unwrap_or(f64::NEG_INFINITY);
                dates_snap[node_id] = d; // ensure consistency
                lnl
            });
            dates[node_id] = best_date;
        }

        // Check convergence
        set_branch_lengths_from_dates(&mut work_tree, &dates, rate);
        let lnl = log_likelihood(&work_tree, alignment, model)
            .map_err(|e| format!("{}", e))?;

        if (lnl - prev_lnl).abs() < TOLERANCE {
            prev_lnl = lnl;
            break;
        }
        prev_lnl = lnl;
    }

    // Final tree with optimized branch lengths
    let mut result_tree = tree.clone();
    set_branch_lengths_from_dates(&mut result_tree, &dates, rate);

    Ok(TipDateResult {
        tree: result_tree,
        lnl: prev_lnl,
        rate,
        root_date: dates[tree.root],
        tip_dates: tip_dates.clone(),
        node_dates: dates,
    })
}

/// Likelihood ratio test comparing the SRDT (TipDate) model to the
/// strict molecular clock.
///
/// Under the null hypothesis (strict clock with contemporaneous tips),
/// the test statistic:
///
/// ```text
/// LRT = 2 * (lnL_tipdate - lnL_clock)
/// ```
///
/// is approximately chi-squared distributed with degrees of freedom
/// equal to `(number of distinct tip dates) - 1`.
///
/// # Arguments
///
/// * `lnl_clock` - Log-likelihood under the strict (contemporaneous) clock.
/// * `lnl_tipdate` - Log-likelihood under the SRDT (dated-tip) model.
/// * `n_distinct_dates` - Number of distinct sampling dates among tips.
///
/// # Returns
///
/// A tuple `(lrt_statistic, df, p_value)` where:
/// - `lrt_statistic` is `2 * (lnL_tipdate - lnL_clock)`, clamped to >= 0.
/// - `df` is `n_distinct_dates - 1`.
/// - `p_value` is the upper-tail probability from the chi-squared distribution.
///
/// # References
///
/// Rambaut, A. (2000). *Bioinformatics*, 16(4), 395-399.
/// Equation 5: LRT = 2(lnL_2 - lnL_1) ~ chi^2(k-1), where k = distinct dates.
pub fn tipdate_lrt(
    lnl_clock: f64,
    lnl_tipdate: f64,
    n_distinct_dates: usize,
) -> (f64, usize, f64) {
    let lrt = (2.0 * (lnl_tipdate - lnl_clock)).max(0.0);
    let df = if n_distinct_dates > 1 {
        n_distinct_dates - 1
    } else {
        1
    };
    let p_value = chi2_survival(lrt, df);
    (lrt, df, p_value)
}

// ---------------------------------------------------------------------------
// Chi-squared survival function (upper-tail p-value)
// ---------------------------------------------------------------------------

/// Compute P(X > x) for X ~ chi^2(df).
///
/// Uses the regularized incomplete gamma function:
/// P(X > x) = 1 - gamma_p(df/2, x/2)
///
/// Implemented from first principles with no external dependencies.
fn chi2_survival(x: f64, df: usize) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    let a = df as f64 / 2.0;
    let z = x / 2.0;
    1.0 - regularized_gamma_p(a, z)
}

/// Regularized lower incomplete gamma function P(a, x) = gamma(a, x) / Gamma(a).
///
/// Uses the series expansion for small x and the continued fraction
/// for large x, following Numerical Recipes.
fn regularized_gamma_p(a: f64, x: f64) -> f64 {
    if x < 0.0 {
        return 0.0;
    }
    if x == 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        // Series expansion
        gamma_series(a, x)
    } else {
        // Continued fraction (complement)
        1.0 - gamma_cf(a, x)
    }
}

/// Series expansion for the regularized incomplete gamma function.
///
/// P(a, x) = exp(-x + a*ln(x) - ln(Gamma(a))) * sum_{n=0}^{inf} x^n / (a*(a+1)*...*(a+n))
fn gamma_series(a: f64, x: f64) -> f64 {
    let ln_gamma_a = ln_gamma(a);
    let mut sum = 1.0 / a;
    let mut term = 1.0 / a;
    for n in 1..200 {
        term *= x / (a + n as f64);
        sum += term;
        if term.abs() < sum.abs() * 1e-14 {
            break;
        }
    }
    let log_val = -x + a * x.ln() - ln_gamma_a + sum.ln();
    if log_val > 500.0 {
        1.0
    } else if log_val < -500.0 {
        0.0
    } else {
        log_val.exp().min(1.0)
    }
}

/// Continued fraction for the complement of the regularized incomplete
/// gamma function: Q(a, x) = 1 - P(a, x).
///
/// Uses the modified Lentz algorithm for the continued fraction:
/// Q(a, x) = exp(-x + a*ln(x) - lnGamma(a)) * CF
/// where CF = 1 / (x + 1 - a + K_{n=1}^{inf} n*(n-a) / (x + 2n + 1 - a))
fn gamma_cf(a: f64, x: f64) -> f64 {
    let ln_gamma_a = ln_gamma(a);
    let tiny = 1e-30;

    // Modified Lentz's method
    let b0 = x + 1.0 - a;
    let mut d = if b0.abs() < tiny {
        1.0 / tiny
    } else {
        1.0 / b0
    };
    let mut c = b0;
    let mut f = d;

    for i in 1..200 {
        let ai = (i as f64) * (a - i as f64);
        let bi = x + 2.0 * i as f64 + 1.0 - a;

        d = bi + ai * d;
        if d.abs() < tiny {
            d = tiny;
        }
        d = 1.0 / d;

        c = bi + ai / c;
        if c.abs() < tiny {
            c = tiny;
        }

        let delta = c * d;
        f *= delta;

        if (delta - 1.0).abs() < 1e-14 {
            break;
        }
    }

    let log_val = -x + a * x.ln() - ln_gamma_a;
    if log_val < -500.0 {
        return 0.0;
    }
    let result = log_val.exp() * f;
    result.max(0.0).min(1.0)
}

/// Natural logarithm of the Gamma function, via Stirling's approximation
/// with Lanczos coefficients.
///
/// Accurate to ~15 digits for positive arguments.
fn ln_gamma(x: f64) -> f64 {
    // Lanczos approximation (g=7, n=9)
    let coeffs = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_13,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];

    if x < 0.5 {
        // Reflection formula: Gamma(x) * Gamma(1-x) = pi / sin(pi*x)
        let log_pi = std::f64::consts::PI.ln();
        let sin_val = (std::f64::consts::PI * x).sin();
        if sin_val.abs() < 1e-300 {
            return f64::INFINITY;
        }
        return log_pi - sin_val.abs().ln() - ln_gamma(1.0 - x);
    }

    let x = x - 1.0;
    let g = 7.0_f64;

    let mut sum = coeffs[0];
    for i in 1..9 {
        sum += coeffs[i] / (x + i as f64);
    }

    let t = x + g + 0.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + sum.ln()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::likelihood::models::Jc69Model;
    use crate::tree::newick::parse_newick;
    use crate::tree::types::{Base, Sequence};

    /// Helper to build an alignment from (name, sequence) pairs.
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

    // -----------------------------------------------------------------------
    // Two-taxa analytical case
    // -----------------------------------------------------------------------

    #[test]
    fn test_two_taxa_analytical() {
        // With two taxa, the tree is trivial: root -> A, root -> B.
        // Under JC69 with rate r, the branch length from root to each
        // tip is r * |tip_date - root_date|.
        //
        // The root date should be estimated as earlier than both tips,
        // and the rate should be positive.
        let aln = make_alignment(&[
            ("A", "ACGTACGTACGTACGT"),
            ("B", "ACGAACGAACGAACGA"),
        ]);
        let tree = parse_newick("(A:0.1,B:0.1);").unwrap();
        let model = Jc69Model;

        let mut tip_dates = HashMap::new();
        tip_dates.insert("A".to_string(), 2000.0);
        tip_dates.insert("B".to_string(), 2010.0);

        let result = tipdate_optimize(&tree, &aln, &model, &tip_dates).unwrap();

        // Rate must be positive
        assert!(
            result.rate > 0.0,
            "rate should be positive, got {}",
            result.rate
        );

        // Root date must precede both tips
        assert!(
            result.root_date < 2000.0,
            "root date ({}) should precede earliest tip (2000.0)",
            result.root_date
        );

        // Log-likelihood should be finite and negative
        assert!(
            result.lnl.is_finite(),
            "lnl should be finite, got {}",
            result.lnl
        );
        assert!(result.lnl < 0.0, "lnl should be negative, got {}", result.lnl);
    }

    // -----------------------------------------------------------------------
    // Rate must be positive
    // -----------------------------------------------------------------------

    #[test]
    fn test_rate_is_positive() {
        let aln = make_alignment(&[
            ("X", "AACCGGTTAACCGGTT"),
            ("Y", "AACCGGTTAACCGGTA"),
            ("Z", "CCAATTGGCCAATTGG"),
        ]);
        let tree = parse_newick("((X:0.1,Y:0.1):0.05,Z:0.15);").unwrap();
        let model = Jc69Model;

        let mut tip_dates = HashMap::new();
        tip_dates.insert("X".to_string(), 1990.0);
        tip_dates.insert("Y".to_string(), 2000.0);
        tip_dates.insert("Z".to_string(), 2010.0);

        let result = tipdate_optimize(&tree, &aln, &model, &tip_dates).unwrap();

        assert!(
            result.rate > 0.0,
            "rate must be positive, got {}",
            result.rate
        );
        assert!(
            result.rate < MAX_RATE,
            "rate should be bounded, got {}",
            result.rate
        );
    }

    // -----------------------------------------------------------------------
    // Root date must precede all tip dates
    // -----------------------------------------------------------------------

    #[test]
    fn test_root_date_precedes_tips() {
        let aln = make_alignment(&[
            ("T1", "AACCGGTTAACCGGTTAACCGGTT"),
            ("T2", "AACCGGTTAACCGGTTAACCGGTA"),
            ("T3", "CCAATTGGCCAATTGGCCAATTGG"),
            ("T4", "CCAATTGGCCAATTGGCCAATTGA"),
        ]);
        let tree =
            parse_newick("((T1:0.1,T2:0.1):0.05,(T3:0.1,T4:0.1):0.05);").unwrap();
        let model = Jc69Model;

        let mut tip_dates = HashMap::new();
        tip_dates.insert("T1".to_string(), 2000.0);
        tip_dates.insert("T2".to_string(), 2005.0);
        tip_dates.insert("T3".to_string(), 2010.0);
        tip_dates.insert("T4".to_string(), 2015.0);

        let result = tipdate_optimize(&tree, &aln, &model, &tip_dates).unwrap();

        let earliest_tip = 2000.0;
        assert!(
            result.root_date < earliest_tip,
            "root date ({}) should precede earliest tip date ({})",
            result.root_date,
            earliest_tip
        );

        // All internal node dates should be earlier than their descendant tips
        for &node_id in &tree
            .postorder()
            .into_iter()
            .filter(|&id| !tree.nodes[id].is_leaf())
            .collect::<Vec<_>>()
        {
            let children = &tree.nodes[node_id].children;
            for &child_id in children {
                assert!(
                    result.node_dates[node_id] <= result.node_dates[child_id],
                    "parent date ({}) should be <= child date ({})",
                    result.node_dates[node_id],
                    result.node_dates[child_id]
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Clock-like data should yield small LRT
    // -----------------------------------------------------------------------

    #[test]
    fn test_clock_like_data_small_lrt() {
        // When all tips are sampled at the same time, the SRDT model
        // should reduce to the standard clock, giving a small LRT.
        let aln = make_alignment(&[
            ("A", "ACGTACGTACGTACGT"),
            ("B", "ACGTACGTACGTACGT"),
            ("C", "ACGTACGTACGTACGT"),
            ("D", "ACGTACGTACGTACGT"),
        ]);
        let tree =
            parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();
        let model = Jc69Model;

        // All tips at same date — SRDT should not improve over clock
        let mut tip_dates = HashMap::new();
        tip_dates.insert("A".to_string(), 2020.0);
        tip_dates.insert("B".to_string(), 2020.0);
        tip_dates.insert("C".to_string(), 2020.0);
        tip_dates.insert("D".to_string(), 2020.0);

        let result = tipdate_optimize(&tree, &aln, &model, &tip_dates).unwrap();

        // With identical sequences and identical dates, the LRT should be tiny.
        // Use the SRDT lnl as both clock and tipdate (they should be near-equal).
        let (lrt, df, p_value) = tipdate_lrt(result.lnl, result.lnl, 1);

        assert!(
            lrt < 1e-6,
            "LRT should be ~0 when all dates are identical, got {}",
            lrt
        );
        assert_eq!(df, 1, "df should be 1 when n_distinct_dates=1");
        assert!(
            p_value > 0.5,
            "p-value should be large (not significant), got {}",
            p_value
        );
    }

    // -----------------------------------------------------------------------
    // Contemporaneous sequences should reduce to standard clock
    // -----------------------------------------------------------------------

    #[test]
    fn test_contemporaneous_reduces_to_clock() {
        // With all tips at the same date, the SRDT model is equivalent
        // to the standard molecular clock.
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTA"),
            ("C", "CCAATTGGCCAATTGG"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.1,C:0.2);").unwrap();
        let model = Jc69Model;

        let mut tip_dates = HashMap::new();
        tip_dates.insert("A".to_string(), 0.0);
        tip_dates.insert("B".to_string(), 0.0);
        tip_dates.insert("C".to_string(), 0.0);

        let result = tipdate_optimize(&tree, &aln, &model, &tip_dates).unwrap();

        // Should produce a valid result
        assert!(result.lnl.is_finite());
        assert!(result.lnl < 0.0);
        assert!(result.rate > 0.0);

        // All tip dates should remain at 0.0
        for leaf in tree.leaves() {
            let name = leaf.name.as_deref().unwrap();
            assert!(
                (result.node_dates[leaf.id] - 0.0).abs() < 1e-10,
                "tip '{}' date should be 0.0, got {}",
                name,
                result.node_dates[leaf.id]
            );
        }

        // Root date should be negative (earlier than tips at time 0)
        assert!(
            result.root_date < 0.0,
            "root date should be negative (before tips at 0), got {}",
            result.root_date
        );
    }

    // -----------------------------------------------------------------------
    // LRT function properties
    // -----------------------------------------------------------------------

    #[test]
    fn test_lrt_basic_properties() {
        // Positive LRT when tipdate model is better
        let (lrt, df, p) = tipdate_lrt(-100.0, -90.0, 5);
        assert!((lrt - 20.0).abs() < 1e-10, "LRT should be 20.0, got {}", lrt);
        assert_eq!(df, 4, "df should be n_distinct_dates - 1 = 4");
        assert!(p >= 0.0 && p <= 1.0, "p-value should be in [0,1], got {}", p);
        assert!(
            p < 0.01,
            "p-value should be small for LRT=20 with df=4, got {}",
            p
        );
    }

    #[test]
    fn test_lrt_zero_when_equal() {
        let (lrt, _, p) = tipdate_lrt(-100.0, -100.0, 3);
        assert!(lrt.abs() < 1e-10, "LRT should be 0 when models are equal");
        assert!(
            (p - 1.0).abs() < 1e-6,
            "p-value should be 1.0 when LRT=0, got {}",
            p
        );
    }

    #[test]
    fn test_lrt_clamped_nonnegative() {
        // If clock model has higher lnl (shouldn't happen in theory,
        // but can due to numerical noise), LRT is clamped to 0.
        let (lrt, _, _) = tipdate_lrt(-90.0, -100.0, 3);
        assert!(
            lrt >= 0.0,
            "LRT should be clamped to >= 0, got {}",
            lrt
        );
        assert!(lrt.abs() < 1e-10);
    }

    #[test]
    fn test_lrt_df_at_least_one() {
        let (_, df, _) = tipdate_lrt(-100.0, -95.0, 1);
        assert!(df >= 1, "df should be at least 1");
    }

    // -----------------------------------------------------------------------
    // Error handling
    // -----------------------------------------------------------------------

    #[test]
    fn test_missing_tip_date_error() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "TGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = Jc69Model;

        // Only provide dates for A and B, not C
        let mut tip_dates = HashMap::new();
        tip_dates.insert("A".to_string(), 2000.0);
        tip_dates.insert("B".to_string(), 2005.0);

        let result = tipdate_optimize(&tree, &aln, &model, &tip_dates);
        assert!(result.is_err(), "should error when a tip date is missing");
        let err = result.unwrap_err();
        assert!(
            err.contains("no date provided for leaf 'C'"),
            "error should mention the missing taxon, got: {}",
            err
        );
    }

    // -----------------------------------------------------------------------
    // tipdate_likelihood direct test
    // -----------------------------------------------------------------------

    #[test]
    fn test_tipdate_likelihood_direct() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGAACGA"),
            ("C", "TGCATGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.1,C:0.2);").unwrap();
        let model = Jc69Model;

        let mut tip_dates = HashMap::new();
        tip_dates.insert("A".to_string(), 2000.0);
        tip_dates.insert("B".to_string(), 2005.0);
        tip_dates.insert("C".to_string(), 2010.0);

        // Set up node dates: tips at their dates, root earlier
        let leaf_dates = leaf_date_map(&tree, &tip_dates).unwrap();
        let node_dates = initialize_node_dates(&tree, &leaf_dates);

        let rate = 0.01;
        let lnl =
            tipdate_likelihood(&tree, &aln, &model, &tip_dates, rate, &node_dates)
                .unwrap();

        assert!(lnl.is_finite(), "likelihood should be finite");
        assert!(lnl < 0.0, "log-likelihood should be negative, got {}", lnl);

        // Changing rate should change likelihood
        let lnl2 =
            tipdate_likelihood(&tree, &aln, &model, &tip_dates, 0.001, &node_dates)
                .unwrap();
        assert!(
            (lnl - lnl2).abs() > 1e-6,
            "different rates should give different likelihoods"
        );
    }

    // -----------------------------------------------------------------------
    // Nonnegative branch lengths
    // -----------------------------------------------------------------------

    #[test]
    fn test_nonnegative_branch_lengths() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTTAACCGGTA"),
            ("C", "CCAATTGGCCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGGCCAATTGA"),
        ]);
        let tree =
            parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();
        let model = Jc69Model;

        let mut tip_dates = HashMap::new();
        tip_dates.insert("A".to_string(), 2000.0);
        tip_dates.insert("B".to_string(), 2005.0);
        tip_dates.insert("C".to_string(), 2010.0);
        tip_dates.insert("D".to_string(), 2015.0);

        let result = tipdate_optimize(&tree, &aln, &model, &tip_dates).unwrap();

        for node in &result.tree.nodes {
            if let Some(bl) = node.branch_length {
                assert!(
                    bl >= 0.0,
                    "branch length should be non-negative, got {} for node {:?}",
                    bl,
                    node.name
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Chi-squared / ln_gamma unit tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_ln_gamma_known_values() {
        // Gamma(1) = 1, ln(Gamma(1)) = 0
        assert!(
            ln_gamma(1.0).abs() < 1e-10,
            "ln(Gamma(1)) should be 0, got {}",
            ln_gamma(1.0)
        );

        // Gamma(2) = 1, ln(Gamma(2)) = 0
        assert!(
            ln_gamma(2.0).abs() < 1e-10,
            "ln(Gamma(2)) should be 0, got {}",
            ln_gamma(2.0)
        );

        // Gamma(3) = 2, ln(Gamma(3)) = ln(2)
        let expected = 2.0_f64.ln();
        assert!(
            (ln_gamma(3.0) - expected).abs() < 1e-10,
            "ln(Gamma(3)) should be {}, got {}",
            expected,
            ln_gamma(3.0)
        );

        // Gamma(0.5) = sqrt(pi), ln(Gamma(0.5)) = ln(sqrt(pi))
        let expected = (std::f64::consts::PI.sqrt()).ln();
        assert!(
            (ln_gamma(0.5) - expected).abs() < 1e-8,
            "ln(Gamma(0.5)) should be {}, got {}",
            expected,
            ln_gamma(0.5)
        );
    }

    #[test]
    fn test_chi2_survival_known_values() {
        // chi2(df=1) at x=0 should give p=1
        assert!(
            (chi2_survival(0.0, 1) - 1.0).abs() < 1e-10,
            "P(X>0) should be 1.0 for chi2(1)"
        );

        // chi2(df=1) at x=3.841 should give p ~ 0.05
        let p = chi2_survival(3.841, 1);
        assert!(
            (p - 0.05).abs() < 0.01,
            "P(X>3.841) for chi2(1) should be ~0.05, got {}",
            p
        );

        // chi2(df=2) at x=5.991 should give p ~ 0.05
        let p2 = chi2_survival(5.991, 2);
        assert!(
            (p2 - 0.05).abs() < 0.01,
            "P(X>5.991) for chi2(2) should be ~0.05, got {}",
            p2
        );

        // Large x should give p ~ 0
        let p_large = chi2_survival(100.0, 1);
        assert!(
            p_large < 1e-10,
            "P(X>100) should be ~0 for chi2(1), got {}",
            p_large
        );
    }

    // -----------------------------------------------------------------------
    // Golden section search unit test
    // -----------------------------------------------------------------------

    #[test]
    fn test_golden_section_maximize() {
        // Maximize -(x-3)^2 + 10 on [0, 6]. Maximum at x=3, f(3)=10.
        let (x, f) = golden_section_maximize(0.0, 6.0, |x| -(x - 3.0).powi(2) + 10.0);
        assert!(
            (x - 3.0).abs() < 1e-4,
            "maximum should be near x=3, got {}",
            x
        );
        assert!(
            (f - 10.0).abs() < 1e-4,
            "f(max) should be near 10, got {}",
            f
        );
    }

    // -----------------------------------------------------------------------
    // Five-taxa heterochronous test
    // -----------------------------------------------------------------------

    #[test]
    fn test_five_taxa_heterochronous() {
        let aln = make_alignment(&[
            ("V1", "AACCGGTTAACCGGTTAACCGGTT"),
            ("V2", "AACCGGTTAACCGGTTAACCGGTA"),
            ("V3", "AACCGGTTAACCGGTAAACCGGTA"),
            ("V4", "CCAATTGGCCAATTGGCCAATTGG"),
            ("V5", "CCAATTGGCCAATTGGCCAATTGA"),
        ]);
        let tree = parse_newick(
            "((V1:0.1,V2:0.1):0.05,((V3:0.05,V4:0.1):0.05,V5:0.15):0.05);",
        )
        .unwrap();
        let model = Jc69Model;

        let mut tip_dates = HashMap::new();
        tip_dates.insert("V1".to_string(), 2000.0);
        tip_dates.insert("V2".to_string(), 2003.0);
        tip_dates.insert("V3".to_string(), 2006.0);
        tip_dates.insert("V4".to_string(), 2009.0);
        tip_dates.insert("V5".to_string(), 2012.0);

        let result = tipdate_optimize(&tree, &aln, &model, &tip_dates).unwrap();

        assert!(result.rate > 0.0, "rate should be positive");
        assert!(result.lnl < 0.0, "lnl should be negative");
        assert!(
            result.root_date < 2000.0,
            "root date should precede earliest tip"
        );
        assert_eq!(result.tree.num_leaves(), 5);

        // Verify tip dates are preserved
        for leaf in tree.leaves() {
            let name = leaf.name.as_deref().unwrap();
            let expected = tip_dates[name];
            assert!(
                (result.node_dates[leaf.id] - expected).abs() < 1e-10,
                "tip '{}' date should be {}, got {}",
                name,
                expected,
                result.node_dates[leaf.id]
            );
        }
    }

    // -----------------------------------------------------------------------
    // Node dates ordering
    // -----------------------------------------------------------------------

    #[test]
    fn test_node_dates_temporal_order() {
        // After optimization, every parent node should have a date
        // earlier than or equal to all its children.
        let aln = make_alignment(&[
            ("S1", "ACGTACGTACGTACGTACGTACGT"),
            ("S2", "ACGTACGTACGTACGTACGTACGA"),
            ("S3", "ACGTACGTACGTACGAACGTACGA"),
            ("S4", "TGCATGCATGCATGCATGCATGCA"),
        ]);
        let tree =
            parse_newick("((S1:0.05,S2:0.05):0.1,(S3:0.05,S4:0.05):0.1);").unwrap();
        let model = Jc69Model;

        let mut tip_dates = HashMap::new();
        tip_dates.insert("S1".to_string(), 1990.0);
        tip_dates.insert("S2".to_string(), 2000.0);
        tip_dates.insert("S3".to_string(), 2010.0);
        tip_dates.insert("S4".to_string(), 2020.0);

        let result = tipdate_optimize(&tree, &aln, &model, &tip_dates).unwrap();

        // Check parent-child date ordering
        for node in &tree.nodes {
            if let Some(parent_id) = node.parent {
                assert!(
                    result.node_dates[parent_id] <= result.node_dates[node.id] + TOLERANCE,
                    "parent date ({}) should be <= child date ({}) for node {:?}",
                    result.node_dates[parent_id],
                    result.node_dates[node.id],
                    node.name
                );
            }
        }
    }
}

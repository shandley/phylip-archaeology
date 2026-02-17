//! Substitution models for likelihood calculations.
//!
//! These models compute transition probability matrices P(t) needed by
//! the pruning algorithm. Each entry P\[i\]\[j\] gives the probability
//! of observing state j at time t given ancestral state i at time 0.
//!
//! States are indexed as: 0=A, 1=C, 2=G, 3=T.
//!
//! Branch lengths t are in units of expected substitutions per site.

use crate::tree::BaseFrequencies;

/// A nucleotide substitution model for the pruning algorithm.
///
/// Provides transition probability matrices P(t) and equilibrium
/// base frequencies needed for computing tree likelihoods.
pub trait SubstitutionModel {
    /// Compute the 4x4 transition probability matrix P(t).
    ///
    /// `P[i][j]` = Pr(state j at descendant | state i at ancestor),
    /// given branch length `t` in expected substitutions per site.
    fn transition_probs(&self, t: f64) -> [[f64; 4]; 4];

    /// Equilibrium base frequencies `[pi_A, pi_C, pi_G, pi_T]`.
    fn frequencies(&self) -> [f64; 4];

    /// Model name.
    fn name(&self) -> &str;
}

/// Jukes-Cantor (1969) substitution model.
///
/// Equal base frequencies (0.25 each) and equal substitution rates.
///
/// P(same; t) = 1/4 + 3/4 * exp(-4t/3)
/// P(diff; t) = 1/4 - 1/4 * exp(-4t/3)
#[derive(Debug, Clone)]
pub struct Jc69Model;

impl SubstitutionModel for Jc69Model {
    fn transition_probs(&self, t: f64) -> [[f64; 4]; 4] {
        let e = (-4.0 * t / 3.0).exp();
        let p_same = 0.25 + 0.75 * e;
        let p_diff = 0.25 - 0.25 * e;
        [
            [p_same, p_diff, p_diff, p_diff],
            [p_diff, p_same, p_diff, p_diff],
            [p_diff, p_diff, p_same, p_diff],
            [p_diff, p_diff, p_diff, p_same],
        ]
    }

    fn frequencies(&self) -> [f64; 4] {
        [0.25, 0.25, 0.25, 0.25]
    }

    fn name(&self) -> &str {
        "JC69"
    }
}

/// Felsenstein 1984 (F84) substitution model.
///
/// Allows unequal base frequencies and distinguishes transitions from
/// transversions. This is the model used in PHYLIP's DNAML program.
///
/// The transition probability is factored into three components
/// (Felsenstein 1981, 1984):
///
/// ```text
/// P(i->j; t) = delta(i,j) * exp1
///            + same_class(i,j) * (pi_j / pi_class) * (exp2 - exp1)
///            + pi_j * (1 - exp2)
/// ```
///
/// where:
/// - `exp1 = exp(-(xi + xv) * t / fracchange)`
/// - `exp2 = exp(-xv * t / fracchange)`
/// - `same_class` is 1 if both purines or both pyrimidines
/// - `pi_class` is `pi_R` for purines, `pi_Y` for pyrimidines
/// - `xi, xv, fracchange` are derived from frequencies and ttratio
#[derive(Debug, Clone)]
pub struct F84Model {
    /// Base frequencies.
    pub freqs: BaseFrequencies,
    /// Transition/transversion ratio.
    pub ttratio: f64,
    // Derived parameters (precomputed)
    xi: f64,
    xv: f64,
    fracchange: f64,
}

impl F84Model {
    /// Create a new F84 model with the given base frequencies and
    /// transition/transversion ratio.
    pub fn new(freqs: BaseFrequencies, ttratio: f64) -> Self {
        let (xi, xv, fracchange) = Self::compute_params(&freqs, ttratio);
        Self {
            freqs,
            ttratio,
            xi,
            xv,
            fracchange,
        }
    }

    /// Create an F84 model with equal frequencies (reduces to K2P-like).
    pub fn equal_freqs(ttratio: f64) -> Self {
        Self::new(BaseFrequencies::equal(), ttratio)
    }

    /// Compute xi, xv, fracchange from frequencies and ttratio.
    ///
    /// This follows PHYLIP's getbasefreqs() logic.
    fn compute_params(freqs: &BaseFrequencies, ttratio: f64) -> (f64, f64, f64) {
        let freq_r = freqs.freq_r();
        let freq_y = freqs.freq_y();

        let freq_ar = freqs.freq_a / freq_r;
        let freq_gr = freqs.freq_g / freq_r;
        let freq_cy = freqs.freq_c / freq_y;
        let freq_ty = freqs.freq_t / freq_y;

        let aa =
            freqs.freq_a * freqs.freq_g / freq_r + freqs.freq_c * freqs.freq_t / freq_y;
        let bb = freqs.freq_a * freqs.freq_g + freqs.freq_c * freqs.freq_t;

        let xv = 1.0;
        let xi_denom = freq_ar * freq_gr * freq_r + freq_cy * freq_ty * freq_y;
        let xi = if xi_denom > 0.0 {
            xv * (ttratio - 1.0) / xi_denom * (freq_r * freq_y)
        } else {
            0.0
        };

        let fracchange = xi * 2.0 * aa + xv * (1.0 - 2.0 * bb);

        (xi, xv, fracchange)
    }

    /// Returns whether state index i is a purine (A=0 or G=2).
    fn is_purine(i: usize) -> bool {
        i == 0 || i == 2
    }
}

impl SubstitutionModel for F84Model {
    fn transition_probs(&self, t: f64) -> [[f64; 4]; 4] {
        if self.fracchange <= 0.0 || t <= 0.0 {
            return [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ];
        }

        let freq = [
            self.freqs.freq_a,
            self.freqs.freq_c,
            self.freqs.freq_g,
            self.freqs.freq_t,
        ];
        let freq_r = self.freqs.freq_r();
        let freq_y = self.freqs.freq_y();

        // Exponential components, normalized so t is in expected subs/site
        let exp1 = (-(self.xi + self.xv) * t / self.fracchange).exp();
        let exp2 = (-self.xv * t / self.fracchange).exp();

        let mut p = [[0.0f64; 4]; 4];

        for i in 0..4 {
            let i_purine = Self::is_purine(i);
            let pi_class = if i_purine { freq_r } else { freq_y };

            for j in 0..4 {
                let j_purine = Self::is_purine(j);
                let same_class = i_purine == j_purine;

                // Three-component factored form
                let mut val = freq[j] * (1.0 - exp2); // transversion component

                if same_class {
                    val += (freq[j] / pi_class) * (exp2 - exp1); // transition component
                }

                if i == j {
                    val += exp1; // identity component
                }

                p[i][j] = val;
            }
        }

        p
    }

    fn frequencies(&self) -> [f64; 4] {
        [
            self.freqs.freq_a,
            self.freqs.freq_c,
            self.freqs.freq_g,
            self.freqs.freq_t,
        ]
    }

    fn name(&self) -> &str {
        "F84"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jc69_identity_at_zero() {
        let model = Jc69Model;
        let p = model.transition_probs(0.0);
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (p[i][j] - expected).abs() < 1e-12,
                    "P[{}][{}] = {}, expected {}",
                    i,
                    j,
                    p[i][j],
                    expected
                );
            }
        }
    }

    #[test]
    fn test_jc69_equilibrium_at_infinity() {
        let model = Jc69Model;
        let p = model.transition_probs(100.0);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (p[i][j] - 0.25).abs() < 1e-10,
                    "P[{}][{}] = {}, expected 0.25",
                    i,
                    j,
                    p[i][j]
                );
            }
        }
    }

    #[test]
    fn test_jc69_rows_sum_to_one() {
        let model = Jc69Model;
        for &t in &[0.0, 0.01, 0.1, 0.5, 1.0, 5.0] {
            let p = model.transition_probs(t);
            for i in 0..4 {
                let row_sum: f64 = p[i].iter().sum();
                assert!(
                    (row_sum - 1.0).abs() < 1e-12,
                    "Row {} sums to {} at t={}",
                    i,
                    row_sum,
                    t
                );
            }
        }
    }

    #[test]
    fn test_jc69_symmetric() {
        let model = Jc69Model;
        let p = model.transition_probs(0.3);
        for i in 0..4 {
            for j in 0..4 {
                // JC69 should have pi_i * P[i][j] = pi_j * P[j][i]
                // With equal freqs: P[i][j] = P[j][i]
                assert!(
                    (p[i][j] - p[j][i]).abs() < 1e-12,
                    "P[{}][{}] != P[{}][{}]",
                    i,
                    j,
                    j,
                    i
                );
            }
        }
    }

    #[test]
    fn test_jc69_known_value() {
        let model = Jc69Model;
        let t = 0.1;
        let p = model.transition_probs(t);
        let e = (-4.0 * t / 3.0).exp();
        let expected_same = 0.25 + 0.75 * e;
        let expected_diff = 0.25 - 0.25 * e;
        assert!((p[0][0] - expected_same).abs() < 1e-12);
        assert!((p[0][1] - expected_diff).abs() < 1e-12);
    }

    #[test]
    fn test_f84_identity_at_zero() {
        let model = F84Model::equal_freqs(2.0);
        let p = model.transition_probs(0.0);
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (p[i][j] - expected).abs() < 1e-12,
                    "P[{}][{}] = {}, expected {}",
                    i,
                    j,
                    p[i][j],
                    expected
                );
            }
        }
    }

    #[test]
    fn test_f84_equilibrium_at_infinity() {
        let freqs = BaseFrequencies::new(0.3, 0.2, 0.2, 0.3);
        let model = F84Model::new(freqs, 2.0);
        let p = model.transition_probs(100.0);
        let pi = model.frequencies();
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (p[i][j] - pi[j]).abs() < 1e-10,
                    "P[{}][{}] = {}, expected pi[{}] = {}",
                    i,
                    j,
                    p[i][j],
                    j,
                    pi[j]
                );
            }
        }
    }

    #[test]
    fn test_f84_rows_sum_to_one() {
        let freqs = BaseFrequencies::new(0.3, 0.2, 0.2, 0.3);
        let model = F84Model::new(freqs, 2.0);
        for &t in &[0.0, 0.01, 0.1, 0.5, 1.0, 5.0] {
            let p = model.transition_probs(t);
            for i in 0..4 {
                let row_sum: f64 = p[i].iter().sum();
                assert!(
                    (row_sum - 1.0).abs() < 1e-12,
                    "Row {} sums to {} at t={}",
                    i,
                    row_sum,
                    t
                );
            }
        }
    }

    #[test]
    fn test_f84_detailed_balance() {
        // Detailed balance: pi_i * P(i->j) = pi_j * P(j->i)
        let freqs = BaseFrequencies::new(0.3, 0.2, 0.2, 0.3);
        let model = F84Model::new(freqs, 2.0);
        let p = model.transition_probs(0.5);
        let pi = model.frequencies();
        for i in 0..4 {
            for j in 0..4 {
                let lhs = pi[i] * p[i][j];
                let rhs = pi[j] * p[j][i];
                assert!(
                    (lhs - rhs).abs() < 1e-12,
                    "Detailed balance violated: pi[{}]*P[{}][{}] = {} != pi[{}]*P[{}][{}] = {}",
                    i,
                    i,
                    j,
                    lhs,
                    j,
                    j,
                    i,
                    rhs
                );
            }
        }
    }

    #[test]
    fn test_f84_reduces_to_jc69() {
        // With equal freqs and ttratio=1, F84 should equal JC69
        let f84 = F84Model::equal_freqs(1.0);
        let jc69 = Jc69Model;
        for &t in &[0.01, 0.1, 0.5, 1.0] {
            let p_f84 = f84.transition_probs(t);
            let p_jc69 = jc69.transition_probs(t);
            for i in 0..4 {
                for j in 0..4 {
                    assert!(
                        (p_f84[i][j] - p_jc69[i][j]).abs() < 1e-10,
                        "F84(ttratio=1) != JC69 at t={}: P[{}][{}] = {} vs {}",
                        t,
                        i,
                        j,
                        p_f84[i][j],
                        p_jc69[i][j]
                    );
                }
            }
        }
    }

    #[test]
    fn test_f84_transitions_faster_than_transversions() {
        // With ttratio > 1, transitions should have higher probability
        // than transversions at moderate t
        let model = F84Model::equal_freqs(4.0);
        let p = model.transition_probs(0.1);
        // A->G is a transition (both purines)
        // A->C is a transversion
        assert!(
            p[0][2] > p[0][1],
            "Transition P(A->G) = {} should be > transversion P(A->C) = {}",
            p[0][2],
            p[0][1]
        );
    }

    #[test]
    fn test_f84_all_positive() {
        let freqs = BaseFrequencies::new(0.4, 0.1, 0.1, 0.4);
        let model = F84Model::new(freqs, 3.0);
        for &t in &[0.001, 0.01, 0.1, 0.5, 1.0, 5.0] {
            let p = model.transition_probs(t);
            for i in 0..4 {
                for j in 0..4 {
                    assert!(
                        p[i][j] >= 0.0,
                        "P[{}][{}] = {} at t={} is negative",
                        i,
                        j,
                        p[i][j],
                        t
                    );
                }
            }
        }
    }
}

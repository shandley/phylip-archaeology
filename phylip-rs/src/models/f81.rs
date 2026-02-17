//! Felsenstein 1981 model (F81).
//!
//! Allows unequal base frequencies while assuming equal substitution
//! rates for all types of substitutions. This model generalizes JC69
//! by relaxing the equal-frequency assumption.
//!
//! # Distance formula
//!
//! Let `p` = fraction of identical sites among comparable positions,
//! and `pi_A, pi_C, pi_G, pi_T` be the base frequencies.
//!
//! Define:
//! ```text
//! b = 1 - (pi_A^2 + pi_C^2 + pi_G^2 + pi_T^2)
//! ```
//!
//! Then:
//! ```text
//! d = -b * ln(1 - p_diff / b)
//! ```
//! where `p_diff = 1 - p` is the fraction of different sites.
//!
//! Equivalently:
//! ```text
//! d = -b * ln((p - (1 - b)) / b)
//! ```
//!
//! The formula requires `p_diff < b`, i.e., the observed proportion of
//! differences must be less than the maximum possible under the model.
//!
//! With gamma rate correction (shape parameter alpha):
//! ```text
//! d = b * alpha * ((1 - p_diff / b)^(-1/alpha) - 1)
//! ```
//!
//! # Note on base frequencies
//!
//! Base frequencies can be:
//! - Specified directly by the user
//! - Estimated from the pairwise comparison (counting bases at comparable sites)
//!
//! When `use_empirical = true` (the default), frequencies are estimated from
//! each pair of sequences being compared. When `false`, the frequencies
//! provided at construction time are used for all comparisons.
//!
//! # References
//!
//! Felsenstein, J. (1981). "Evolutionary trees from DNA sequences: a maximum
//! likelihood approach." *Journal of Molecular Evolution*, 17:368-376.

use crate::tree::{Base, BaseFrequencies};

use super::{count_bases_pairwise, count_differences, DistanceError, DistanceModel};

/// Felsenstein 1981 substitution model.
///
/// Allows unequal base frequencies while assuming equal rates of
/// substitution between all pairs of nucleotides. This generalizes
/// JC69 by relaxing the equal-frequency constraint.
#[derive(Debug, Clone)]
pub struct F81 {
    /// Base frequencies to use. If `use_empirical` is true, these serve
    /// as defaults only and are overridden by empirical estimates for
    /// each pairwise comparison.
    pub freqs: BaseFrequencies,

    /// If true, estimate base frequencies from each pair of sequences
    /// being compared. If false, use the fixed frequencies in `freqs`.
    pub use_empirical: bool,

    /// Optional gamma distribution shape parameter (alpha) for
    /// among-site rate variation.
    pub alpha: Option<f64>,
}

impl F81 {
    /// Create a new F81 model with the given base frequencies.
    ///
    /// By default, uses fixed frequencies (not empirical estimation).
    pub fn new(freqs: BaseFrequencies) -> Self {
        Self {
            freqs,
            use_empirical: false,
            alpha: None,
        }
    }

    /// Create a new F81 model that estimates base frequencies empirically
    /// from each pairwise comparison.
    pub fn empirical() -> Self {
        Self {
            freqs: BaseFrequencies::equal(),
            use_empirical: true,
            alpha: None,
        }
    }

    /// Create a new F81 model with gamma-distributed rate correction.
    pub fn with_gamma(freqs: BaseFrequencies, alpha: f64) -> Self {
        Self {
            freqs,
            use_empirical: false,
            alpha: Some(alpha),
        }
    }

    /// Compute `b = 1 - sum(pi_i^2)` from base frequencies.
    ///
    /// This is the maximum fraction of sites that can differ under the model
    /// as time goes to infinity (the equilibrium dissimilarity).
    fn compute_b(freqs: &BaseFrequencies) -> f64 {
        1.0 - (freqs.freq_a * freqs.freq_a
            + freqs.freq_c * freqs.freq_c
            + freqs.freq_g * freqs.freq_g
            + freqs.freq_t * freqs.freq_t)
    }

    /// Estimate base frequencies from a pairwise comparison.
    fn estimate_freqs(seq1: &[Base], seq2: &[Base]) -> BaseFrequencies {
        let counts = count_bases_pairwise(seq1, seq2);
        let total: usize = counts.iter().sum();
        if total == 0 {
            return BaseFrequencies::equal();
        }
        let total_f = total as f64;
        BaseFrequencies::new(
            counts[0] as f64 / total_f,
            counts[1] as f64 / total_f,
            counts[2] as f64 / total_f,
            counts[3] as f64 / total_f,
        )
    }
}

impl Default for F81 {
    fn default() -> Self {
        Self::new(BaseFrequencies::equal())
    }
}

impl DistanceModel for F81 {
    fn distance(&self, seq1: &[Base], seq2: &[Base]) -> Result<f64, DistanceError> {
        if seq1.len() != seq2.len() {
            return Err(DistanceError::UnequalLengths);
        }

        let (identical, total) = count_differences(seq1, seq2);

        if total == 0 {
            return Err(DistanceError::NoOverlap { seq1: 0, seq2: 0 });
        }

        // Determine base frequencies to use
        let freqs = if self.use_empirical {
            Self::estimate_freqs(seq1, seq2)
        } else {
            self.freqs
        };

        // b = 1 - sum(pi_i^2)
        let b = Self::compute_b(&freqs);

        // Fraction of different sites
        let p_diff = 1.0 - (identical as f64 / total as f64);

        // The argument to the log: 1 - p_diff / b
        // Must be positive, i.e., p_diff < b
        let inner = 1.0 - p_diff / b;

        if inner <= 0.0 {
            return Err(DistanceError::InfiniteDistance {
                seq1: 0,
                seq2: 0,
                reason: format!(
                    "fraction of differences ({:.4}) exceeds maximum ({:.4}) under F81 model",
                    p_diff, b
                ),
            });
        }

        let d = match self.alpha {
            None => {
                // Standard F81: d = -b * ln(1 - p_diff / b)
                -b * inner.ln()
            }
            Some(alpha) => {
                // Gamma-corrected: d = b * alpha * ((1 - p_diff / b)^(-1/alpha) - 1)
                b * alpha * (inner.powf(-1.0 / alpha) - 1.0)
            }
        };

        if !d.is_finite() || d < 0.0 {
            return Err(DistanceError::InfiniteDistance {
                seq1: 0,
                seq2: 0,
                reason: "computed distance is not finite".to_string(),
            });
        }

        Ok(d)
    }

    fn name(&self) -> &str {
        "F81"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build a pair of sequences from strings.
    fn make_seqs(s1: &str, s2: &str) -> (Vec<Base>, Vec<Base>) {
        let seq1: Vec<Base> = s1.chars().map(|c| Base::from_char(c).unwrap()).collect();
        let seq2: Vec<Base> = s2.chars().map(|c| Base::from_char(c).unwrap()).collect();
        (seq1, seq2)
    }

    #[test]
    fn test_identical_sequences() {
        let model = F81::new(BaseFrequencies::equal());
        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACGT");
        let d = model.distance(&s1, &s2).unwrap();
        assert!(
            d.abs() < 1e-10,
            "identical sequences should have distance 0, got {}",
            d
        );
    }

    #[test]
    fn test_equal_freqs_matches_jc69() {
        // With equal base frequencies, F81 should give the same result as JC69.
        // b = 1 - 4 * 0.25^2 = 1 - 0.25 = 0.75
        // F81: d = -0.75 * ln(1 - p_diff / 0.75)
        // JC69: d = -0.75 * ln((4p - 1) / 3)
        // p_diff = 1 - p, so 1 - p_diff / 0.75 = 1 - (1-p)/0.75 = (0.75 - 1 + p)/0.75 = (p - 0.25)/0.75
        // And (4p - 1)/3 = (4p - 1)/3
        // (p - 0.25)/0.75 = (p - 1/4) / (3/4) = (4p - 1)/3 -- same!
        let f81 = F81::new(BaseFrequencies::equal());
        let jc69 = super::super::jc69::JC69::new();

        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACGA");
        let d_f81 = f81.distance(&s1, &s2).unwrap();
        let d_jc69 = jc69.distance(&s1, &s2).unwrap();

        assert!(
            (d_f81 - d_jc69).abs() < 1e-10,
            "F81 with equal freqs ({}) should match JC69 ({})",
            d_f81,
            d_jc69
        );
    }

    #[test]
    fn test_unequal_frequencies() {
        // With unequal frequencies, the distance should differ from JC69
        let freqs = BaseFrequencies::new(0.4, 0.1, 0.1, 0.4);
        let model = F81::new(freqs);

        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACGA");
        let d = model.distance(&s1, &s2).unwrap();

        // b = 1 - (0.16 + 0.01 + 0.01 + 0.16) = 1 - 0.34 = 0.66
        let b: f64 = 1.0 - (0.4 * 0.4 + 0.1 * 0.1 + 0.1 * 0.1 + 0.4 * 0.4);
        let p_diff: f64 = 1.0 / 8.0;
        let inner: f64 = 1.0 - p_diff / b;
        let expected: f64 = -b * inner.ln();
        assert!(
            (d - expected).abs() < 1e-10,
            "expected {}, got {}",
            expected,
            d
        );
    }

    #[test]
    fn test_known_value() {
        // Specific numerical test:
        // freqs = (0.3, 0.2, 0.2, 0.3)
        // b = 1 - (0.09 + 0.04 + 0.04 + 0.09) = 1 - 0.26 = 0.74
        // 8 sites, 2 differences -> p_diff = 0.25
        // d = -0.74 * ln(1 - 0.25/0.74) = -0.74 * ln(0.662162...)
        let freqs = BaseFrequencies::new(0.3, 0.2, 0.2, 0.3);
        let model = F81::new(freqs);

        let (s1, s2) = make_seqs("ACGTACGT", "ACGAACGA");
        let d = model.distance(&s1, &s2).unwrap();

        let b: f64 = 0.74;
        let p_diff: f64 = 2.0 / 8.0;
        let inner: f64 = 1.0 - p_diff / b;
        let expected: f64 = -b * inner.ln();
        assert!(
            (d - expected).abs() < 1e-10,
            "expected {}, got {}",
            expected,
            d
        );
    }

    #[test]
    fn test_too_divergent() {
        // If all sites are different, p_diff = 1.0 > b for any valid frequencies
        let model = F81::new(BaseFrequencies::equal());
        let (s1, s2) = make_seqs("AAAA", "CCCC");
        let result = model.distance(&s1, &s2);
        assert!(result.is_err());
    }

    #[test]
    fn test_gaps_skipped() {
        let model = F81::new(BaseFrequencies::equal());
        let (s1, s2) = make_seqs("A-GT", "ACGT");
        let d = model.distance(&s1, &s2).unwrap();
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn test_ambiguous_bases_skipped() {
        let model = F81::new(BaseFrequencies::equal());
        let (s1, s2) = make_seqs("ANGT", "ACGT");
        let d = model.distance(&s1, &s2).unwrap();
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn test_no_overlap() {
        let model = F81::new(BaseFrequencies::equal());
        let (s1, s2) = make_seqs("----", "ACGT");
        let result = model.distance(&s1, &s2);
        assert!(result.is_err());
    }

    #[test]
    fn test_unequal_lengths() {
        let model = F81::new(BaseFrequencies::equal());
        let s1: Vec<Base> = "ACGT".chars().map(|c| Base::from_char(c).unwrap()).collect();
        let s2: Vec<Base> = "AC".chars().map(|c| Base::from_char(c).unwrap()).collect();
        let result = model.distance(&s1, &s2);
        assert!(result.is_err());
        match result.unwrap_err() {
            DistanceError::UnequalLengths => {}
            e => panic!("expected UnequalLengths, got {:?}", e),
        }
    }

    #[test]
    fn test_empirical_frequencies() {
        let model = F81::empirical();
        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACGA");
        let d = model.distance(&s1, &s2).unwrap();
        // With empirical frequencies estimated from these sequences,
        // the result should still be positive
        assert!(d > 0.0);
    }

    #[test]
    fn test_gamma_correction() {
        let alpha = 1.0;
        let freqs = BaseFrequencies::equal();
        let model = F81::with_gamma(freqs, alpha);

        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACGA");
        let d = model.distance(&s1, &s2).unwrap();

        let b = F81::compute_b(&freqs);
        let p_diff = 1.0 / 8.0;
        let inner = 1.0 - p_diff / b;
        let expected = b * alpha * (inner.powf(-1.0 / alpha) - 1.0);
        assert!(
            (d - expected).abs() < 1e-10,
            "expected {}, got {}",
            expected,
            d
        );
    }

    #[test]
    fn test_gamma_correction_equal_freqs_matches_jc69_gamma() {
        // With equal frequencies and gamma correction, F81 should match JC69 with gamma
        let alpha = 0.5;
        let f81 = F81::with_gamma(BaseFrequencies::equal(), alpha);
        let jc69 = super::super::jc69::JC69::with_gamma(alpha);

        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACGA");
        let d_f81 = f81.distance(&s1, &s2).unwrap();
        let d_jc69 = jc69.distance(&s1, &s2).unwrap();

        assert!(
            (d_f81 - d_jc69).abs() < 1e-10,
            "F81+gamma with equal freqs ({}) should match JC69+gamma ({})",
            d_f81,
            d_jc69
        );
    }

    #[test]
    fn test_distance_increases_with_differences() {
        let model = F81::new(BaseFrequencies::equal());

        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACGT");
        let d0 = model.distance(&s1, &s2).unwrap();

        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACGA");
        let d1 = model.distance(&s1, &s2).unwrap();

        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACAA");
        let d2 = model.distance(&s1, &s2).unwrap();

        assert!(d0 < d1);
        assert!(d1 < d2);
    }

    #[test]
    fn test_b_value_equal_freqs() {
        let freqs = BaseFrequencies::equal();
        let b = F81::compute_b(&freqs);
        assert!(
            (b - 0.75).abs() < 1e-10,
            "b with equal freqs should be 0.75, got {}",
            b
        );
    }

    #[test]
    fn test_b_value_unequal_freqs() {
        let freqs = BaseFrequencies::new(0.4, 0.1, 0.1, 0.4);
        let b = F81::compute_b(&freqs);
        let expected = 1.0 - (0.16 + 0.01 + 0.01 + 0.16);
        assert!(
            (b - expected).abs() < 1e-10,
            "expected b = {}, got {}",
            expected,
            b
        );
    }

    #[test]
    fn test_higher_b_gives_higher_distance() {
        // Higher b (more uniform frequencies) should give a smaller distance
        // correction than lower b (more skewed frequencies) for the same
        // observed proportion of differences.
        let freqs_equal = BaseFrequencies::equal(); // b = 0.75
        let freqs_skewed = BaseFrequencies::new(0.7, 0.1, 0.1, 0.1); // b = 1 - (0.49+0.01+0.01+0.01) = 0.48

        let model_equal = F81::new(freqs_equal);
        let model_skewed = F81::new(freqs_skewed);

        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACGA");
        let d_equal = model_equal.distance(&s1, &s2).unwrap();
        let d_skewed = model_skewed.distance(&s1, &s2).unwrap();

        // With skewed frequencies, the same observed difference represents more
        // evolutionary change, so the distance should be larger
        assert!(
            d_skewed > d_equal,
            "skewed ({}) should give larger distance than equal ({})",
            d_skewed,
            d_equal
        );
    }
}

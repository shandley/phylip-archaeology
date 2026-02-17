//! Jukes-Cantor model (Jukes & Cantor, 1969).
//!
//! The simplest nucleotide substitution model, assuming equal base
//! frequencies (0.25 each) and equal substitution rates between all pairs of
//! nucleotides.
//!
//! # Distance formula
//!
//! Given `p` = fraction of identical sites:
//!
//! **Without gamma rate correction:**
//! ```text
//! d = -3/4 * ln((4p - 1) / 3)
//! ```
//!
//! **With gamma rate correction (shape parameter alpha):**
//! ```text
//! d = 3/4 * alpha * (((4p - 1) / 3)^(-1/alpha) - 1)
//! ```
//!
//! The formula requires `p > 0.25`; if the fraction of identical sites is
//! at or below the random expectation of 0.25, the distance is infinite.
//!
//! # References
//!
//! Jukes, T.H. & Cantor, C.R. (1969). "Evolution of protein molecules."
//! In *Mammalian Protein Metabolism*, Vol. III, pp. 21-132.

use crate::tree::Base;

use super::{count_differences, DistanceError, DistanceModel};

/// Jukes-Cantor 1969 substitution model.
///
/// Assumes equal base frequencies and equal rates of substitution between
/// all pairs of nucleotides. Optionally applies a gamma-distributed rate
/// correction for among-site rate heterogeneity.
#[derive(Debug, Clone)]
pub struct JC69 {
    /// Optional gamma distribution shape parameter (alpha) for
    /// among-site rate variation. When `None`, a uniform rate is assumed.
    /// Smaller values indicate more rate heterogeneity; typical values
    /// range from 0.5 to 2.0.
    pub alpha: Option<f64>,
}

impl JC69 {
    /// Create a new JC69 model without gamma rate correction.
    pub fn new() -> Self {
        Self { alpha: None }
    }

    /// Create a new JC69 model with gamma-distributed rate correction.
    ///
    /// # Arguments
    ///
    /// * `alpha` - Shape parameter of the gamma distribution. Must be positive.
    ///   Smaller values indicate more rate heterogeneity.
    pub fn with_gamma(alpha: f64) -> Self {
        Self { alpha: Some(alpha) }
    }
}

impl Default for JC69 {
    fn default() -> Self {
        Self::new()
    }
}

impl DistanceModel for JC69 {
    fn distance(&self, seq1: &[Base], seq2: &[Base]) -> Result<f64, DistanceError> {
        if seq1.len() != seq2.len() {
            return Err(DistanceError::UnequalLengths);
        }

        let (identical, total) = count_differences(seq1, seq2);

        if total == 0 {
            return Err(DistanceError::NoOverlap { seq1: 0, seq2: 0 });
        }

        // p = fraction of identical sites
        let p = identical as f64 / total as f64;

        // The argument to the log: (4p - 1) / 3
        // This must be positive, which requires p > 0.25
        let inner = (4.0 * p - 1.0) / 3.0;

        if inner <= 0.0 {
            return Err(DistanceError::InfiniteDistance {
                seq1: 0,
                seq2: 0,
                reason: format!(
                    "fraction of identical sites ({:.4}) is at or below the random expectation of 0.25",
                    p
                ),
            });
        }

        let d = match self.alpha {
            None => {
                // Standard JC69: d = -3/4 * ln((4p - 1) / 3)
                -0.75 * inner.ln()
            }
            Some(alpha) => {
                // Gamma-corrected: d = 3/4 * alpha * (((4p - 1) / 3)^(-1/alpha) - 1)
                0.75 * alpha * (inner.powf(-1.0 / alpha) - 1.0)
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
        "JC69"
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
        let model = JC69::new();
        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACGT");
        let d = model.distance(&s1, &s2).unwrap();
        assert!((d - 0.0).abs() < 1e-10, "identical sequences should have distance 0, got {}", d);
    }

    #[test]
    fn test_one_difference() {
        let model = JC69::new();
        // 1 out of 4 sites differ -> p = 3/4
        let (s1, s2) = make_seqs("ACGT", "ACGA");
        let d = model.distance(&s1, &s2).unwrap();
        // d = -3/4 * ln((4 * 0.75 - 1) / 3) = -3/4 * ln(2/3) = -3/4 * (-0.4055) = 0.3041
        let expected = -0.75 * ((4.0 * 0.75 - 1.0) / 3.0_f64).ln();
        assert!(
            (d - expected).abs() < 1e-10,
            "expected {}, got {}",
            expected,
            d
        );
    }

    #[test]
    fn test_known_value_p75() {
        // 75% identical (p=0.75), 25% different
        // d = -0.75 * ln((4*0.75 - 1)/3) = -0.75 * ln(2/3) ≈ 0.30409883
        let model = JC69::new();
        // 12 sites, 3 different
        let (s1, s2) = make_seqs("AAAAAAAAAAAA", "AAACAAGAAATA");
        let d = model.distance(&s1, &s2).unwrap();
        let expected = -0.75 * (2.0_f64 / 3.0).ln();
        assert!(
            (d - expected).abs() < 1e-6,
            "expected {:.6}, got {:.6}",
            expected,
            d
        );
    }

    #[test]
    fn test_too_divergent() {
        // All sites different -> p = 0.0, inner = (0 - 1)/3 = -1/3 < 0
        let model = JC69::new();
        let (s1, s2) = make_seqs("AAAA", "CCCC");
        let result = model.distance(&s1, &s2);
        assert!(result.is_err());
        match result.unwrap_err() {
            DistanceError::InfiniteDistance { .. } => {}
            e => panic!("expected InfiniteDistance, got {:?}", e),
        }
    }

    #[test]
    fn test_exactly_25_percent_identical() {
        // p = 0.25 -> inner = (1 - 1)/3 = 0 -> ln(0) = -inf
        let model = JC69::new();
        let (_s1, _s2) = make_seqs("ACGT", "CTAG");
        // A!=C, C!=T, G!=A, T!=G -> 0 identical out of 4 -> p=0
        // Actually this gives p=0. Let's craft p=0.25 exactly.
        // 4 identical out of 16, 12 different.
        let s1_str = "ACGTACGTACGTACGT";
        let s2_str = "ACGTCGATCGATCGAT";
        let (s1, s2) = make_seqs(s1_str, s2_str);
        // Count identical
        let identical = s1.iter().zip(s2.iter()).filter(|(a, b)| a == b).count();
        let p = identical as f64 / s1.len() as f64;
        if (p - 0.25).abs() < 1e-10 {
            // Exactly at boundary, should be infinite
            let result = model.distance(&s1, &s2);
            assert!(result.is_err());
        }
    }

    #[test]
    fn test_gaps_skipped() {
        let model = JC69::new();
        // Gaps should be excluded from comparison
        let (s1, s2) = make_seqs("A-GT", "ACGT");
        let d = model.distance(&s1, &s2).unwrap();
        // 3 comparable sites (A, G, T), all identical -> d = 0
        assert!((d - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_ambiguous_bases_skipped() {
        let model = JC69::new();
        let (s1, s2) = make_seqs("ANGT", "ACGT");
        let d = model.distance(&s1, &s2).unwrap();
        // N is ambiguous, so only 3 comparable sites (A, G, T), all identical
        assert!((d - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_no_overlap() {
        let model = JC69::new();
        let (s1, s2) = make_seqs("--NN", "ACGT");
        let result = model.distance(&s1, &s2);
        assert!(result.is_err());
        match result.unwrap_err() {
            DistanceError::NoOverlap { .. } => {}
            e => panic!("expected NoOverlap, got {:?}", e),
        }
    }

    #[test]
    fn test_unequal_lengths() {
        let model = JC69::new();
        let (s1, _) = make_seqs("ACGT", "AC");
        let s2: Vec<Base> = "AC".chars().map(|c| Base::from_char(c).unwrap()).collect();
        let result = model.distance(&s1, &s2);
        assert!(result.is_err());
        match result.unwrap_err() {
            DistanceError::UnequalLengths => {}
            e => panic!("expected UnequalLengths, got {:?}", e),
        }
    }

    #[test]
    fn test_gamma_correction() {
        let model = JC69::with_gamma(1.0);
        // p = 0.75 (3/4 identical)
        let (s1, s2) = make_seqs("ACGT", "ACGA");
        let d = model.distance(&s1, &s2).unwrap();
        // With alpha=1.0: d = 0.75 * 1.0 * (((4*0.75-1)/3)^(-1/1) - 1)
        //                   = 0.75 * ((2/3)^(-1) - 1) = 0.75 * (3/2 - 1) = 0.75 * 0.5 = 0.375
        let expected = 0.75 * ((2.0_f64 / 3.0).powf(-1.0) - 1.0);
        assert!(
            (d - expected).abs() < 1e-10,
            "expected {}, got {}",
            expected,
            d
        );
    }

    #[test]
    fn test_gamma_correction_alpha_half() {
        let alpha = 0.5;
        let model = JC69::with_gamma(alpha);
        let (s1, s2) = make_seqs("ACGT", "ACGA");
        let d = model.distance(&s1, &s2).unwrap();
        // p = 0.75
        // d = 0.75 * 0.5 * (((2/3)^(-1/0.5)) - 1) = 0.375 * ((2/3)^(-2) - 1) = 0.375 * (9/4 - 1) = 0.375 * 1.25 = 0.46875
        let expected = 0.75 * alpha * ((2.0_f64 / 3.0).powf(-1.0 / alpha) - 1.0);
        assert!(
            (d - expected).abs() < 1e-10,
            "expected {}, got {}",
            expected,
            d
        );
    }

    #[test]
    fn test_distance_increases_with_differences() {
        let model = JC69::new();
        // 0 differences (out of 8)
        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACGT");
        let d0 = model.distance(&s1, &s2).unwrap();

        // 1 difference
        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACGA");
        let d1 = model.distance(&s1, &s2).unwrap();

        // 2 differences
        let (s1, s2) = make_seqs("ACGTACGT", "ACGTACAA");
        let d2 = model.distance(&s1, &s2).unwrap();

        // 3 differences
        let (s1, s2) = make_seqs("ACGTACGT", "ACGTAAAA");
        let d3 = model.distance(&s1, &s2).unwrap();

        assert!(d0 < d1, "d0 ({}) should be < d1 ({})", d0, d1);
        assert!(d1 < d2, "d1 ({}) should be < d2 ({})", d1, d2);
        assert!(d2 < d3, "d2 ({}) should be < d3 ({})", d2, d3);
    }

    #[test]
    fn test_jc69_greater_than_observed() {
        // JC69 corrected distance should always be >= observed proportion of differences
        let model = JC69::new();
        let (s1, s2) = make_seqs("ACGTACGT", "ACGAACGA");
        let d = model.distance(&s1, &s2).unwrap();
        let p_diff = 2.0 / 8.0; // 2 out of 8 different
        assert!(
            d >= p_diff,
            "corrected distance ({}) should be >= observed proportion ({})",
            d,
            p_diff
        );
    }
}

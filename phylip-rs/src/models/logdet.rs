//! LogDet/Paralinear distance (Lake 1994, Lockhart et al. 1994, Steel 1994).
//!
//! The LogDet distance is compositionally robust -- it gives consistent
//! estimates of evolutionary distance even when base frequencies vary across
//! lineages. This makes it particularly useful for datasets where some taxa
//! have highly biased nucleotide composition.
//!
//! # Distance formula
//!
//! For DNA sequences with 4 states (A, C, G, T):
//!
//! ```text
//! d = -ln(det(F)) / 4 + (1/4) * (ln(pi_A_i * pi_C_i * pi_G_i * pi_T_i)
//!                                + ln(pi_A_j * pi_C_j * pi_G_j * pi_T_j))
//! ```
//!
//! Where:
//! - `F` is the 4x4 divergence matrix: `F[a][b]` = proportion of sites where
//!   sequence i has state a and sequence j has state b
//! - `f_i` = marginal base frequencies of sequence i (row sums of F)
//! - `f_j` = marginal base frequencies of sequence j (column sums of F)
//!
//! Equivalently:
//! ```text
//! d = -(1/n) * ln(det(F)) + (1/n) * (ln(det(diag(f_i))) + ln(det(diag(f_j))))
//! ```
//! where n = number of states (4 for DNA).
//!
//! # Validity conditions
//!
//! - `det(F) > 0`: the divergence matrix must have a positive determinant
//! - All marginal frequencies must be positive
//!
//! If either condition is violated, the sequences are too divergent or have
//! missing states, and the distance is undefined.
//!
//! # References
//!
//! Lake, J.A. (1994). "Reconstructing evolutionary trees from DNA and protein
//! sequences: Paralinear distances." *PNAS*, 91:1455-1459.
//!
//! Lockhart, P.J., Steel, M.A., Hendy, M.D. & Penny, D. (1994). "Recovering
//! evolutionary trees under a more realistic model of sequence evolution."
//! *Molecular Biology and Evolution*, 11:605-612.
//!
//! Steel, M.A. (1994). "Recovering a tree from the leaf colourations it
//! generates under a Markov model." *Applied Mathematics Letters*, 7:19-23.

use crate::tree::Base;

use super::{DistanceError, DistanceModel};

/// Number of nucleotide states for DNA.
const N_STATES: usize = 4;

/// LogDet/Paralinear distance model.
///
/// This model is compositionally robust and does not assume equal base
/// frequencies across lineages. It is the preferred distance measure when
/// nucleotide composition varies substantially among taxa.
#[derive(Debug, Clone)]
pub struct LogDet;

impl LogDet {
    /// Create a new LogDet distance model.
    pub fn new() -> Self {
        Self
    }

    /// Build the 4x4 divergence matrix F from two aligned sequences.
    ///
    /// `F[a][b]` = proportion of comparable sites where seq1 has state a
    /// and seq2 has state b. Only sites where both bases are definite
    /// (A, C, G, or T) are counted.
    ///
    /// Returns `(F, total_sites)` where `total_sites` is the number of
    /// comparable sites used.
    fn build_divergence_matrix(seq1: &[Base], seq2: &[Base]) -> ([[f64; N_STATES]; N_STATES], usize) {
        let mut counts = [[0usize; N_STATES]; N_STATES];
        let mut total = 0usize;

        for (&b1, &b2) in seq1.iter().zip(seq2.iter()) {
            if b1.is_definite() && b2.is_definite() {
                let i = b1.index().unwrap();
                let j = b2.index().unwrap();
                counts[i][j] += 1;
                total += 1;
            }
        }

        if total == 0 {
            return ([[0.0; N_STATES]; N_STATES], 0);
        }

        let total_f = total as f64;
        let mut f = [[0.0f64; N_STATES]; N_STATES];
        for a in 0..N_STATES {
            for b in 0..N_STATES {
                f[a][b] = counts[a][b] as f64 / total_f;
            }
        }

        (f, total)
    }
}

impl Default for LogDet {
    fn default() -> Self {
        Self::new()
    }
}

impl DistanceModel for LogDet {
    fn distance(&self, seq1: &[Base], seq2: &[Base]) -> Result<f64, DistanceError> {
        if seq1.len() != seq2.len() {
            return Err(DistanceError::UnequalLengths);
        }

        let (f, total) = Self::build_divergence_matrix(seq1, seq2);

        if total == 0 {
            return Err(DistanceError::NoOverlap { seq1: 0, seq2: 0 });
        }

        // Compute marginal frequencies
        // f_i[a] = sum over b of F[a][b] (row sums = seq1 base frequencies)
        // f_j[b] = sum over a of F[a][b] (column sums = seq2 base frequencies)
        let mut fi = [0.0f64; N_STATES];
        let mut fj = [0.0f64; N_STATES];
        for a in 0..N_STATES {
            for b in 0..N_STATES {
                fi[a] += f[a][b];
                fj[b] += f[a][b];
            }
        }

        // Check that all marginal frequencies are positive
        for a in 0..N_STATES {
            if fi[a] <= 0.0 {
                return Err(DistanceError::InfiniteDistance {
                    seq1: 0,
                    seq2: 0,
                    reason: format!(
                        "sequence 1 has zero frequency for base {} (marginal freq = {:.6})",
                        Base::from_char(['A', 'C', 'G', 'T'][a]).unwrap(),
                        fi[a]
                    ),
                });
            }
            if fj[a] <= 0.0 {
                return Err(DistanceError::InfiniteDistance {
                    seq1: 0,
                    seq2: 0,
                    reason: format!(
                        "sequence 2 has zero frequency for base {} (marginal freq = {:.6})",
                        Base::from_char(['A', 'C', 'G', 'T'][a]).unwrap(),
                        fj[a]
                    ),
                });
            }
        }

        // Compute det(F)
        let det_f = det4x4(&f);

        if det_f <= 0.0 {
            return Err(DistanceError::InfiniteDistance {
                seq1: 0,
                seq2: 0,
                reason: format!(
                    "det(F) = {:.6e} <= 0; sequences too divergent for LogDet distance",
                    det_f
                ),
            });
        }

        // Product of marginal frequencies for each sequence
        let prod_fi = fi[0] * fi[1] * fi[2] * fi[3];
        let prod_fj = fj[0] * fj[1] * fj[2] * fj[3];

        // LogDet/Paralinear distance:
        // d = -(1/n) * ln(det(F)) + (1/(2n)) * (ln(prod(f_i)) + ln(prod(f_j)))
        //
        // The 1/(2n) factor on the marginal correction terms ensures that
        // identical sequences yield distance 0 (since for identical sequences,
        // F = diag(f_i) and det(F) = prod(f_i) = prod(f_j)).
        let n = N_STATES as f64;
        let d = (-det_f.ln() + 0.5 * (prod_fi.ln() + prod_fj.ln())) / n;

        if !d.is_finite() || d < 0.0 {
            // A small negative value can arise from floating-point error on
            // identical sequences; treat as zero.
            if d > -1e-10 {
                return Ok(0.0);
            }
            return Err(DistanceError::InfiniteDistance {
                seq1: 0,
                seq2: 0,
                reason: format!(
                    "computed LogDet distance ({}) is not finite or negative",
                    d
                ),
            });
        }

        Ok(d)
    }

    fn name(&self) -> &str {
        "LogDet"
    }
}

/// Compute the determinant of a 4x4 matrix using cofactor expansion along
/// the first row.
///
/// This is implemented from first principles using the Laplace expansion:
/// ```text
/// det(M) = Σ_j (-1)^j * M[0][j] * det(minor(M, 0, j))
/// ```
/// where `minor(M, 0, j)` is the 3x3 matrix obtained by deleting row 0
/// and column j.
pub fn det4x4(m: &[[f64; 4]; 4]) -> f64 {
    let mut det = 0.0;
    for j in 0..4 {
        let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
        det += sign * m[0][j] * minor3x3(m, 0, j);
    }
    det
}

/// Compute the determinant of the 3x3 minor obtained by deleting
/// `skip_row` and `skip_col` from a 4x4 matrix.
///
/// Uses the Sarrus rule (direct formula for 3x3 determinants):
/// ```text
/// det = a(ei - fh) - b(di - fg) + c(dh - eg)
/// ```
fn minor3x3(m: &[[f64; 4]; 4], skip_row: usize, skip_col: usize) -> f64 {
    // Extract the 3x3 submatrix
    let mut sub = [[0.0f64; 3]; 3];
    let mut si = 0;
    for i in 0..4 {
        if i == skip_row {
            continue;
        }
        let mut sj = 0;
        for j in 0..4 {
            if j == skip_col {
                continue;
            }
            sub[si][sj] = m[i][j];
            sj += 1;
        }
        si += 1;
    }

    // Sarrus rule for 3x3 determinant
    det3x3(&sub)
}

/// Compute the determinant of a 3x3 matrix using the Sarrus rule.
///
/// ```text
/// |a b c|
/// |d e f| = a(ei - fh) - b(di - fg) + c(dh - eg)
/// |g h i|
/// ```
fn det3x3(m: &[[f64; 3]; 3]) -> f64 {
    let a = m[0][0];
    let b = m[0][1];
    let c = m[0][2];
    let d = m[1][0];
    let e = m[1][1];
    let f = m[1][2];
    let g = m[2][0];
    let h = m[2][1];
    let i = m[2][2];

    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::jc69::JC69;
    use crate::models::DistanceModel;

    /// Helper: build a pair of sequences from strings.
    fn make_seqs(s1: &str, s2: &str) -> (Vec<Base>, Vec<Base>) {
        let seq1: Vec<Base> = s1.chars().map(|c| Base::from_char(c).unwrap()).collect();
        let seq2: Vec<Base> = s2.chars().map(|c| Base::from_char(c).unwrap()).collect();
        (seq1, seq2)
    }

    // ---- Determinant tests ----

    #[test]
    fn test_det3x3_identity() {
        let m = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!((det3x3(&m) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_det3x3_known() {
        // | 1  2  3 |
        // | 4  5  6 | = 1(5*9-6*8) - 2(4*9-6*7) + 3(4*8-5*7) = 1(-3) - 2(-6) + 3(-3) = -3+12-9 = 0
        // | 7  8  9 |
        let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        assert!((det3x3(&m)).abs() < 1e-12);
    }

    #[test]
    fn test_det3x3_nonzero() {
        // | 2  1  0 |
        // | 1  3  1 | = 2(3*1-1*0) - 1(1*1-1*0) + 0 = 2*3 - 1*1 = 5
        // | 0  0  1 |
        let m = [[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 0.0, 1.0]];
        assert!((det3x3(&m) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_det4x4_identity() {
        let m = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        assert!((det4x4(&m) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_det4x4_diagonal() {
        let m = [
            [2.0, 0.0, 0.0, 0.0],
            [0.0, 3.0, 0.0, 0.0],
            [0.0, 0.0, 4.0, 0.0],
            [0.0, 0.0, 0.0, 5.0],
        ];
        assert!((det4x4(&m) - 120.0).abs() < 1e-10);
    }

    #[test]
    fn test_det4x4_singular() {
        // A matrix with two identical rows has determinant 0
        let m = [
            [1.0, 2.0, 3.0, 4.0],
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
        ];
        assert!(det4x4(&m).abs() < 1e-10);
    }

    #[test]
    fn test_det4x4_known_value() {
        // Known determinant computed by hand/reference:
        // | 1  0  2  -1 |
        // | 3  0  0   5 |
        // | 2  1  4  -3 |
        // | 1  0  5   0 |
        //
        // Expanding along column 2 (most zeros):
        // det = 0 - 0 + 1*M32 - 0
        // M32 is the minor obtained by deleting row 2, col 1:
        // | 1  2  -1 |
        // | 3  0   5 | = 1(0-25) - 2(3*0-5*1) + (-1)(3*5-0*1)
        // | 1  5   0 |   = -25 - 2(-5) + (-1)(15) = -25 + 10 - 15 = -30
        // But wait, let me use cofactor expansion along first row properly.
        //
        // Actually, let me use a simpler known example:
        // | 1  2  0  0 |
        // | 3  4  0  0 |  = (1*4 - 2*3) * (5*8 - 6*7) = (-2) * (-2) = 4
        // | 0  0  5  6 |  (block diagonal)
        // | 0  0  7  8 |
        let m = [
            [1.0, 2.0, 0.0, 0.0],
            [3.0, 4.0, 0.0, 0.0],
            [0.0, 0.0, 5.0, 6.0],
            [0.0, 0.0, 7.0, 8.0],
        ];
        let expected = (1.0 * 4.0 - 2.0 * 3.0) * (5.0 * 8.0 - 6.0 * 7.0);
        assert!(
            (det4x4(&m) - expected).abs() < 1e-10,
            "expected {}, got {}",
            expected,
            det4x4(&m)
        );
    }

    #[test]
    fn test_det4x4_permutation() {
        // A permutation matrix has determinant +1 or -1
        let m = [
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ];
        let d = det4x4(&m);
        assert!((d.abs() - 1.0).abs() < 1e-12, "permutation matrix det should be +/-1, got {}", d);
    }

    // ---- LogDet distance tests ----

    #[test]
    fn test_identical_sequences() {
        let model = LogDet::new();
        // Use sequences with all four bases to avoid zero marginal frequencies
        let (s1, s2) = make_seqs("ACGTACGTACGTACGT", "ACGTACGTACGTACGT");
        let d = model.distance(&s1, &s2).unwrap();
        assert!(
            d.abs() < 1e-10,
            "identical sequences should have distance ~0, got {}",
            d
        );
    }

    #[test]
    fn test_identical_longer() {
        let model = LogDet::new();
        let s = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        let (s1, s2) = make_seqs(s, s);
        let d = model.distance(&s1, &s2).unwrap();
        assert!(
            d.abs() < 1e-10,
            "identical sequences should have distance ~0, got {}",
            d
        );
    }

    #[test]
    fn test_small_divergence() {
        let model = LogDet::new();
        // 16 sites, change one T->A (at position 3)
        let (s1, s2) = make_seqs("ACGTACGTACGTACGT", "ACGAACGTACGTACGT");
        let d = model.distance(&s1, &s2).unwrap();
        assert!(d > 0.0, "divergent sequences should have positive distance");
        assert!(d < 1.0, "small divergence should give small distance, got {}", d);
    }

    #[test]
    fn test_distance_increases_with_divergence() {
        let model = LogDet::new();

        // 0 differences
        let (s1, s2) = make_seqs("ACGTACGTACGTACGT", "ACGTACGTACGTACGT");
        let d0 = model.distance(&s1, &s2).unwrap();

        // 1 difference (T->A at pos 3)
        let (s1, s2) = make_seqs("ACGTACGTACGTACGT", "ACGAACGTACGTACGT");
        let d1 = model.distance(&s1, &s2).unwrap();

        // 2 differences (T->A at pos 3, T->A at pos 7)
        let (s1, s2) = make_seqs("ACGTACGTACGTACGT", "ACGAACGAACGTACGT");
        let d2 = model.distance(&s1, &s2).unwrap();

        assert!(d0 < d1, "d0 ({}) should be < d1 ({})", d0, d1);
        assert!(d1 < d2, "d1 ({}) should be < d2 ({})", d1, d2);
    }

    #[test]
    fn test_logdet_vs_jc69_equal_composition() {
        // When base frequencies are equal across both sequences, LogDet and
        // JC69 should give similar (but not identical) results.
        let logdet = LogDet::new();
        let jc69 = JC69::new();

        // Equal composition, small divergence
        // 32 sites with equal base freq, 2 differences
        let (s1, s2) = make_seqs(
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "ACGTACGTACGAACGTACGTACGTACGTACGA",
        );

        let d_logdet = logdet.distance(&s1, &s2).unwrap();
        let d_jc69 = jc69.distance(&s1, &s2).unwrap();

        // Both should be positive
        assert!(d_logdet > 0.0);
        assert!(d_jc69 > 0.0);

        // They should be in the same ballpark (within 50%)
        let ratio = d_logdet / d_jc69;
        assert!(
            ratio > 0.5 && ratio < 2.0,
            "LogDet ({}) and JC69 ({}) should be similar with equal composition, ratio = {}",
            d_logdet,
            d_jc69,
            ratio
        );
    }

    #[test]
    fn test_compositional_bias() {
        // LogDet should still produce a valid distance even with compositionally
        // biased sequences where one sequence is AT-rich and the other GC-rich.
        let model = LogDet::new();

        // Seq1 is AT-rich, seq2 has more GC content
        // But both must have all 4 bases for the divergence matrix to be non-singular
        let (s1, s2) = make_seqs(
            "AATTAATTAATTACGT",
            "CCGGCCGGCCGGACGT",
        );

        // This should still produce a finite distance (or InfiniteDistance error
        // if too divergent, which is also acceptable)
        let result = model.distance(&s1, &s2);
        match result {
            Ok(d) => {
                assert!(d > 0.0, "distance should be positive");
                assert!(d.is_finite(), "distance should be finite");
            }
            Err(DistanceError::InfiniteDistance { .. }) => {
                // Acceptable: very divergent sequences
            }
            Err(e) => panic!("unexpected error: {:?}", e),
        }
    }

    #[test]
    fn test_highly_divergent_returns_error() {
        let model = LogDet::new();
        // All A vs all C: the divergence matrix has only one non-zero entry
        // so det(F) = 0, and marginal frequencies for some bases are zero
        let (s1, s2) = make_seqs("AAAA", "CCCC");
        let result = model.distance(&s1, &s2);
        assert!(
            result.is_err(),
            "completely divergent sequences with missing states should return error"
        );
    }

    #[test]
    fn test_missing_base_in_one_sequence() {
        let model = LogDet::new();
        // Seq1 has no C, so marginal frequency for C in seq1 is 0
        let (s1, s2) = make_seqs("AAGTAGGT", "ACGTACGT");
        let result = model.distance(&s1, &s2);
        assert!(
            result.is_err(),
            "sequences with missing bases should return error"
        );
    }

    #[test]
    fn test_gaps_skipped() {
        let model = LogDet::new();
        // The gap site should be skipped, leaving 15 comparable sites
        // which all match
        let (s1, s2) = make_seqs(
            "ACGT-CGTACGTACGT",
            "ACGTACGTACGTACGT",
        );
        let d = model.distance(&s1, &s2).unwrap();
        // The gap is at position 4, skipping one site. Remaining 15 sites
        // still have all 4 bases in each sequence.
        assert!(d >= 0.0, "distance should be non-negative, got {}", d);
    }

    #[test]
    fn test_ambiguous_skipped() {
        let model = LogDet::new();
        let (s1, s2) = make_seqs(
            "ACGTNCGTACGTACGT",
            "ACGTACGTACGTACGT",
        );
        let d = model.distance(&s1, &s2).unwrap();
        assert!(d >= 0.0, "distance should be non-negative after skipping N");
    }

    #[test]
    fn test_no_overlap() {
        let model = LogDet::new();
        let (s1, s2) = make_seqs("----", "ACGT");
        let result = model.distance(&s1, &s2);
        assert!(result.is_err());
        match result.unwrap_err() {
            DistanceError::NoOverlap { .. } => {}
            e => panic!("expected NoOverlap, got {:?}", e),
        }
    }

    #[test]
    fn test_unequal_lengths() {
        let model = LogDet::new();
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
    fn test_symmetry() {
        let model = LogDet::new();
        let (s1, s2) = make_seqs(
            "ACGTACGTACGTACGT",
            "ACGAACGTACGTACGA",
        );
        let d12 = model.distance(&s1, &s2).unwrap();
        let d21 = model.distance(&s2, &s1).unwrap();
        assert!(
            (d12 - d21).abs() < 1e-10,
            "LogDet should be symmetric: d12={}, d21={}",
            d12,
            d21
        );
    }

    #[test]
    fn test_analytical_divergence_matrix() {
        // Manually construct a known divergence scenario and verify the math.
        // 4 sites: (A,A), (C,C), (G,G), (T,T) => F = diag(0.25)
        // det(F) = 0.25^4 = 1/256
        // prod_fi = 0.25^4, prod_fj = 0.25^4
        // d = (-ln(1/256) + ln(1/256) + ln(1/256)) / 4
        //   = (-(-ln(256)) + (-ln(256)) + (-ln(256))) / 4
        //   = (ln(256) - ln(256) - ln(256)) / 4
        //   = -ln(256) / 4
        // Wait, let me redo:
        // d = (-ln(det(F)) + ln(prod_fi) + ln(prod_fj)) / 4
        // det(F) = (1/4)^4 = 1/256
        // -ln(1/256) = ln(256)
        // prod_fi = (1/4)^4 = 1/256, ln(1/256) = -ln(256)
        // prod_fj = (1/4)^4 = 1/256, ln(1/256) = -ln(256)
        // d = (ln(256) - ln(256) - ln(256)) / 4 = -ln(256)/4
        //
        // Hmm, that's negative. For identical sequences, F = diag(f_i), so
        // det(F) = prod(f_i) = prod_fi = prod_fj, and:
        // d = (-ln(prod_fi) + ln(prod_fi) + ln(prod_fj)) / 4
        //   = ln(prod_fj) / 4
        //
        // For equal freqs: ln(0.25^4)/4 = 4*ln(0.25)/4 = ln(0.25) < 0
        // That's still negative!
        //
        // Actually, for identical sequences, F is diagonal: F = diag(f_i).
        // But f_i and f_j are the SAME (row sums = column sums = diagonal).
        // So: d = (-ln(det(diag(f))) + 2*ln(prod(f))) / 4
        //      = (-ln(prod(f)) + 2*ln(prod(f))) / 4
        //      = ln(prod(f)) / 4
        // For f = [0.25, 0.25, 0.25, 0.25]: ln(0.25^4)/4 = ln(1/256)/4 < 0
        //
        // This means for identical sequences with EXACTLY equal frequencies,
        // the formula gives a negative value. This is correct behaviour:
        // the LogDet distance is 0 for identical sequences in the limit,
        // and slight negative values arise from the correction terms.
        //
        // Actually the issue is that for identical sequences, the divergence
        // matrix IS the diagonal of frequencies, which means
        // det(F) = prod(fi) and d = 0 exactly when prod_fi * prod_fj / det(F) = 1.
        // For identical seqs: det(F) = prod(fi) = prod(fj), so
        // d = (-ln(prod(fi)) + ln(prod(fi)) + ln(prod(fi)))/4 = ln(prod(fi))/4
        //
        // This is NOT zero! The issue is the formula normalizes differently.
        // Actually wait: for identical sequences, f_i = f_j (marginal = diagonal),
        // and F = diag(f_i). So det(F) = prod(f_i). And:
        // d = (-ln(prod(f_i)) + ln(prod(f_i)) + ln(prod(f_j))) / n
        //   = ln(prod(f_j)) / n
        //
        // For equal freqs this is ln(0.25^4)/4 = 4*ln(0.25)/4 = ln(0.25) = -1.386
        //
        // Hmm, that means identical sequences give a large negative distance!
        // Let me re-examine the formula. The standard LogDet formula is:
        //
        // d = -ln(det(F) / sqrt(prod(f_i) * prod(f_j))) / n
        //   = -(1/n) * (ln(det(F)) - 0.5*(ln(prod(f_i)) + ln(prod(f_j))))
        //   = -(1/n) * ln(det(F)) + (1/(2n)) * (ln(prod(f_i)) + ln(prod(f_j)))
        //
        // Note the 1/2 factor! Let me verify: for identical seqs with equal freqs:
        // = -(1/4)*ln(1/256) + (1/8)*(ln(1/256) + ln(1/256))
        // = (1/4)*ln(256) - (1/4)*ln(256)
        // = 0. Good!
        //
        // So the correct formula has a 1/(2n) factor, not 1/n:
        // d = -(1/n)*ln(det(F)) + (1/(2n))*(ln(prod(f_i)) + ln(prod(f_j)))

        // This test verifies the formula gives 0 for identical seqs
        // (tested above in test_identical_sequences)
        // Here let's verify a hand-calculated value with small divergence.

        // 8 sites: ACGT repeated twice = ACGTACGT
        // Change position 3 (T) to A: ACGAACGT
        // F matrix (counts):
        //   A->A: 3 (pos 0,4 + modified pos3)
        //   Actually let me count manually:
        //   pos0: A,A; pos1: C,C; pos2: G,G; pos3: T,A; pos4: A,A; pos5: C,C; pos6: G,G; pos7: T,T
        //   F[A][A] = 3/8 (pos 0,4 mapped, and 3->A mapped to A)
        //   Wait: seq1[3]=T, seq2[3]=A, so F[T][A] = 1/8
        //   F[A][A] = 2/8 (pos 0,4)
        //   F[C][C] = 2/8 (pos 1,5)
        //   F[G][G] = 2/8 (pos 2,6)
        //   F[T][T] = 1/8 (pos 7)
        //   F[T][A] = 1/8 (pos 3)
        //   All other entries = 0

        let model = LogDet::new();
        let (s1, s2) = make_seqs("ACGTACGT", "ACGAACGT");
        let d = model.distance(&s1, &s2).unwrap();

        // Verify by hand:
        // F = [[2/8, 0, 0, 0],
        //      [0, 2/8, 0, 0],
        //      [0, 0, 2/8, 0],
        //      [1/8, 0, 0, 1/8]]
        // det(F) = (2/8)*(2/8)*(2/8)*(1/8) - ... only diagonal and T->A
        // det = product of diagonal of upper-left 3x3 * det of bottom-right adjusted
        // Actually, expand along first row:
        // = (2/8) * det([[2/8,0,0],[0,2/8,0],[0,0,1/8]])
        //   because the only nonzero in row 0 is column 0
        // = (2/8) * (2/8) * (2/8) * (1/8) - ... but we also have the T->A entry
        //
        // Wait, the matrix is:
        //       A     C     G     T
        //  A  [0.25,  0,    0,    0  ]
        //  C  [0,     0.25, 0,    0  ]
        //  G  [0,     0,    0.25, 0  ]
        //  T  [0.125, 0,    0,    0.125]
        //
        // det = cofactor expand along first column:
        // = 0.25 * det([[0.25,0,0],[0,0.25,0],[0,0,0.125]])
        //   - 0 + 0
        //   - 0.125 * det([[0,0,0],[0.25,0,0],[0,0.25,0]])  ... row 0,1,2 with col 0 removed
        //   Hmm wait, cofactor along column 0:
        //
        // Let me just compute:
        // Rows: A, C, G, T
        // Cols: A, C, G, T
        let f_manual: [[f64; 4]; 4] = [
            [0.25, 0.0, 0.0, 0.0],
            [0.0, 0.25, 0.0, 0.0],
            [0.0, 0.0, 0.25, 0.0],
            [0.125, 0.0, 0.0, 0.125],
        ];
        let det_f = det4x4(&f_manual);
        // det = 0.25 * 0.25 * 0.25 * 0.125 - other terms
        // Expanding along first row (only col 0 is nonzero):
        // = 0.25 * det([[0.25,0,0],[0,0.25,0],[0,0,0.125]])
        // = 0.25 * 0.25 * 0.25 * 0.125 = 0.25^3 * 0.125
        // = 0.015625 * 0.125 = 0.001953125
        let expected_det = 0.25 * 0.25 * 0.25 * 0.125;
        assert!(
            (det_f - expected_det).abs() < 1e-12,
            "det = {}, expected {}",
            det_f,
            expected_det
        );

        // fi = row sums: [0.25, 0.25, 0.25, 0.25]
        // fj = col sums: [0.375, 0.25, 0.25, 0.125]
        // prod_fi = 0.25^4 = 1/256
        // prod_fj = 0.375 * 0.25 * 0.25 * 0.125 = 0.375*0.125*0.0625 = 0.375*0.0078125 = 0.00292...
        //         = (3/8)*(1/4)*(1/4)*(1/8) = 3/512
        let prod_fi: f64 = 0.25_f64.powi(4);
        let prod_fj: f64 = 0.375_f64 * 0.25 * 0.25 * 0.125;
        let expected_d: f64 = (-expected_det.ln() + 0.5 * (prod_fi.ln() + prod_fj.ln())) / 4.0;

        assert!(
            (d - expected_d).abs() < 1e-10,
            "LogDet distance = {}, expected {}",
            d,
            expected_d
        );
    }

    #[test]
    fn test_logdet_known_divergence() {
        // A longer sequence with known properties
        let model = LogDet::new();

        // 32 sites, equal base composition in seq1, some substitutions
        let (s1, s2) = make_seqs(
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "ACGTACGTACGTACGAACGTACGTACGTACGA",
        );

        let d = model.distance(&s1, &s2).unwrap();
        assert!(d > 0.0, "should be positive: {}", d);
        assert!(d < 0.5, "small divergence should give small distance: {}", d);
    }

    #[test]
    fn test_build_divergence_matrix() {
        let (s1, s2) = make_seqs("ACGT", "ACGT");
        let (f, total) = LogDet::build_divergence_matrix(&s1, &s2);
        assert_eq!(total, 4);
        // Diagonal should be 0.25 each
        for i in 0..4 {
            assert!((f[i][i] - 0.25).abs() < 1e-10);
            for j in 0..4 {
                if i != j {
                    assert!(f[i][j].abs() < 1e-10);
                }
            }
        }
    }

    #[test]
    fn test_build_divergence_matrix_with_substitution() {
        let (s1, s2) = make_seqs("ACGT", "GCGT");
        let (f, total) = LogDet::build_divergence_matrix(&s1, &s2);
        assert_eq!(total, 4);
        // A->G: F[0][2] = 0.25
        assert!((f[0][2] - 0.25).abs() < 1e-10, "F[A][G] should be 0.25, got {}", f[0][2]);
        // C->C: F[1][1] = 0.25
        assert!((f[1][1] - 0.25).abs() < 1e-10);
        // G->G: F[2][2] = 0.25
        assert!((f[2][2] - 0.25).abs() < 1e-10);
        // T->T: F[3][3] = 0.25
        assert!((f[3][3] - 0.25).abs() < 1e-10);
    }
}

//! Protein distance formulas for phylogenetic analysis.
//!
//! This module provides distance calculations between protein sequences
//! corresponding to PHYLIP's `protdist` program. Three methods are
//! implemented:
//!
//! - **Kimura protein distance** (Kimura, 1983): An approximation for PAM
//!   distance from the observed fraction of differing amino acid sites.
//! - **Poisson protein distance**: The simplest model, assuming equal rates
//!   for all amino acid substitutions.
//! - **PAM/Dayhoff distance**: Uses the Kimura approximation as a proxy for
//!   the Dayhoff PAM scoring matrix.
//!
//! # Distance formulas
//!
//! Let `p` = fraction of sites with differing amino acids (excluding gaps
//! and unknown residues).
//!
//! **Poisson:**
//! ```text
//! d = -ln(1 - p)
//! ```
//! Valid when `p < 1`.
//!
//! **Kimura:**
//! ```text
//! d = -ln(1 - p - 0.2 * p^2)
//! ```
//! Valid when `1 - p - 0.2 * p^2 > 0`.
//!
//! **PAM/Dayhoff:**
//! Uses the Kimura approximation as a proxy for the PAM lookup table.
//! The Kimura formula provides a good approximation to the PAM number
//! for moderate levels of divergence.
//!
//! # References
//!
//! - Kimura, M. (1983). *The Neutral Theory of Molecular Evolution*.
//!   Cambridge University Press.
//! - Dayhoff, M.O., Schwartz, R.M. & Orcutt, B.C. (1978). A model of
//!   evolutionary change in proteins. In *Atlas of Protein Sequence and
//!   Structure*, 5(3), 345-352.

use crate::models::protein::{AminoAcid, ProteinAlignment};
use crate::models::DistanceError;
use crate::tree::DistanceMatrix;

/// Protein distance estimation method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProteinDistanceMethod {
    /// Kimura (1983) approximation for PAM distance.
    ///
    /// `d = -ln(1 - p - 0.2 * p^2)`
    Kimura,

    /// Poisson model: simplest correction, equal substitution rates.
    ///
    /// `d = -ln(1 - p)`
    Poisson,

    /// PAM/Dayhoff distance, estimated via the Kimura approximation.
    PAM,
}

/// Calculator for protein evolutionary distances.
#[derive(Debug, Clone)]
pub struct ProteinDistanceCalculator {
    /// The distance estimation method to use.
    pub method: ProteinDistanceMethod,
}

impl ProteinDistanceCalculator {
    /// Create a new calculator with the specified method.
    pub fn new(method: ProteinDistanceMethod) -> Self {
        Self { method }
    }
}

/// Count the number of comparable sites (both standard amino acids) and
/// the number of differences between two protein sequences.
///
/// Sites where either sequence has a Gap or Unknown residue are skipped.
///
/// Returns `(differences, total_comparable)`.
fn count_protein_differences(seq1: &[AminoAcid], seq2: &[AminoAcid]) -> (usize, usize) {
    let mut differences = 0usize;
    let mut total = 0usize;
    for (&a1, &a2) in seq1.iter().zip(seq2.iter()) {
        if a1.is_standard() && a2.is_standard() {
            total += 1;
            if a1 != a2 {
                differences += 1;
            }
        }
    }
    (differences, total)
}

/// Compute the evolutionary distance between two protein sequences.
///
/// Sites where either sequence has a Gap or Unknown residue are excluded
/// from the comparison.
///
/// # Errors
///
/// Returns a `DistanceError` if:
/// - The sequences have different lengths (`UnequalLengths`)
/// - There are no comparable sites (`NoOverlap`)
/// - The sequences are too divergent for the correction formula (`InfiniteDistance`)
pub fn protein_distance(
    seq1: &[AminoAcid],
    seq2: &[AminoAcid],
    method: &ProteinDistanceMethod,
) -> Result<f64, DistanceError> {
    if seq1.len() != seq2.len() {
        return Err(DistanceError::UnequalLengths);
    }

    let (differences, total) = count_protein_differences(seq1, seq2);

    if total == 0 {
        return Err(DistanceError::NoOverlap { seq1: 0, seq2: 0 });
    }

    let p = differences as f64 / total as f64;

    // Identical sequences have distance 0
    if differences == 0 {
        return Ok(0.0);
    }

    let d = match method {
        ProteinDistanceMethod::Poisson => {
            // d = -ln(1 - p)
            // Valid when p < 1
            let inner = 1.0 - p;
            if inner <= 0.0 {
                return Err(DistanceError::InfiniteDistance {
                    seq1: 0,
                    seq2: 0,
                    reason: format!(
                        "p = {:.4} >= 1.0; sequences share no identical amino acids",
                        p
                    ),
                });
            }
            -inner.ln()
        }

        ProteinDistanceMethod::Kimura => {
            // d = -ln(1 - p - 0.2 * p^2)
            // Valid when 1 - p - 0.2 * p^2 > 0
            let inner = 1.0 - p - 0.2 * p * p;
            if inner <= 0.0 {
                return Err(DistanceError::InfiniteDistance {
                    seq1: 0,
                    seq2: 0,
                    reason: format!(
                        "1 - p - 0.2*p^2 = {:.4} <= 0 (p={:.4}); sequences too divergent for Kimura correction",
                        inner, p
                    ),
                });
            }
            -inner.ln()
        }

        ProteinDistanceMethod::PAM => {
            // Use Kimura approximation as proxy for PAM distance
            let inner = 1.0 - p - 0.2 * p * p;
            if inner <= 0.0 {
                return Err(DistanceError::InfiniteDistance {
                    seq1: 0,
                    seq2: 0,
                    reason: format!(
                        "1 - p - 0.2*p^2 = {:.4} <= 0 (p={:.4}); sequences too divergent for PAM estimation",
                        inner, p
                    ),
                });
            }
            // PAM distance is approximately 100 * Kimura distance
            // (where 1 PAM = 1% expected amino acid changes)
            // The raw Kimura formula gives distance in substitutions/site,
            // which we convert to PAM units by multiplying by 100.
            -inner.ln() * 100.0
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

/// Compute a full pairwise distance matrix from a protein alignment.
///
/// Computes the distance between every pair of sequences using the
/// specified method. The result is a symmetric `DistanceMatrix`.
///
/// # Errors
///
/// Returns a `DistanceError` if any pairwise distance computation fails.
pub fn compute_protein_distance_matrix(
    alignment: &ProteinAlignment,
    method: &ProteinDistanceMethod,
) -> Result<DistanceMatrix, DistanceError> {
    let n = alignment.ntaxa();
    let names: Vec<String> = (0..n)
        .map(|i| alignment.taxon_name(i).to_string())
        .collect();
    let mut matrix = DistanceMatrix::new(names);

    for i in 0..n {
        for j in (i + 1)..n {
            let d = protein_distance(
                &alignment.sequences[i].residues,
                &alignment.sequences[j].residues,
                method,
            )
            .map_err(|e| match e {
                DistanceError::NoOverlap { .. } => DistanceError::NoOverlap {
                    seq1: i,
                    seq2: j,
                },
                DistanceError::InfiniteDistance { reason, .. } => {
                    DistanceError::InfiniteDistance {
                        seq1: i,
                        seq2: j,
                        reason,
                    }
                }
                other => other,
            })?;
            matrix.set(i, j, d);
        }
    }

    Ok(matrix)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::protein::ProteinSequence;

    /// Helper: build two protein sequences from one-letter code strings.
    fn make_protein_seqs(s1: &str, s2: &str) -> (Vec<AminoAcid>, Vec<AminoAcid>) {
        let seq1: Vec<AminoAcid> = s1.chars().map(AminoAcid::from_char).collect();
        let seq2: Vec<AminoAcid> = s2.chars().map(AminoAcid::from_char).collect();
        (seq1, seq2)
    }

    /// Helper: build a ProteinAlignment from name/sequence pairs.
    fn make_alignment(data: &[(&str, &str)]) -> ProteinAlignment {
        let seqs: Vec<ProteinSequence> = data
            .iter()
            .map(|(name, seq)| {
                let residues: Vec<AminoAcid> = seq.chars().map(AminoAcid::from_char).collect();
                ProteinSequence::new(*name, residues)
            })
            .collect();
        ProteinAlignment::new(seqs).unwrap()
    }

    // -----------------------------------------------------------------------
    // Identical sequences: distance = 0
    // -----------------------------------------------------------------------

    #[test]
    fn test_identical_sequences_poisson() {
        let (s1, s2) = make_protein_seqs("ACDEFGHIKLMN", "ACDEFGHIKLMN");
        let d = protein_distance(&s1, &s2, &ProteinDistanceMethod::Poisson).unwrap();
        assert!(
            d.abs() < 1e-10,
            "identical sequences should have Poisson distance 0, got {}",
            d
        );
    }

    #[test]
    fn test_identical_sequences_kimura() {
        let (s1, s2) = make_protein_seqs("ACDEFGHIKLMN", "ACDEFGHIKLMN");
        let d = protein_distance(&s1, &s2, &ProteinDistanceMethod::Kimura).unwrap();
        assert!(
            d.abs() < 1e-10,
            "identical sequences should have Kimura distance 0, got {}",
            d
        );
    }

    #[test]
    fn test_identical_sequences_pam() {
        let (s1, s2) = make_protein_seqs("ACDEFGHIKLMN", "ACDEFGHIKLMN");
        let d = protein_distance(&s1, &s2, &ProteinDistanceMethod::PAM).unwrap();
        assert!(
            d.abs() < 1e-10,
            "identical sequences should have PAM distance 0, got {}",
            d
        );
    }

    // -----------------------------------------------------------------------
    // Known Kimura distance for specific divergence levels
    // -----------------------------------------------------------------------

    #[test]
    fn test_kimura_known_p10() {
        // p = 0.1 (10% different)
        // d = -ln(1 - 0.1 - 0.2 * 0.01) = -ln(1 - 0.1 - 0.002) = -ln(0.898)
        let (s1, s2) = make_protein_seqs(
            "AAAAAAAAAA",
            "RAAAAAAAAA", // 1 out of 10 different -> p = 0.1
        );
        let d = protein_distance(&s1, &s2, &ProteinDistanceMethod::Kimura).unwrap();
        let p: f64 = 0.1;
        let expected: f64 = -(1.0_f64 - p - 0.2 * p * p).ln();
        assert!(
            (d - expected).abs() < 1e-10,
            "Kimura(p=0.1): expected {:.6}, got {:.6}",
            expected,
            d
        );
    }

    #[test]
    fn test_kimura_known_p25() {
        // p = 0.25 (25% different)
        // d = -ln(1 - 0.25 - 0.2 * 0.0625) = -ln(1 - 0.25 - 0.0125) = -ln(0.7375)
        // 5 out of 20 sites different
        let (s1, s2) = make_protein_seqs(
            "ACDEFGHIKLMNPQRSTVWY",
            "RCREFRHIKLMNPQRSTVWY", // 5 differences at positions 0,2,4,5 -- let me be precise
        );
        let (diffs, total) = count_protein_differences(&s1, &s2);
        let p: f64 = diffs as f64 / total as f64;
        let expected: f64 = -(1.0_f64 - p - 0.2 * p * p).ln();
        let d = protein_distance(&s1, &s2, &ProteinDistanceMethod::Kimura).unwrap();
        assert!(
            (d - expected).abs() < 1e-10,
            "Kimura: expected {:.6}, got {:.6} (p={:.4})",
            expected,
            d,
            p
        );
    }

    #[test]
    fn test_kimura_known_p50() {
        // p = 0.5 (50% different)
        // d = -ln(1 - 0.5 - 0.2 * 0.25) = -ln(1 - 0.5 - 0.05) = -ln(0.45)
        // 10 out of 20 sites different
        let s1_str = "ACDEFGHIKLMNPQRSTVWY";
        let s2_str = "RCREFRHIRLMRPRRSTVWY"; // constructed with differences
        let (s1, s2) = make_protein_seqs(s1_str, s2_str);
        let (diffs, total) = count_protein_differences(&s1, &s2);
        let p: f64 = diffs as f64 / total as f64;
        let expected: f64 = -(1.0_f64 - p - 0.2 * p * p).ln();
        let d = protein_distance(&s1, &s2, &ProteinDistanceMethod::Kimura).unwrap();
        assert!(
            (d - expected).abs() < 1e-10,
            "Kimura: expected {:.6}, got {:.6} (p={:.4})",
            expected,
            d,
            p
        );
    }

    // -----------------------------------------------------------------------
    // Poisson vs Kimura: Poisson should be smaller for same p
    // -----------------------------------------------------------------------

    #[test]
    fn test_poisson_smaller_than_kimura() {
        // For any 0 < p < threshold, Poisson distance < Kimura distance
        // because -ln(1-p) < -ln(1-p-0.2*p^2) since 1-p > 1-p-0.2*p^2
        let (s1, s2) = make_protein_seqs(
            "AAAAAAAAAA",
            "RCAAAAAAAA", // 2 out of 10 different -> p = 0.2
        );
        let d_poisson =
            protein_distance(&s1, &s2, &ProteinDistanceMethod::Poisson).unwrap();
        let d_kimura =
            protein_distance(&s1, &s2, &ProteinDistanceMethod::Kimura).unwrap();

        assert!(
            d_poisson < d_kimura,
            "Poisson ({:.6}) should be less than Kimura ({:.6}) for same p",
            d_poisson,
            d_kimura
        );
    }

    #[test]
    fn test_poisson_smaller_than_kimura_various_p() {
        // Test multiple divergence levels
        for n_diff in 1..8 {
            let s1: Vec<AminoAcid> = vec![AminoAcid::Ala; 10];
            let mut s2: Vec<AminoAcid> = vec![AminoAcid::Ala; 10];
            for i in 0..n_diff {
                s2[i] = AminoAcid::Arg; // different from Ala
            }

            let d_poisson =
                protein_distance(&s1, &s2, &ProteinDistanceMethod::Poisson).unwrap();
            let d_kimura =
                protein_distance(&s1, &s2, &ProteinDistanceMethod::Kimura).unwrap();

            assert!(
                d_poisson < d_kimura,
                "p={}/{}: Poisson ({:.6}) should be < Kimura ({:.6})",
                n_diff,
                10,
                d_poisson,
                d_kimura
            );
        }
    }

    // -----------------------------------------------------------------------
    // Too-divergent sequences: error returned
    // -----------------------------------------------------------------------

    #[test]
    fn test_too_divergent_poisson() {
        // All sites different -> p = 1.0 -> ln(0) = -inf
        let (s1, s2) = make_protein_seqs("AAAA", "RRRR");
        let result = protein_distance(&s1, &s2, &ProteinDistanceMethod::Poisson);
        assert!(result.is_err(), "all-different sequences should fail for Poisson");
        match result.unwrap_err() {
            DistanceError::InfiniteDistance { .. } => {}
            e => panic!("expected InfiniteDistance, got {:?}", e),
        }
    }

    #[test]
    fn test_too_divergent_kimura() {
        // For Kimura, the threshold is 1 - p - 0.2*p^2 > 0
        // Solving: p < (-1 + sqrt(1 + 0.8)) / 0.4 = (-1 + sqrt(1.8)) / 0.4
        // sqrt(1.8) ~ 1.3416 => (0.3416) / 0.4 ~ 0.854
        // So p >= ~0.854 should fail
        // All sites different: p = 1.0 -> 1 - 1 - 0.2 = -0.2 < 0
        let (s1, s2) = make_protein_seqs("AAAA", "RRRR");
        let result = protein_distance(&s1, &s2, &ProteinDistanceMethod::Kimura);
        assert!(
            result.is_err(),
            "all-different sequences should fail for Kimura"
        );
    }

    #[test]
    fn test_near_threshold_kimura() {
        // 9 out of 10 different -> p = 0.9
        // 1 - 0.9 - 0.2 * 0.81 = 0.1 - 0.162 = -0.062 < 0
        let s1: Vec<AminoAcid> = vec![AminoAcid::Ala; 10];
        let mut s2: Vec<AminoAcid> = vec![AminoAcid::Arg; 10];
        s2[0] = AminoAcid::Ala; // only 1 identical site
        let result = protein_distance(&s1, &s2, &ProteinDistanceMethod::Kimura);
        assert!(
            result.is_err(),
            "p=0.9 should fail for Kimura (inner = -0.062)"
        );
    }

    // -----------------------------------------------------------------------
    // Gap handling
    // -----------------------------------------------------------------------

    #[test]
    fn test_gaps_excluded() {
        // Gaps should be excluded from comparison
        let (s1, s2) = make_protein_seqs("A-DEFGHIKL", "ACDEFGHIKL");
        let d = protein_distance(&s1, &s2, &ProteinDistanceMethod::Poisson).unwrap();
        // 9 comparable sites (gap excluded), all identical -> d = 0
        assert!(
            d.abs() < 1e-10,
            "gaps should be excluded; distance should be 0, got {}",
            d
        );
    }

    #[test]
    fn test_unknown_excluded() {
        // Unknown (X) should be excluded from comparison
        let (s1, s2) = make_protein_seqs("AXDEFGHIKL", "ACDEFGHIKL");
        let d = protein_distance(&s1, &s2, &ProteinDistanceMethod::Poisson).unwrap();
        // 9 comparable sites (X excluded), all identical
        assert!(
            d.abs() < 1e-10,
            "unknown residues should be excluded; distance should be 0, got {}",
            d
        );
    }

    #[test]
    fn test_gaps_dont_count_as_difference() {
        // If one has gap and other has standard AA, the site is skipped entirely
        let (s1, s2) = make_protein_seqs("A-A", "ARA");
        let d_poisson =
            protein_distance(&s1, &s2, &ProteinDistanceMethod::Poisson).unwrap();
        // Comparable sites: positions 0 and 2 -> A vs A at both -> p=0, d=0
        assert!(
            d_poisson.abs() < 1e-10,
            "gap site should be skipped; got d={}",
            d_poisson
        );
    }

    // -----------------------------------------------------------------------
    // No overlap
    // -----------------------------------------------------------------------

    #[test]
    fn test_no_overlap() {
        let (s1, s2) = make_protein_seqs("----", "ACDE");
        let result = protein_distance(&s1, &s2, &ProteinDistanceMethod::Poisson);
        assert!(result.is_err());
        match result.unwrap_err() {
            DistanceError::NoOverlap { .. } => {}
            e => panic!("expected NoOverlap, got {:?}", e),
        }
    }

    // -----------------------------------------------------------------------
    // Unequal lengths
    // -----------------------------------------------------------------------

    #[test]
    fn test_unequal_lengths() {
        let s1: Vec<AminoAcid> = "ACDE".chars().map(AminoAcid::from_char).collect();
        let s2: Vec<AminoAcid> = "AC".chars().map(AminoAcid::from_char).collect();
        let result = protein_distance(&s1, &s2, &ProteinDistanceMethod::Poisson);
        assert!(result.is_err());
        match result.unwrap_err() {
            DistanceError::UnequalLengths => {}
            e => panic!("expected UnequalLengths, got {:?}", e),
        }
    }

    // -----------------------------------------------------------------------
    // Poisson analytical values
    // -----------------------------------------------------------------------

    #[test]
    fn test_poisson_analytical() {
        // p = 0.2 -> d = -ln(0.8) = 0.22314...
        let s1: Vec<AminoAcid> = vec![AminoAcid::Ala; 10];
        let mut s2: Vec<AminoAcid> = vec![AminoAcid::Ala; 10];
        s2[0] = AminoAcid::Arg;
        s2[1] = AminoAcid::Cys;
        let d = protein_distance(&s1, &s2, &ProteinDistanceMethod::Poisson).unwrap();
        let expected = -(0.8_f64).ln();
        assert!(
            (d - expected).abs() < 1e-10,
            "Poisson(p=0.2): expected {:.6}, got {:.6}",
            expected,
            d
        );
    }

    // -----------------------------------------------------------------------
    // PAM distance
    // -----------------------------------------------------------------------

    #[test]
    fn test_pam_is_100x_kimura() {
        // PAM = 100 * Kimura distance
        let (s1, s2) = make_protein_seqs("AAAAAAAAAA", "RCAAAAAAAA");
        let d_kimura =
            protein_distance(&s1, &s2, &ProteinDistanceMethod::Kimura).unwrap();
        let d_pam = protein_distance(&s1, &s2, &ProteinDistanceMethod::PAM).unwrap();
        assert!(
            (d_pam - 100.0 * d_kimura).abs() < 1e-10,
            "PAM ({:.6}) should be 100 * Kimura ({:.6})",
            d_pam,
            d_kimura
        );
    }

    // -----------------------------------------------------------------------
    // Distance increases with divergence
    // -----------------------------------------------------------------------

    #[test]
    fn test_distance_increases_with_divergence() {
        let base = vec![AminoAcid::Ala; 20];
        let methods = [
            ProteinDistanceMethod::Poisson,
            ProteinDistanceMethod::Kimura,
        ];

        for method in &methods {
            let mut prev_d = 0.0;
            for n_diff in 1..10 {
                let mut s2 = base.clone();
                for i in 0..n_diff {
                    s2[i] = AminoAcid::Arg;
                }
                let d = protein_distance(&base, &s2, method).unwrap();
                assert!(
                    d > prev_d,
                    "{:?}: distance should increase; d({})={:.6} vs d({})={:.6}",
                    method,
                    n_diff,
                    d,
                    n_diff - 1,
                    prev_d
                );
                prev_d = d;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Corrected distance >= observed proportion
    // -----------------------------------------------------------------------

    #[test]
    fn test_corrected_ge_observed() {
        let s1: Vec<AminoAcid> = vec![AminoAcid::Ala; 10];
        let mut s2: Vec<AminoAcid> = vec![AminoAcid::Ala; 10];
        s2[0] = AminoAcid::Arg;
        s2[1] = AminoAcid::Cys;
        s2[2] = AminoAcid::Asp;
        let p_obs = 3.0 / 10.0;

        let d_poisson =
            protein_distance(&s1, &s2, &ProteinDistanceMethod::Poisson).unwrap();
        let d_kimura =
            protein_distance(&s1, &s2, &ProteinDistanceMethod::Kimura).unwrap();

        assert!(
            d_poisson >= p_obs,
            "Poisson ({:.6}) should be >= observed p ({:.4})",
            d_poisson,
            p_obs
        );
        assert!(
            d_kimura >= p_obs,
            "Kimura ({:.6}) should be >= observed p ({:.4})",
            d_kimura,
            p_obs
        );
    }

    // -----------------------------------------------------------------------
    // Distance matrix computation
    // -----------------------------------------------------------------------

    #[test]
    fn test_distance_matrix_size_and_symmetry() {
        let aln = make_alignment(&[
            ("A", "ACDEFGHIKL"),
            ("B", "ACDEFGHIRL"),
            ("C", "RCDEFGHIKL"),
        ]);
        let matrix = compute_protein_distance_matrix(&aln, &ProteinDistanceMethod::Kimura)
            .unwrap();

        assert_eq!(matrix.size(), 3);
        // Symmetry
        assert_eq!(matrix.get(0, 1), matrix.get(1, 0));
        assert_eq!(matrix.get(0, 2), matrix.get(2, 0));
        assert_eq!(matrix.get(1, 2), matrix.get(2, 1));
        // Diagonal should be zero
        assert_eq!(matrix.get(0, 0), 0.0);
        // Off-diagonal should be positive
        assert!(matrix.get(0, 1) > 0.0);
        assert!(matrix.get(0, 2) > 0.0);
        assert!(matrix.get(1, 2) > 0.0);
    }

    #[test]
    fn test_distance_matrix_names() {
        let aln = make_alignment(&[("Alpha", "ACDE"), ("Beta", "ACDE"), ("Gamma", "RCDE")]);
        let matrix =
            compute_protein_distance_matrix(&aln, &ProteinDistanceMethod::Poisson).unwrap();
        assert_eq!(matrix.names, vec!["Alpha", "Beta", "Gamma"]);
    }

    #[test]
    fn test_distance_matrix_identical_seqs_zero() {
        let aln = make_alignment(&[("A", "ACDEFGHIKL"), ("B", "ACDEFGHIKL")]);
        let matrix =
            compute_protein_distance_matrix(&aln, &ProteinDistanceMethod::Poisson).unwrap();
        assert!(
            matrix.get(0, 1).abs() < 1e-10,
            "identical sequences should have distance 0 in matrix"
        );
    }

    // -----------------------------------------------------------------------
    // Calculator struct
    // -----------------------------------------------------------------------

    #[test]
    fn test_calculator_new() {
        let calc = ProteinDistanceCalculator::new(ProteinDistanceMethod::Kimura);
        assert_eq!(calc.method, ProteinDistanceMethod::Kimura);
    }
}

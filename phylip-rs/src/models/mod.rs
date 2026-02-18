//! Nucleotide substitution models and distance measures for phylogenetics.
//!
//! This module provides implementations of classical nucleotide substitution
//! models used in molecular phylogenetics, plus gene frequency distance
//! measures for population genetics:
//!
//! **Nucleotide substitution models** (implement [`DistanceModel`]):
//!
//! - [`jc69::JC69`] — Jukes-Cantor (1969): equal rates, equal frequencies
//! - [`k2p::K2P`] — Kimura 2-parameter (1980): transitions vs transversions
//! - [`f81::F81`] — Felsenstein (1981): unequal frequencies, equal rates
//! - [`f84::F84`] — Felsenstein (1984): unequal frequencies + ts/tv distinction
//! - [`logdet::LogDet`] — LogDet/Paralinear (1994): compositionally robust
//!
//! **Gene frequency distances** (see [`gene_freq`]):
//!
//! - Nei's genetic distance (Nei 1972)
//! - Cavalli-Sforza chord distance (Cavalli-Sforza & Edwards 1967)
//! - Reynolds coancestry distance (Reynolds et al. 1983)
//!
//! All nucleotide substitution models implement the [`DistanceModel`] trait,
//! which provides a uniform interface for computing pairwise evolutionary
//! distances between DNA sequences.
//!
//! # Computing a distance matrix
//!
//! Use [`compute_distance_matrix`] to compute all pairwise distances from an
//! alignment:
//!
//! ```
//! use phylip_rs::tree::{Alignment, Sequence, Base};
//! use phylip_rs::models::{compute_distance_matrix, DistanceModel};
//! use phylip_rs::models::jc69::JC69;
//!
//! let seq1 = Sequence::new("A", vec![Base::A, Base::C, Base::G, Base::T]);
//! let seq2 = Sequence::new("B", vec![Base::A, Base::C, Base::G, Base::A]);
//! let alignment = Alignment::new(vec![seq1, seq2]).unwrap();
//! let model = JC69::new();
//! let matrix = compute_distance_matrix(&alignment, &model).unwrap();
//! ```

pub mod f81;
pub mod f84;
pub mod gene_freq;
pub mod jc69;
pub mod k2p;
pub mod logdet;
pub mod protein;
pub mod protein_distances;
pub mod restriction;

use crate::tree::{Alignment, Base, DistanceMatrix};
use std::fmt;

/// Errors that can occur during evolutionary distance computation.
#[derive(Debug, Clone)]
pub enum DistanceError {
    /// The two sequences share no sites where both have definite (non-gap,
    /// non-ambiguous) bases, so no distance can be estimated.
    NoOverlap {
        /// Index of the first sequence.
        seq1: usize,
        /// Index of the second sequence.
        seq2: usize,
    },
    /// The computed distance is infinite or undefined, typically because
    /// the sequences are too divergent for the model's correction formula.
    /// For example, JC69 requires the fraction of identical sites > 0.25.
    InfiniteDistance {
        /// Index of the first sequence.
        seq1: usize,
        /// Index of the second sequence.
        seq2: usize,
        /// A description of why the distance is infinite.
        reason: String,
    },
    /// The iterative distance estimation procedure did not converge within
    /// the maximum number of iterations.
    ConvergenceFailure {
        /// Index of the first sequence.
        seq1: usize,
        /// Index of the second sequence.
        seq2: usize,
    },
    /// The sequences have different lengths and cannot be compared.
    UnequalLengths,
}

impl fmt::Display for DistanceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DistanceError::NoOverlap { seq1, seq2 } => {
                write!(
                    f,
                    "sequences {} and {} share no sites with definite bases",
                    seq1, seq2
                )
            }
            DistanceError::InfiniteDistance {
                seq1,
                seq2,
                reason,
            } => {
                write!(
                    f,
                    "infinite distance between sequences {} and {}: {}",
                    seq1, seq2, reason
                )
            }
            DistanceError::ConvergenceFailure { seq1, seq2 } => {
                write!(
                    f,
                    "distance estimation did not converge for sequences {} and {}",
                    seq1, seq2
                )
            }
            DistanceError::UnequalLengths => {
                write!(f, "sequences have different lengths")
            }
        }
    }
}

impl std::error::Error for DistanceError {}

/// Trait for nucleotide substitution models that compute pairwise
/// evolutionary distances between DNA sequences.
///
/// Implementors must provide a method to compute the distance between
/// two aligned sequences (represented as slices of [`Base`]) and a
/// human-readable model name.
pub trait DistanceModel {
    /// Compute the evolutionary distance between two aligned DNA sequences.
    ///
    /// Sites where either sequence has a gap or ambiguous base are skipped.
    /// Returns an error if the sequences have no overlapping definite sites,
    /// if the distance is infinite (sequences too divergent), or if the
    /// sequences have different lengths.
    fn distance(&self, seq1: &[Base], seq2: &[Base]) -> Result<f64, DistanceError>;

    /// Returns the name of this substitution model (e.g., "JC69", "K2P").
    fn name(&self) -> &str;
}

/// Count the number of sites where both sequences have definite bases
/// and determine how many are identical vs different.
///
/// Returns `(identical, total)` where `total` is the number of comparable
/// sites (both definite) and `identical` is the number of those that match.
pub(crate) fn count_differences(seq1: &[Base], seq2: &[Base]) -> (usize, usize) {
    let mut identical = 0usize;
    let mut total = 0usize;
    for (&b1, &b2) in seq1.iter().zip(seq2.iter()) {
        if b1.is_definite() && b2.is_definite() {
            total += 1;
            if b1 == b2 {
                identical += 1;
            }
        }
    }
    (identical, total)
}

/// Classify differences between two sequences into transitions and transversions.
///
/// A transition is a purine-to-purine (A<->G) or pyrimidine-to-pyrimidine (C<->T)
/// substitution. A transversion is any other substitution (purine <-> pyrimidine).
///
/// Returns `(transitions, transversions, total)` where `total` is the number
/// of comparable sites with both bases definite.
pub(crate) fn count_transitions_transversions(
    seq1: &[Base],
    seq2: &[Base],
) -> (usize, usize, usize) {
    let mut transitions = 0usize;
    let mut transversions = 0usize;
    let mut total = 0usize;
    for (&b1, &b2) in seq1.iter().zip(seq2.iter()) {
        if b1.is_definite() && b2.is_definite() {
            total += 1;
            if b1 != b2 {
                if is_transition(b1, b2) {
                    transitions += 1;
                } else {
                    transversions += 1;
                }
            }
        }
    }
    (transitions, transversions, total)
}

/// Returns true if the substitution between two definite bases is a transition
/// (purine <-> purine or pyrimidine <-> pyrimidine).
fn is_transition(b1: Base, b2: Base) -> bool {
    matches!(
        (b1, b2),
        (Base::A, Base::G)
            | (Base::G, Base::A)
            | (Base::C, Base::T)
            | (Base::T, Base::C)
    )
}

/// Count bases at comparable (both-definite) sites for estimating base frequencies
/// from a pairwise comparison.
///
/// Returns counts `[A, C, G, T]` summed over both sequences at sites where
/// both have definite bases.
pub(crate) fn count_bases_pairwise(seq1: &[Base], seq2: &[Base]) -> [usize; 4] {
    let mut counts = [0usize; 4];
    for (&b1, &b2) in seq1.iter().zip(seq2.iter()) {
        if b1.is_definite() && b2.is_definite() {
            counts[b1.index().unwrap()] += 1;
            counts[b2.index().unwrap()] += 1;
        }
    }
    counts
}

/// Compute a full pairwise distance matrix from an alignment using the given
/// substitution model.
///
/// This function computes the distance between every pair of sequences in the
/// alignment and stores the results in a symmetric [`DistanceMatrix`].
///
/// # Errors
///
/// Returns a [`DistanceError`] if any pairwise distance computation fails.
/// The error will contain the indices of the problematic pair.
pub fn compute_distance_matrix(
    alignment: &Alignment,
    model: &dyn DistanceModel,
) -> Result<DistanceMatrix, DistanceError> {
    let n = alignment.ntaxa();
    let names: Vec<String> = (0..n)
        .map(|i| alignment.taxon_name(i).to_string())
        .collect();
    let mut matrix = DistanceMatrix::new(names);

    for i in 0..n {
        for j in (i + 1)..n {
            let d = model
                .distance(&alignment.sequences[i].bases, &alignment.sequences[j].bases)
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
                    DistanceError::ConvergenceFailure { .. } => {
                        DistanceError::ConvergenceFailure {
                            seq1: i,
                            seq2: j,
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
    use crate::tree::{Alignment, Base, Sequence};

    #[test]
    fn test_count_differences_identical() {
        let seq = vec![Base::A, Base::C, Base::G, Base::T];
        let (identical, total) = count_differences(&seq, &seq);
        assert_eq!(identical, 4);
        assert_eq!(total, 4);
    }

    #[test]
    fn test_count_differences_all_different() {
        let seq1 = vec![Base::A, Base::C, Base::G, Base::T];
        let seq2 = vec![Base::T, Base::G, Base::C, Base::A];
        let (identical, total) = count_differences(&seq1, &seq2);
        assert_eq!(identical, 0);
        assert_eq!(total, 4);
    }

    #[test]
    fn test_count_differences_skips_gaps() {
        let seq1 = vec![Base::A, Base::Gap, Base::G, Base::T];
        let seq2 = vec![Base::A, Base::C, Base::G, Base::A];
        let (identical, total) = count_differences(&seq1, &seq2);
        assert_eq!(identical, 2); // A=A, G=G
        assert_eq!(total, 3); // gap site skipped
    }

    #[test]
    fn test_count_differences_skips_ambiguous() {
        let seq1 = vec![Base::A, Base::N, Base::G, Base::T];
        let seq2 = vec![Base::A, Base::C, Base::R, Base::T];
        let (identical, total) = count_differences(&seq1, &seq2);
        assert_eq!(identical, 2); // A=A, T=T
        assert_eq!(total, 2); // N and R sites skipped
    }

    #[test]
    fn test_count_transitions_transversions() {
        // A->G is transition, A->C is transversion, G=G is identical, T->C is transition
        let seq1 = vec![Base::A, Base::A, Base::G, Base::T];
        let seq2 = vec![Base::G, Base::C, Base::G, Base::C];
        let (ts, tv, total) = count_transitions_transversions(&seq1, &seq2);
        assert_eq!(ts, 2); // A->G, T->C
        assert_eq!(tv, 1); // A->C
        assert_eq!(total, 4);
    }

    #[test]
    fn test_compute_distance_matrix() {
        let seq1 = Sequence::new("A", vec![Base::A, Base::C, Base::G, Base::T]);
        let seq2 = Sequence::new("B", vec![Base::A, Base::C, Base::G, Base::A]);
        let seq3 = Sequence::new("C", vec![Base::A, Base::T, Base::G, Base::T]);
        let alignment = Alignment::new(vec![seq1, seq2, seq3]).unwrap();
        let model = jc69::JC69::new();
        let matrix = compute_distance_matrix(&alignment, &model).unwrap();
        assert_eq!(matrix.size(), 3);
        assert!(matrix.get(0, 1) > 0.0);
        assert!(matrix.get(0, 2) > 0.0);
        assert!(matrix.get(1, 2) > 0.0);
        // Symmetry
        assert_eq!(matrix.get(0, 1), matrix.get(1, 0));
        // d(A,B) should be less than d(B,C) since B differs from C at 2 sites
        // and B differs from A at 1 site
        assert!(matrix.get(0, 1) < matrix.get(1, 2));
    }
}

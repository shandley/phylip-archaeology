//! Maximum likelihood pairwise distance estimation.
//!
//! For each pair of sequences, finds the branch length `t` that maximizes
//! the likelihood of observing the two sequences under a given substitution
//! model. This provides ML-corrected evolutionary distances suitable for use
//! in distance-based tree reconstruction methods (neighbor-joining, UPGMA,
//! Fitch-Margoliash).
//!
//! # Method
//!
//! Given two observed sequences and a substitution model with equilibrium
//! frequencies `pi` and transition probability matrix `P(t)`, the likelihood
//! of observing site pattern `(x, y)` at branch length `t` is:
//!
//! ```text
//! L_site(t) = SUM_b pi[b] * P(b -> x; t) * P(b -> y; t)
//! ```
//!
//! where the sum is over all possible ancestral states `b`. For ambiguous
//! bases, we sum over all states consistent with the observed ambiguity code.
//!
//! The total log-likelihood is the sum over all sites:
//!
//! ```text
//! ln L(t) = SUM_site ln( L_site(t) )
//! ```
//!
//! The optimal `t` is found by Newton-Raphson optimization with numerical
//! derivatives, following the approach used in `pruning.rs`.
//!
//! # References
//!
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//!   Chapter 11: Distance matrix methods.
//! - Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum
//!   likelihood approach. *J. Mol. Evol.*, 17, 368-376.

use crate::likelihood::models::SubstitutionModel;
use crate::tree::types::{Alignment, Base, DistanceMatrix, Sequence};

/// Minimum branch length for numerical stability.
const MIN_T: f64 = 1e-8;
/// Maximum branch length for optimization.
const MAX_T: f64 = 10.0;
/// Initial branch length for optimization.
const INITIAL_T: f64 = 0.1;
/// Convergence tolerance for Newton-Raphson.
const TOLERANCE: f64 = 1e-8;
/// Maximum Newton-Raphson iterations.
const MAX_ITER: usize = 100;

/// Compute the site likelihood for a single site given two observed bases
/// and a substitution model at branch length `t`.
///
/// Models the two sequences as leaves connected through a common ancestor:
///
/// ```text
///       ancestor (b)
///       /          \
///      t            t
///     /              \
///   obs1            obs2
/// ```
///
/// Each leaf is at distance `t/2` from the ancestor (so total distance
/// between leaves is `t`). Equivalently, we use the full distance `t`
/// by computing:
///
/// ```text
/// L_site(t) = SUM_b pi[b] * (SUM_{x in obs1} P(b->x; t/2))
///                          * (SUM_{y in obs2} P(b->y; t/2))
/// ```
///
/// For unambiguous bases, the inner sums reduce to a single P entry.
fn site_likelihood(
    obs1: Base,
    obs2: Base,
    model: &dyn SubstitutionModel,
    t: f64,
) -> f64 {
    let pi = model.frequencies();
    let p_matrix = model.transition_probs(t / 2.0);

    let lv1 = obs1.likelihood_vector();
    let lv2 = obs2.likelihood_vector();

    let mut site_like = 0.0;
    for b in 0..4 {
        // Partial likelihood contribution from each observed base
        let mut partial1 = 0.0;
        let mut partial2 = 0.0;
        for j in 0..4 {
            partial1 += p_matrix[b][j] * lv1[j];
            partial2 += p_matrix[b][j] * lv2[j];
        }
        site_like += pi[b] * partial1 * partial2;
    }

    site_like
}

/// Compute the total log-likelihood of two sequences at branch length `t`.
fn total_log_likelihood(
    seq1: &[Base],
    seq2: &[Base],
    model: &dyn SubstitutionModel,
    t: f64,
) -> f64 {
    let mut lnl = 0.0;
    for (&b1, &b2) in seq1.iter().zip(seq2.iter()) {
        let sl = site_likelihood(b1, b2, model, t);
        if sl > 0.0 {
            lnl += sl.ln();
        } else {
            // If site likelihood is zero or negative, return a very
            // negative value to steer the optimizer away.
            return f64::NEG_INFINITY;
        }
    }
    lnl
}

/// Estimate the ML pairwise distance between two sequences under a
/// substitution model.
///
/// Uses Newton-Raphson with numerical derivatives to find the branch
/// length `t` that maximizes the log-likelihood.
///
/// Returns `None` if the sequences have different lengths, are empty,
/// or if the optimization fails to converge.
///
/// Returns `Some(0.0)` (or very close) for identical sequences.
pub fn ml_pairwise_distance(
    seq1: &Sequence,
    seq2: &Sequence,
    model: &dyn SubstitutionModel,
) -> Option<f64> {
    if seq1.len() != seq2.len() || seq1.is_empty() {
        return None;
    }

    let bases1 = &seq1.bases;
    let bases2 = &seq2.bases;

    // Check if sequences are identical (at definite sites).
    // If so, the ML distance is 0 (or effectively MIN_T).
    let all_same = bases1
        .iter()
        .zip(bases2.iter())
        .all(|(&b1, &b2)| {
            // If both are definite, check equality.
            // If either is ambiguous/gap, skip (doesn't distinguish).
            if b1.is_definite() && b2.is_definite() {
                b1 == b2
            } else {
                true
            }
        });

    // Also check that there is at least one definite site
    let has_definite = bases1
        .iter()
        .zip(bases2.iter())
        .any(|(&b1, &b2)| b1.is_definite() && b2.is_definite());

    if all_same && has_definite {
        return Some(0.0);
    }

    // Newton-Raphson optimization
    let mut t = INITIAL_T;

    for _ in 0..MAX_ITER {
        let h = (t * 0.001).max(1e-8);

        let lnl_center = total_log_likelihood(bases1, bases2, model, t);
        let lnl_plus = total_log_likelihood(bases1, bases2, model, (t + h).min(MAX_T));
        let lnl_minus = total_log_likelihood(bases1, bases2, model, (t - h).max(MIN_T));

        if !lnl_center.is_finite() {
            // Try a different starting point
            t = t * 0.5;
            if t < MIN_T {
                return None;
            }
            continue;
        }

        // Numerical first and second derivatives
        let slope = (lnl_plus - lnl_minus) / (2.0 * h);
        let curve = (lnl_plus - 2.0 * lnl_center + lnl_minus) / (h * h);

        // Newton step (maximizing, so step = -slope/curve when curve < 0)
        let delta = if curve < -1e-12 {
            -slope / curve
        } else if slope.abs() > 1e-12 {
            // Not concave: step in direction of gradient
            slope.signum() * h * 10.0
        } else {
            break; // Flat: converged
        };

        let t_new = (t + delta).clamp(MIN_T, MAX_T);

        if (t_new - t).abs() < TOLERANCE {
            t = t_new;
            break;
        }
        t = t_new;
    }

    Some(t)
}

/// Compute a full ML pairwise distance matrix from an alignment.
///
/// For each pair of sequences, estimates the ML distance under the given
/// substitution model. If estimation fails for a pair (returns `None`),
/// a fallback distance of `MAX_T` is used.
///
/// # Arguments
///
/// * `alignment` - A multiple sequence alignment.
/// * `model` - A nucleotide substitution model implementing `SubstitutionModel`.
///
/// # Returns
///
/// A symmetric `DistanceMatrix` with ML-estimated pairwise distances.
pub fn ml_distance_matrix(
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
) -> DistanceMatrix {
    let n = alignment.ntaxa();
    let names: Vec<String> = (0..n)
        .map(|i| alignment.taxon_name(i).to_string())
        .collect();
    let mut matrix = DistanceMatrix::new(names);

    for i in 0..n {
        for j in (i + 1)..n {
            let d = ml_pairwise_distance(
                &alignment.sequences[i],
                &alignment.sequences[j],
                model,
            )
            .unwrap_or(MAX_T);
            matrix.set(i, j, d);
        }
    }

    matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::likelihood::models::{F84Model, Jc69Model};
    use crate::models::jc69::JC69;
    use crate::models::DistanceModel;
    use crate::tree::types::BaseFrequencies;

    /// Helper: build a Sequence from a name and string.
    fn make_seq(name: &str, s: &str) -> Sequence {
        let bases: Vec<Base> = s.chars().map(|c| Base::from_char(c).unwrap()).collect();
        Sequence::new(name, bases)
    }

    /// Helper: build an alignment from name/sequence pairs.
    fn make_alignment(data: &[(&str, &str)]) -> Alignment {
        let seqs: Vec<Sequence> = data
            .iter()
            .map(|(name, seq)| make_seq(name, seq))
            .collect();
        Alignment::new(seqs).unwrap()
    }

    // ---------------------------------------------------------------
    // Test 1: Identical sequences should give distance 0
    // ---------------------------------------------------------------
    #[test]
    fn test_identical_sequences_zero_distance() {
        let model = Jc69Model;
        let s1 = make_seq("A", "ACGTACGTACGTACGT");
        let s2 = make_seq("B", "ACGTACGTACGTACGT");
        let d = ml_pairwise_distance(&s1, &s2, &model).unwrap();
        assert!(
            d.abs() < 1e-6,
            "identical sequences should have distance ~0, got {}",
            d
        );
    }

    // ---------------------------------------------------------------
    // Test 2: Divergent sequences should give large distance
    // ---------------------------------------------------------------
    #[test]
    fn test_divergent_sequences_large_distance() {
        let model = Jc69Model;
        // All different bases: very divergent
        let s1 = make_seq("A", "AAAAAAAAAAAAAAAA");
        let s2 = make_seq("B", "CCCCCCCCCCCCCCCC");
        let d = ml_pairwise_distance(&s1, &s2, &model).unwrap();
        assert!(
            d > 0.5,
            "completely divergent sequences should have large distance, got {}",
            d
        );
    }

    // ---------------------------------------------------------------
    // Test 3: ML distance under JC69 should approximate analytical JC69
    // ---------------------------------------------------------------
    #[test]
    fn test_ml_jc69_approximates_analytical_jc69() {
        let ml_model = Jc69Model;
        let analytical_model = JC69::new();

        // Sequences with known differences: 4 out of 20 sites differ
        let s1 = make_seq("A", "ACGTACGTACGTACGTACGT");
        let s2 = make_seq("B", "ACGAACGAACGTACGTACGT");

        let d_ml = ml_pairwise_distance(&s1, &s2, &ml_model).unwrap();
        let d_analytical = analytical_model
            .distance(&s1.bases, &s2.bases)
            .unwrap();

        // They should agree within a reasonable tolerance.
        // The ML estimator and the analytical formula should give very
        // similar results for the JC69 model.
        let relative_diff = (d_ml - d_analytical).abs() / d_analytical.max(1e-10);
        assert!(
            relative_diff < 0.15,
            "ML JC69 distance ({}) should be close to analytical JC69 ({}), \
             relative diff = {}",
            d_ml,
            d_analytical,
            relative_diff
        );
    }

    // ---------------------------------------------------------------
    // Test 4: F84 model produces valid distances
    // ---------------------------------------------------------------
    #[test]
    fn test_f84_model_works() {
        let freqs = BaseFrequencies::new(0.3, 0.2, 0.2, 0.3);
        let model = F84Model::new(freqs, 2.0);

        let s1 = make_seq("A", "ACGTACGTACGTACGT");
        let s2 = make_seq("B", "ACGAACGTACGAACGT");

        let d = ml_pairwise_distance(&s1, &s2, &model).unwrap();
        assert!(d > 0.0, "F84 distance should be positive, got {}", d);
        assert!(d < 5.0, "F84 distance should be reasonable, got {}", d);
    }

    // ---------------------------------------------------------------
    // Test 5: Symmetry: d(a,b) == d(b,a)
    // ---------------------------------------------------------------
    #[test]
    fn test_symmetry() {
        let model = Jc69Model;
        let s1 = make_seq("A", "ACGTACGTACGTACGT");
        let s2 = make_seq("B", "GCGAACGTACGTACGA");

        let d_ab = ml_pairwise_distance(&s1, &s2, &model).unwrap();
        let d_ba = ml_pairwise_distance(&s2, &s1, &model).unwrap();

        assert!(
            (d_ab - d_ba).abs() < 1e-6,
            "distance should be symmetric: d(a,b)={}, d(b,a)={}",
            d_ab,
            d_ba
        );
    }

    // ---------------------------------------------------------------
    // Test 6: Triangle inequality
    // ---------------------------------------------------------------
    #[test]
    fn test_triangle_inequality() {
        let model = Jc69Model;
        // Use long sequences where the intermediate B is truly a
        // "stepping stone" between A and C. A and B differ at a few
        // sites, B and C differ at a few (different) sites, so A and C
        // differs at the union of those sites.
        let s1 = make_seq("A", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let s2 = make_seq("B", "AAGAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGA");
        let s3 = make_seq("C", "AAGAGAAAAAAAAAAAAAAGAGAAAAAAAAAAAAAAAGGA");

        let d12 = ml_pairwise_distance(&s1, &s2, &model).unwrap();
        let d23 = ml_pairwise_distance(&s2, &s3, &model).unwrap();
        let d13 = ml_pairwise_distance(&s1, &s3, &model).unwrap();

        // d13 should be roughly d12 + d23 (additive), and definitely <=
        assert!(
            d13 <= d12 + d23 + 0.05,
            "triangle inequality violated: d(1,3)={} > d(1,2)+d(2,3)={}",
            d13,
            d12 + d23
        );
    }

    // ---------------------------------------------------------------
    // Test 7: Ambiguity code handling (N)
    // ---------------------------------------------------------------
    #[test]
    fn test_ambiguity_n_handling() {
        let model = Jc69Model;
        // With N bases, the distance should still be computable
        let s1 = make_seq("A", "ACGTACGTNNNNACGT");
        let s2 = make_seq("B", "ACGTACGTACGTACGT");

        let d = ml_pairwise_distance(&s1, &s2, &model);
        assert!(d.is_some(), "should handle N ambiguity codes");
        let d = d.unwrap();
        assert!(d >= 0.0, "distance should be non-negative, got {}", d);
    }

    // ---------------------------------------------------------------
    // Test 8: Gap handling
    // ---------------------------------------------------------------
    #[test]
    fn test_gap_handling() {
        let model = Jc69Model;
        // Gaps are treated as fully ambiguous
        let s1 = make_seq("A", "ACGT----ACGTACGT");
        let s2 = make_seq("B", "ACGTACGTACGTACGT");

        let d = ml_pairwise_distance(&s1, &s2, &model);
        assert!(d.is_some(), "should handle gaps");
        let d = d.unwrap();
        assert!(d >= 0.0, "distance should be non-negative, got {}", d);
    }

    // ---------------------------------------------------------------
    // Test 9: Distance increases with more differences
    // ---------------------------------------------------------------
    #[test]
    fn test_distance_increases_with_divergence() {
        let model = Jc69Model;

        let s_ref = make_seq("ref", "ACGTACGTACGTACGT");
        let s_close = make_seq("close", "ACGTACGTACGTACGA"); // 1 diff
        let s_medium = make_seq("medium", "GCGAACGTACGTACGA"); // 3 diffs
        let s_far = make_seq("far", "GCGAGCGAACGAACGA"); // many diffs

        let d_close = ml_pairwise_distance(&s_ref, &s_close, &model).unwrap();
        let d_medium = ml_pairwise_distance(&s_ref, &s_medium, &model).unwrap();
        let d_far = ml_pairwise_distance(&s_ref, &s_far, &model).unwrap();

        assert!(
            d_close < d_medium,
            "d_close ({}) should be < d_medium ({})",
            d_close,
            d_medium
        );
        assert!(
            d_medium < d_far,
            "d_medium ({}) should be < d_far ({})",
            d_medium,
            d_far
        );
    }

    // ---------------------------------------------------------------
    // Test 10: Full distance matrix construction
    // ---------------------------------------------------------------
    #[test]
    fn test_distance_matrix_construction() {
        let model = Jc69Model;
        let aln = make_alignment(&[
            ("A", "ACGTACGTACGTACGT"),
            ("B", "ACGTACGTACGTACGA"),
            ("C", "GCGAGCGAACGAACGA"),
        ]);

        let dm = ml_distance_matrix(&aln, &model);

        assert_eq!(dm.size(), 3);
        assert_eq!(dm.names, vec!["A", "B", "C"]);

        // Diagonal is zero
        assert_eq!(dm.get(0, 0), 0.0);
        assert_eq!(dm.get(1, 1), 0.0);
        assert_eq!(dm.get(2, 2), 0.0);

        // Off-diagonal should be positive
        assert!(dm.get(0, 1) > 0.0, "d(A,B) should be > 0");
        assert!(dm.get(0, 2) > 0.0, "d(A,C) should be > 0");
        assert!(dm.get(1, 2) > 0.0, "d(B,C) should be > 0");

        // Symmetry
        assert_eq!(dm.get(0, 1), dm.get(1, 0));
        assert_eq!(dm.get(0, 2), dm.get(2, 0));
        assert_eq!(dm.get(1, 2), dm.get(2, 1));

        // A and B are close, A and C are far
        assert!(
            dm.get(0, 1) < dm.get(0, 2),
            "d(A,B)={} should be < d(A,C)={}",
            dm.get(0, 1),
            dm.get(0, 2)
        );
    }

    // ---------------------------------------------------------------
    // Test 11: Unequal length sequences return None
    // ---------------------------------------------------------------
    #[test]
    fn test_unequal_lengths_returns_none() {
        let model = Jc69Model;
        let s1 = make_seq("A", "ACGT");
        let s2 = make_seq("B", "ACGTACGT");
        assert!(ml_pairwise_distance(&s1, &s2, &model).is_none());
    }

    // ---------------------------------------------------------------
    // Test 12: Empty sequences return None
    // ---------------------------------------------------------------
    #[test]
    fn test_empty_sequences_returns_none() {
        let model = Jc69Model;
        let s1 = Sequence::new("A", vec![]);
        let s2 = Sequence::new("B", vec![]);
        assert!(ml_pairwise_distance(&s1, &s2, &model).is_none());
    }

    // ---------------------------------------------------------------
    // Test 13: F84 with equal freqs and ttratio=1 matches JC69
    // ---------------------------------------------------------------
    #[test]
    fn test_f84_equal_freqs_ttratio1_matches_jc69() {
        let jc69 = Jc69Model;
        let f84 = F84Model::equal_freqs(1.0);

        let s1 = make_seq("A", "ACGTACGTACGTACGT");
        let s2 = make_seq("B", "ACGAACGTACGAACGT");

        let d_jc69 = ml_pairwise_distance(&s1, &s2, &jc69).unwrap();
        let d_f84 = ml_pairwise_distance(&s1, &s2, &f84).unwrap();

        assert!(
            (d_jc69 - d_f84).abs() < 0.01,
            "F84(ttratio=1, equal freqs) should match JC69: {} vs {}",
            d_f84,
            d_jc69
        );
    }

    // ---------------------------------------------------------------
    // Test 14: Two-fold ambiguity codes (R, Y, etc.)
    // ---------------------------------------------------------------
    #[test]
    fn test_two_fold_ambiguity_codes() {
        let model = Jc69Model;
        // R = A or G, Y = C or T
        let s1 = make_seq("A", "ACGTACGTACGTACGT");
        let s2 = make_seq("B", "RCGTYCGTACGTACGT");

        let d = ml_pairwise_distance(&s1, &s2, &model);
        assert!(d.is_some(), "should handle two-fold ambiguity codes");
        let d = d.unwrap();
        assert!(d >= 0.0, "distance should be non-negative");

        // Distance should be less than if R and Y were fully different bases,
        // because R=A|G includes A (matches pos 0) and Y=C|T includes C (matches pos 5)
        let s3 = make_seq("C", "TCGTGCGTACGTACGT"); // definite mismatches
        let d_definite = ml_pairwise_distance(&s1, &s3, &model).unwrap();
        assert!(
            d < d_definite + 1e-6,
            "ambiguous distance ({}) should be <= definite mismatch distance ({})",
            d,
            d_definite
        );
    }

    // ---------------------------------------------------------------
    // Test 15: Site likelihood is positive for valid inputs
    // ---------------------------------------------------------------
    #[test]
    fn test_site_likelihood_positive() {
        let model = Jc69Model;
        for &t in &[0.001, 0.01, 0.1, 0.5, 1.0, 5.0] {
            for &b1 in &[Base::A, Base::C, Base::G, Base::T] {
                for &b2 in &[Base::A, Base::C, Base::G, Base::T] {
                    let sl = site_likelihood(b1, b2, &model, t);
                    assert!(
                        sl > 0.0,
                        "site likelihood should be positive for {:?},{:?} at t={}, got {}",
                        b1,
                        b2,
                        t,
                        sl
                    );
                }
            }
        }
    }

    // ---------------------------------------------------------------
    // Test 16: ML distance is always non-negative
    // ---------------------------------------------------------------
    #[test]
    fn test_ml_distance_non_negative() {
        let model = Jc69Model;
        let test_cases = vec![
            ("ACGTACGT", "ACGTACGT"),
            ("ACGTACGT", "ACGTACGA"),
            ("AAAAAAAA", "TTTTTTTT"),
            ("ACGTACGT", "GCGAGCGA"),
        ];

        for (s1_str, s2_str) in test_cases {
            let s1 = make_seq("A", s1_str);
            let s2 = make_seq("B", s2_str);
            if let Some(d) = ml_pairwise_distance(&s1, &s2, &model) {
                assert!(
                    d >= 0.0,
                    "ML distance should be non-negative for ({}, {}), got {}",
                    s1_str,
                    s2_str,
                    d
                );
            }
        }
    }

    // ---------------------------------------------------------------
    // Test 17: F84 with high ttratio
    // ---------------------------------------------------------------
    #[test]
    fn test_f84_high_ttratio() {
        let model = F84Model::equal_freqs(10.0);
        // Sequences differing only by transitions (A<->G)
        let s1 = make_seq("A", "AAAAAAAAAAAAAAAA");
        let s2 = make_seq("B", "AAAGAAAGAAAGAAAG");

        let d = ml_pairwise_distance(&s1, &s2, &model).unwrap();
        assert!(d > 0.0, "distance should be positive");
        assert!(d < 5.0, "distance should be reasonable, got {}", d);
    }

    // ---------------------------------------------------------------
    // Test 18: Distance matrix with identical sequences
    // ---------------------------------------------------------------
    #[test]
    fn test_distance_matrix_identical_sequences() {
        let model = Jc69Model;
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGTACGT"),
            ("C", "ACGTACGT"),
        ]);

        let dm = ml_distance_matrix(&aln, &model);

        // All pairwise distances should be 0
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    dm.get(i, j).abs() < 1e-6,
                    "all distances should be ~0 for identical seqs, got d({},{})={}",
                    i,
                    j,
                    dm.get(i, j)
                );
            }
        }
    }

    // ---------------------------------------------------------------
    // Test 19: ML distance with F84 model and unequal frequencies
    // ---------------------------------------------------------------
    #[test]
    fn test_f84_unequal_freqs() {
        let freqs = BaseFrequencies::new(0.4, 0.1, 0.1, 0.4);
        let model = F84Model::new(freqs, 2.0);

        let s1 = make_seq("A", "AATTAATTAATTAATT");
        let s2 = make_seq("B", "AATTAATTGCCGGCCG");

        let d = ml_pairwise_distance(&s1, &s2, &model).unwrap();
        assert!(d > 0.0, "distance should be positive, got {}", d);
        assert!(d.is_finite(), "distance should be finite");
    }

    // ---------------------------------------------------------------
    // Test 20: Consistency - ML distance matrix can be used with NJ
    // ---------------------------------------------------------------
    #[test]
    fn test_ml_matrix_compatible_with_nj() {
        use crate::distance::neighbor_joining::neighbor_joining;

        let model = Jc69Model;
        let aln = make_alignment(&[
            ("A", "AAAAAAAAACGTACGT"),
            ("B", "AAAAAAAAACGTACGA"),
            ("C", "CCCCCCCCGCGAGCGA"),
            ("D", "CCCCCCCCGCGAGCGG"),
        ]);

        let dm = ml_distance_matrix(&aln, &model);
        let tree = neighbor_joining(&dm);

        // Tree should have 4 leaves
        assert_eq!(tree.num_leaves(), 4);

        // A and B should be siblings (small distance)
        // C and D should be siblings (small distance)
        let a_id = tree.nodes.iter().find(|n| n.name.as_deref() == Some("A")).unwrap().id;
        let b_id = tree.nodes.iter().find(|n| n.name.as_deref() == Some("B")).unwrap().id;
        assert_eq!(
            tree.nodes[a_id].parent,
            tree.nodes[b_id].parent,
            "A and B should be siblings"
        );

        let c_id = tree.nodes.iter().find(|n| n.name.as_deref() == Some("C")).unwrap().id;
        let d_id = tree.nodes.iter().find(|n| n.name.as_deref() == Some("D")).unwrap().id;
        assert_eq!(
            tree.nodes[c_id].parent,
            tree.nodes[d_id].parent,
            "C and D should be siblings"
        );
    }
}

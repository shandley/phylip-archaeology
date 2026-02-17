//! Transition/transversion ratio (ts/tv) estimation.
//!
//! Provides methods for estimating the transition/transversion ratio from
//! DNA sequence data, a key parameter for substitution models like F84 that
//! distinguish between transitions (purine-purine or pyrimidine-pyrimidine
//! substitutions: A<->G, C<->T) and transversions (all other substitutions).
//!
//! Three estimation approaches are provided:
//!
//! 1. **Counting method**: Count observed transitions and transversions from
//!    pairwise sequence comparisons and compute the ratio directly.
//!
//! 2. **Empirical base frequencies**: Count nucleotide frequencies from alignment
//!    data, useful as a companion to ttratio estimation.
//!
//! 3. **Maximum likelihood estimation**: Optimize ttratio by maximizing the
//!    log-likelihood of the data given a tree, using golden section search.
//!
//! # Example
//!
//! ```
//! use phylip_rs::tree::{Alignment, Base, BaseFrequencies, Sequence};
//! use phylip_rs::tree::newick::parse_newick;
//! use phylip_rs::likelihood::ttratio::{
//!     count_ts_tv, estimate_ttratio_counting, estimate_base_frequencies,
//!     estimate_ttratio_ml,
//! };
//!
//! // Create a simple alignment with transitions and transversions
//! let alignment = Alignment::new(vec![
//!     Sequence::new("A", vec![Base::A, Base::C, Base::G, Base::T, Base::A, Base::C]),
//!     Sequence::new("B", vec![Base::G, Base::T, Base::G, Base::T, Base::C, Base::A]),
//!     Sequence::new("C", vec![Base::A, Base::C, Base::A, Base::C, Base::G, Base::T]),
//! ]).unwrap();
//!
//! // Count transitions and transversions between two sequences
//! let (ts, tv) = count_ts_tv(&alignment.sequences[0], &alignment.sequences[1]);
//! assert!(ts > 0 || tv > 0);
//!
//! // Estimate ttratio from all pairwise comparisons
//! let ratio = estimate_ttratio_counting(&alignment);
//! assert!(ratio > 0.0);
//!
//! // Estimate base frequencies
//! let freqs = estimate_base_frequencies(&alignment);
//! assert!((freqs.freq_a + freqs.freq_c + freqs.freq_g + freqs.freq_t - 1.0).abs() < 1e-10);
//!
//! // ML estimation given a tree
//! let tree = parse_newick("((A:0.1,B:0.2):0.05,C:0.15);").unwrap();
//! let ml_ratio = estimate_ttratio_ml(&tree, &alignment, &freqs).unwrap();
//! assert!(ml_ratio > 0.0);
//! ```
//!
//! # References
//!
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//! - Felsenstein, J. (1984). Distance methods for inferring phylogenies:
//!   a justification. *Evolution*, 38(1), 16-24.

use crate::tree::types::{Alignment, Base, BaseFrequencies, Sequence, Tree};

use super::models::F84Model;
use super::pruning::{log_likelihood, LikelihoodError};

/// Default ttratio returned when estimation is not possible
/// (e.g., zero transversions observed).
const DEFAULT_TTRATIO: f64 = 2.0;

/// Minimum ttratio for ML search bracket.
const MIN_TTRATIO: f64 = 0.5;

/// Maximum ttratio for ML search bracket.
const MAX_TTRATIO: f64 = 20.0;

/// Golden ratio constant used in golden section search.
const GOLDEN_RATIO: f64 = 1.618033988749895;

/// Convergence tolerance for golden section search (bracket width).
const TOLERANCE: f64 = 1e-4;

/// Maximum iterations for golden section search.
const MAX_ITERATIONS: usize = 50;

/// Returns true if the substitution from base `a` to base `b` is a transition.
///
/// Transitions are purine-purine (A<->G) or pyrimidine-pyrimidine (C<->T)
/// substitutions. All other substitutions are transversions.
fn is_transition(a: Base, b: Base) -> bool {
    matches!(
        (a, b),
        (Base::A, Base::G)
            | (Base::G, Base::A)
            | (Base::C, Base::T)
            | (Base::T, Base::C)
    )
}

/// Count transitions and transversions between two sequences.
///
/// Only sites where both sequences have definite (unambiguous) bases
/// are considered. Sites with gaps, unknown, or ambiguous bases are skipped.
/// Sites where both sequences have the same base are also skipped (no
/// substitution occurred).
///
/// Returns a tuple `(transitions, transversions)`.
///
/// # Panics
///
/// Does not panic. If sequences have different lengths, only the overlapping
/// portion is compared.
pub fn count_ts_tv(seq1: &Sequence, seq2: &Sequence) -> (usize, usize) {
    let len = seq1.len().min(seq2.len());
    let mut transitions = 0usize;
    let mut transversions = 0usize;

    for i in 0..len {
        let b1 = seq1.bases[i];
        let b2 = seq2.bases[i];

        // Skip non-definite bases (gaps, ambiguities, unknowns)
        if !b1.is_definite() || !b2.is_definite() {
            continue;
        }

        // Skip identical sites
        if b1 == b2 {
            continue;
        }

        if is_transition(b1, b2) {
            transitions += 1;
        } else {
            transversions += 1;
        }
    }

    (transitions, transversions)
}

/// Estimate the transition/transversion ratio from pairwise counts across
/// all sequence pairs in an alignment.
///
/// Sums the transition and transversion counts over all unique pairs of
/// sequences in the alignment. The ratio is computed as:
///
/// ```text
/// ttratio = total_transitions / total_transversions
/// ```
///
/// If no transversions are observed (e.g., all differences are transitions,
/// or sequences are identical), returns the default value of 2.0.
///
/// If no differences are observed at all (all sequences identical at
/// informative sites), returns the default value of 2.0.
pub fn estimate_ttratio_counting(alignment: &Alignment) -> f64 {
    let ntaxa = alignment.ntaxa();
    let mut total_ts = 0usize;
    let mut total_tv = 0usize;

    for i in 0..ntaxa {
        for j in (i + 1)..ntaxa {
            let (ts, tv) = count_ts_tv(&alignment.sequences[i], &alignment.sequences[j]);
            total_ts += ts;
            total_tv += tv;
        }
    }

    if total_tv == 0 {
        // No transversions observed; cannot compute ratio from counts
        return DEFAULT_TTRATIO;
    }

    let ratio = total_ts as f64 / total_tv as f64;

    // Clamp to a reasonable range
    if ratio < MIN_TTRATIO {
        MIN_TTRATIO
    } else if ratio > MAX_TTRATIO {
        MAX_TTRATIO
    } else {
        ratio
    }
}

/// Estimate base frequencies from alignment data.
///
/// Counts the occurrences of each definite base (A, C, G, T) across all
/// sequences and sites, then normalizes to produce frequency estimates.
///
/// Sites with gaps, unknown, or ambiguous bases are excluded from the count.
///
/// If the alignment contains no definite bases, returns equal frequencies
/// (0.25 each).
pub fn estimate_base_frequencies(alignment: &Alignment) -> BaseFrequencies {
    BaseFrequencies::from_alignment(alignment)
}

/// Estimate the transition/transversion ratio by maximum likelihood.
///
/// Given a tree topology (with branch lengths) and alignment, finds the
/// ttratio that maximizes the log-likelihood under the F84 substitution
/// model, using golden section search over the bracket
/// [`MIN_TTRATIO`, `MAX_TTRATIO`] (i.e., [0.5, 20.0]).
///
/// The search terminates when the bracket width is less than [`TOLERANCE`]
/// (1e-4) or after [`MAX_ITERATIONS`] (50) iterations.
///
/// # Errors
///
/// Returns a [`LikelihoodError`] if the likelihood computation fails
/// (e.g., taxa mismatch between tree and alignment).
pub fn estimate_ttratio_ml(
    tree: &Tree,
    alignment: &Alignment,
    freqs: &BaseFrequencies,
) -> Result<f64, LikelihoodError> {
    // Objective function: log-likelihood as a function of ttratio
    let lnl_at = |ratio: f64| -> Result<f64, LikelihoodError> {
        let model = F84Model::new(*freqs, ratio);
        log_likelihood(tree, alignment, &model)
    };

    // Golden section search to maximize log-likelihood.
    //
    // We maintain a bracket [a, b] with two interior probe points x1 < x2:
    //   x1 = a + resphi * (b - a)     (closer to a)
    //   x2 = b - resphi * (b - a)     (closer to b)
    // where resphi = 2 - phi ≈ 0.382.
    //
    // At each step we discard the portion of the bracket on the side
    // of the worse probe point.
    let resphi = 2.0 - GOLDEN_RATIO; // ≈ 0.382

    let mut a = MIN_TTRATIO;
    let mut b = MAX_TTRATIO;

    // x1 is the left interior point, x2 is the right interior point
    let mut x1 = a + resphi * (b - a);
    let mut x2 = b - resphi * (b - a);

    let mut f1 = lnl_at(x1)?;
    let mut f2 = lnl_at(x2)?;

    for _ in 0..MAX_ITERATIONS {
        if (b - a).abs() < TOLERANCE {
            break;
        }

        if f1 > f2 {
            // Left probe is better => maximum is in [a, x2].
            // Discard the right portion: set b = x2.
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + resphi * (b - a);
            f1 = lnl_at(x1)?;
        } else {
            // Right probe is better => maximum is in [x1, b].
            // Discard the left portion: set a = x1.
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = b - resphi * (b - a);
            f2 = lnl_at(x2)?;
        }
    }

    Ok((a + b) / 2.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::newick::parse_newick;

    /// Helper to create a Sequence from a string.
    fn make_seq(name: &str, s: &str) -> Sequence {
        let bases: Vec<Base> = s.chars().map(|c| Base::from_char(c).unwrap()).collect();
        Sequence::new(name, bases)
    }

    /// Helper to create an Alignment from (name, seq_str) pairs.
    fn make_alignment(data: &[(&str, &str)]) -> Alignment {
        let seqs: Vec<Sequence> = data.iter().map(|(n, s)| make_seq(n, s)).collect();
        Alignment::new(seqs).unwrap()
    }

    // ---------------------------------------------------------------
    // Tests for is_transition
    // ---------------------------------------------------------------

    #[test]
    fn test_is_transition_ag() {
        assert!(is_transition(Base::A, Base::G));
        assert!(is_transition(Base::G, Base::A));
    }

    #[test]
    fn test_is_transition_ct() {
        assert!(is_transition(Base::C, Base::T));
        assert!(is_transition(Base::T, Base::C));
    }

    #[test]
    fn test_transversions_are_not_transitions() {
        // A<->C, A<->T, G<->C, G<->T are transversions
        assert!(!is_transition(Base::A, Base::C));
        assert!(!is_transition(Base::C, Base::A));
        assert!(!is_transition(Base::A, Base::T));
        assert!(!is_transition(Base::T, Base::A));
        assert!(!is_transition(Base::G, Base::C));
        assert!(!is_transition(Base::C, Base::G));
        assert!(!is_transition(Base::G, Base::T));
        assert!(!is_transition(Base::T, Base::G));
    }

    #[test]
    fn test_same_base_not_transition() {
        assert!(!is_transition(Base::A, Base::A));
        assert!(!is_transition(Base::C, Base::C));
        assert!(!is_transition(Base::G, Base::G));
        assert!(!is_transition(Base::T, Base::T));
    }

    // ---------------------------------------------------------------
    // Tests for count_ts_tv
    // ---------------------------------------------------------------

    #[test]
    fn test_count_ts_tv_basic() {
        // A->G is transition, A->C is transversion
        let seq1 = make_seq("s1", "AACG");
        let seq2 = make_seq("s2", "GCGT");
        // site 0: A->G = transition
        // site 1: A->C = transversion
        // site 2: C->G = transversion
        // site 3: G->T = transversion
        let (ts, tv) = count_ts_tv(&seq1, &seq2);
        assert_eq!(ts, 1, "expected 1 transition");
        assert_eq!(tv, 3, "expected 3 transversions");
    }

    #[test]
    fn test_count_ts_tv_all_transitions() {
        // A<->G and C<->T
        let seq1 = make_seq("s1", "AACCGG");
        let seq2 = make_seq("s2", "GGTTAA");
        // site 0: A->G ts
        // site 1: A->G ts
        // site 2: C->T ts
        // site 3: C->T ts
        // site 4: G->A ts
        // site 5: G->A ts
        let (ts, tv) = count_ts_tv(&seq1, &seq2);
        assert_eq!(ts, 6);
        assert_eq!(tv, 0);
    }

    #[test]
    fn test_count_ts_tv_all_transversions() {
        let seq1 = make_seq("s1", "AACCGGTT");
        let seq2 = make_seq("s2", "CCAATTGG");
        // A->C tv, A->C tv, C->A tv, C->A tv, G->T tv, G->T tv, T->G tv, T->G tv
        let (ts, tv) = count_ts_tv(&seq1, &seq2);
        assert_eq!(ts, 0);
        assert_eq!(tv, 8);
    }

    #[test]
    fn test_count_ts_tv_identical_sequences() {
        let seq1 = make_seq("s1", "ACGTACGT");
        let seq2 = make_seq("s2", "ACGTACGT");
        let (ts, tv) = count_ts_tv(&seq1, &seq2);
        assert_eq!(ts, 0);
        assert_eq!(tv, 0);
    }

    #[test]
    fn test_count_ts_tv_with_gaps() {
        // Gaps should be skipped
        let seq1 = make_seq("s1", "A-GT");
        let seq2 = make_seq("s2", "G-CA");
        // site 0: A->G ts
        // site 1: gap, skipped
        // site 2: G->C tv
        // site 3: T->A tv
        let (ts, tv) = count_ts_tv(&seq1, &seq2);
        assert_eq!(ts, 1);
        assert_eq!(tv, 2);
    }

    #[test]
    fn test_count_ts_tv_with_unknown() {
        // N (unknown) should be skipped
        let seq1 = make_seq("s1", "ANGT");
        let seq2 = make_seq("s2", "GNCA");
        // site 0: A->G ts
        // site 1: N, skipped
        // site 2: G->C tv
        // site 3: T->A tv
        let (ts, tv) = count_ts_tv(&seq1, &seq2);
        assert_eq!(ts, 1);
        assert_eq!(tv, 2);
    }

    #[test]
    fn test_count_ts_tv_both_gaps() {
        let seq1 = make_seq("s1", "---");
        let seq2 = make_seq("s2", "---");
        let (ts, tv) = count_ts_tv(&seq1, &seq2);
        assert_eq!(ts, 0);
        assert_eq!(tv, 0);
    }

    #[test]
    fn test_count_ts_tv_one_gap_one_base() {
        let seq1 = make_seq("s1", "A-C");
        let seq2 = make_seq("s2", "-GT");
        // site 0: A vs Gap -> skipped
        // site 1: Gap vs G -> skipped
        // site 2: C vs T -> transition
        let (ts, tv) = count_ts_tv(&seq1, &seq2);
        assert_eq!(ts, 1);
        assert_eq!(tv, 0);
    }

    #[test]
    fn test_count_ts_tv_different_length_sequences() {
        // Only overlapping portion compared
        let seq1 = make_seq("s1", "ACGT");
        let seq2 = make_seq("s2", "GC");
        // site 0: A->G ts
        // site 1: C->C identical
        let (ts, tv) = count_ts_tv(&seq1, &seq2);
        assert_eq!(ts, 1);
        assert_eq!(tv, 0);
    }

    #[test]
    fn test_count_ts_tv_empty_sequences() {
        let seq1 = Sequence::new("s1", vec![]);
        let seq2 = Sequence::new("s2", vec![]);
        let (ts, tv) = count_ts_tv(&seq1, &seq2);
        assert_eq!(ts, 0);
        assert_eq!(tv, 0);
    }

    // ---------------------------------------------------------------
    // Tests for estimate_ttratio_counting
    // ---------------------------------------------------------------

    #[test]
    fn test_estimate_ttratio_counting_basic() {
        // Create data with known ts/tv counts
        // seq1 vs seq2: A->G (ts), C->T (ts), G->C (tv), T->A (tv) => 2ts, 2tv
        // seq1 vs seq3: A->A (same), C->C (same), G->A (ts), T->C (ts) => 2ts, 0tv
        // seq2 vs seq3: G->A (ts), T->C (ts), C->A (tv), A->C (tv) => 2ts, 2tv (wait, need to check)
        // Let me use a cleaner example:
        let aln = make_alignment(&[
            ("s1", "AG"),
            ("s2", "GA"),
        ]);
        // s1 vs s2: A->G (ts), G->A (ts) => 2ts, 0tv
        let ratio = estimate_ttratio_counting(&aln);
        // No transversions => default 2.0
        assert!((ratio - DEFAULT_TTRATIO).abs() < 1e-10);
    }

    #[test]
    fn test_estimate_ttratio_counting_mixed() {
        // Build alignment with known ratio
        // We want 4 transitions and 2 transversions => ratio = 2.0
        let aln = make_alignment(&[
            ("s1", "AACCGG"),
            ("s2", "GGTTCA"),
            // s1 vs s2:
            //   A->G ts, A->G ts, C->T ts, C->T ts, G->C tv, G->A ts
            // That's 5ts, 1tv => ratio = 5.0
        ]);
        let ratio = estimate_ttratio_counting(&aln);
        assert!(
            (ratio - 5.0).abs() < 1e-10,
            "expected ratio 5.0, got {}",
            ratio
        );
    }

    #[test]
    fn test_estimate_ttratio_counting_identical_sequences() {
        let aln = make_alignment(&[
            ("s1", "ACGTACGT"),
            ("s2", "ACGTACGT"),
            ("s3", "ACGTACGT"),
        ]);
        // No differences at all => default
        let ratio = estimate_ttratio_counting(&aln);
        assert!((ratio - DEFAULT_TTRATIO).abs() < 1e-10);
    }

    #[test]
    fn test_estimate_ttratio_counting_all_transversions() {
        // Only transversions: A<->C, A<->T, G<->C, G<->T
        let aln = make_alignment(&[
            ("s1", "AAGG"),
            ("s2", "CCTT"),
        ]);
        // A->C tv, A->C tv, G->T tv, G->T tv => 0ts, 4tv
        let ratio = estimate_ttratio_counting(&aln);
        // 0/4 = 0.0, but clamped to MIN_TTRATIO = 0.5
        assert!(
            (ratio - MIN_TTRATIO).abs() < 1e-10,
            "expected clamped to MIN_TTRATIO, got {}",
            ratio
        );
    }

    #[test]
    fn test_estimate_ttratio_counting_with_gaps() {
        // Gaps should be excluded
        let aln = make_alignment(&[
            ("s1", "A-GT"),
            ("s2", "G-CA"),
        ]);
        // site 0: A->G ts
        // site 1: skipped (gap)
        // site 2: G->C tv
        // site 3: T->A tv
        // => 1ts / 2tv = 0.5
        let ratio = estimate_ttratio_counting(&aln);
        assert!(
            (ratio - 0.5).abs() < 1e-10,
            "expected 0.5, got {}",
            ratio
        );
    }

    #[test]
    fn test_estimate_ttratio_counting_single_sequence() {
        let aln = Alignment::new(vec![make_seq("s1", "ACGT")]).unwrap();
        // Only one sequence, no pairs
        let ratio = estimate_ttratio_counting(&aln);
        assert!((ratio - DEFAULT_TTRATIO).abs() < 1e-10);
    }

    // ---------------------------------------------------------------
    // Tests for estimate_base_frequencies
    // ---------------------------------------------------------------

    #[test]
    fn test_estimate_base_frequencies_equal() {
        let aln = make_alignment(&[("s1", "ACGT")]);
        let freqs = estimate_base_frequencies(&aln);
        assert!((freqs.freq_a - 0.25).abs() < 1e-10);
        assert!((freqs.freq_c - 0.25).abs() < 1e-10);
        assert!((freqs.freq_g - 0.25).abs() < 1e-10);
        assert!((freqs.freq_t - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_estimate_base_frequencies_biased() {
        let aln = make_alignment(&[("s1", "AAAACGT")]);
        // 4A, 1C, 1G, 1T out of 7
        let freqs = estimate_base_frequencies(&aln);
        assert!((freqs.freq_a - 4.0 / 7.0).abs() < 1e-10);
        assert!((freqs.freq_c - 1.0 / 7.0).abs() < 1e-10);
        assert!((freqs.freq_g - 1.0 / 7.0).abs() < 1e-10);
        assert!((freqs.freq_t - 1.0 / 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_estimate_base_frequencies_sum_to_one() {
        let aln = make_alignment(&[
            ("s1", "AACCGGTTAACCGGTT"),
            ("s2", "AAAAAACCCCCCGGGG"),
        ]);
        let freqs = estimate_base_frequencies(&aln);
        let total = freqs.freq_a + freqs.freq_c + freqs.freq_g + freqs.freq_t;
        assert!(
            (total - 1.0).abs() < 1e-10,
            "frequencies should sum to 1.0, got {}",
            total
        );
    }

    #[test]
    fn test_estimate_base_frequencies_with_gaps() {
        // Gaps are excluded from counting
        let aln = make_alignment(&[("s1", "AA--CCGG")]);
        // 2A, 2C, 2G out of 6 definite bases (no T)
        let freqs = estimate_base_frequencies(&aln);
        let total = freqs.freq_a + freqs.freq_c + freqs.freq_g + freqs.freq_t;
        assert!((total - 1.0).abs() < 1e-10);
        assert!(
            (freqs.freq_a - 2.0 / 6.0).abs() < 1e-10,
            "expected freq_a=2/6, got {}",
            freqs.freq_a
        );
        assert!(
            (freqs.freq_t).abs() < 1e-10,
            "expected freq_t=0, got {}",
            freqs.freq_t
        );
    }

    #[test]
    fn test_estimate_base_frequencies_all_gaps() {
        // All gaps => fall back to equal frequencies
        let aln = make_alignment(&[("s1", "----")]);
        let freqs = estimate_base_frequencies(&aln);
        assert!((freqs.freq_a - 0.25).abs() < 1e-10);
        assert!((freqs.freq_c - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_estimate_base_frequencies_multiple_sequences() {
        let aln = make_alignment(&[
            ("s1", "AAAA"),
            ("s2", "CCCC"),
            ("s3", "GGGG"),
            ("s4", "TTTT"),
        ]);
        // 16 total: 4 each
        let freqs = estimate_base_frequencies(&aln);
        assert!((freqs.freq_a - 0.25).abs() < 1e-10);
        assert!((freqs.freq_c - 0.25).abs() < 1e-10);
        assert!((freqs.freq_g - 0.25).abs() < 1e-10);
        assert!((freqs.freq_t - 0.25).abs() < 1e-10);
    }

    // ---------------------------------------------------------------
    // Tests for estimate_ttratio_ml
    // ---------------------------------------------------------------

    #[test]
    fn test_estimate_ttratio_ml_basic() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let freqs = estimate_base_frequencies(&aln);
        let ratio = estimate_ttratio_ml(&tree, &aln, &freqs).unwrap();
        assert!(
            ratio >= MIN_TTRATIO && ratio <= MAX_TTRATIO,
            "ML ttratio {} should be in [{}, {}]",
            ratio,
            MIN_TTRATIO,
            MAX_TTRATIO
        );
    }

    #[test]
    fn test_estimate_ttratio_ml_with_transitions() {
        // Data enriched with transitions (A<->G, C<->T) should yield
        // higher ttratio than data enriched with transversions.
        //
        // We include a mix of both types in each dataset so that the
        // ML optimum falls in the interior of the search range rather
        // than saturating at a boundary.

        // Transition-heavy: 8 transitions, 2 transversions per pair
        let aln_ts = make_alignment(&[
            ("A", "AAAAGGGGCCAAAAGGGGCC"),
            ("B", "GGGGAAAATTGGGGAAAATT"), // A<->G ts, A<->G ts, ..., C<->T ts, C<->T ts, + 2 tv
            ("C", "AAAAGGGGCCAAAAGGGGCC"),
        ]);
        // Transversion-heavy: 2 transitions, 8 transversions per pair
        let aln_tv = make_alignment(&[
            ("A", "AACCAAGGCCAACCAAGGCC"),
            ("B", "CCAACCTTAACCAACCTTAA"), // A<->C tv, A<->C tv, ..., + some ts
            ("C", "AACCAAGGCCAACCAAGGCC"),
        ]);
        let tree = parse_newick("((A:0.3,B:0.3):0.1,C:0.2);").unwrap();

        let freqs_ts = estimate_base_frequencies(&aln_ts);
        let freqs_tv = estimate_base_frequencies(&aln_tv);

        let ratio_ts = estimate_ttratio_ml(&tree, &aln_ts, &freqs_ts).unwrap();
        let ratio_tv = estimate_ttratio_ml(&tree, &aln_tv, &freqs_tv).unwrap();

        assert!(
            ratio_ts > ratio_tv,
            "transition-enriched data should yield higher ttratio ({}) than transversion-enriched ({})",
            ratio_ts,
            ratio_tv
        );
    }

    #[test]
    fn test_estimate_ttratio_ml_taxa_mismatch() {
        let aln = make_alignment(&[
            ("X", "ACGT"),
            ("Y", "ACGT"),
            ("Z", "ACGT"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let freqs = BaseFrequencies::equal();
        let result = estimate_ttratio_ml(&tree, &aln, &freqs);
        assert!(result.is_err(), "should fail on taxa mismatch");
    }

    #[test]
    fn test_estimate_ttratio_ml_returns_finite() {
        let aln = make_alignment(&[
            ("A", "ACGTACGTACGT"),
            ("B", "ACGAACGAACGA"),
            ("C", "TGCATGCATGCA"),
        ]);
        let tree = parse_newick("((A:0.15,B:0.15):0.1,C:0.25);").unwrap();
        let freqs = estimate_base_frequencies(&aln);
        let ratio = estimate_ttratio_ml(&tree, &aln, &freqs).unwrap();
        assert!(ratio.is_finite(), "ML ttratio should be finite, got {}", ratio);
        assert!(ratio > 0.0, "ML ttratio should be positive, got {}", ratio);
    }

    #[test]
    fn test_estimate_ttratio_ml_with_equal_freqs() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "ACGAACGAACGAACGA"),
            ("C", "TGCATGCATGCATGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let freqs = BaseFrequencies::equal();
        let ratio = estimate_ttratio_ml(&tree, &aln, &freqs).unwrap();
        assert!(
            ratio >= MIN_TTRATIO && ratio <= MAX_TTRATIO,
            "ML ttratio with equal freqs should be in range, got {}",
            ratio
        );
    }

    #[test]
    fn test_estimate_ttratio_ml_four_taxa() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGG"),
        ]);
        let tree =
            parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();
        let freqs = estimate_base_frequencies(&aln);
        let ratio = estimate_ttratio_ml(&tree, &aln, &freqs).unwrap();
        assert!(ratio.is_finite());
        assert!(ratio > 0.0);
    }

    #[test]
    fn test_estimate_ttratio_ml_improves_likelihood() {
        // The ML ttratio should yield a likelihood at least as good as
        // several fixed ratios across the interior of the search range.
        // (Boundary effects at the bracket edges may cause tiny discrepancies.)
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let freqs = estimate_base_frequencies(&aln);

        let ml_ratio = estimate_ttratio_ml(&tree, &aln, &freqs).unwrap();
        let model_ml = F84Model::new(freqs, ml_ratio);
        let lnl_ml = log_likelihood(&tree, &aln, &model_ml).unwrap();

        // ML estimate should be at least as good as interior probe points
        // (not at the bracket boundaries where numerical precision matters)
        for &probe in &[1.0, 2.0, 5.0, 10.0, 15.0] {
            let model_probe = F84Model::new(freqs, probe);
            let lnl_probe = log_likelihood(&tree, &aln, &model_probe).unwrap();
            assert!(
                lnl_ml >= lnl_probe - 1e-3,
                "ML ttratio ({}, lnl={}) should be >= probe ttratio ({}, lnl={})",
                ml_ratio,
                lnl_ml,
                probe,
                lnl_probe
            );
        }
    }

    #[test]
    fn test_estimate_ttratio_ml_with_gaps_in_alignment() {
        let aln = make_alignment(&[
            ("A", "A-CGTACGT"),
            ("B", "ACGAAC-GT"),
            ("C", "TGCATG-CA"),
        ]);
        let tree = parse_newick("((A:0.15,B:0.15):0.1,C:0.25);").unwrap();
        let freqs = estimate_base_frequencies(&aln);
        let ratio = estimate_ttratio_ml(&tree, &aln, &freqs).unwrap();
        assert!(ratio.is_finite());
        assert!(ratio > 0.0);
    }

    // ---------------------------------------------------------------
    // Integration / consistency tests
    // ---------------------------------------------------------------

    #[test]
    fn test_counting_and_ml_agree_direction() {
        // When data is heavily transition-biased, both methods should
        // give a ratio > 2.0
        let aln = make_alignment(&[
            ("A", "AAAAAAAACCCCCCCCAAAAAAAACCCCCCCC"),
            ("B", "GGGGGGGGTTTTTTTTUUUUUUUUTTTTTTTT"),
            ("C", "AAAAAAAACCCCCCCCAAAAAAAACCCCCCCC"),
        ]);
        let tree = parse_newick("((A:0.3,B:0.3):0.1,C:0.2);").unwrap();
        let freqs = estimate_base_frequencies(&aln);

        let counting_ratio = estimate_ttratio_counting(&aln);
        let ml_ratio = estimate_ttratio_ml(&tree, &aln, &freqs).unwrap();

        // Both should be high (transition-biased data)
        assert!(
            counting_ratio > 2.0,
            "counting ratio should be > 2 for transition-biased data, got {}",
            counting_ratio
        );
        // ML may differ but should also reflect high ts/tv
        assert!(
            ml_ratio > 1.0,
            "ML ratio should be > 1 for transition-biased data, got {}",
            ml_ratio
        );
    }

    #[test]
    fn test_base_frequencies_consistent_with_counting() {
        let aln = make_alignment(&[
            ("A", "AAACGT"),
            ("B", "AAACGT"),
            ("C", "AAACGT"),
        ]);
        let freqs = estimate_base_frequencies(&aln);
        // 3 sequences x 6 sites = 18 bases total
        // A: 3*3=9, C: 3*1=3, G: 3*1=3, T: 3*1=3
        assert!((freqs.freq_a - 9.0 / 18.0).abs() < 1e-10);
        assert!((freqs.freq_c - 3.0 / 18.0).abs() < 1e-10);
        assert!((freqs.freq_g - 3.0 / 18.0).abs() < 1e-10);
        assert!((freqs.freq_t - 3.0 / 18.0).abs() < 1e-10);
    }
}

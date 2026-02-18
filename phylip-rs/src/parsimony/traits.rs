//! Pluggable parsimony scoring criteria.
//!
//! Provides the [`ParsimonyScorer`] trait that abstracts the state-combination
//! rules used during the postorder (upward) pass of parsimony scoring.
//! Different parsimony criteria—Fitch (unordered), Dollo (irreversible loss),
//! Camin-Sokal (irreversible gain)—implement this trait to supply their own
//! rules for combining child state sets at internal nodes.
//!
//! # Example
//!
//! ```
//! use phylip_rs::parsimony::traits::{ParsimonyScorer, FitchScorer};
//! use phylip_rs::parsimony::wagner::StateSet;
//! use phylip_rs::tree::types::Base;
//!
//! let scorer = FitchScorer;
//! let left = scorer.leaf_state(Base::A);
//! let right = scorer.leaf_state(Base::C);
//! let (combined, cost) = scorer.combine(left, right);
//! // A and C have empty intersection, so Fitch takes union with cost 1
//! assert_eq!(cost, 1);
//! assert!(!combined.intersection(StateSet::A).is_empty());
//! assert!(!combined.intersection(StateSet::C).is_empty());
//! ```
//!
//! # References
//!
//! - Fitch, W.M. (1971). Toward defining the course of evolution: minimum
//!   change for a specific tree topology. *Systematic Zoology*, 20, 406-416.
//! - Le Quesne, W.J. (1974). The uniquely evolved character concept and its
//!   cladistic application. *Systematic Zoology*, 23, 513-517.
//! - Camin, J.H. & Sokal, R.R. (1965). A method for deducing branching
//!   sequences in phylogeny. *Evolution*, 19, 311-326.

use crate::parsimony::wagner::StateSet;
use crate::tree::types::Base;

/// Trait for pluggable parsimony scoring criteria.
///
/// Different criteria (Fitch, Dollo, Camin-Sokal) implement this trait
/// to provide their own state combination rules for the postorder
/// (upward) pass of the Fitch-style dynamic programming algorithm.
pub trait ParsimonyScorer {
    /// Initialize a state set at a leaf from an observed base.
    fn leaf_state(&self, base: Base) -> StateSet;

    /// Combine two child state sets at an internal node.
    ///
    /// Returns `(combined_state_set, additional_steps_cost)`.
    /// The cost is the number of extra evolutionary steps implied
    /// by the combination under this parsimony criterion.
    fn combine(&self, left: StateSet, right: StateSet) -> (StateSet, usize);

    /// Name of the parsimony criterion (e.g., "Fitch", "Dollo").
    fn name(&self) -> &str;
}

/// Fitch (1971) unordered parsimony scorer.
///
/// The standard maximum parsimony criterion for unordered characters:
/// - Intersection of children's state sets if non-empty (0 extra steps).
/// - Union of children's state sets if intersection is empty (+1 step).
///
/// This is symmetric: gains and losses are equally likely.
pub struct FitchScorer;

impl ParsimonyScorer for FitchScorer {
    fn leaf_state(&self, base: Base) -> StateSet {
        StateSet::from_base(base)
    }

    fn combine(&self, left: StateSet, right: StateSet) -> (StateSet, usize) {
        let isect = left.intersection(right);
        if isect.is_empty() {
            (left.union(right), 1)
        } else {
            (isect, 0)
        }
    }

    fn name(&self) -> &str {
        "Fitch"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parsimony::wagner::{parsimony_score, parsimony_score_with};
    use crate::tree::newick::parse_newick;
    use crate::tree::types::{Alignment, Sequence};

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

    // --- FitchScorer unit tests ---

    #[test]
    fn test_fitch_leaf_state() {
        let scorer = FitchScorer;
        assert_eq!(scorer.leaf_state(Base::A), StateSet::A);
        assert_eq!(scorer.leaf_state(Base::C), StateSet::C);
        assert_eq!(scorer.leaf_state(Base::G), StateSet::G);
        assert_eq!(scorer.leaf_state(Base::T), StateSet::T);
        assert_eq!(scorer.leaf_state(Base::Gap), StateSet::GAP);
    }

    #[test]
    fn test_fitch_combine_intersection() {
        let scorer = FitchScorer;
        // Same state: intersection is non-empty, cost 0
        let (combined, cost) = scorer.combine(StateSet::A, StateSet::A);
        assert_eq!(combined, StateSet::A);
        assert_eq!(cost, 0);
    }

    #[test]
    fn test_fitch_combine_union() {
        let scorer = FitchScorer;
        // Different states: intersection is empty, take union, cost 1
        let (combined, cost) = scorer.combine(StateSet::A, StateSet::C);
        assert_eq!(combined, StateSet::A.union(StateSet::C));
        assert_eq!(cost, 1);
    }

    #[test]
    fn test_fitch_combine_with_ambiguity() {
        let scorer = FitchScorer;
        // {A,C} and {C,G}: intersection = {C}, cost 0
        let ac = StateSet::A.union(StateSet::C);
        let cg = StateSet::C.union(StateSet::G);
        let (combined, cost) = scorer.combine(ac, cg);
        assert_eq!(combined, StateSet::C);
        assert_eq!(cost, 0);
    }

    #[test]
    fn test_fitch_name() {
        let scorer = FitchScorer;
        assert_eq!(scorer.name(), "Fitch");
    }

    // --- FitchScorer produces same results as original parsimony_score ---

    #[test]
    fn test_fitch_scorer_matches_original_identical_seqs() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "ACGT"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let (orig_score, orig_steps) = parsimony_score(&tree, &aln).unwrap();
        let (new_score, new_steps) =
            parsimony_score_with(&tree, &aln, &FitchScorer).unwrap();
        assert_eq!(orig_score, new_score);
        assert_eq!(orig_steps, new_steps);
    }

    #[test]
    fn test_fitch_scorer_matches_original_one_change() {
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "A"),
            ("C", "T"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let (orig_score, orig_steps) = parsimony_score(&tree, &aln).unwrap();
        let (new_score, new_steps) =
            parsimony_score_with(&tree, &aln, &FitchScorer).unwrap();
        assert_eq!(orig_score, new_score);
        assert_eq!(orig_steps, new_steps);
    }

    #[test]
    fn test_fitch_scorer_matches_original_four_taxa() {
        let aln = make_alignment(&[
            ("A", "AACC"),
            ("B", "AACC"),
            ("C", "CCAA"),
            ("D", "CCAA"),
        ]);
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let (orig_score, orig_steps) = parsimony_score(&tree, &aln).unwrap();
        let (new_score, new_steps) =
            parsimony_score_with(&tree, &aln, &FitchScorer).unwrap();
        assert_eq!(orig_score, new_score);
        assert_eq!(orig_steps, new_steps);
    }

    #[test]
    fn test_fitch_scorer_matches_original_all_different() {
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "C"),
            ("C", "G"),
            ("D", "T"),
        ]);
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let (orig_score, _) = parsimony_score(&tree, &aln).unwrap();
        let (new_score, _) =
            parsimony_score_with(&tree, &aln, &FitchScorer).unwrap();
        assert_eq!(orig_score, new_score);
    }

    #[test]
    fn test_fitch_scorer_matches_original_multisite() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGTACGT"),
            ("C", "TGCATGCA"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let (orig_score, orig_steps) = parsimony_score(&tree, &aln).unwrap();
        let (new_score, new_steps) =
            parsimony_score_with(&tree, &aln, &FitchScorer).unwrap();
        assert_eq!(orig_score, new_score);
        assert_eq!(orig_steps, new_steps);
    }
}

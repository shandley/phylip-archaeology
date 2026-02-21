//! Dollo parsimony for binary characters.
//!
//! Implements Dollo parsimony (Le Quesne 1974, Farris 1977), in which the
//! derived state (1) can arise only **once** on the tree but can be lost
//! (reversed to 0) any number of times.  The parsimony score under Dollo's
//! criterion is the total number of loss (reversal) events across all
//! characters.
//!
//! ## Binary encoding
//!
//! Binary character data is represented using DNA bases:
//! - `Base::A` encodes state **0** (ancestral / absent)
//! - `Base::T` encodes state **1** (derived / present)
//!
//! ## Algorithm
//!
//! During the postorder (upward) pass the Dollo rule is:
//! - If **either** child has state 1, the parent must be 1 (because the
//!   derived state cannot be re-gained once lost).
//! - The cost at a node equals the number of children that have state 0
//!   when the parent is forced to 1 (each such child implies one loss).
//! - If both children are 0, the parent is 0 with cost 0.
//!
//! # Example
//!
//! ```
//! use phylip_rs::parsimony::dollo::{dollo_parsimony, DolloParsimonyResult};
//! use phylip_rs::tree::{Alignment, Base, Sequence, Tree};
//! use phylip_rs::tree::newick::parse_newick;
//!
//! // Binary data: A=0, T=1
//! let alignment = Alignment::new(vec![
//!     Sequence::new("A", vec![Base::T]),
//!     Sequence::new("B", vec![Base::T]),
//!     Sequence::new("C", vec![Base::A]),
//! ]).unwrap();
//!
//! let tree = parse_newick("((A,B),C);").unwrap();
//! let result = dollo_parsimony(&tree, &alignment);
//! // One gain at MRCA of A,B; one loss on branch to C? No — C is already 0.
//! // The gain is at the (A,B) ancestor. C has 0. The root must be 1
//! // because child (A,B) has 1. That forces a loss on the C branch.
//! assert_eq!(result.score, 1); // one loss (C's branch)
//! assert_eq!(result.gains, 1); // one character with derived state
//! ```
//!
//! # References
//!
//! - Le Quesne, W.J. (1974). The uniquely evolved character concept and
//!   its cladistic application. *Systematic Zoology*, 23, 513-517.
//! - Farris, J.S. (1977). Phylogenetic analysis under Dollo's law.
//!   *Systematic Zoology*, 26, 77-88.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//!   Chapter 8.

use crate::parsimony::traits::ParsimonyScorer;
use crate::parsimony::wagner::{parsimony_score_with, search_with, ParsimonyResult, StateSet};
use crate::tree::types::{Alignment, Base, Tree};

/// Dollo parsimony scorer.
///
/// Under Dollo's assumption of irreversible loss, the derived state (1)
/// can originate only once.  Any descendant of the gaining ancestor that
/// lacks the trait must have lost it independently.
///
/// In terms of state-set combination:
/// - If **either** child is {T} (=1), the parent must be {T}.
/// - Cost = number of children that are {A} (=0) when the parent is
///   forced to {T} (each is one reversal / loss event).
/// - If **both** children are {A} (=0), the parent is {A} with 0 cost.
///
/// **Limitation:** This scorer implements only the upward (postorder) pass,
/// which always propagates the derived state to the root.  PHYLIP's Dollo
/// algorithm uses a two-pass approach: an upward pass followed by a downward
/// correction that places the gain at the MRCA of derived-state taxa rather
/// than the root.  The upward-only approach can overcount losses when the
/// optimal gain placement is below the root.  A future improvement would
/// add the downward correction pass to match PHYLIP's exact scoring.
pub struct DolloScorer;

impl ParsimonyScorer for DolloScorer {
    fn leaf_state(&self, base: Base) -> StateSet {
        StateSet::from_base(base)
    }

    fn combine(&self, left: StateSet, right: StateSet) -> (StateSet, usize) {
        let left_has_one = !left.intersection(StateSet::T).is_empty();
        let right_has_one = !right.intersection(StateSet::T).is_empty();

        if left_has_one || right_has_one {
            // Parent must be 1 (derived state cannot re-arise).
            // Each child that is 0 represents a loss event.
            let mut cost = 0;
            if !left_has_one {
                cost += 1;
            }
            if !right_has_one {
                cost += 1;
            }
            (StateSet::T, cost)
        } else {
            // Both children are 0; parent is 0.
            (StateSet::A, 0)
        }
    }

    fn name(&self) -> &str {
        "Dollo"
    }
}

/// Result of a Dollo parsimony analysis.
#[derive(Debug, Clone)]
pub struct DolloParsimonyResult {
    /// The best tree found (or the input tree, for evaluation).
    pub tree: Tree,
    /// The Dollo parsimony score (total number of loss events).
    pub score: usize,
    /// Number of characters that have at least one taxon with the
    /// derived state (each such character has exactly one gain under
    /// Dollo's assumption).
    pub gains: usize,
    /// Per-site loss counts.
    pub site_steps: Vec<usize>,
}

/// Evaluate a given tree under Dollo parsimony.
///
/// Returns the total number of loss events, the number of characters
/// with at least one derived-state taxon (= number of gains), and
/// per-site loss counts.
pub fn dollo_parsimony(tree: &Tree, alignment: &Alignment) -> DolloParsimonyResult {
    let scorer = DolloScorer;
    let (score, site_steps) = parsimony_score_with(tree, alignment, &scorer)
        .expect("Dollo parsimony scoring failed");

    // Count gains: each character (site) that has at least one taxon
    // with state T (=1) has exactly one gain under Dollo.
    let nsites = alignment.nsites();
    let mut gains = 0;
    for site in 0..nsites {
        let has_derived = (0..alignment.ntaxa())
            .any(|taxon| alignment.base(taxon, site) == Base::T);
        if has_derived {
            gains += 1;
        }
    }

    DolloParsimonyResult {
        tree: tree.clone(),
        score,
        gains,
        site_steps,
    }
}

/// Search for the most parsimonious tree under Dollo parsimony.
///
/// Uses stepwise addition + SPR rearrangement (same heuristic as
/// Wagner/Fitch search) but scores trees under the Dollo criterion.
///
/// # Arguments
/// * `alignment` — the binary-encoded alignment (A=0, T=1)
/// * `seed` — optional random seed for taxon addition order
pub fn dollo_search(alignment: &Alignment, seed: Option<u64>) -> DolloParsimonyResult {
    let scorer = DolloScorer;
    let result: ParsimonyResult = search_with(alignment, &scorer, seed);

    // Count gains
    let nsites = alignment.nsites();
    let mut gains = 0;
    for site in 0..nsites {
        let has_derived = (0..alignment.ntaxa())
            .any(|taxon| alignment.base(taxon, site) == Base::T);
        if has_derived {
            gains += 1;
        }
    }

    DolloParsimonyResult {
        tree: result.tree,
        score: result.score,
        gains,
        site_steps: result.site_steps,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parsimony::wagner::parsimony_score;
    use crate::tree::newick::parse_newick;
    use crate::tree::types::Sequence;

    /// Helper: build an alignment from (name, sequence_str) pairs.
    /// Uses A=0, T=1 for binary data.
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

    // --- DolloScorer unit tests ---

    #[test]
    fn test_dollo_combine_both_zero() {
        let scorer = DolloScorer;
        let (state, cost) = scorer.combine(StateSet::A, StateSet::A);
        assert_eq!(state, StateSet::A);
        assert_eq!(cost, 0);
    }

    #[test]
    fn test_dollo_combine_both_one() {
        let scorer = DolloScorer;
        let (state, cost) = scorer.combine(StateSet::T, StateSet::T);
        assert_eq!(state, StateSet::T);
        assert_eq!(cost, 0);
    }

    #[test]
    fn test_dollo_combine_mixed() {
        let scorer = DolloScorer;
        // Left=1, right=0: parent must be 1, one loss on right
        let (state, cost) = scorer.combine(StateSet::T, StateSet::A);
        assert_eq!(state, StateSet::T);
        assert_eq!(cost, 1);

        // Left=0, right=1: parent must be 1, one loss on left
        let (state, cost) = scorer.combine(StateSet::A, StateSet::T);
        assert_eq!(state, StateSet::T);
        assert_eq!(cost, 1);
    }

    // --- Single character with one gain ---

    #[test]
    fn test_dollo_single_char_one_gain_no_losses() {
        // All taxa have "1" -> one gain at root, no losses
        let aln = make_alignment(&[
            ("A", "T"),
            ("B", "T"),
            ("C", "T"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let result = dollo_parsimony(&tree, &aln);
        assert_eq!(result.score, 0); // no losses
        assert_eq!(result.gains, 1); // one character with derived state
    }

    #[test]
    fn test_dollo_one_taxon_has_derived() {
        // Only A has "1", B and C have "0"
        // On tree ((A,B),C):
        //   Internal (A,B): A=T, B=A -> parent=T, cost=1 (loss on B)
        //   Root: (A,B)=T, C=A -> parent=T, cost=1 (loss on C)
        //   Total losses = 2
        let aln = make_alignment(&[
            ("A", "T"),
            ("B", "A"),
            ("C", "A"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let result = dollo_parsimony(&tree, &aln);
        assert_eq!(result.score, 2); // losses on B and C branches
        assert_eq!(result.gains, 1);
    }

    #[test]
    fn test_dollo_two_of_three_have_derived() {
        // A=1, B=1, C=0
        // On tree ((A,B),C):
        //   Internal (A,B): T & T -> T, cost=0
        //   Root: T & A -> T, cost=1 (loss on C)
        //   Total losses = 1
        let aln = make_alignment(&[
            ("A", "T"),
            ("B", "T"),
            ("C", "A"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let result = dollo_parsimony(&tree, &aln);
        assert_eq!(result.score, 1);
        assert_eq!(result.gains, 1);
    }

    #[test]
    fn test_dollo_all_ancestral() {
        // All taxa have "0" -> no gains, no losses
        let aln = make_alignment(&[
            ("A", "A"),
            ("B", "A"),
            ("C", "A"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let result = dollo_parsimony(&tree, &aln);
        assert_eq!(result.score, 0);
        assert_eq!(result.gains, 0);
    }

    // --- Dollo vs Fitch comparison ---

    #[test]
    fn test_dollo_vs_fitch_different_scores() {
        // A=1, B=0, C=0, D=1 on ((A,B),(C,D))
        // Fitch: {T}&{A}=empty -> union {A,T}, +1; {A}&{T}=empty -> union, +1;
        //        {A,T}&{A,T}={A,T}, +0 -> total Fitch = 2
        // Dollo: Internal(A,B): T&A -> T, cost=1; Internal(C,D): A&T -> T, cost=1
        //        Root: T&T -> T, cost=0 -> total Dollo = 2
        // Actually same here. Let's pick a case where they differ.
        //
        // A=1, B=1, C=0, D=0 on ((A,C),(B,D))
        // Fitch: {T}&{A}=empty -> {A,T}, +1; {T}&{A}=empty -> {A,T}, +1;
        //        {A,T}&{A,T}={A,T}, +0 -> Fitch = 2
        // Dollo: (A,C): T&A -> T, cost=1; (B,D): T&A -> T, cost=1;
        //        Root: T&T -> T, cost=0 -> Dollo = 2
        //
        // Let's use a topology where Dollo is worse:
        // A=1, B=0, C=1, D=0 on ((A,B),(C,D))
        // Fitch: {T}&{A}=empty -> union, +1; {T}&{A}=empty -> union, +1;
        //        {A,T}&{A,T}={A,T}, +0 -> Fitch = 2
        // Dollo: (A,B): T & A -> T, cost=1; (C,D): T & A -> T, cost=1
        //        Root: T & T -> T, cost=0 -> Dollo = 2
        //
        // For Dollo > Fitch, we need scattered 1s in a topology that Fitch
        // resolves more cheaply. Consider:
        // One taxon with 1, three with 0, on ((A,B),(C,D)) where only A=1
        // Fitch: {T}&{A}=empty -> {A,T}, +1; {A}&{A}={A}, +0;
        //        {A,T}&{A}={A}, +0 -> Fitch = 1
        // Dollo: (A,B): T&A -> T, cost=1; (C,D): A&A -> A, cost=0;
        //        Root: T&A -> T, cost=1 -> Dollo = 2
        let aln = make_alignment(&[
            ("A", "T"),
            ("B", "A"),
            ("C", "A"),
            ("D", "A"),
        ]);
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let dollo_result = dollo_parsimony(&tree, &aln);
        let (fitch_score, _) = parsimony_score(&tree, &aln).unwrap();

        // Dollo = 2 (loss on B, loss on (C,D) subtree)
        // Fitch = 1 (one change somewhere)
        assert_eq!(dollo_result.score, 2);
        assert_eq!(fitch_score, 1);
        assert!(dollo_result.score > fitch_score);
    }

    // --- Multi-character tests ---

    #[test]
    fn test_dollo_multiple_characters() {
        // Two characters:
        //   Site 0: A=1, B=1, C=0 -> 1 loss (C)
        //   Site 1: A=0, B=0, C=0 -> 0 losses, 0 gains
        let aln = make_alignment(&[
            ("A", "TA"),
            ("B", "TA"),
            ("C", "AA"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let result = dollo_parsimony(&tree, &aln);
        assert_eq!(result.score, 1);
        assert_eq!(result.gains, 1); // only site 0 has derived state
        assert_eq!(result.site_steps, vec![1, 0]);
    }

    // --- Known textbook example ---

    #[test]
    fn test_dollo_four_taxa_known() {
        // ((A,B),(C,D)) with A=1, B=1, C=1, D=0
        // (A,B): T&T -> T, cost=0
        // (C,D): T&A -> T, cost=1 (loss on D)
        // Root: T&T -> T, cost=0
        // Total Dollo = 1
        let aln = make_alignment(&[
            ("A", "T"),
            ("B", "T"),
            ("C", "T"),
            ("D", "A"),
        ]);
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let result = dollo_parsimony(&tree, &aln);
        assert_eq!(result.score, 1);
        assert_eq!(result.gains, 1);
    }

    // --- Search test ---

    #[test]
    fn test_dollo_search_basic() {
        // Clear signal: A,B have derived state; C,D are ancestral
        let aln = make_alignment(&[
            ("A", "TTTT"),
            ("B", "TTTT"),
            ("C", "AAAA"),
            ("D", "AAAA"),
        ]);
        let result = dollo_search(&aln, Some(42));
        assert_eq!(result.tree.num_leaves(), 4);
        assert_eq!(result.gains, 4); // all 4 sites have derived state
        // Optimal tree groups A,B together -> 1 loss per site on the
        // (C,D) side = 4 total. But if A,B are grouped and C,D are grouped,
        // the root must be T (because A,B child is T), and C,D child is A,
        // so cost = 1 per site * 4 sites = 4.
        assert_eq!(result.score, 4);
    }

    #[test]
    fn test_dollo_search_all_derived() {
        let aln = make_alignment(&[
            ("A", "TT"),
            ("B", "TT"),
            ("C", "TT"),
        ]);
        let result = dollo_search(&aln, Some(1));
        assert_eq!(result.score, 0); // no losses
        assert_eq!(result.gains, 2);
    }

    #[test]
    fn test_dollo_scorer_name() {
        assert_eq!(DolloScorer.name(), "Dollo");
    }
}

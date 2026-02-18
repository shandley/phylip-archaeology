//! Phylogenetic invariants for topology inference.
//!
//! Phylogenetic invariants are algebraic functions of site-pattern frequencies
//! that equal zero for one tree topology but not others. They provide a
//! statistically consistent method for inferring tree topology that is
//! independent of branch-length estimation.
//!
//! # Methods
//!
//! - **Lake's invariants**: Uses transversion-type parsimony-informative
//!   site patterns to choose among the three unrooted 4-taxon topologies.
//!
//! - **Cavender's invariants**: Uses linear functions of RY-coded site
//!   pattern frequencies under the symmetric (CFN) model. The true
//!   topology's invariant should be zero.
//!
//! # References
//!
//! - Lake, J.A. (1987). A rate-independent technique for analysis of
//!   nucleic acid sequences. *Mol. Biol. Evol.*, 4:167-191.
//! - Cavender, J.A. (1978). Taxonomy with confidence. *Math. Biosci.*,
//!   40:271-280.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*, Chapter 29.

pub mod lake;

pub use lake::{
    cavender_invariants, lake_invariants, CavenderResult, InvariantResult, PatternCounts,
};

//! Maximum likelihood methods for phylogenetic inference.
//!
//! Includes Felsenstein's pruning algorithm (1981), substitution models,
//! and heuristic tree search via stepwise addition with SPR rearrangement.

pub mod models;
pub mod nni;
pub mod pruning;
pub mod search;

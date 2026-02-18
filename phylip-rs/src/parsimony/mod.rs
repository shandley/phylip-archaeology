//! Maximum parsimony methods for phylogenetic inference.
//!
//! This module provides multiple parsimony criteria for phylogenetic
//! analysis:
//!
//! - **Wagner / Fitch** (`wagner`) — standard unordered parsimony where
//!   any character-state change costs one step.
//! - **Dollo** (`dollo`) — the derived state can originate only once;
//!   losses (reversals) are freely allowed.
//! - **Camin-Sokal** (`camin_sokal`) — characters can only gain the
//!   derived state (irreversible gain); parallel gains explain homoplasy.
//! - **Branch-and-bound** (`branch_and_bound`) — exact search that
//!   guarantees finding all globally optimal parsimony trees.
//! - **Multistate** (`multistate`) — parsimony for arbitrary discrete
//!   character alphabets (up to 32 states), with unordered or weighted
//!   step matrices.
//!
//! The [`traits`] module defines the [`ParsimonyScorer`](traits::ParsimonyScorer)
//! trait, which allows the same tree-scoring and tree-search machinery
//! to be reused across all criteria.

pub mod branch_and_bound;
pub mod camin_sokal;
pub mod dollo;
pub mod multistate;
pub mod protein_parsimony;
pub mod traits;
pub mod wagner;

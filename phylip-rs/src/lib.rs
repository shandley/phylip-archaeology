//! # phylip-rs
//!
//! Modern Rust reimplementations of classic phylogenetic algorithms from
//! Joe Felsenstein's PHYLIP (PHYLogeny Inference Package).
//!
//! This crate preserves and modernizes the algorithmic legacy of PHYLIP,
//! one of the most influential bioinformatics software packages ever created.
//! Every implementation is validated against original PHYLIP output to ensure
//! fidelity to the original algorithms.
//!
//! ## Modules
//!
//! - `tree` - Phylogenetic tree data structures and Newick format I/O
//! - `distance` - Distance matrix methods (neighbor-joining, UPGMA, Fitch-Margoliash)
//! - `parsimony` - Maximum parsimony methods
//! - `likelihood` - Maximum likelihood methods (Felsenstein pruning algorithm)
//! - `bootstrap` - Bootstrap resampling for confidence estimation
//! - `compatibility` - Character compatibility analysis and clique finding
//! - `consensus` - Consensus tree construction (strict, majority-rule, extended)
//! - `io` - File format I/O (PHYLIP format, FASTA)
//! - `models` - Nucleotide substitution models (JC69, K2P, F81, F84)
//! - `comparative` - Continuous character comparative methods (independent contrasts, contml)
//! - `biogeography` - Dispersal-Vicariance Analysis (DIVA, Ronquist 1997)
//! - `reconciliation` - Host-parasite tree reconciliation (TREEMAP, Page 1994)

pub mod tree;
pub mod distance;
pub mod parsimony;
pub mod likelihood;
pub mod bootstrap;
pub mod compatibility;
pub mod comparative;
pub mod consensus;
pub mod invariants;
pub mod io;
pub mod models;
pub mod biogeography;
pub mod reconciliation;

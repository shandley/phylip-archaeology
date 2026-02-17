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
//! - `consensus` - Consensus tree construction (strict, majority-rule, extended)
//! - `io` - File format I/O (PHYLIP format, FASTA)
//! - `models` - Nucleotide substitution models (JC69, K2P, F81, F84)

pub mod tree;
pub mod distance;
pub mod parsimony;
pub mod likelihood;
pub mod bootstrap;
pub mod consensus;
pub mod io;
pub mod models;

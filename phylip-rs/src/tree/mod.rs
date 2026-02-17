//! Phylogenetic tree data structures and Newick format I/O.
//!
//! Core types for representing rooted and unrooted phylogenetic trees,
//! with support for reading and writing Newick format strings.

pub mod newick;
pub mod types;

pub use types::{
    Alignment, AlignmentError, Base, BaseFrequencies, DistanceMatrix, Sequence, Tree, TreeNode,
};

//! Distance matrix methods for phylogenetic inference.
//!
//! Implementations of neighbor-joining (Saitou & Nei, 1987),
//! UPGMA (Sokal & Michener, 1958), and Fitch-Margoliash (1967).

pub mod neighbor_joining;
pub mod upgma;
pub mod fitch_margoliash;
pub mod kitsch;
pub mod ml_distances;

pub use neighbor_joining::neighbor_joining;
pub use upgma::upgma;
pub use fitch_margoliash::fitch_margoliash;
pub use kitsch::kitsch;

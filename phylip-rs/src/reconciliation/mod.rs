//! Host-parasite tree reconciliation: mapping parasite evolutionary history
//! onto host phylogenies.
//!
//! This module implements the TREEMAP algorithm of Page (1994) for
//! reconciling a parasite (or symbiont) tree with a host tree, given a
//! known association between parasite and host leaf taxa.
//!
//! # Background
//!
//! When two lineages share an intimate ecological relationship (e.g.,
//! parasites and their hosts, mutualists and their partners), their
//! phylogenies may exhibit congruent branching patterns produced by
//! *cospeciation* -- simultaneous speciation driven by the association.
//! However, several other processes can produce incongruence:
//!
//! - **Cospeciation** -- host and parasite speciate together (cost 0)
//! - **Duplication** -- the parasite speciates within a single host lineage
//!   without a corresponding host speciation event (cost 1)
//! - **Host-switching** -- a parasite lineage colonises a new, distantly
//!   related host (cost 2)
//! - **Sorting (lineage loss)** -- a parasite lineage fails to persist in
//!   one of the descendant host lineages after a host speciation event
//!   (cost 1)
//!
//! The reconciliation problem is: given trees H and P and a leaf-level
//! mapping phi, find the mapping of internal parasite nodes to host nodes
//! that minimises total event cost.
//!
//! # Algorithm
//!
//! The core algorithm proceeds bottom-up (postorder) on the parasite tree.
//! For each internal parasite node, the Lowest Common Ancestor (LCA) of the
//! images of its children in the host tree determines the optimal placement.
//! The event type (cospeciation vs. duplication) and the number of sorting
//! events are then read off from the relative positions of the mapped nodes
//! in the host tree.
//!
//! # References
//!
//! - Page, R.D.M. (1994). Maps between trees and cladistic analysis of
//!   historical associations among genes, organisms, and areas.
//!   *Systematic Biology*, 43:58-77.
//! - Page, R.D.M. (1994). Parallel phylogenies: reconstructing the history
//!   of host-parasite assemblages. *Cladistics*, 10:155-173.
//! - Charleston, M.A. (1998). Jungles: a new solution to the host/parasite
//!   phylogeny reconciliation problem. *Mathematical Biosciences*, 149:191-223.

pub mod treemap;

pub use treemap::{
    reconcile, lca, Event, EventCosts, ReconciliationResult,
};

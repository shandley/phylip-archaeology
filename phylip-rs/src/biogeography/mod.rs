//! Historical biogeography: reconstructing ancestral geographic distributions.
//!
//! This module implements methods for inferring the geographic history of
//! lineages on a phylogenetic tree, given the present-day distributions of
//! terminal taxa across discrete geographic areas.
//!
//! # Methods
//!
//! - **DIVA (Dispersal-Vicariance Analysis)** — Ronquist's (1997) event-based
//!   method that reconstructs ancestral areas by minimizing dispersal and
//!   extinction events under the assumption that speciation occurs primarily
//!   by vicariance (geographic splitting of ranges).
//!
//! # Representation
//!
//! Geographic areas are numbered 0 through `max_areas - 1`. Distributions
//! (the set of areas occupied by a taxon or ancestor) are represented as
//! bit-vectors using the [`AreaSet`] type, which wraps a `u32` and supports
//! up to 32 areas (practical analyses typically use 2-10).
//!
//! # References
//!
//! - Ronquist, F. (1997). Dispersal-vicariance analysis: a new approach to the
//!   quantification of historical biogeography. *Systematic Biology*, 46:195-203.
//! - Ronquist, F. (1996). DIVA version 1.1. Computer program and manual available
//!   from Uppsala University.

pub mod diva;

pub use diva::{
    diva_optimize, AreaSet, DivaInput, DivaResult, EventCosts, NodeReconstruction,
};

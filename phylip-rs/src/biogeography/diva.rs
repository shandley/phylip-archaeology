//! Dispersal-Vicariance Analysis (DIVA).
//!
//! Implements the DIVA algorithm of Ronquist (1997) for reconstructing
//! ancestral geographic distributions on a phylogenetic tree. The method
//! minimizes the number of dispersal and extinction events required to
//! explain present-day distributions, assuming that speciation proceeds
//! by vicariance (geographic splitting of ancestral ranges).
//!
//! # Algorithm
//!
//! DIVA uses bottom-up dynamic programming on a rooted tree:
//!
//! 1. **Leaves**: The cost of the observed distribution is 0; all other
//!    distributions have infinite cost.
//!
//! 2. **Internal nodes**: For each candidate ancestral distribution *A*,
//!    the algorithm considers all ways to partition *A* between the two
//!    daughter lineages via vicariance (disjoint split of *A*) or
//!    within-area speciation (duplication, both daughters inherit *A*).
//!    Dispersal and extinction events bridge any mismatch between the
//!    inherited distribution and the optimal distribution reconstructed
//!    for each daughter subtree.
//!
//! 3. **Root**: The distribution with minimum total cost is the optimal
//!    ancestral reconstruction. Multiple equally optimal distributions
//!    may exist.
//!
//! # Event costs (Ronquist 1997 defaults)
//!
//! | Event | Cost | Description |
//! |-------|------|-------------|
//! | Vicariance | 0 | Range splits at speciation |
//! | Duplication | 0 | Within-area speciation (no range change) |
//! | Dispersal | 1 | Gain of one area along a branch |
//! | Extinction | 1 | Loss of one area along a branch |
//!
//! # Complexity
//!
//! For *k* areas, there are 2^k - 1 possible non-empty distributions.
//! The algorithm enumerates all bipartitions of each distribution, giving
//! O(3^k) work per node (by the binomial theorem). Practical for k <= 10;
//! the implementation supports up to k = 20 but values above ~12 may be
//! slow.
//!
//! # References
//!
//! - Ronquist, F. (1997). Dispersal-vicariance analysis: a new approach to
//!   the quantification of historical biogeography. *Systematic Biology*,
//!   46:195-203.
//! - Ronquist, F. (1996). DIVA version 1.1. Computer program for
//!   dispersal-vicariance analysis.

use std::collections::HashMap;
use std::fmt;

use crate::tree::types::Tree;

// ---------------------------------------------------------------------------
// AreaSet: bit-vector representation of geographic distributions
// ---------------------------------------------------------------------------

/// A set of geographic areas represented as a bit-vector.
///
/// Each bit position corresponds to a numbered area (0-indexed). Bit *i*
/// being set means the taxon (or ancestor) is present in area *i*. This
/// representation allows efficient set operations using bitwise operators.
///
/// # Examples
///
/// ```
/// use phylip_rs::biogeography::AreaSet;
///
/// let a = AreaSet::singleton(0); // area 0 only
/// let b = AreaSet::singleton(1); // area 1 only
/// let ab = a.union(b);          // areas 0 and 1
/// assert_eq!(ab.count(), 2);
/// assert!(ab.contains(0));
/// assert!(ab.contains(1));
/// assert!(!ab.contains(2));
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct AreaSet(u32);

impl AreaSet {
    /// The empty set (no areas).
    pub const EMPTY: AreaSet = AreaSet(0);

    /// Create an area set from a raw bitmask.
    #[inline]
    pub fn from_bits(bits: u32) -> Self {
        AreaSet(bits)
    }

    /// Return the raw bitmask.
    #[inline]
    pub fn bits(self) -> u32 {
        self.0
    }

    /// Create a set containing a single area.
    ///
    /// # Panics
    ///
    /// Panics if `area >= 32`.
    #[inline]
    pub fn singleton(area: usize) -> Self {
        assert!(area < 32, "area index must be < 32");
        AreaSet(1u32 << area)
    }

    /// Create a set from a slice of area indices.
    pub fn from_areas(areas: &[usize]) -> Self {
        let mut bits = 0u32;
        for &a in areas {
            assert!(a < 32, "area index must be < 32");
            bits |= 1u32 << a;
        }
        AreaSet(bits)
    }

    /// True if this set contains no areas.
    #[inline]
    pub fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// Number of areas in this set.
    #[inline]
    pub fn count(self) -> u32 {
        self.0.count_ones()
    }

    /// True if this set contains the given area.
    #[inline]
    pub fn contains(self, area: usize) -> bool {
        area < 32 && (self.0 & (1u32 << area)) != 0
    }

    /// Set union (areas in either set).
    #[inline]
    pub fn union(self, other: AreaSet) -> AreaSet {
        AreaSet(self.0 | other.0)
    }

    /// Set intersection (areas in both sets).
    #[inline]
    pub fn intersection(self, other: AreaSet) -> AreaSet {
        AreaSet(self.0 & other.0)
    }

    /// Set difference (areas in self but not in other).
    #[inline]
    pub fn difference(self, other: AreaSet) -> AreaSet {
        AreaSet(self.0 & !other.0)
    }

    /// True if self is a subset of other.
    #[inline]
    pub fn is_subset_of(self, other: AreaSet) -> bool {
        self.0 & !other.0 == 0
    }

    /// True if self and other share no areas.
    #[inline]
    pub fn is_disjoint(self, other: AreaSet) -> bool {
        self.0 & other.0 == 0
    }

    /// Iterate over the area indices in this set (ascending order).
    pub fn iter(self) -> AreaSetIter {
        AreaSetIter { remaining: self.0 }
    }

    /// Return the full area set for `k` areas: {0, 1, ..., k-1}.
    ///
    /// # Panics
    ///
    /// Panics if `k > 32`.
    pub fn full(k: usize) -> Self {
        assert!(k <= 32, "max_areas must be <= 32");
        if k == 32 {
            AreaSet(u32::MAX)
        } else {
            AreaSet((1u32 << k) - 1)
        }
    }
}

impl fmt::Debug for AreaSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "AreaSet{{")?;
        let mut first = true;
        for a in self.iter() {
            if !first {
                write!(f, ",")?;
            }
            write!(f, "{}", a)?;
            first = false;
        }
        write!(f, "}}")
    }
}

impl fmt::Display for AreaSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{{")?;
        let mut first = true;
        for a in self.iter() {
            if !first {
                write!(f, ", ")?;
            }
            write!(f, "{}", a)?;
            first = false;
        }
        write!(f, "}}")
    }
}

/// Iterator over the area indices in an [`AreaSet`].
pub struct AreaSetIter {
    remaining: u32,
}

impl Iterator for AreaSetIter {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        if self.remaining == 0 {
            None
        } else {
            let bit = self.remaining.trailing_zeros() as usize;
            self.remaining &= self.remaining - 1; // clear lowest set bit
            Some(bit)
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.remaining.count_ones() as usize;
        (n, Some(n))
    }
}

impl ExactSizeIterator for AreaSetIter {}

// ---------------------------------------------------------------------------
// Subset enumeration helpers
// ---------------------------------------------------------------------------

/// Enumerate all non-empty subsets of a given area set.
///
/// Uses the classic bitmask subset enumeration technique:
/// starting from the full set, repeatedly subtract 1 and mask
/// to get the next smaller subset.
///
/// Returns subsets in descending order by bitmask value.
pub fn nonempty_subsets(set: AreaSet) -> Vec<AreaSet> {
    let mask = set.bits();
    if mask == 0 {
        return Vec::new();
    }
    let mut result = Vec::new();
    let mut sub = mask;
    loop {
        result.push(AreaSet::from_bits(sub));
        if sub == 0 {
            break;
        }
        sub = (sub - 1) & mask;
        if sub == 0 {
            // We have enumerated down to empty; include nothing more.
            break;
        }
    }
    result
}

/// Enumerate all disjoint bipartitions (A, B) of a set S such that
/// A | B = S, A & B = 0, A != 0, B != 0.
///
/// Each bipartition is returned once with A <= B (by bitmask value).
/// This avoids double-counting (A,B) and (B,A).
pub fn disjoint_bipartitions(set: AreaSet) -> Vec<(AreaSet, AreaSet)> {
    let mask = set.bits();
    if mask == 0 || mask.count_ones() < 2 {
        return Vec::new();
    }
    let mut result = Vec::new();
    // Enumerate all non-empty proper subsets
    let mut sub = (mask - 1) & mask;
    while sub != 0 {
        let complement = mask & !sub;
        // Only include each pair once: sub < complement by bitmask value
        if sub < complement {
            result.push((AreaSet::from_bits(sub), AreaSet::from_bits(complement)));
        }
        sub = (sub - 1) & mask;
    }
    result
}

// ---------------------------------------------------------------------------
// Event costs
// ---------------------------------------------------------------------------

/// Cost parameters for biogeographic events.
///
/// Following Ronquist (1997), the default costs are:
/// - Vicariance (range splitting at speciation): 0
/// - Duplication (within-area speciation): 0
/// - Dispersal (gain of one area): 1
/// - Extinction (loss of one area): 1
///
/// These defaults encode the assumption that vicariance is the primary mode
/// of speciation, while dispersal and extinction are secondary processes
/// whose frequency we wish to minimize.
#[derive(Debug, Clone, Copy)]
pub struct EventCosts {
    /// Cost per area gained (dispersal). Default: 1.
    pub dispersal: u32,
    /// Cost per area lost (extinction). Default: 1.
    pub extinction: u32,
}

impl EventCosts {
    /// Create custom event costs.
    pub fn new(dispersal: u32, extinction: u32) -> Self {
        Self {
            dispersal,
            extinction,
        }
    }

    /// Compute the cost of going from an inherited distribution to an
    /// observed (or optimal) distribution along a branch.
    ///
    /// - **Dispersal cost**: number of areas in `observed` but not in
    ///   `inherited` (areas gained), times `self.dispersal`.
    /// - **Extinction cost**: number of areas in `inherited` but not in
    ///   `observed` (areas lost), times `self.extinction`.
    #[inline]
    pub fn branch_cost(self, inherited: AreaSet, observed: AreaSet) -> u32 {
        let gained = observed.difference(inherited).count(); // dispersal
        let lost = inherited.difference(observed).count(); // extinction
        gained * self.dispersal + lost * self.extinction
    }
}

impl Default for EventCosts {
    fn default() -> Self {
        Self {
            dispersal: 1,
            extinction: 1,
        }
    }
}

// ---------------------------------------------------------------------------
// Input / output types
// ---------------------------------------------------------------------------

/// Input for DIVA analysis.
///
/// Combines a rooted phylogenetic tree with the known geographic
/// distributions of terminal taxa and analysis parameters.
#[derive(Debug, Clone)]
pub struct DivaInput {
    /// The rooted phylogenetic tree.
    pub tree: Tree,
    /// Geographic distribution of each terminal taxon, keyed by taxon name.
    pub leaf_areas: HashMap<String, AreaSet>,
    /// Number of geographic areas in the analysis.
    pub max_areas: usize,
    /// Event cost parameters.
    pub costs: EventCosts,
    /// Optional constraint on the maximum number of areas in any ancestral
    /// distribution. If `None`, no constraint is applied (any non-empty
    /// subset of all areas is allowed). Ronquist (1997) recommended
    /// constraining this to reduce the solution space and produce more
    /// biologically realistic reconstructions.
    pub max_range_size: Option<usize>,
}

/// Per-node reconstruction result.
#[derive(Debug, Clone)]
pub struct NodeReconstruction {
    /// Node index in the tree.
    pub node_id: usize,
    /// The optimal ancestral distribution(s) (may be multiple if tied).
    pub optimal_areas: Vec<AreaSet>,
    /// The minimum cost achievable at this node.
    pub min_cost: u32,
}

/// Result of a DIVA optimization.
///
/// Contains the optimal ancestral distribution for each node and the
/// total minimum cost (number of dispersal + extinction events).
#[derive(Debug, Clone)]
pub struct DivaResult {
    /// Per-node reconstructions, indexed by node ID.
    pub reconstructions: Vec<NodeReconstruction>,
    /// Total minimum cost at the root (dispersal + extinction events).
    pub total_cost: u32,
    /// Optimal distribution(s) at the root.
    pub root_areas: Vec<AreaSet>,
    /// Number of areas used in the analysis.
    pub max_areas: usize,
}

impl DivaResult {
    /// Get the reconstruction for a specific node.
    pub fn node(&self, node_id: usize) -> &NodeReconstruction {
        &self.reconstructions[node_id]
    }

    /// Get the names of the optimal areas at a node, given area labels.
    pub fn node_area_names(&self, node_id: usize, area_names: &[&str]) -> Vec<Vec<String>> {
        self.reconstructions[node_id]
            .optimal_areas
            .iter()
            .map(|areas| {
                areas
                    .iter()
                    .map(|i| {
                        if i < area_names.len() {
                            area_names[i].to_string()
                        } else {
                            format!("Area{}", i)
                        }
                    })
                    .collect()
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// DIVA error type
// ---------------------------------------------------------------------------

/// Errors that can occur during DIVA analysis.
#[derive(Debug, Clone)]
pub enum DivaError {
    /// A leaf taxon in the tree has no area assignment.
    MissingLeafArea(String),
    /// The tree is not rooted or has no nodes.
    InvalidTree(String),
    /// A leaf has an empty area set.
    EmptyLeafArea(String),
    /// max_areas is out of the supported range.
    InvalidMaxAreas(usize),
    /// An internal node does not have exactly two children (non-binary tree).
    NonBinaryNode(usize),
}

impl fmt::Display for DivaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DivaError::MissingLeafArea(name) => {
                write!(f, "no area assignment for leaf taxon '{}'", name)
            }
            DivaError::InvalidTree(msg) => write!(f, "invalid tree: {}", msg),
            DivaError::EmptyLeafArea(name) => {
                write!(f, "leaf taxon '{}' has empty area set", name)
            }
            DivaError::InvalidMaxAreas(k) => {
                write!(f, "max_areas = {} is out of supported range 1..=20", k)
            }
            DivaError::NonBinaryNode(id) => {
                write!(f, "node {} is not binary (DIVA requires a fully bifurcating tree)", id)
            }
        }
    }
}

impl std::error::Error for DivaError {}

// ---------------------------------------------------------------------------
// Sentinel cost value
// ---------------------------------------------------------------------------

/// Sentinel for impossible / infinite cost.
const INF_COST: u32 = u32::MAX / 2;

// ---------------------------------------------------------------------------
// Core DIVA algorithm
// ---------------------------------------------------------------------------

/// Run Dispersal-Vicariance Analysis (DIVA) on a rooted phylogenetic tree.
///
/// Given a rooted tree and the geographic distributions of terminal taxa,
/// this function reconstructs ancestral distributions by minimizing the
/// total number of dispersal and extinction events, following Ronquist (1997).
///
/// # Arguments
///
/// * `tree` - A rooted, fully bifurcating phylogenetic tree.
/// * `leaf_areas` - A map from leaf taxon name to its geographic distribution.
/// * `max_areas` - The number of geographic areas (areas are numbered 0..max_areas-1).
///
/// # Returns
///
/// A [`DivaResult`] containing the optimal ancestral distribution for each
/// node and the total minimum cost.
///
/// # Errors
///
/// Returns a [`DivaError`] if the tree is invalid, any leaf lacks an area
/// assignment, or the tree is not fully bifurcating.
///
/// # Example
///
/// ```
/// use phylip_rs::biogeography::{diva_optimize, AreaSet};
/// use phylip_rs::tree::types::Tree;
/// use std::collections::HashMap;
///
/// // Build a simple 3-taxon tree: ((A,B),C)
/// let mut tree = Tree::new();
/// let root = tree.add_node(None, None, None);
/// let ab = tree.add_node(None, None, Some(root));
/// let c = tree.add_node(Some("C".to_string()), Some(1.0), Some(root));
/// let a = tree.add_node(Some("A".to_string()), Some(1.0), Some(ab));
/// let b = tree.add_node(Some("B".to_string()), Some(1.0), Some(ab));
///
/// let mut leaf_areas = HashMap::new();
/// leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
/// leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
/// leaf_areas.insert("C".to_string(), AreaSet::singleton(0));
///
/// let result = diva_optimize(&tree, &leaf_areas, 2).unwrap();
/// assert!(result.total_cost <= 1);
/// ```
///
/// # References
///
/// Ronquist, F. (1997). Dispersal-vicariance analysis: a new approach to the
/// quantification of historical biogeography. *Systematic Biology*, 46:195-203.
pub fn diva_optimize(
    tree: &Tree,
    leaf_areas: &HashMap<String, AreaSet>,
    max_areas: usize,
) -> Result<DivaResult, DivaError> {
    diva_optimize_with_costs(tree, leaf_areas, max_areas, EventCosts::default(), None)
}

/// Run DIVA with custom event costs and optional range-size constraint.
///
/// This is the full-featured entry point. See [`diva_optimize`] for a
/// simplified interface using default costs and no range constraint.
///
/// # Arguments
///
/// * `tree` - A rooted, fully bifurcating phylogenetic tree.
/// * `leaf_areas` - Map from leaf taxon name to geographic distribution.
/// * `max_areas` - Number of geographic areas.
/// * `costs` - Event cost parameters (dispersal and extinction costs).
/// * `max_range_size` - Optional cap on the number of areas in any ancestral
///   distribution. `None` means no constraint.
pub fn diva_optimize_with_costs(
    tree: &Tree,
    leaf_areas: &HashMap<String, AreaSet>,
    max_areas: usize,
    costs: EventCosts,
    max_range_size: Option<usize>,
) -> Result<DivaResult, DivaError> {
    // Validate inputs
    if max_areas == 0 || max_areas > 20 {
        return Err(DivaError::InvalidMaxAreas(max_areas));
    }
    if tree.num_nodes() == 0 {
        return Err(DivaError::InvalidTree("tree has no nodes".to_string()));
    }

    let n_nodes = tree.num_nodes();
    let n_dists = (1usize << max_areas) - 1; // number of non-empty distributions

    // Validate that all leaves have area assignments
    for leaf in tree.leaves() {
        let name = leaf
            .name
            .as_ref()
            .ok_or_else(|| DivaError::InvalidTree("leaf node has no name".to_string()))?;
        let areas = leaf_areas
            .get(name)
            .ok_or_else(|| DivaError::MissingLeafArea(name.clone()))?;
        if areas.is_empty() {
            return Err(DivaError::EmptyLeafArea(name.clone()));
        }
    }

    // Validate that all internal nodes are binary
    for node in &tree.nodes {
        if !node.is_leaf() && node.children.len() != 2 {
            return Err(DivaError::NonBinaryNode(node.id));
        }
    }

    // cost_table[node_id][dist_index] = minimum cost for node to have distribution dist_index+1
    // Distribution index d corresponds to AreaSet::from_bits(d+1), so index 0 = bits 1, etc.
    // We use dist = bits value (1..=n_dists), stored at index dist-1.
    let mut cost_table: Vec<Vec<u32>> = vec![vec![INF_COST; n_dists]; n_nodes];

    // Determine which distributions are allowed (range-size constraint)
    let allowed: Vec<bool> = (1..=n_dists)
        .map(|bits| {
            let a = AreaSet::from_bits(bits as u32);
            // Must be a subset of the full area set
            let within_range = a.is_subset_of(AreaSet::full(max_areas));
            let within_size = match max_range_size {
                Some(max_sz) => (a.count() as usize) <= max_sz,
                None => true,
            };
            within_range && within_size
        })
        .collect();

    // Bottom-up (postorder) pass
    let postorder = tree.postorder();
    for &node_id in &postorder {
        let node = &tree.nodes[node_id];

        if node.is_leaf() {
            // Leaf: cost is 0 for observed distribution, INF for everything else
            let name = node.name.as_ref().unwrap();
            let obs = leaf_areas[name];
            let idx = obs.bits() as usize - 1;
            if idx < n_dists {
                cost_table[node_id][idx] = 0;
            }
        } else {
            // Internal node with two children
            let left_id = node.children[0];
            let right_id = node.children[1];

            // For each possible ancestral distribution A at this node:
            for a_bits in 1..=n_dists {
                if !allowed[a_bits - 1] {
                    continue;
                }
                let a = AreaSet::from_bits(a_bits as u32);

                let mut best_cost = INF_COST;

                // Strategy 1: Vicariance -- split A into two disjoint, non-empty
                // subsets (left_part, right_part) where left_part | right_part = A.
                // Cost = 0 for the split itself, plus the cost of each child
                // reaching its optimal distribution from the inherited part.
                if a.count() >= 2 {
                    let biparts = disjoint_bipartitions(a);
                    for &(part1, part2) in &biparts {
                        // Try both assignments: (left=part1, right=part2) and vice versa
                        for &(left_part, right_part) in &[(part1, part2), (part2, part1)] {
                            let cost = vicariance_cost(
                                left_part,
                                right_part,
                                left_id,
                                right_id,
                                &cost_table,
                                n_dists,
                                &allowed,
                                costs,
                            );
                            if cost < best_cost {
                                best_cost = cost;
                            }
                        }
                    }
                }

                // Strategy 2: Duplication (within-area speciation) -- both children
                // inherit the full distribution A. Cost = 0 for the event itself.
                {
                    let cost = duplication_cost(
                        a, left_id, right_id, &cost_table, n_dists, &allowed, costs,
                    );
                    if cost < best_cost {
                        best_cost = cost;
                    }
                }

                cost_table[node_id][a_bits - 1] = best_cost;
            }
        }
    }

    // Build reconstructions: for each node, find all distributions with minimum cost
    let mut reconstructions: Vec<NodeReconstruction> = Vec::with_capacity(n_nodes);
    for node_id in 0..n_nodes {
        let costs_for_node = &cost_table[node_id];
        let min_cost = *costs_for_node.iter().min().unwrap_or(&INF_COST);
        let optimal: Vec<AreaSet> = if min_cost >= INF_COST {
            Vec::new()
        } else {
            (1..=n_dists)
                .filter(|&d| costs_for_node[d - 1] == min_cost)
                .map(|d| AreaSet::from_bits(d as u32))
                .collect()
        };
        reconstructions.push(NodeReconstruction {
            node_id,
            optimal_areas: optimal,
            min_cost,
        });
    }

    let root_recon = &reconstructions[tree.root];
    let total_cost = root_recon.min_cost;
    let root_areas = root_recon.optimal_areas.clone();

    Ok(DivaResult {
        reconstructions,
        total_cost,
        root_areas,
        max_areas,
    })
}

/// Compute the cost of a vicariance split at an internal node.
///
/// The ancestral distribution A is split into (left_part, right_part).
/// Each child then needs dispersal/extinction events to reach its optimal
/// distribution from the inherited part.
fn vicariance_cost(
    left_part: AreaSet,
    right_part: AreaSet,
    left_id: usize,
    right_id: usize,
    cost_table: &[Vec<u32>],
    n_dists: usize,
    allowed: &[bool],
    costs: EventCosts,
) -> u32 {
    let left_cost = child_transition_cost(left_part, left_id, cost_table, n_dists, allowed, costs);
    if left_cost >= INF_COST {
        return INF_COST;
    }
    let right_cost =
        child_transition_cost(right_part, right_id, cost_table, n_dists, allowed, costs);
    if right_cost >= INF_COST {
        return INF_COST;
    }
    left_cost.saturating_add(right_cost)
}

/// Compute the cost of a duplication (within-area speciation) at an internal node.
///
/// Both children inherit the full ancestral distribution A.
fn duplication_cost(
    ancestor: AreaSet,
    left_id: usize,
    right_id: usize,
    cost_table: &[Vec<u32>],
    n_dists: usize,
    allowed: &[bool],
    costs: EventCosts,
) -> u32 {
    let left_cost =
        child_transition_cost(ancestor, left_id, cost_table, n_dists, allowed, costs);
    if left_cost >= INF_COST {
        return INF_COST;
    }
    let right_cost =
        child_transition_cost(ancestor, right_id, cost_table, n_dists, allowed, costs);
    if right_cost >= INF_COST {
        return INF_COST;
    }
    left_cost.saturating_add(right_cost)
}

/// Compute the minimum cost for a child to transition from an inherited
/// distribution to its optimal distribution.
///
/// The child inherits `inherited` and must reach some distribution `d` for
/// which `cost_table[child_id][d-1]` is defined. The transition cost is:
///   dispersal_cost(areas gained) + extinction_cost(areas lost) + cost_table[child_id][d-1]
///
/// We search over all allowed distributions `d` to find the minimum.
fn child_transition_cost(
    inherited: AreaSet,
    child_id: usize,
    cost_table: &[Vec<u32>],
    n_dists: usize,
    allowed: &[bool],
    costs: EventCosts,
) -> u32 {
    let child_costs = &cost_table[child_id];
    let mut best = INF_COST;

    for d_bits in 1..=n_dists {
        if !allowed[d_bits - 1] {
            continue;
        }
        let child_opt_cost = child_costs[d_bits - 1];
        if child_opt_cost >= INF_COST {
            continue;
        }
        let d = AreaSet::from_bits(d_bits as u32);
        let transition = costs.branch_cost(inherited, d);
        let total = transition.saturating_add(child_opt_cost);
        if total < best {
            best = total;
        }
    }

    best
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // AreaSet unit tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_area_set_singleton() {
        let a = AreaSet::singleton(0);
        assert!(a.contains(0));
        assert!(!a.contains(1));
        assert_eq!(a.count(), 1);
    }

    #[test]
    fn test_area_set_from_areas() {
        let a = AreaSet::from_areas(&[0, 2, 4]);
        assert!(a.contains(0));
        assert!(!a.contains(1));
        assert!(a.contains(2));
        assert!(!a.contains(3));
        assert!(a.contains(4));
        assert_eq!(a.count(), 3);
    }

    #[test]
    fn test_area_set_union() {
        let a = AreaSet::from_areas(&[0, 1]);
        let b = AreaSet::from_areas(&[1, 2]);
        let u = a.union(b);
        assert_eq!(u.count(), 3);
        assert!(u.contains(0));
        assert!(u.contains(1));
        assert!(u.contains(2));
    }

    #[test]
    fn test_area_set_intersection() {
        let a = AreaSet::from_areas(&[0, 1, 2]);
        let b = AreaSet::from_areas(&[1, 2, 3]);
        let i = a.intersection(b);
        assert_eq!(i.count(), 2);
        assert!(i.contains(1));
        assert!(i.contains(2));
        assert!(!i.contains(0));
        assert!(!i.contains(3));
    }

    #[test]
    fn test_area_set_difference() {
        let a = AreaSet::from_areas(&[0, 1, 2]);
        let b = AreaSet::from_areas(&[1, 2, 3]);
        let d = a.difference(b);
        assert_eq!(d.count(), 1);
        assert!(d.contains(0));
        assert!(!d.contains(1));
    }

    #[test]
    fn test_area_set_subset() {
        let a = AreaSet::from_areas(&[1, 2]);
        let b = AreaSet::from_areas(&[0, 1, 2, 3]);
        assert!(a.is_subset_of(b));
        assert!(!b.is_subset_of(a));
        assert!(a.is_subset_of(a));
    }

    #[test]
    fn test_area_set_disjoint() {
        let a = AreaSet::from_areas(&[0, 1]);
        let b = AreaSet::from_areas(&[2, 3]);
        assert!(a.is_disjoint(b));
        let c = AreaSet::from_areas(&[1, 2]);
        assert!(!a.is_disjoint(c));
    }

    #[test]
    fn test_area_set_iter() {
        let a = AreaSet::from_areas(&[0, 3, 5]);
        let indices: Vec<usize> = a.iter().collect();
        assert_eq!(indices, vec![0, 3, 5]);
    }

    #[test]
    fn test_area_set_empty() {
        let a = AreaSet::EMPTY;
        assert!(a.is_empty());
        assert_eq!(a.count(), 0);
        assert_eq!(a.iter().count(), 0);
    }

    #[test]
    fn test_area_set_full() {
        let f = AreaSet::full(4);
        assert_eq!(f.count(), 4);
        assert!(f.contains(0));
        assert!(f.contains(1));
        assert!(f.contains(2));
        assert!(f.contains(3));
        assert!(!f.contains(4));
    }

    #[test]
    fn test_area_set_display() {
        let a = AreaSet::from_areas(&[0, 2]);
        let s = format!("{}", a);
        assert_eq!(s, "{0, 2}");
    }

    // -----------------------------------------------------------------------
    // Subset enumeration tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_nonempty_subsets() {
        // {0, 1} has 3 non-empty subsets: {0,1}, {0}, {1}
        let s = AreaSet::from_areas(&[0, 1]);
        let subs = nonempty_subsets(s);
        assert_eq!(subs.len(), 3);
        // Check they are all non-empty subsets of s
        for sub in &subs {
            assert!(!sub.is_empty());
            assert!(sub.is_subset_of(s));
        }
    }

    #[test]
    fn test_nonempty_subsets_three_areas() {
        // {0, 1, 2} has 2^3 - 1 = 7 non-empty subsets
        let s = AreaSet::from_areas(&[0, 1, 2]);
        let subs = nonempty_subsets(s);
        assert_eq!(subs.len(), 7);
    }

    #[test]
    fn test_nonempty_subsets_empty() {
        let subs = nonempty_subsets(AreaSet::EMPTY);
        assert!(subs.is_empty());
    }

    #[test]
    fn test_disjoint_bipartitions() {
        // {0, 1} has exactly 1 disjoint bipartition: ({0}, {1})
        let s = AreaSet::from_areas(&[0, 1]);
        let parts = disjoint_bipartitions(s);
        assert_eq!(parts.len(), 1);
        let (a, b) = parts[0];
        assert!(a.is_disjoint(b));
        assert_eq!(a.union(b), s);
    }

    #[test]
    fn test_disjoint_bipartitions_three() {
        // {0, 1, 2} has 3 disjoint bipartitions:
        // ({0}, {1,2}), ({1}, {0,2}), ({2}, {0,1})
        let s = AreaSet::from_areas(&[0, 1, 2]);
        let parts = disjoint_bipartitions(s);
        assert_eq!(parts.len(), 3);
        for &(a, b) in &parts {
            assert!(a.is_disjoint(b));
            assert_eq!(a.union(b), s);
            assert!(!a.is_empty());
            assert!(!b.is_empty());
        }
    }

    #[test]
    fn test_disjoint_bipartitions_four() {
        // {0,1,2,3} has 7 bipartitions: C(4,1)/2 wouldn't work;
        // actually 2^(4-1) - 1 = 7 bipartitions
        let s = AreaSet::from_areas(&[0, 1, 2, 3]);
        let parts = disjoint_bipartitions(s);
        assert_eq!(parts.len(), 7);
    }

    #[test]
    fn test_disjoint_bipartitions_singleton() {
        // A singleton set has no bipartitions
        let s = AreaSet::singleton(0);
        assert!(disjoint_bipartitions(s).is_empty());
    }

    // -----------------------------------------------------------------------
    // Event cost tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_branch_cost_no_change() {
        let costs = EventCosts::default();
        let a = AreaSet::from_areas(&[0, 1]);
        assert_eq!(costs.branch_cost(a, a), 0);
    }

    #[test]
    fn test_branch_cost_dispersal() {
        let costs = EventCosts::default();
        let inherited = AreaSet::singleton(0);
        let observed = AreaSet::from_areas(&[0, 1, 2]);
        // Gained areas 1 and 2: cost = 2 * dispersal = 2
        assert_eq!(costs.branch_cost(inherited, observed), 2);
    }

    #[test]
    fn test_branch_cost_extinction() {
        let costs = EventCosts::default();
        let inherited = AreaSet::from_areas(&[0, 1, 2]);
        let observed = AreaSet::singleton(0);
        // Lost areas 1 and 2: cost = 2 * extinction = 2
        assert_eq!(costs.branch_cost(inherited, observed), 2);
    }

    #[test]
    fn test_branch_cost_mixed() {
        let costs = EventCosts::default();
        let inherited = AreaSet::from_areas(&[0, 1]);
        let observed = AreaSet::from_areas(&[1, 2]);
        // Gained area 2 (dispersal = 1), lost area 0 (extinction = 1): total = 2
        assert_eq!(costs.branch_cost(inherited, observed), 2);
    }

    #[test]
    fn test_branch_cost_custom() {
        let costs = EventCosts::new(2, 3);
        let inherited = AreaSet::from_areas(&[0, 1]);
        let observed = AreaSet::from_areas(&[1, 2, 3]);
        // Gained {2,3} = 2 areas * cost 2 = 4, lost {0} = 1 area * cost 3 = 3: total = 7
        assert_eq!(costs.branch_cost(inherited, observed), 7);
    }

    // -----------------------------------------------------------------------
    // Helper: build a simple tree
    // -----------------------------------------------------------------------

    /// Build a tree: ((A, B), C)
    ///
    /// ```text
    ///       root (0)
    ///      /        \
    ///   ab (1)     C (2)
    ///   /    \
    ///  A (3)  B (4)
    /// ```
    fn make_tree_abc() -> Tree {
        let mut tree = Tree::new();
        let root = tree.add_node(None, None, None);
        let ab = tree.add_node(None, None, Some(root));
        let _c = tree.add_node(Some("C".to_string()), Some(1.0), Some(root));
        let _a = tree.add_node(Some("A".to_string()), Some(1.0), Some(ab));
        let _b = tree.add_node(Some("B".to_string()), Some(1.0), Some(ab));
        tree
    }

    /// Build a tree: ((A, B), (C, D))
    ///
    /// ```text
    ///        root (0)
    ///       /        \
    ///    ab (1)     cd (2)
    ///    /   \      /   \
    ///  A (3) B (4) C (5) D (6)
    /// ```
    fn make_tree_abcd() -> Tree {
        let mut tree = Tree::new();
        let root = tree.add_node(None, None, None);
        let ab = tree.add_node(None, None, Some(root));
        let cd = tree.add_node(None, None, Some(root));
        let _a = tree.add_node(Some("A".to_string()), Some(1.0), Some(ab));
        let _b = tree.add_node(Some("B".to_string()), Some(1.0), Some(ab));
        let _c = tree.add_node(Some("C".to_string()), Some(1.0), Some(cd));
        let _d = tree.add_node(Some("D".to_string()), Some(1.0), Some(cd));
        tree
    }

    // -----------------------------------------------------------------------
    // DIVA optimization tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_pure_vicariance_two_areas() {
        // ((A, B), C) where A={0}, B={1}, C={0}
        // Optimal: ancestor of (A,B) = {0,1} via vicariance, root = {0,1} via vicariance
        // or root = {0} with extinction of area 1 somewhere
        // Vicariance: root={0,1}, split to ab={0,1} and C={0}.
        //   C inherits {1} from vicariance split -> needs dispersal/extinction to become {0}.
        //   But better: root = {0,1}, vicariance split: left={0,1}, right depends.
        //
        // Actually with ((A,B),C), A={0}, B={1}, C={0}:
        //   Node ab: vicariance {0,1} -> A gets {0}, B gets {1}: cost = 0+0 = 0
        //   Root: tries ancestor = {0,1}:
        //     vicariance split {0},{1}: ab inherits {0} needs to reach {0,1} = 1 dispersal + 0 = 1
        //     vicariance split {1},{0}: ab inherits {1} needs to reach {0,1} = 1 dispersal, C inherits {0} cost 0 = 1
        //     duplication: both get {0,1}: ab cost 0, C needs to go from {0,1} to {0} = 1 extinction = 1
        //   Root with ancestor = {0}: duplication: ab gets {0}, needs {0,1} = 1 disp; C gets {0} = 0. Total = 1
        //   So minimum cost = 1 (one dispersal or one extinction event)
        let tree = make_tree_abc();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(0));

        let result = diva_optimize(&tree, &leaf_areas, 2).unwrap();
        assert_eq!(result.total_cost, 1);
    }

    #[test]
    fn test_perfect_vicariance_zero_cost() {
        // ((A, B), (C, D)) where A={0}, B={1}, C={0}, D={1}
        // Perfect vicariance at every node:
        //   ab ancestor = {0,1}, split: A gets {0}, B gets {1}: cost 0
        //   cd ancestor = {0,1}, split: C gets {0}, D gets {1}: cost 0
        //   root = {0,1}, duplication: both ab and cd inherit {0,1}: cost 0+0 = 0
        //
        // Total cost = 0
        let tree = make_tree_abcd();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("D".to_string(), AreaSet::singleton(1));

        let result = diva_optimize(&tree, &leaf_areas, 2).unwrap();
        assert_eq!(result.total_cost, 0);

        // The root optimal areas should include {0,1}
        assert!(result.root_areas.contains(&AreaSet::from_areas(&[0, 1])));
    }

    #[test]
    fn test_perfect_vicariance_four_areas() {
        // ((A, B), (C, D)) where each taxon in a unique area
        // A={0}, B={1}, C={2}, D={3}
        //   ab ancestor = {0,1}, vicariance: cost 0
        //   cd ancestor = {2,3}, vicariance: cost 0
        //   root = {0,1,2,3}, vicariance {0,1}|{2,3}: cost 0
        // Total = 0
        let tree = make_tree_abcd();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(2));
        leaf_areas.insert("D".to_string(), AreaSet::singleton(3));

        let result = diva_optimize(&tree, &leaf_areas, 4).unwrap();
        assert_eq!(result.total_cost, 0);
        assert!(result
            .root_areas
            .contains(&AreaSet::from_areas(&[0, 1, 2, 3])));
    }

    #[test]
    fn test_all_same_area_zero_cost() {
        // ((A, B), C) where all taxa in area 0
        // Duplication everywhere: root={0}, ab={0}: cost 0
        let tree = make_tree_abc();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(0));

        let result = diva_optimize(&tree, &leaf_areas, 1).unwrap();
        assert_eq!(result.total_cost, 0);
        assert_eq!(result.root_areas, vec![AreaSet::singleton(0)]);
    }

    #[test]
    fn test_dispersal_needed() {
        // ((A, B), C) with A={0}, B={0}, C={1}
        // If ancestor of all = {0}: C needs dispersal to {1}, cost = 1+1 = 2 (lose 0, gain 1)
        //   Actually: branch_cost({0}, {1}) = 1 dispersal + 1 extinction = 2
        // If ancestor = {0,1}: duplication to ab={0,1}, C={0,1}
        //   ab: transition cost to {0} = 1 extinction
        //   C: transition cost to {1} = 1 extinction
        //   total = 2
        // If ancestor = {0,1}: vicariance {0}|{1}: ab gets {0} cost 0, C gets {1} cost 0 = 0!
        //   Wait, ab has optimal cost for dist {0} which is 0 (since both A and B are {0}).
        //   So total cost at root = 0.
        // Total = 0. This is a perfect vicariance scenario!
        let tree = make_tree_abc();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(1));

        let result = diva_optimize(&tree, &leaf_areas, 2).unwrap();
        assert_eq!(result.total_cost, 0);
    }

    #[test]
    fn test_three_areas_requires_dispersal() {
        // ((A, B), C) where A={0}, B={1}, C={2}
        // All different areas. At ab: {0,1} vicariance, cost 0.
        // Root: {0,1,2} vicariance: {0,1}|{2}: ab cost 0, C cost 0 => total 0!
        // Actually this is also perfect vicariance.
        let tree = make_tree_abc();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(2));

        let result = diva_optimize(&tree, &leaf_areas, 3).unwrap();
        assert_eq!(result.total_cost, 0);
    }

    #[test]
    fn test_widespread_taxa() {
        // ((A, B), C) where A={0,1}, B={0,1}, C={0}
        // Duplication at ab: ancestor {0,1}, both children inherit {0,1}: cost 0
        // Root: {0,1} duplication: ab gets {0,1} cost 0, C transitions from {0,1} to {0} = 1 extinction
        //   or: root = {0} duplication: ab transitions {0} -> {0,1} = 1 dispersal, C = 0. Total = 1
        //   or: root = {0,1} vicariance {0}|{1}: ab gets {0}, needs {0,1} = 1 dispersal.
        //     C gets {1}, needs {0} = 1 dispersal + 1 extinction = 2. Total = 3.
        //   Best: root = {0,1}, duplication, total = 1 or root = {0}, duplication, total = 1.
        let tree = make_tree_abc();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::from_areas(&[0, 1]));
        leaf_areas.insert("B".to_string(), AreaSet::from_areas(&[0, 1]));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(0));

        let result = diva_optimize(&tree, &leaf_areas, 2).unwrap();
        assert_eq!(result.total_cost, 1);
    }

    #[test]
    fn test_single_area_parsimony() {
        // ((A, B), (C, D)) where A={0}, B={0}, C={1}, D={1}
        // ab: {0} duplication, cost 0
        // cd: {1} duplication, cost 0
        // root: {0,1} vicariance {0}|{1}: ab gets {0} cost 0, cd gets {1} cost 0
        // Total = 0
        let tree = make_tree_abcd();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("D".to_string(), AreaSet::singleton(1));

        let result = diva_optimize(&tree, &leaf_areas, 2).unwrap();
        assert_eq!(result.total_cost, 0);

        // Check internal node reconstructions
        // Node 1 (ab): should have {0} with cost 0
        let ab_recon = result.node(1);
        assert_eq!(ab_recon.min_cost, 0);
        assert!(ab_recon.optimal_areas.contains(&AreaSet::singleton(0)));

        // Node 2 (cd): should have {1} with cost 0
        let cd_recon = result.node(2);
        assert_eq!(cd_recon.min_cost, 0);
        assert!(cd_recon.optimal_areas.contains(&AreaSet::singleton(1)));
    }

    #[test]
    fn test_max_range_constraint() {
        // ((A, B), (C, D)) with A={0}, B={1}, C={2}, D={3}
        // Without constraint: root = {0,1,2,3}, cost = 0
        // With max_range_size = 2: root cannot be {0,1,2,3}, must be at most 2 areas.
        //   Root = {0,1}: vicariance {0}|{1}: ab gets {0} cost 0, cd gets {1}.
        //     cd cost for dist {1}: must split {1}. But {1} is singleton, can only duplicate.
        //     C needs {2}: branch_cost({1},{2}) = 2. D needs {3}: branch_cost({1},{3}) = 2. Total = 4.
        //   Might be cheaper with other root. In any case, cost > 0.
        let tree = make_tree_abcd();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(2));
        leaf_areas.insert("D".to_string(), AreaSet::singleton(3));

        let result_unconstrained = diva_optimize(&tree, &leaf_areas, 4).unwrap();
        assert_eq!(result_unconstrained.total_cost, 0);

        let result_constrained =
            diva_optimize_with_costs(&tree, &leaf_areas, 4, EventCosts::default(), Some(2))
                .unwrap();
        // With max range size 2, cost must be > 0
        assert!(
            result_constrained.total_cost > 0,
            "expected cost > 0 with range constraint, got {}",
            result_constrained.total_cost
        );
        // All optimal distributions should have at most 2 areas
        for recon in &result_constrained.reconstructions {
            for areas in &recon.optimal_areas {
                assert!(
                    areas.count() <= 2,
                    "distribution {} exceeds max range size 2",
                    areas
                );
            }
        }
    }

    #[test]
    fn test_increasing_areas_expands_solution_space() {
        // Same tree and leaf distributions, but increasing max_areas
        // should never increase the cost (more areas = more solution space)
        let tree = make_tree_abc();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(0));

        let result2 = diva_optimize(&tree, &leaf_areas, 2).unwrap();
        let result3 = diva_optimize(&tree, &leaf_areas, 3).unwrap();
        let result4 = diva_optimize(&tree, &leaf_areas, 4).unwrap();

        assert!(
            result3.total_cost <= result2.total_cost,
            "more areas should not increase cost"
        );
        assert!(
            result4.total_cost <= result3.total_cost,
            "more areas should not increase cost"
        );
    }

    #[test]
    fn test_custom_event_costs() {
        // With high dispersal cost and low extinction cost, the algorithm
        // should prefer explanations with more extinctions over dispersals
        let tree = make_tree_abc();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(0));

        let default_result = diva_optimize(&tree, &leaf_areas, 2).unwrap();

        let high_disp = EventCosts::new(10, 1);
        let high_disp_result =
            diva_optimize_with_costs(&tree, &leaf_areas, 2, high_disp, None).unwrap();

        // With higher dispersal cost, total cost should be >= default
        assert!(high_disp_result.total_cost >= default_result.total_cost);
    }

    #[test]
    fn test_symmetric_tree_symmetric_areas() {
        // ((A, B), (C, D)) with A={0}, B={1}, C={1}, D={0}
        // This is the same as the 4-taxon perfect vicariance but areas swapped for C,D.
        // ab: {0,1} vicariance, cost 0
        // cd: {0,1} vicariance ({1}|{0}), cost 0
        // root: {0,1} duplication, cost 0
        let tree = make_tree_abcd();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("D".to_string(), AreaSet::singleton(0));

        let result = diva_optimize(&tree, &leaf_areas, 2).unwrap();
        assert_eq!(result.total_cost, 0);
    }

    #[test]
    fn test_cost_with_required_dispersal() {
        // ((A, B), (C, D)) where A={0}, B={0}, C={0}, D={1}
        // ab: {0} duplication, cost 0
        // cd: {0,1} vicariance {0}|{1}, cost 0
        // root: {0} duplication: ab=0, cd needs {0,1}, branch_cost({0},{0,1}) = 1 dispersal. Total=1
        // root: {0,1} vicariance {0}|{1}: ab gets {0}, cost 0; cd gets {1}, cd optimal for {1}?
        //   cd dist {1}: duplication: C inherits {1}, needs {0}: cost 2. D inherits {1}, cost 0. Total=2.
        //   worse.
        // root: {0,1} duplication: ab gets {0,1}, needs {0}: cost 1 ext. cd gets {0,1}, cost 0. Total=1.
        // root: {0} duplication: ab=0, cd transition {0}->{0,1}=1 disp. Total=1.
        // Best = 1
        let tree = make_tree_abcd();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("D".to_string(), AreaSet::singleton(1));

        let result = diva_optimize(&tree, &leaf_areas, 2).unwrap();
        assert_eq!(result.total_cost, 1);
    }

    #[test]
    fn test_leaf_reconstructions() {
        // Leaf nodes should have their observed distribution as optimal
        let tree = make_tree_abc();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("C".to_string(), AreaSet::from_areas(&[0, 1]));

        let result = diva_optimize(&tree, &leaf_areas, 2).unwrap();

        // Node 3 = A, should have {0}
        assert_eq!(result.node(3).optimal_areas, vec![AreaSet::singleton(0)]);
        assert_eq!(result.node(3).min_cost, 0);

        // Node 4 = B, should have {1}
        assert_eq!(result.node(4).optimal_areas, vec![AreaSet::singleton(1)]);
        assert_eq!(result.node(4).min_cost, 0);

        // Node 2 = C, should have {0,1}
        assert_eq!(
            result.node(2).optimal_areas,
            vec![AreaSet::from_areas(&[0, 1])]
        );
        assert_eq!(result.node(2).min_cost, 0);
    }

    #[test]
    fn test_node_area_names() {
        let tree = make_tree_abc();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(0));

        let result = diva_optimize(&tree, &leaf_areas, 2).unwrap();
        let area_names = ["North", "South"];
        let root_names = result.node_area_names(0, &area_names);
        // Root should be {0} = "North"
        assert!(root_names.iter().any(|names| names == &["North"]));
    }

    // -----------------------------------------------------------------------
    // Error handling tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_error_missing_leaf_area() {
        let tree = make_tree_abc();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        // Missing C

        let result = diva_optimize(&tree, &leaf_areas, 2);
        assert!(result.is_err());
        match result.unwrap_err() {
            DivaError::MissingLeafArea(name) => assert_eq!(name, "C"),
            e => panic!("expected MissingLeafArea, got {:?}", e),
        }
    }

    #[test]
    fn test_error_empty_leaf_area() {
        let tree = make_tree_abc();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("C".to_string(), AreaSet::EMPTY);

        let result = diva_optimize(&tree, &leaf_areas, 2);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), DivaError::EmptyLeafArea(_)));
    }

    #[test]
    fn test_error_invalid_max_areas() {
        let tree = make_tree_abc();
        let leaf_areas = HashMap::new();
        let result = diva_optimize(&tree, &leaf_areas, 0);
        assert!(matches!(result.unwrap_err(), DivaError::InvalidMaxAreas(0)));

        let result = diva_optimize(&tree, &leaf_areas, 25);
        assert!(matches!(
            result.unwrap_err(),
            DivaError::InvalidMaxAreas(25)
        ));
    }

    #[test]
    fn test_error_non_binary_tree() {
        // Build a tree with a trifurcation
        let mut tree = Tree::new();
        let root = tree.add_node(None, None, None);
        let _a = tree.add_node(Some("A".to_string()), Some(1.0), Some(root));
        let _b = tree.add_node(Some("B".to_string()), Some(1.0), Some(root));
        let _c = tree.add_node(Some("C".to_string()), Some(1.0), Some(root));

        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(0));

        let result = diva_optimize(&tree, &leaf_areas, 2);
        assert!(matches!(result.unwrap_err(), DivaError::NonBinaryNode(0)));
    }

    // -----------------------------------------------------------------------
    // Biological scenario tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_gondwanan_vicariance_scenario() {
        // Model: a Gondwanan group split by continental breakup
        // ((Africa_sp, S.America_sp), (Australia_sp, Antarctica_sp))
        // Areas: 0=Africa, 1=S.America, 2=Australia, 3=Antarctica
        // Perfect vicariance: cost should be 0
        let tree = make_tree_abcd();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0)); // Africa
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1)); // S.America
        leaf_areas.insert("C".to_string(), AreaSet::singleton(2)); // Australia
        leaf_areas.insert("D".to_string(), AreaSet::singleton(3)); // Antarctica

        let result = diva_optimize(&tree, &leaf_areas, 4).unwrap();
        assert_eq!(result.total_cost, 0);

        // Root should span all 4 areas
        assert!(result
            .root_areas
            .contains(&AreaSet::from_areas(&[0, 1, 2, 3])));

        // ab ancestor: Africa + S.America
        let ab = result.node(1);
        assert!(ab
            .optimal_areas
            .contains(&AreaSet::from_areas(&[0, 1])));

        // cd ancestor: Australia + Antarctica
        let cd = result.node(2);
        assert!(cd
            .optimal_areas
            .contains(&AreaSet::from_areas(&[2, 3])));
    }

    #[test]
    fn test_island_colonization_scenario() {
        // Model: mainland ancestor colonizes islands via dispersal
        // ((island1_sp, island2_sp), mainland_sp)
        // Areas: 0=mainland, 1=island1, 2=island2
        //
        // If ancestor was mainland {0}:
        //   island clade: ancestor must have dispersed to get to {1,2}
        //   island clade ancestor = {0}: duplication, then children each need
        //   to transition. But actual: vicariance of {1,2}: cost 0 at that node,
        //   but ancestor of island clade needs to transition from inherited to {1,2}.
        //
        // Best: root = {0,1,2}, vicariance {0}|{1,2}: mainland=0 (cost 0),
        //   island clade gets {1,2}, vicariance {1}|{2}: cost 0. Total = 0.
        //
        // But that says ancestor was in all 3 areas, which is geologically unrealistic.
        // With max_range_size = 1, cost will be > 0 (dispersal needed).
        let tree = make_tree_abc();
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(1)); // island1
        leaf_areas.insert("B".to_string(), AreaSet::singleton(2)); // island2
        leaf_areas.insert("C".to_string(), AreaSet::singleton(0)); // mainland

        let result_unconstrained = diva_optimize(&tree, &leaf_areas, 3).unwrap();
        assert_eq!(result_unconstrained.total_cost, 0);

        let result_constrained = diva_optimize_with_costs(
            &tree,
            &leaf_areas,
            3,
            EventCosts::default(),
            Some(1),
        )
        .unwrap();
        // With max range size 1, must have dispersal events
        assert!(result_constrained.total_cost > 0);
    }

    #[test]
    fn test_two_taxon_tree() {
        // Simplest possible tree: (A, B)
        let mut tree = Tree::new();
        let root = tree.add_node(None, None, None);
        let _a = tree.add_node(Some("A".to_string()), Some(1.0), Some(root));
        let _b = tree.add_node(Some("B".to_string()), Some(1.0), Some(root));

        // A={0}, B={1}: vicariance at root, cost 0
        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));

        let result = diva_optimize(&tree, &leaf_areas, 2).unwrap();
        assert_eq!(result.total_cost, 0);
        assert!(result.root_areas.contains(&AreaSet::from_areas(&[0, 1])));
    }

    #[test]
    fn test_two_taxon_same_area() {
        // (A, B) both in area 0: duplication, cost 0
        let mut tree = Tree::new();
        let root = tree.add_node(None, None, None);
        let _a = tree.add_node(Some("A".to_string()), Some(1.0), Some(root));
        let _b = tree.add_node(Some("B".to_string()), Some(1.0), Some(root));

        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(0));

        let result = diva_optimize(&tree, &leaf_areas, 2).unwrap();
        assert_eq!(result.total_cost, 0);
        assert!(result.root_areas.contains(&AreaSet::singleton(0)));
    }

    #[test]
    fn test_five_taxon_asymmetric() {
        // (((A, B), C), (D, E))
        // Areas: A={0}, B={1}, C={0}, D={2}, E={2}
        //
        // Build tree:
        //         root (0)
        //        /        \
        //     abc (1)    de (2)
        //     /    \     /   \
        //   ab (3)  C(4) D(5) E(6)
        //   / \
        //  A(7) B(8)
        let mut tree = Tree::new();
        let root = tree.add_node(None, None, None);
        let abc = tree.add_node(None, None, Some(root));
        let de = tree.add_node(None, None, Some(root));
        let ab = tree.add_node(None, None, Some(abc));
        let _c = tree.add_node(Some("C".to_string()), Some(1.0), Some(abc));
        let _d = tree.add_node(Some("D".to_string()), Some(1.0), Some(de));
        let _e = tree.add_node(Some("E".to_string()), Some(1.0), Some(de));
        let _a = tree.add_node(Some("A".to_string()), Some(1.0), Some(ab));
        let _b = tree.add_node(Some("B".to_string()), Some(1.0), Some(ab));

        let mut leaf_areas = HashMap::new();
        leaf_areas.insert("A".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("B".to_string(), AreaSet::singleton(1));
        leaf_areas.insert("C".to_string(), AreaSet::singleton(0));
        leaf_areas.insert("D".to_string(), AreaSet::singleton(2));
        leaf_areas.insert("E".to_string(), AreaSet::singleton(2));

        let result = diva_optimize(&tree, &leaf_areas, 3).unwrap();

        // de ancestor should be {2} (duplication), cost 0
        let de_recon = result.node(de);
        assert_eq!(de_recon.min_cost, 0);
        assert!(de_recon.optimal_areas.contains(&AreaSet::singleton(2)));

        // ab ancestor should be {0,1} (vicariance), cost 0
        let ab_recon = result.node(ab);
        assert_eq!(ab_recon.min_cost, 0);
        assert!(ab_recon
            .optimal_areas
            .contains(&AreaSet::from_areas(&[0, 1])));

        // Total cost should be small (abc is {0,1} or {0}, root adds vicariance)
        // abc: dist {0,1}: vicariance {0}|{1}: ab gets {0}, cost=1 (needs {0,1}); C gets {1}, cost=2.
        //   or duplication {0,1}: ab gets {0,1} cost 0; C transitions {0,1}->{0} = 1 ext. Total=1.
        //   or abc = {0}: duplication: ab transitions {0}->{0,1} = 1 disp. C=0. Total=1.
        // root: {0,2} vicariance: abc gets {0}, total=1; de gets {2}, cost=0. Total=1.
        //   or root {0,1,2}: vicariance {0,1}|{2}: abc gets {0,1}, cost=1; de gets {2}, cost=0. Total=1.
        assert!(result.total_cost <= 1);
    }
}

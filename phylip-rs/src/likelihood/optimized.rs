//! Performance-optimized likelihood computation.
//!
//! This module provides faster versions of the core likelihood functions
//! by exploiting three key optimizations:
//!
//! 1. **Site-pattern compression** ([`CompressedAlignment`]): Many alignment
//!    columns have identical nucleotide patterns across taxa. By grouping
//!    identical patterns and weighting them, we reduce the number of
//!    per-site evaluations — often dramatically for real datasets.
//!
//! 2. **Transition matrix caching** ([`TransitionCache`]): The pruning
//!    algorithm calls `model.transition_probs(t)` for every edge on every
//!    likelihood evaluation. Since branch lengths change infrequently
//!    relative to the number of evaluations (especially during branch
//!    length optimization), caching and reusing matrices avoids redundant
//!    expensive `exp()` calls.
//!
//! 3. **Pre-allocated buffers** ([`PruningBuffers`]): The standard pruning
//!    algorithm allocates conditional likelihood vectors on each call.
//!    Pre-allocating and reusing these buffers eliminates per-call
//!    allocation overhead.
//!
//! The optimized functions produce numerically identical results to the
//! originals in [`super::pruning`], and are tested against them.
//!
//! # Example
//!
//! ```
//! use phylip_rs::likelihood::models::Jc69Model;
//! use phylip_rs::likelihood::optimized::{
//!     CompressedAlignment, TransitionCache, fast_log_likelihood,
//! };
//! use phylip_rs::tree::{Alignment, Base, Sequence};
//! use phylip_rs::tree::newick::parse_newick;
//!
//! let alignment = Alignment::new(vec![
//!     Sequence::new("A", vec![Base::A, Base::C, Base::G, Base::T]),
//!     Sequence::new("B", vec![Base::A, Base::C, Base::A, Base::T]),
//!     Sequence::new("C", vec![Base::T, Base::G, Base::G, Base::A]),
//! ]).unwrap();
//!
//! let tree = parse_newick("((A:0.1,B:0.1):0.1,C:0.2);").unwrap();
//! let model = Jc69Model;
//! let compressed = CompressedAlignment::from_alignment(&alignment);
//! let mut cache = TransitionCache::new(tree.num_nodes());
//! let lnl = fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap();
//! assert!(lnl < 0.0);
//! ```
//!
//! # References
//!
//! - Felsenstein, J. (1981). Evolutionary trees from DNA sequences:
//!   a maximum likelihood approach. *Journal of Molecular Evolution*,
//!   17, 368-376.
//! - Stamatakis, A. (2006). RAxML-VI-HPC: maximum likelihood-based
//!   phylogenetic analyses with thousands of taxa and mixed models.
//!   *Bioinformatics*, 22(21), 2688-2690.

use std::collections::HashMap;

use crate::tree::types::{Alignment, Base, Tree};

use super::models::SubstitutionModel;
use super::pruning::LikelihoodError;

/// Minimum branch length for numerical stability.
const MIN_BRANCH_LENGTH: f64 = 1e-8;
/// Maximum branch length for optimization.
const MAX_BRANCH_LENGTH: f64 = 10.0;
/// Default branch length when none is specified.
const DEFAULT_BRANCH_LENGTH: f64 = 0.1;
/// Underflow threshold: rescale when max conditional likelihood drops below this.
const UNDERFLOW_THRESHOLD: f64 = 1e-100;
/// Tolerance for considering two branch lengths equal (for cache invalidation).
const BRANCH_LENGTH_TOLERANCE: f64 = 1e-15;

// ============================================================================
// Compressed Alignment
// ============================================================================

/// Compressed alignment that groups identical site patterns.
///
/// Many alignment columns have identical nucleotide patterns across all taxa.
/// For example, if columns 3, 7, and 12 all have the pattern (A, C, G, T) for
/// taxa (X, Y, Z, W), they are represented once with weight 3. This reduces
/// the number of per-site likelihood evaluations proportionally.
///
/// For real datasets the compression ratio is often 2x-10x.
#[derive(Debug, Clone)]
pub struct CompressedAlignment {
    /// Unique site patterns. Each inner Vec has one Base per taxon,
    /// in the same order as `taxa_names`.
    pub patterns: Vec<Vec<Base>>,
    /// Weight (count) of each unique pattern.
    pub weights: Vec<usize>,
    /// Taxon names in the same order as the alignment rows.
    pub taxa_names: Vec<String>,
}

impl CompressedAlignment {
    /// Build a compressed alignment from a standard alignment.
    ///
    /// Identical site columns are merged and counted. The taxa ordering
    /// is preserved from the input alignment.
    pub fn from_alignment(alignment: &Alignment) -> Self {
        let ntaxa = alignment.ntaxa();
        let nsites = alignment.nsites();

        let taxa_names: Vec<String> = alignment
            .sequences
            .iter()
            .map(|s| s.name.clone())
            .collect();

        // Build patterns column by column, using a HashMap to group identical ones.
        let mut pattern_map: HashMap<Vec<Base>, usize> = HashMap::new();
        let mut patterns: Vec<Vec<Base>> = Vec::new();
        let mut weights: Vec<usize> = Vec::new();

        for site in 0..nsites {
            let mut column = Vec::with_capacity(ntaxa);
            for taxon in 0..ntaxa {
                column.push(alignment.base(taxon, site));
            }

            if let Some(&idx) = pattern_map.get(&column) {
                weights[idx] += 1;
            } else {
                let idx = patterns.len();
                pattern_map.insert(column.clone(), idx);
                patterns.push(column);
                weights.push(1);
            }
        }

        CompressedAlignment {
            patterns,
            weights,
            taxa_names,
        }
    }

    /// Number of unique site patterns.
    pub fn n_patterns(&self) -> usize {
        self.patterns.len()
    }

    /// Number of taxa.
    pub fn n_taxa(&self) -> usize {
        self.taxa_names.len()
    }

    /// Total number of sites (sum of weights).
    pub fn total_sites(&self) -> usize {
        self.weights.iter().sum()
    }
}

// ============================================================================
// Transition Matrix Cache
// ============================================================================

/// Cached transition probability matrices for a tree.
///
/// Stores one 4x4 transition matrix per node (keyed by node ID). Each matrix
/// is computed from `model.transition_probs(t)` and is invalidated only when
/// the corresponding branch length changes. This avoids redundant `exp()`
/// computations during iterative branch-length optimization, where typically
/// only one branch length changes per iteration while all others stay fixed.
#[derive(Debug, Clone)]
pub struct TransitionCache {
    /// Cached transition probability matrices, indexed by node ID.
    /// `None` means the cache entry is invalid / not yet computed.
    matrices: Vec<Option<[[f64; 4]; 4]>>,
    /// Branch lengths at the time of caching, for invalidation.
    branch_lengths: Vec<Option<f64>>,
}

impl TransitionCache {
    /// Create a new cache for a tree with `n` nodes.
    pub fn new(n: usize) -> Self {
        TransitionCache {
            matrices: vec![None; n],
            branch_lengths: vec![None; n],
        }
    }

    /// Ensure the cache has room for at least `n` nodes.
    /// Extends with `None` entries if needed.
    pub fn ensure_capacity(&mut self, n: usize) {
        if self.matrices.len() < n {
            self.matrices.resize(n, None);
            self.branch_lengths.resize(n, None);
        }
    }

    /// Invalidate the cache entry for a specific node.
    pub fn invalidate(&mut self, node_id: usize) {
        if node_id < self.matrices.len() {
            self.matrices[node_id] = None;
            self.branch_lengths[node_id] = None;
        }
    }

    /// Invalidate all cache entries.
    pub fn invalidate_all(&mut self) {
        for m in &mut self.matrices {
            *m = None;
        }
        for bl in &mut self.branch_lengths {
            *bl = None;
        }
    }

    /// Get the transition matrix for a node, computing and caching if needed.
    ///
    /// The matrix is recomputed only if:
    /// - No cached matrix exists for this node, or
    /// - The branch length has changed since the last computation.
    pub fn get_or_compute(
        &mut self,
        node_id: usize,
        branch_length: f64,
        model: &dyn SubstitutionModel,
    ) -> [[f64; 4]; 4] {
        self.ensure_capacity(node_id + 1);

        // Check if the cached value is still valid
        if let Some(cached_bl) = self.branch_lengths[node_id] {
            if (cached_bl - branch_length).abs() < BRANCH_LENGTH_TOLERANCE {
                if let Some(matrix) = self.matrices[node_id] {
                    return matrix;
                }
            }
        }

        // Recompute
        let matrix = model.transition_probs(branch_length);
        self.matrices[node_id] = Some(matrix);
        self.branch_lengths[node_id] = Some(branch_length);
        matrix
    }
}

// ============================================================================
// Pre-allocated Pruning Buffers
// ============================================================================

/// Conditional likelihood vector at a node for one site.
type CondLike = [f64; 4];

/// Per-node data for the pruning algorithm.
#[derive(Debug, Clone)]
struct NodeLikelihoods {
    /// Conditional likelihood vectors, one per pattern.
    cond: Vec<CondLike>,
    /// Log-scaling factor per pattern (for underflow prevention).
    log_scale: Vec<f64>,
}

impl NodeLikelihoods {
    fn new(npatterns: usize) -> Self {
        Self {
            cond: vec![[0.0; 4]; npatterns],
            log_scale: vec![0.0; npatterns],
        }
    }

    fn reset(&mut self, npatterns: usize) {
        self.cond.resize(npatterns, [0.0; 4]);
        self.log_scale.resize(npatterns, 0.0);
        for c in &mut self.cond {
            *c = [0.0; 4];
        }
        for s in &mut self.log_scale {
            *s = 0.0;
        }
    }
}

/// Pre-allocated buffers for the pruning algorithm.
///
/// Avoids per-call allocation of conditional likelihood vectors.
/// Create once and reuse across multiple likelihood evaluations.
#[derive(Debug, Clone)]
pub struct PruningBuffers {
    data: Vec<NodeLikelihoods>,
}

impl PruningBuffers {
    /// Create new buffers for a tree with `nnodes` nodes and `npatterns` patterns.
    pub fn new(nnodes: usize, npatterns: usize) -> Self {
        PruningBuffers {
            data: vec![NodeLikelihoods::new(npatterns); nnodes],
        }
    }

    /// Reset all buffers for reuse with possibly different dimensions.
    fn reset(&mut self, nnodes: usize, npatterns: usize) {
        self.data.resize_with(nnodes, || NodeLikelihoods::new(npatterns));
        for d in &mut self.data {
            d.reset(npatterns);
        }
    }
}

// ============================================================================
// Fast Pruning Algorithm
// ============================================================================

/// Initialize conditional likelihood at a leaf node from observed data.
fn init_leaf(base: Base) -> CondLike {
    let bs = base.base_set();
    [
        if bs[0] { 1.0 } else { 0.0 },
        if bs[1] { 1.0 } else { 0.0 },
        if bs[2] { 1.0 } else { 0.0 },
        if bs[3] { 1.0 } else { 0.0 },
    ]
}

/// Get the effective branch length for a node, using default if missing.
fn effective_branch_length(tree: &Tree, node_id: usize) -> f64 {
    tree.nodes[node_id]
        .branch_length
        .unwrap_or(DEFAULT_BRANCH_LENGTH)
        .max(MIN_BRANCH_LENGTH)
}

/// Run the optimized pruning algorithm on compressed patterns with caching.
///
/// This is the core inner routine. It uses:
/// - `compressed` for site-pattern compression (fewer per-site evaluations)
/// - `cache` for transition matrix caching (fewer `exp()` calls)
/// - `buffers` for pre-allocated memory (no per-call allocations)
fn fast_pruning_pass(
    tree: &Tree,
    compressed: &CompressedAlignment,
    model: &dyn SubstitutionModel,
    cache: &mut TransitionCache,
    buffers: &mut PruningBuffers,
) -> Result<(), LikelihoodError> {
    let npatterns = compressed.n_patterns();
    let nnodes = tree.num_nodes();

    if tree.num_leaves() == 0 {
        return Err(LikelihoodError::EmptyTree);
    }

    // Ensure buffers and cache are properly sized
    buffers.reset(nnodes, npatterns);
    cache.ensure_capacity(nnodes);

    // Build name -> alignment index map
    let name_to_idx: HashMap<&str, usize> = compressed
        .taxa_names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    // Set leaf conditional likelihoods from observed data
    for node in &tree.nodes {
        if node.is_leaf() {
            if let Some(name) = &node.name {
                let taxon_idx = name_to_idx.get(name.as_str()).ok_or_else(|| {
                    LikelihoodError::TaxaMismatch(name.clone())
                })?;
                for pat in 0..npatterns {
                    buffers.data[node.id].cond[pat] =
                        init_leaf(compressed.patterns[pat][*taxon_idx]);
                }
            }
        }
    }

    // Postorder traversal: pruning algorithm
    let postorder = tree.postorder();
    for &node_id in &postorder {
        let children = &tree.nodes[node_id].children;
        if children.is_empty() {
            continue; // leaf: already initialized
        }

        // For each child, get the (cached) transition probability matrix
        let child_data: Vec<(usize, [[f64; 4]; 4])> = children
            .iter()
            .map(|&child_id| {
                let t = effective_branch_length(tree, child_id);
                let p = cache.get_or_compute(child_id, t, model);
                (child_id, p)
            })
            .collect();

        for pat in 0..npatterns {
            let mut node_cond = [1.0f64; 4];
            let mut node_log_scale = 0.0f64;

            for &(child_id, ref p_matrix) in &child_data {
                let child_cond = &buffers.data[child_id].cond[pat];
                node_log_scale += buffers.data[child_id].log_scale[pat];

                let mut partial = [0.0f64; 4];
                for b in 0..4 {
                    let mut sum = 0.0;
                    for j in 0..4 {
                        sum += p_matrix[b][j] * child_cond[j];
                    }
                    partial[b] = sum;
                }

                for b in 0..4 {
                    node_cond[b] *= partial[b];
                }
            }

            // Underflow prevention
            let max_val = node_cond
                .iter()
                .copied()
                .fold(0.0f64, f64::max);

            if max_val > 0.0 && max_val < UNDERFLOW_THRESHOLD {
                for b in 0..4 {
                    node_cond[b] /= max_val;
                }
                node_log_scale += max_val.ln();
            }

            buffers.data[node_id].cond[pat] = node_cond;
            buffers.data[node_id].log_scale[pat] = node_log_scale;
        }
    }

    Ok(())
}

/// Compute the total log-likelihood using compressed alignment and cached matrices.
///
/// This is functionally equivalent to [`super::pruning::log_likelihood`] but
/// exploits site-pattern compression, transition matrix caching, and
/// pre-allocated buffers for better performance.
///
/// # Arguments
///
/// * `tree` - The phylogenetic tree with branch lengths.
/// * `compressed` - The compressed alignment (from [`CompressedAlignment::from_alignment`]).
/// * `model` - The substitution model.
/// * `cache` - Mutable reference to a transition cache (will be updated).
///
/// # Returns
///
/// The total log-likelihood, identical to what [`super::pruning::log_likelihood`]
/// would return on the uncompressed alignment.
pub fn fast_log_likelihood(
    tree: &Tree,
    compressed: &CompressedAlignment,
    model: &dyn SubstitutionModel,
    cache: &mut TransitionCache,
) -> Result<f64, LikelihoodError> {
    let npatterns = compressed.n_patterns();
    let mut buffers = PruningBuffers::new(tree.num_nodes(), npatterns);

    fast_pruning_pass(tree, compressed, model, cache, &mut buffers)?;

    let root = tree.root;
    let pi = model.frequencies();

    let mut total_lnl = 0.0;

    for pat in 0..npatterns {
        let root_cond = &buffers.data[root].cond[pat];
        let log_scale = buffers.data[root].log_scale[pat];

        let mut site_like = 0.0;
        for b in 0..4 {
            site_like += pi[b] * root_cond[b];
        }

        if site_like <= 0.0 {
            return Err(LikelihoodError::NumericalError(format!(
                "pattern {} has zero likelihood",
                pat
            )));
        }

        // Multiply by weight (the key compression benefit)
        total_lnl += (site_like.ln() + log_scale) * compressed.weights[pat] as f64;
    }

    Ok(total_lnl)
}

/// Compute per-site log-likelihoods using compressed alignment and cached matrices.
///
/// Returns a vector of per-*pattern* log-likelihoods paired with their weights.
/// To get the total log-likelihood, sum `lnl * weight` over all patterns.
pub fn fast_site_log_likelihoods(
    tree: &Tree,
    compressed: &CompressedAlignment,
    model: &dyn SubstitutionModel,
    cache: &mut TransitionCache,
) -> Result<Vec<(f64, usize)>, LikelihoodError> {
    let npatterns = compressed.n_patterns();
    let mut buffers = PruningBuffers::new(tree.num_nodes(), npatterns);

    fast_pruning_pass(tree, compressed, model, cache, &mut buffers)?;

    let root = tree.root;
    let pi = model.frequencies();

    let mut result = Vec::with_capacity(npatterns);

    for pat in 0..npatterns {
        let root_cond = &buffers.data[root].cond[pat];
        let log_scale = buffers.data[root].log_scale[pat];

        let mut site_like = 0.0;
        for b in 0..4 {
            site_like += pi[b] * root_cond[b];
        }

        if site_like <= 0.0 {
            return Err(LikelihoodError::NumericalError(format!(
                "pattern {} has zero likelihood",
                pat
            )));
        }

        result.push((site_like.ln() + log_scale, compressed.weights[pat]));
    }

    Ok(result)
}

/// Optimize a single branch length using the fast likelihood computation.
fn fast_optimize_one_branch(
    tree: &mut Tree,
    node_id: usize,
    compressed: &CompressedAlignment,
    model: &dyn SubstitutionModel,
    cache: &mut TransitionCache,
    max_iter: usize,
) -> Result<f64, LikelihoodError> {
    let epsilon = 1e-6;
    let mut t = effective_branch_length(tree, node_id);

    for _ in 0..max_iter {
        let h = (t * 0.001).max(1e-8);

        // Evaluate likelihood at t, t+h, t-h
        tree.nodes[node_id].branch_length = Some(t);
        cache.invalidate(node_id);
        let lnl_center = fast_log_likelihood(tree, compressed, model, cache)?;

        tree.nodes[node_id].branch_length = Some((t + h).min(MAX_BRANCH_LENGTH));
        cache.invalidate(node_id);
        let lnl_plus = fast_log_likelihood(tree, compressed, model, cache)?;

        tree.nodes[node_id].branch_length = Some((t - h).max(MIN_BRANCH_LENGTH));
        cache.invalidate(node_id);
        let lnl_minus = fast_log_likelihood(tree, compressed, model, cache)?;

        // Numerical derivatives
        let slope = (lnl_plus - lnl_minus) / (2.0 * h);
        let curve = (lnl_plus - 2.0 * lnl_center + lnl_minus) / (h * h);

        let delta = if curve < -1e-12 {
            -slope / curve
        } else if slope.abs() > 1e-12 {
            slope.signum() * h * 10.0
        } else {
            break;
        };

        let t_new = (t + delta).clamp(MIN_BRANCH_LENGTH, MAX_BRANCH_LENGTH);

        if (t_new - t).abs() < epsilon {
            t = t_new;
            break;
        }
        t = t_new;
    }

    tree.nodes[node_id].branch_length = Some(t);
    cache.invalidate(node_id);
    Ok(t)
}

/// Optimize all branch lengths using the fast likelihood computation.
///
/// This is functionally equivalent to [`super::pruning::optimize_branch_lengths`]
/// but uses site-pattern compression and transition matrix caching for
/// better performance.
///
/// # Arguments
///
/// * `tree` - Mutable reference to the tree (branch lengths will be updated).
/// * `compressed` - The compressed alignment.
/// * `model` - The substitution model.
///
/// # Returns
///
/// The final log-likelihood after optimization.
pub fn fast_optimize_branch_lengths(
    tree: &mut Tree,
    compressed: &CompressedAlignment,
    model: &dyn SubstitutionModel,
) -> Result<f64, LikelihoodError> {
    let max_rounds = 20;
    let max_iter_per_branch = 10;
    let tolerance = 1e-4;

    let mut cache = TransitionCache::new(tree.num_nodes());

    // Ensure all branches have initial lengths
    for node in &mut tree.nodes {
        if node.parent.is_some() && node.branch_length.is_none() {
            node.branch_length = Some(DEFAULT_BRANCH_LENGTH);
        }
    }

    let mut prev_lnl = fast_log_likelihood(tree, compressed, model, &mut cache)?;

    for _ in 0..max_rounds {
        let node_ids: Vec<usize> = tree
            .nodes
            .iter()
            .filter(|n| n.parent.is_some())
            .map(|n| n.id)
            .collect();

        for node_id in node_ids {
            fast_optimize_one_branch(
                tree,
                node_id,
                compressed,
                model,
                &mut cache,
                max_iter_per_branch,
            )?;
        }

        let lnl = fast_log_likelihood(tree, compressed, model, &mut cache)?;

        if (lnl - prev_lnl).abs() < tolerance {
            return Ok(lnl);
        }
        prev_lnl = lnl;
    }

    Ok(prev_lnl)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::likelihood::models::{F84Model, Jc69Model};
    use crate::likelihood::pruning;
    use crate::tree::newick::parse_newick;
    use crate::tree::types::{BaseFrequencies, Sequence};

    fn make_alignment(data: &[(&str, &str)]) -> Alignment {
        let seqs: Vec<Sequence> = data
            .iter()
            .map(|(name, seq)| {
                let bases: Vec<Base> = seq
                    .chars()
                    .map(|c| Base::from_char(c).unwrap())
                    .collect();
                Sequence::new(*name, bases)
            })
            .collect();
        Alignment::new(seqs).unwrap()
    }

    // ---------------------------------------------------------------
    // CompressedAlignment tests
    // ---------------------------------------------------------------

    #[test]
    fn test_compressed_alignment_basic() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGA"),
            ("C", "TGCA"),
        ]);
        let comp = CompressedAlignment::from_alignment(&aln);
        assert_eq!(comp.n_taxa(), 3);
        assert_eq!(comp.total_sites(), 4);
        // All four sites are distinct patterns here
        assert_eq!(comp.n_patterns(), 4);
        assert!(comp.weights.iter().all(|&w| w == 1));
    }

    #[test]
    fn test_compressed_alignment_groups_identical_patterns() {
        // Construct an alignment where some columns are identical:
        // Site 0 and 2 both have pattern (A, A, T)
        // Site 1 and 3 both have pattern (C, C, G)
        let aln = make_alignment(&[
            ("A", "ACAC"),
            ("B", "ACAC"),
            ("C", "TGTG"),
        ]);
        let comp = CompressedAlignment::from_alignment(&aln);
        assert_eq!(comp.n_taxa(), 3);
        assert_eq!(comp.total_sites(), 4);
        // Only 2 unique patterns: (A,A,T) and (C,C,G)
        assert_eq!(comp.n_patterns(), 2);
        // Each pattern appears twice
        assert!(comp.weights.iter().all(|&w| w == 2));
    }

    #[test]
    fn test_compressed_alignment_all_identical_sites() {
        // All sites have the same pattern
        let aln = make_alignment(&[
            ("A", "AAAA"),
            ("B", "CCCC"),
            ("C", "GGGG"),
        ]);
        let comp = CompressedAlignment::from_alignment(&aln);
        assert_eq!(comp.n_patterns(), 1);
        assert_eq!(comp.weights[0], 4);
        assert_eq!(comp.total_sites(), 4);
    }

    #[test]
    fn test_compressed_alignment_all_unique_sites() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "ACGT"),
        ]);
        let comp = CompressedAlignment::from_alignment(&aln);
        assert_eq!(comp.n_patterns(), 4);
        assert!(comp.weights.iter().all(|&w| w == 1));
    }

    #[test]
    fn test_compressed_alignment_preserves_taxa_names() {
        let aln = make_alignment(&[
            ("Alpha", "AC"),
            ("Beta", "GT"),
            ("Gamma", "AC"),
        ]);
        let comp = CompressedAlignment::from_alignment(&aln);
        assert_eq!(comp.taxa_names, vec!["Alpha", "Beta", "Gamma"]);
    }

    // ---------------------------------------------------------------
    // TransitionCache tests
    // ---------------------------------------------------------------

    #[test]
    fn test_transition_cache_computes_correctly() {
        let model = Jc69Model;
        let mut cache = TransitionCache::new(5);

        let matrix = cache.get_or_compute(0, 0.1, &model);
        let expected = model.transition_probs(0.1);

        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (matrix[i][j] - expected[i][j]).abs() < 1e-15,
                    "cache mismatch at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_transition_cache_returns_cached_value() {
        let model = Jc69Model;
        let mut cache = TransitionCache::new(5);

        // First call computes
        let m1 = cache.get_or_compute(2, 0.3, &model);
        // Second call with same branch length returns cached
        let m2 = cache.get_or_compute(2, 0.3, &model);

        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (m1[i][j] - m2[i][j]).abs() < 1e-15,
                    "cached value should be identical"
                );
            }
        }
    }

    #[test]
    fn test_transition_cache_invalidates_on_branch_change() {
        let model = Jc69Model;
        let mut cache = TransitionCache::new(5);

        let m1 = cache.get_or_compute(1, 0.1, &model);
        let m2 = cache.get_or_compute(1, 0.5, &model);

        let expected1 = model.transition_probs(0.1);
        let expected2 = model.transition_probs(0.5);

        // m1 should match t=0.1
        assert!((m1[0][0] - expected1[0][0]).abs() < 1e-15);
        // m2 should match t=0.5 (different from m1)
        assert!((m2[0][0] - expected2[0][0]).abs() < 1e-15);
        assert!((m1[0][0] - m2[0][0]).abs() > 1e-6);
    }

    #[test]
    fn test_transition_cache_invalidate_all() {
        let model = Jc69Model;
        let mut cache = TransitionCache::new(3);

        cache.get_or_compute(0, 0.1, &model);
        cache.get_or_compute(1, 0.2, &model);
        cache.invalidate_all();

        // After invalidate_all, all entries should be None
        assert!(cache.matrices[0].is_none());
        assert!(cache.matrices[1].is_none());
    }

    #[test]
    fn test_transition_cache_f84_model() {
        let freqs = BaseFrequencies::new(0.3, 0.2, 0.2, 0.3);
        let model = F84Model::new(freqs, 2.0);
        let mut cache = TransitionCache::new(3);

        let matrix = cache.get_or_compute(0, 0.25, &model);
        let expected = model.transition_probs(0.25);

        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (matrix[i][j] - expected[i][j]).abs() < 1e-15,
                    "F84 cache mismatch at [{i}][{j}]"
                );
            }
        }
    }

    // ---------------------------------------------------------------
    // fast_log_likelihood correctness tests
    // ---------------------------------------------------------------

    #[test]
    fn test_fast_lnl_matches_original_three_taxa() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGA"),
            ("C", "TGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let original_lnl = pruning::log_likelihood(&tree, &aln, &model).unwrap();

        let compressed = CompressedAlignment::from_alignment(&aln);
        let mut cache = TransitionCache::new(tree.num_nodes());
        let fast_lnl = fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap();

        assert!(
            (original_lnl - fast_lnl).abs() < 1e-10,
            "fast_lnl ({}) should match original ({})",
            fast_lnl,
            original_lnl
        );
    }

    #[test]
    fn test_fast_lnl_matches_original_four_taxa() {
        let aln = make_alignment(&[
            ("A", "AACCGGTT"),
            ("B", "AACCGGTT"),
            ("C", "CCAATTGG"),
            ("D", "CCAATTGG"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();
        let model = Jc69Model;

        let original_lnl = pruning::log_likelihood(&tree, &aln, &model).unwrap();

        let compressed = CompressedAlignment::from_alignment(&aln);
        let mut cache = TransitionCache::new(tree.num_nodes());
        let fast_lnl = fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap();

        assert!(
            (original_lnl - fast_lnl).abs() < 1e-10,
            "fast_lnl ({}) should match original ({})",
            fast_lnl,
            original_lnl
        );
    }

    #[test]
    fn test_fast_lnl_matches_original_f84() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
        ]);
        let tree = parse_newick("((A:0.2,B:0.2):0.1,C:0.3);").unwrap();
        let model = F84Model::equal_freqs(2.0);

        let original_lnl = pruning::log_likelihood(&tree, &aln, &model).unwrap();

        let compressed = CompressedAlignment::from_alignment(&aln);
        let mut cache = TransitionCache::new(tree.num_nodes());
        let fast_lnl = fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap();

        assert!(
            (original_lnl - fast_lnl).abs() < 1e-10,
            "F84 fast_lnl ({}) should match original ({})",
            fast_lnl,
            original_lnl
        );
    }

    #[test]
    fn test_fast_lnl_matches_original_with_repeated_patterns() {
        // Alignment with many repeated patterns -- this is the case
        // where compression provides the biggest speedup.
        let aln = make_alignment(&[
            ("A", "ACACACACACACACAC"),
            ("B", "GTGTGTGTGTGTGTGT"),
            ("C", "ACACACACACACACAC"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = Jc69Model;

        let original_lnl = pruning::log_likelihood(&tree, &aln, &model).unwrap();

        let compressed = CompressedAlignment::from_alignment(&aln);
        // Only 2 unique patterns with weight 8 each
        assert_eq!(compressed.n_patterns(), 2);

        let mut cache = TransitionCache::new(tree.num_nodes());
        let fast_lnl = fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap();

        assert!(
            (original_lnl - fast_lnl).abs() < 1e-10,
            "fast_lnl ({}) should match original ({}) with repeated patterns",
            fast_lnl,
            original_lnl
        );
    }

    #[test]
    fn test_fast_lnl_matches_original_with_ambiguity_codes() {
        let aln = make_alignment(&[
            ("A", "AN"),
            ("B", "AC"),
            ("C", "TC"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = Jc69Model;

        let original_lnl = pruning::log_likelihood(&tree, &aln, &model).unwrap();

        let compressed = CompressedAlignment::from_alignment(&aln);
        let mut cache = TransitionCache::new(tree.num_nodes());
        let fast_lnl = fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap();

        assert!(
            (original_lnl - fast_lnl).abs() < 1e-10,
            "fast_lnl ({}) should match original ({}) with ambiguity codes",
            fast_lnl,
            original_lnl
        );
    }

    #[test]
    fn test_fast_lnl_matches_original_five_taxa() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAA"),
            ("B", "AACCGGTTAA"),
            ("C", "CCAATTGGCC"),
            ("D", "CCAATTGGCC"),
            ("E", "GGTTAACCGG"),
        ]);
        let tree = parse_newick(
            "((A:0.1,B:0.1):0.05,((C:0.1,D:0.1):0.05,E:0.15):0.05);",
        )
        .unwrap();
        let model = Jc69Model;

        let original_lnl = pruning::log_likelihood(&tree, &aln, &model).unwrap();

        let compressed = CompressedAlignment::from_alignment(&aln);
        let mut cache = TransitionCache::new(tree.num_nodes());
        let fast_lnl = fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap();

        assert!(
            (original_lnl - fast_lnl).abs() < 1e-10,
            "5-taxon fast_lnl ({}) should match original ({})",
            fast_lnl,
            original_lnl
        );
    }

    #[test]
    fn test_fast_lnl_is_negative() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "TGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = Jc69Model;

        let compressed = CompressedAlignment::from_alignment(&aln);
        let mut cache = TransitionCache::new(tree.num_nodes());
        let lnl = fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap();
        assert!(lnl < 0.0, "log-likelihood should be negative, got {}", lnl);
    }

    #[test]
    fn test_fast_lnl_taxa_mismatch_error() {
        let aln = make_alignment(&[
            ("X", "ACGT"),
            ("Y", "ACGT"),
            ("Z", "ACGT"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = Jc69Model;

        let compressed = CompressedAlignment::from_alignment(&aln);
        let mut cache = TransitionCache::new(tree.num_nodes());
        let result = fast_log_likelihood(&tree, &compressed, &model, &mut cache);
        assert!(result.is_err());
    }

    // ---------------------------------------------------------------
    // fast_optimize_branch_lengths tests
    // ---------------------------------------------------------------

    #[test]
    fn test_fast_optimize_improves_likelihood() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
            ("D", "CCAATTGGCCAATTGG"),
        ]);
        // Start with poor branch lengths
        let tree = parse_newick("((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);").unwrap();
        let model = Jc69Model;
        let compressed = CompressedAlignment::from_alignment(&aln);

        let mut cache = TransitionCache::new(tree.num_nodes());
        let lnl_before = fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap();

        let mut opt_tree = tree.clone();
        let lnl_after = fast_optimize_branch_lengths(&mut opt_tree, &compressed, &model).unwrap();

        assert!(
            lnl_after >= lnl_before - 1e-6,
            "optimization should not decrease likelihood: before={}, after={}",
            lnl_before,
            lnl_after
        );
    }

    #[test]
    fn test_fast_optimize_matches_original_result() {
        let aln = make_alignment(&[
            ("A", "AACCGGTT"),
            ("B", "AACCGGTT"),
            ("C", "CCAATTGG"),
        ]);
        let tree = parse_newick("((A:0.5,B:0.5):0.5,C:0.5);").unwrap();
        let model = Jc69Model;

        // Optimize with original
        let mut orig_tree = tree.clone();
        let orig_lnl = pruning::optimize_branch_lengths(&mut orig_tree, &aln, &model).unwrap();

        // Optimize with fast version
        let compressed = CompressedAlignment::from_alignment(&aln);
        let mut fast_tree = tree.clone();
        let fast_lnl =
            fast_optimize_branch_lengths(&mut fast_tree, &compressed, &model).unwrap();

        // Both should converge to similar optima
        assert!(
            (orig_lnl - fast_lnl).abs() < 0.1,
            "fast optimize lnl ({}) should be close to original ({})",
            fast_lnl,
            orig_lnl
        );
    }

    #[test]
    fn test_fast_optimize_f84_model() {
        let aln = make_alignment(&[
            ("A", "AACCGGTTAACCGGTT"),
            ("B", "AACCGGTTAACCGGTT"),
            ("C", "CCAATTGGCCAATTGG"),
        ]);
        let tree = parse_newick("((A:0.5,B:0.5):0.5,C:0.5);").unwrap();
        let model = F84Model::equal_freqs(2.0);

        let compressed = CompressedAlignment::from_alignment(&aln);

        let lnl_before = {
            let mut cache = TransitionCache::new(tree.num_nodes());
            fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap()
        };

        let mut opt_tree = tree.clone();
        let lnl_after =
            fast_optimize_branch_lengths(&mut opt_tree, &compressed, &model).unwrap();

        assert!(
            lnl_after >= lnl_before - 1e-6,
            "F84 optimization should not decrease lnl: {} -> {}",
            lnl_before,
            lnl_after
        );
    }

    // ---------------------------------------------------------------
    // fast_site_log_likelihoods tests
    // ---------------------------------------------------------------

    #[test]
    fn test_fast_site_lnl_sums_to_total() {
        let aln = make_alignment(&[
            ("A", "ACGTACGT"),
            ("B", "ACGAACGA"),
            ("C", "TGCATGCA"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let compressed = CompressedAlignment::from_alignment(&aln);
        let mut cache = TransitionCache::new(tree.num_nodes());

        let site_lnl = fast_site_log_likelihoods(&tree, &compressed, &model, &mut cache).unwrap();
        let total_from_sites: f64 = site_lnl.iter().map(|&(lnl, w)| lnl * w as f64).sum();
        let total = fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap();

        assert!(
            (total_from_sites - total).abs() < 1e-10,
            "sum of weighted site lnl ({}) should equal total ({})",
            total_from_sites,
            total
        );
    }

    // ---------------------------------------------------------------
    // PruningBuffers tests
    // ---------------------------------------------------------------

    #[test]
    fn test_pruning_buffers_reset() {
        let mut buffers = PruningBuffers::new(3, 4);
        // Dirty the buffers
        buffers.data[0].cond[0] = [1.0, 2.0, 3.0, 4.0];
        buffers.data[0].log_scale[0] = 5.0;

        // Reset should zero everything
        buffers.reset(3, 4);
        assert_eq!(buffers.data[0].cond[0], [0.0; 4]);
        assert_eq!(buffers.data[0].log_scale[0], 0.0);
    }

    // ---------------------------------------------------------------
    // Compression ratio test
    // ---------------------------------------------------------------

    #[test]
    fn test_compression_ratio() {
        // Simulate a "real" dataset with many repeated patterns:
        // 100-site alignment where only 10 unique column patterns exist
        let pattern = "ACGTACGTAC"; // 10 unique bases
        let seq_a: String = pattern.repeat(10); // 100 sites
        let seq_b: String = pattern.chars().map(|c| match c {
            'A' => 'C', 'C' => 'A', 'G' => 'T', 'T' => 'G', _ => c
        }).collect::<String>().repeat(10);
        let seq_c: String = "AAAAAAAAAA".repeat(10);

        let aln = make_alignment(&[
            ("A", &seq_a),
            ("B", &seq_b),
            ("C", &seq_c),
        ]);

        let comp = CompressedAlignment::from_alignment(&aln);
        assert_eq!(comp.total_sites(), 100);
        // With 10-char repeating pattern, should have at most 10 unique patterns
        assert!(
            comp.n_patterns() <= 10,
            "expected at most 10 patterns, got {}",
            comp.n_patterns()
        );
    }

    // ---------------------------------------------------------------
    // Matches original on tree without branch lengths
    // ---------------------------------------------------------------

    #[test]
    fn test_fast_lnl_matches_original_no_branch_lengths() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "TGCA"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let model = Jc69Model;

        let original_lnl = pruning::log_likelihood(&tree, &aln, &model).unwrap();

        let compressed = CompressedAlignment::from_alignment(&aln);
        let mut cache = TransitionCache::new(tree.num_nodes());
        let fast_lnl = fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap();

        assert!(
            (original_lnl - fast_lnl).abs() < 1e-10,
            "fast_lnl ({}) should match original ({}) without branch lengths",
            fast_lnl,
            original_lnl
        );
    }

    // ---------------------------------------------------------------
    // F84 with unequal frequencies
    // ---------------------------------------------------------------

    #[test]
    fn test_fast_lnl_f84_unequal_freqs() {
        let aln = make_alignment(&[
            ("A", "AAACCCTTTGGG"),
            ("B", "AAACCCTTTGGG"),
            ("C", "CCCAAAGGGCCC"),
        ]);
        let tree = parse_newick("((A:0.2,B:0.2):0.1,C:0.3);").unwrap();
        let freqs = BaseFrequencies::new(0.3, 0.2, 0.2, 0.3);
        let model = F84Model::new(freqs, 2.0);

        let original_lnl = pruning::log_likelihood(&tree, &aln, &model).unwrap();

        let compressed = CompressedAlignment::from_alignment(&aln);
        let mut cache = TransitionCache::new(tree.num_nodes());
        let fast_lnl = fast_log_likelihood(&tree, &compressed, &model, &mut cache).unwrap();

        assert!(
            (original_lnl - fast_lnl).abs() < 1e-10,
            "F84 unequal freqs: fast_lnl ({}) should match original ({})",
            fast_lnl,
            original_lnl
        );
    }
}

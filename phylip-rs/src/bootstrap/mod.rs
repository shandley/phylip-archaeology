//! Bootstrap resampling for phylogenetic confidence estimation.
//!
//! Implementation of Felsenstein's bootstrap method (1985) for
//! assessing confidence in phylogenetic tree topology.
//!
//! The bootstrap resamples alignment columns with replacement, producing
//! pseudoreplicate datasets. Trees inferred from each replicate are then
//! compared to estimate support for each branch (clade) in the original tree.
//!
//! # Supported resampling methods
//!
//! - **Bootstrap**: Sample sites with replacement (Felsenstein 1985)
//! - **Block bootstrap**: Sample contiguous blocks of sites (for correlated data)
//! - **Jackknife**: Sample a fraction of sites without replacement
//! - **Permutation**: Permute site assignments across taxa (for ILD tests)
//!
//! # Example
//!
//! ```
//! use phylip_rs::bootstrap::{ResamplingMethod, bootstrap_weights, resample_alignment};
//! use phylip_rs::tree::{Alignment, Sequence, Base};
//!
//! let alignment = Alignment::new(vec![
//!     Sequence::new("A", vec![Base::A, Base::C, Base::G, Base::T]),
//!     Sequence::new("B", vec![Base::T, Base::G, Base::C, Base::A]),
//! ]).unwrap();
//!
//! let mut rng = phylip_rs::bootstrap::SimpleRng::seed(42);
//! let weights = bootstrap_weights(alignment.nsites(), &ResamplingMethod::Bootstrap, &mut rng);
//! let resampled = resample_alignment(&alignment, &weights);
//! assert_eq!(resampled.nsites(), alignment.nsites());
//! ```
//!
//! Reference: Felsenstein, J. (1985). Confidence limits on phylogenies:
//! an approach using the bootstrap. Evolution, 39, 783-791.

pub mod ml;
pub mod support;

use crate::tree::{Alignment, Sequence};
use std::fmt;

// ============================================================
// Random number generator (no external dependencies)
// ============================================================

/// A simple linear congruential generator for reproducible resampling.
///
/// Uses the same constants as glibc's `rand()` for wide compatibility.
/// For production use with security requirements, use a cryptographic RNG;
/// for phylogenetic bootstrap, reproducibility matters more than randomness quality.
pub struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    /// Create a new RNG with the given seed.
    pub fn seed(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Generate the next pseudo-random u64.
    fn next_u64(&mut self) -> u64 {
        // LCG parameters from Knuth MMIX
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.state
    }

    /// Generate a random usize in [0, n).
    pub fn gen_range(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    /// Generate a random f64 in [0, 1).
    pub fn gen_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

// ============================================================
// Resampling method configuration
// ============================================================

/// Resampling method for generating pseudoreplicate datasets.
#[derive(Debug, Clone)]
pub enum ResamplingMethod {
    /// Standard bootstrap: sample n sites with replacement.
    Bootstrap,

    /// Block bootstrap: sample contiguous blocks of sites with replacement.
    /// The block_size parameter controls the length of each block.
    /// Blocks wrap around at sequence boundaries (circular).
    BlockBootstrap { block_size: usize },

    /// Partial bootstrap: sample a fraction of n sites with replacement.
    PartialBootstrap { fraction: f64 },

    /// Delete-fraction jackknife: sample a fraction of sites without replacement.
    /// Default fraction is 0.5 (delete-half jackknife).
    Jackknife { fraction: f64 },
}

/// Errors that can occur during bootstrap operations.
#[derive(Debug, Clone)]
pub enum BootstrapError {
    /// Alignment has zero sites.
    EmptySites,
    /// Invalid parameter value.
    InvalidParameter(String),
    /// No replicate trees provided.
    NoReplicates,
    /// Tree has no leaves.
    EmptyTree,
}

impl fmt::Display for BootstrapError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BootstrapError::EmptySites => write!(f, "alignment has zero sites"),
            BootstrapError::InvalidParameter(msg) => write!(f, "invalid parameter: {}", msg),
            BootstrapError::NoReplicates => write!(f, "no replicate trees provided"),
            BootstrapError::EmptyTree => write!(f, "tree has no leaves"),
        }
    }
}

impl std::error::Error for BootstrapError {}

// ============================================================
// Weight generation
// ============================================================

/// Generate a weight vector for one bootstrap/jackknife replicate.
///
/// Each element `weights[i]` indicates how many times site `i` appears
/// in the resampled dataset. For bootstrap, weights follow a multinomial
/// distribution; for jackknife, each weight is 0 or 1.
///
/// The sum of weights equals the resampled dataset size:
/// - Bootstrap: sum = n (original number of sites)
/// - Partial bootstrap: sum = floor(fraction * n)
/// - Jackknife: sum = floor(fraction * n)
pub fn bootstrap_weights(
    nsites: usize,
    method: &ResamplingMethod,
    rng: &mut SimpleRng,
) -> Vec<usize> {
    match method {
        ResamplingMethod::Bootstrap => {
            standard_bootstrap_weights(nsites, nsites, 1, rng)
        }
        ResamplingMethod::BlockBootstrap { block_size } => {
            let bs = (*block_size).max(1);
            standard_bootstrap_weights(nsites, nsites, bs, rng)
        }
        ResamplingMethod::PartialBootstrap { fraction } => {
            let frac = fraction.clamp(0.0, 1.0);
            let count = (nsites as f64 * frac + 0.5) as usize;
            standard_bootstrap_weights(nsites, count, 1, rng)
        }
        ResamplingMethod::Jackknife { fraction } => {
            jackknife_weights(nsites, *fraction, rng)
        }
    }
}

/// Standard (block) bootstrap: sample `count` sites in blocks of `block_size`.
fn standard_bootstrap_weights(
    nsites: usize,
    count: usize,
    block_size: usize,
    rng: &mut SimpleRng,
) -> Vec<usize> {
    let mut weights = vec![0usize; nsites];
    if nsites == 0 {
        return weights;
    }

    let blocks = count / block_size;
    for _ in 0..blocks {
        let start = rng.gen_range(nsites);
        for k in 0..block_size {
            let idx = (start + k) % nsites; // circular wrap
            weights[idx] += 1;
        }
    }
    weights
}

/// Jackknife: select exactly `floor(fraction * n)` sites without replacement.
///
/// Uses a streaming selection algorithm (like PHYLIP's seqboot):
/// each site is included with a probability that adapts to maintain
/// the exact target count.
fn jackknife_weights(nsites: usize, fraction: f64, rng: &mut SimpleRng) -> Vec<usize> {
    let mut weights = vec![0usize; nsites];
    if nsites == 0 {
        return weights;
    }

    let frac = fraction.clamp(0.0, 1.0);
    let mut remaining_to_select = (nsites as f64 * frac + 0.5) as usize;
    let mut remaining_sites = nsites;

    for i in 0..nsites {
        if remaining_sites == 0 {
            break;
        }
        let p = remaining_to_select as f64 / remaining_sites as f64;
        if rng.gen_f64() < p {
            weights[i] = 1;
            remaining_to_select -= 1;
        }
        remaining_sites -= 1;
    }
    weights
}

// ============================================================
// Alignment resampling
// ============================================================

/// Produce a resampled alignment from weights.
///
/// Each site `i` in the original alignment is repeated `weights[i]` times
/// in the output. The total length of the resampled alignment is the
/// sum of all weights.
pub fn resample_alignment(alignment: &Alignment, weights: &[usize]) -> Alignment {
    let nsites = alignment.nsites();
    let new_len: usize = weights.iter().sum();

    let sequences: Vec<Sequence> = alignment
        .sequences
        .iter()
        .map(|seq| {
            let mut new_bases = Vec::with_capacity(new_len);
            for (site, &w) in weights.iter().enumerate().take(nsites) {
                for _ in 0..w {
                    new_bases.push(seq.bases[site]);
                }
            }
            Sequence::new(seq.name.clone(), new_bases)
        })
        .collect();

    // Safety: all sequences have the same length (new_len)
    Alignment::new(sequences).expect("resampled alignment should be valid")
}

/// Generate multiple bootstrap replicate alignments.
///
/// Returns an iterator over `nreps` resampled alignments.
pub fn bootstrap_replicates<'a>(
    alignment: &'a Alignment,
    nreps: usize,
    method: &'a ResamplingMethod,
    rng: &'a mut SimpleRng,
) -> impl Iterator<Item = Alignment> + 'a {
    (0..nreps).map(move |_| {
        let weights = bootstrap_weights(alignment.nsites(), method, rng);
        resample_alignment(alignment, &weights)
    })
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::{Alignment, Sequence, Base};

    fn make_test_alignment() -> Alignment {
        Alignment::new(vec![
            Sequence::new("TaxonA", vec![Base::A, Base::C, Base::G, Base::T, Base::A, Base::C, Base::G, Base::T]),
            Sequence::new("TaxonB", vec![Base::T, Base::G, Base::C, Base::A, Base::T, Base::G, Base::C, Base::A]),
            Sequence::new("TaxonC", vec![Base::A, Base::A, Base::G, Base::G, Base::T, Base::T, Base::C, Base::C]),
        ])
        .unwrap()
    }

    #[test]
    fn test_bootstrap_weights_sum() {
        let mut rng = SimpleRng::seed(12345);
        for _ in 0..100 {
            let weights = bootstrap_weights(20, &ResamplingMethod::Bootstrap, &mut rng);
            assert_eq!(weights.len(), 20);
            assert_eq!(weights.iter().sum::<usize>(), 20, "bootstrap weights must sum to n");
        }
    }

    #[test]
    fn test_bootstrap_weights_nonempty() {
        let mut rng = SimpleRng::seed(42);
        let weights = bootstrap_weights(10, &ResamplingMethod::Bootstrap, &mut rng);
        // At least one site must have been sampled
        assert!(weights.iter().any(|&w| w > 0));
    }

    #[test]
    fn test_bootstrap_empty_sites() {
        let mut rng = SimpleRng::seed(1);
        let weights = bootstrap_weights(0, &ResamplingMethod::Bootstrap, &mut rng);
        assert_eq!(weights.len(), 0);
    }

    #[test]
    fn test_block_bootstrap_weights_sum() {
        let mut rng = SimpleRng::seed(999);
        let method = ResamplingMethod::BlockBootstrap { block_size: 3 };
        for _ in 0..100 {
            let weights = bootstrap_weights(12, &method, &mut rng);
            assert_eq!(weights.len(), 12);
            // blocks = 12 / 3 = 4, each block adds 3, total = 12
            assert_eq!(weights.iter().sum::<usize>(), 12);
        }
    }

    #[test]
    fn test_block_bootstrap_contiguity() {
        // With block_size equal to nsites, we should get one contiguous block
        let mut rng = SimpleRng::seed(77);
        let method = ResamplingMethod::BlockBootstrap { block_size: 5 };
        let weights = bootstrap_weights(5, &method, &mut rng);
        // One block of 5 from nsites=5 -> all weights should be 1
        assert_eq!(weights.iter().sum::<usize>(), 5);
        for &w in &weights {
            assert_eq!(w, 1, "single full-size block should sample each site once");
        }
    }

    #[test]
    fn test_jackknife_weights_sum() {
        let mut rng = SimpleRng::seed(54321);
        let method = ResamplingMethod::Jackknife { fraction: 0.5 };
        for _ in 0..100 {
            let weights = bootstrap_weights(20, &method, &mut rng);
            assert_eq!(weights.len(), 20);
            let sum: usize = weights.iter().sum();
            assert_eq!(sum, 10, "delete-half jackknife should select exactly n/2 sites");
            // All weights are 0 or 1
            assert!(weights.iter().all(|&w| w <= 1));
        }
    }

    #[test]
    fn test_jackknife_fraction() {
        let mut rng = SimpleRng::seed(11111);
        let method = ResamplingMethod::Jackknife { fraction: 0.75 };
        for _ in 0..100 {
            let weights = bootstrap_weights(20, &method, &mut rng);
            let sum: usize = weights.iter().sum();
            assert_eq!(sum, 15, "75% jackknife of 20 sites should select 15");
        }
    }

    #[test]
    fn test_partial_bootstrap_weights_sum() {
        let mut rng = SimpleRng::seed(7777);
        let method = ResamplingMethod::PartialBootstrap { fraction: 0.5 };
        for _ in 0..100 {
            let weights = bootstrap_weights(20, &method, &mut rng);
            assert_eq!(weights.len(), 20);
            assert_eq!(weights.iter().sum::<usize>(), 10);
        }
    }

    #[test]
    fn test_resample_alignment_preserves_taxa() {
        let aln = make_test_alignment();
        let mut rng = SimpleRng::seed(42);
        let weights = bootstrap_weights(aln.nsites(), &ResamplingMethod::Bootstrap, &mut rng);
        let resampled = resample_alignment(&aln, &weights);

        assert_eq!(resampled.ntaxa(), 3);
        assert_eq!(resampled.taxon_name(0), "TaxonA");
        assert_eq!(resampled.taxon_name(1), "TaxonB");
        assert_eq!(resampled.taxon_name(2), "TaxonC");
    }

    #[test]
    fn test_resample_alignment_length() {
        let aln = make_test_alignment();
        let mut rng = SimpleRng::seed(42);
        let weights = bootstrap_weights(aln.nsites(), &ResamplingMethod::Bootstrap, &mut rng);
        let resampled = resample_alignment(&aln, &weights);

        assert_eq!(resampled.nsites(), aln.nsites(),
                   "resampled alignment should have same length as original");
    }

    #[test]
    fn test_resample_alignment_content() {
        // With known weights, verify exact content
        let aln = make_test_alignment();
        let weights = vec![2, 0, 1, 0, 0, 0, 1, 0]; // sites 0,0,2,6
        let resampled = resample_alignment(&aln, &weights);

        assert_eq!(resampled.nsites(), 4);
        // TaxonA: site0=A, site0=A, site2=G, site6=G -> AAGG
        assert_eq!(resampled.base(0, 0), Base::A);
        assert_eq!(resampled.base(0, 1), Base::A);
        assert_eq!(resampled.base(0, 2), Base::G);
        assert_eq!(resampled.base(0, 3), Base::G);
    }

    #[test]
    fn test_reproducibility() {
        let aln = make_test_alignment();
        let method = ResamplingMethod::Bootstrap;

        let mut rng1 = SimpleRng::seed(42);
        let w1 = bootstrap_weights(aln.nsites(), &method, &mut rng1);

        let mut rng2 = SimpleRng::seed(42);
        let w2 = bootstrap_weights(aln.nsites(), &method, &mut rng2);

        assert_eq!(w1, w2, "same seed should produce identical weights");
    }

    #[test]
    fn test_different_seeds_differ() {
        let method = ResamplingMethod::Bootstrap;

        let mut rng1 = SimpleRng::seed(1);
        let w1 = bootstrap_weights(100, &method, &mut rng1);

        let mut rng2 = SimpleRng::seed(2);
        let w2 = bootstrap_weights(100, &method, &mut rng2);

        assert_ne!(w1, w2, "different seeds should produce different weights");
    }

    #[test]
    fn test_bootstrap_replicates_count() {
        let aln = make_test_alignment();
        let method = ResamplingMethod::Bootstrap;
        let mut rng = SimpleRng::seed(42);
        let reps: Vec<_> = bootstrap_replicates(&aln, 10, &method, &mut rng).collect();
        assert_eq!(reps.len(), 10);
        for rep in &reps {
            assert_eq!(rep.ntaxa(), aln.ntaxa());
            assert_eq!(rep.nsites(), aln.nsites());
        }
    }

    #[test]
    fn test_bootstrap_statistical_properties() {
        // Over many replicates, each site should be sampled ~n times on average
        let nsites = 100;
        let nreps = 10000;
        let method = ResamplingMethod::Bootstrap;
        let mut rng = SimpleRng::seed(42);

        let mut total_weights = vec![0usize; nsites];
        for _ in 0..nreps {
            let weights = bootstrap_weights(nsites, &method, &mut rng);
            for (i, &w) in weights.iter().enumerate() {
                total_weights[i] += w;
            }
        }

        // Each site should be sampled ~nreps times on average
        let expected = nreps as f64;
        for &total in &total_weights {
            let ratio = total as f64 / expected;
            assert!(
                (0.9..1.1).contains(&ratio),
                "site sampled {} times, expected ~{} (ratio {})",
                total, expected, ratio
            );
        }
    }

    #[test]
    fn test_rng_basic() {
        let mut rng = SimpleRng::seed(0);
        // Should produce different values
        let a = rng.gen_range(1000);
        let b = rng.gen_range(1000);
        let c = rng.gen_range(1000);
        // Extremely unlikely all three are equal
        assert!(!(a == b && b == c), "RNG should produce varying values");
    }

    #[test]
    fn test_rng_range() {
        let mut rng = SimpleRng::seed(42);
        for _ in 0..1000 {
            let v = rng.gen_range(10);
            assert!(v < 10, "gen_range(10) produced {}", v);
        }
    }

    #[test]
    fn test_rng_f64_range() {
        let mut rng = SimpleRng::seed(42);
        for _ in 0..1000 {
            let v = rng.gen_f64();
            assert!((0.0..1.0).contains(&v), "gen_f64() produced {}", v);
        }
    }
}

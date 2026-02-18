//! Lake's and Cavender's phylogenetic invariants for 4-taxon trees.
//!
//! Phylogenetic invariants are functions of site-pattern frequencies that
//! equal zero under a particular tree topology but are non-zero under
//! alternative topologies. They can be used to infer which of the three
//! possible unrooted 4-taxon topologies best explains an alignment.
//!
//! # Lake's invariants
//!
//! Lake (1987) proposed using transversion-type site patterns to choose
//! among the three unrooted 4-taxon topologies. For taxa (1,2,3,4),
//! the three possible unrooted topologies are:
//!
//! - T1: (1,2)|(3,4)
//! - T2: (1,3)|(2,4)
//! - T3: (1,4)|(2,3)
//!
//! A site pattern is "informative" under Lake's method if exactly two
//! of the four taxa share one state and the other two share a different
//! state, and the two states are related by transversion (purine vs
//! pyrimidine). These "parsimony-informative transversion" patterns
//! support one of the three topologies.
//!
//! # Cavender's invariants
//!
//! Cavender (1978) showed that for the symmetric (CFN/Cavender-Farris-Neyman)
//! model, certain linear functions of site pattern frequencies vanish for
//! the true topology. For each topology, one can compute an invariant
//! value; the topology whose invariant is closest to zero is favored.
//!
//! # References
//!
//! - Lake, J.A. (1987). A rate-independent technique for analysis of
//!   nucleic acid sequences: evolutionary parsimony. *Molecular Biology
//!   and Evolution*, 4:167-191.
//! - Cavender, J.A. (1978). Taxonomy with confidence. *Mathematical
//!   Biosciences*, 40:271-280.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//!   Chapter 29: Invariants.

use crate::tree::Alignment;

/// Site pattern counts for 4-taxon analysis.
///
/// Counts of each of the 256 possible 4-taxon site patterns (4 bases
/// at 4 positions = 4^4 = 256 patterns).
#[derive(Debug, Clone)]
pub struct PatternCounts {
    /// Counts indexed by (b1*64 + b2*16 + b3*4 + b4) where bi in {0=A,1=C,2=G,3=T}.
    pub counts: [usize; 256],
    /// Total number of sites counted (those where all 4 taxa have definite bases).
    pub total_sites: usize,
}

impl PatternCounts {
    /// Create a new empty pattern count array.
    fn new() -> Self {
        Self {
            counts: [0; 256],
            total_sites: 0,
        }
    }

    /// Compute the pattern index for four base indices (0-3).
    fn index(b1: usize, b2: usize, b3: usize, b4: usize) -> usize {
        b1 * 64 + b2 * 16 + b3 * 4 + b4
    }

    /// Get count for a particular base combination.
    pub fn get(&self, b1: usize, b2: usize, b3: usize, b4: usize) -> usize {
        self.counts[Self::index(b1, b2, b3, b4)]
    }
}

/// Result of Lake's invariants analysis.
#[derive(Debug, Clone)]
pub struct InvariantResult {
    /// Support counts for each of the 3 unrooted topologies.
    /// Index 0: (1,2)|(3,4), Index 1: (1,3)|(2,4), Index 2: (1,4)|(2,3).
    pub topology_support: [usize; 3],
    /// Index of the best-supported topology (0, 1, or 2).
    pub best_topology_index: usize,
    /// Detailed site pattern frequencies.
    pub pattern_counts: PatternCounts,
    /// Chi-squared statistic testing whether the support values differ
    /// significantly from equal.
    pub chi_squared: f64,
}

/// Result of Cavender's invariants analysis.
#[derive(Debug, Clone)]
pub struct CavenderResult {
    /// Invariant values for each of the 3 unrooted topologies.
    /// The topology whose invariant is closest to zero is favored.
    pub invariant_values: [f64; 3],
    /// Index of the best-supported topology (0, 1, or 2).
    pub best_topology_index: usize,
}

/// Returns true if a base is a purine (A or G).
fn is_purine(b: usize) -> bool {
    b == 0 || b == 2 // A=0, G=2
}

/// Returns true if a base is a pyrimidine (C or T).
fn is_pyrimidine(b: usize) -> bool {
    b == 1 || b == 3 // C=1, T=3
}

/// Returns true if two bases differ by a transversion (purine <-> pyrimidine).
fn is_transversion(b1: usize, b2: usize) -> bool {
    (is_purine(b1) && is_pyrimidine(b2)) || (is_pyrimidine(b1) && is_purine(b2))
}

/// Count all 256 site patterns from an alignment of exactly 4 taxa.
///
/// Only sites where all 4 taxa have definite (unambiguous) bases are counted.
fn count_patterns(alignment: &Alignment) -> Result<PatternCounts, String> {
    if alignment.ntaxa() != 4 {
        return Err(format!(
            "Lake's invariants require exactly 4 taxa, got {}",
            alignment.ntaxa()
        ));
    }

    let mut pc = PatternCounts::new();
    let nsites = alignment.nsites();

    for site in 0..nsites {
        let b1 = alignment.base(0, site);
        let b2 = alignment.base(1, site);
        let b3 = alignment.base(2, site);
        let b4 = alignment.base(3, site);

        // Only count sites where all 4 taxa have definite bases
        if let (Some(i1), Some(i2), Some(i3), Some(i4)) =
            (b1.index(), b2.index(), b3.index(), b4.index())
        {
            pc.counts[PatternCounts::index(i1, i2, i3, i4)] += 1;
            pc.total_sites += 1;
        }
    }

    Ok(pc)
}

/// Compute Lake's phylogenetic invariants for a 4-taxon alignment.
///
/// For 4 taxa labeled 1, 2, 3, 4, the three unrooted topologies are:
/// - T1: (1,2)|(3,4) — taxa 1 and 2 are grouped
/// - T2: (1,3)|(2,4) — taxa 1 and 3 are grouped
/// - T3: (1,4)|(2,3) — taxa 1 and 4 are grouped
///
/// Lake's method counts "parsimony-informative transversion" patterns:
/// site patterns where exactly two taxa share one state and the other
/// two share a different state, and the two states are transversions
/// of each other (purine vs pyrimidine).
///
/// Pattern xxyy supports T1: (1,2)|(3,4)
/// Pattern xyxy supports T2: (1,3)|(2,4)
/// Pattern xyyx supports T3: (1,4)|(2,3)
///
/// The topology with the most supporting patterns is favored.
///
/// # Errors
///
/// Returns an error if the alignment does not have exactly 4 taxa.
pub fn lake_invariants(alignment: &Alignment) -> Result<InvariantResult, String> {
    let pc = count_patterns(alignment)?;

    let mut support = [0usize; 3];

    // Enumerate all 256 patterns and classify informative transversion patterns
    for b1 in 0..4 {
        for b2 in 0..4 {
            for b3 in 0..4 {
                for b4 in 0..4 {
                    let count = pc.get(b1, b2, b3, b4);
                    if count == 0 {
                        continue;
                    }

                    // Pattern xxyy: b1==b2, b3==b4, b1!=b3, transversion
                    if b1 == b2 && b3 == b4 && b1 != b3 && is_transversion(b1, b3) {
                        support[0] += count; // supports T1: (1,2)|(3,4)
                    }
                    // Pattern xyxy: b1==b3, b2==b4, b1!=b2, transversion
                    if b1 == b3 && b2 == b4 && b1 != b2 && is_transversion(b1, b2) {
                        support[1] += count; // supports T2: (1,3)|(2,4)
                    }
                    // Pattern xyyx: b1==b4, b2==b3, b1!=b2, transversion
                    if b1 == b4 && b2 == b3 && b1 != b2 && is_transversion(b1, b2) {
                        support[2] += count; // supports T3: (1,4)|(2,3)
                    }
                }
            }
        }
    }

    // Find best topology
    let best = if support[0] >= support[1] && support[0] >= support[2] {
        0
    } else if support[1] >= support[2] {
        1
    } else {
        2
    };

    // Compute chi-squared statistic: tests if support counts deviate from equal
    let total_informative: f64 = (support[0] + support[1] + support[2]) as f64;
    let chi_sq = if total_informative > 0.0 {
        let expected = total_informative / 3.0;
        support
            .iter()
            .map(|&obs| {
                let diff = obs as f64 - expected;
                diff * diff / expected
            })
            .sum()
    } else {
        0.0
    };

    Ok(InvariantResult {
        topology_support: support,
        best_topology_index: best,
        pattern_counts: pc,
        chi_squared: chi_sq,
    })
}

/// Compute Cavender's phylogenetic invariants for a 4-taxon alignment.
///
/// Under the symmetric (CFN/Cavender-Farris-Neyman) model, certain linear
/// functions of site pattern frequencies vanish for the true topology.
///
/// First, we reduce each nucleotide to R (purine: A,G) or Y (pyrimidine: C,T),
/// giving 2^4 = 16 binary patterns. Then we count the informative patterns:
///
/// - n1 = count(RRYY) + count(YYRR): patterns grouping taxa 1,2 vs 3,4
/// - n2 = count(RYRY) + count(YRYR): patterns grouping taxa 1,3 vs 2,4
/// - n3 = count(RYYR) + count(YRRY): patterns grouping taxa 1,4 vs 2,3
///
/// Under the CFN model:
/// - If T1: (1,2)|(3,4) is true, then n2 = n3 (invariant: n2 - n3 = 0)
/// - If T2: (1,3)|(2,4) is true, then n1 = n3 (invariant: n1 - n3 = 0)
/// - If T3: (1,4)|(2,3) is true, then n1 = n2 (invariant: n1 - n2 = 0)
///
/// The topology whose invariant is closest to zero is favored.
///
/// # Errors
///
/// Returns an error if the alignment does not have exactly 4 taxa.
pub fn cavender_invariants(alignment: &Alignment) -> Result<CavenderResult, String> {
    if alignment.ntaxa() != 4 {
        return Err(format!(
            "Cavender's invariants require exactly 4 taxa, got {}",
            alignment.ntaxa()
        ));
    }

    // Count RY patterns (reduce to binary: R=0, Y=1)
    // 2^4 = 16 patterns
    let mut ry_counts = [0usize; 16];
    let mut total = 0usize;

    let nsites = alignment.nsites();
    for site in 0..nsites {
        let bases: Vec<Option<usize>> = (0..4)
            .map(|t| alignment.base(t, site).index())
            .collect();

        // Skip if any taxon has an ambiguous base
        if bases.iter().any(|b| b.is_none()) {
            continue;
        }

        let indices: Vec<usize> = bases.iter().map(|b| b.unwrap()).collect();

        // Map to R=0, Y=1
        let ry: Vec<usize> = indices
            .iter()
            .map(|&b| if is_purine(b) { 0 } else { 1 })
            .collect();

        let idx = ry[0] * 8 + ry[1] * 4 + ry[2] * 2 + ry[3];
        ry_counts[idx] += 1;
        total += 1;
    }

    if total == 0 {
        return Ok(CavenderResult {
            invariant_values: [0.0; 3],
            best_topology_index: 0,
        });
    }

    // RY pattern indices:
    // R=0, Y=1 for each position (taxon)
    // RRYY = 0011 = idx 3,  YYRR = 1100 = idx 12
    // RYRY = 0101 = idx 5,  YRYR = 1010 = idx 10
    // RYYR = 0110 = idx 6,  YRRY = 1001 = idx 9

    // Informative pattern group counts:
    // n1: patterns grouping taxa 1,2 vs 3,4 (supports T1)
    let n1 = (ry_counts[3] + ry_counts[12]) as f64;
    // n2: patterns grouping taxa 1,3 vs 2,4 (supports T2)
    let n2 = (ry_counts[5] + ry_counts[10]) as f64;
    // n3: patterns grouping taxa 1,4 vs 2,3 (supports T3)
    let n3 = (ry_counts[6] + ry_counts[9]) as f64;

    // Cavender's linear invariants:
    // For T1 to be true: n2 - n3 = 0
    // For T2 to be true: n1 - n3 = 0
    // For T3 to be true: n1 - n2 = 0
    let inv1 = n2 - n3; // Should be zero if T1 is true
    let inv2 = n1 - n3; // Should be zero if T2 is true
    let inv3 = n1 - n2; // Should be zero if T3 is true

    let invariant_values = [inv1, inv2, inv3];

    // Find topology with invariant closest to zero
    let best = invariant_values
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);

    Ok(CavenderResult {
        invariant_values,
        best_topology_index: best,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::{Base, Sequence};

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

    // --- Pattern counting tests ---

    #[test]
    fn test_count_patterns_identical() {
        let aln = make_alignment(&[
            ("A", "AAAA"),
            ("B", "AAAA"),
            ("C", "AAAA"),
            ("D", "AAAA"),
        ]);
        let pc = count_patterns(&aln).unwrap();
        assert_eq!(pc.total_sites, 4);
        assert_eq!(pc.get(0, 0, 0, 0), 4); // AAAA
    }

    #[test]
    fn test_count_patterns_wrong_ntaxa() {
        let aln = make_alignment(&[
            ("A", "AAAA"),
            ("B", "AAAA"),
            ("C", "AAAA"),
        ]);
        assert!(count_patterns(&aln).is_err());
    }

    #[test]
    fn test_count_patterns_skips_ambiguous() {
        let aln = make_alignment(&[
            ("A", "AN"),
            ("B", "AA"),
            ("C", "AA"),
            ("D", "AA"),
        ]);
        let pc = count_patterns(&aln).unwrap();
        assert_eq!(pc.total_sites, 1); // Only site 0, site 1 has N
    }

    // --- Lake's invariants tests ---

    #[test]
    fn test_lake_clear_signal_t1() {
        // Create alignment with strong signal for T1: (1,2)|(3,4)
        // Use many xxyy transversion patterns: AACC, AATT, GGCC, GGTT, CCAA, CCGG, TTAA, TTGG
        let mut seq1 = String::new();
        let mut seq2 = String::new();
        let mut seq3 = String::new();
        let mut seq4 = String::new();

        // Add many xxyy transversion sites supporting T1
        for _ in 0..50 {
            seq1.push('A'); seq2.push('A'); seq3.push('C'); seq4.push('C');
            seq1.push('G'); seq2.push('G'); seq3.push('T'); seq4.push('T');
        }
        // Add some constant sites (noise)
        for _ in 0..20 {
            seq1.push('A'); seq2.push('A'); seq3.push('A'); seq4.push('A');
        }

        let aln = make_alignment(&[
            ("T1", &seq1),
            ("T2", &seq2),
            ("T3", &seq3),
            ("T4", &seq4),
        ]);

        let result = lake_invariants(&aln).unwrap();
        assert_eq!(result.best_topology_index, 0, "T1 should be favored");
        assert!(
            result.topology_support[0] > result.topology_support[1],
            "T1 support should exceed T2 support"
        );
        assert!(
            result.topology_support[0] > result.topology_support[2],
            "T1 support should exceed T3 support"
        );
    }

    #[test]
    fn test_lake_clear_signal_t2() {
        // Create alignment with strong signal for T2: (1,3)|(2,4)
        // Use xyxy transversion patterns
        let mut seq1 = String::new();
        let mut seq2 = String::new();
        let mut seq3 = String::new();
        let mut seq4 = String::new();

        for _ in 0..50 {
            seq1.push('A'); seq2.push('C'); seq3.push('A'); seq4.push('C');
            seq1.push('G'); seq2.push('T'); seq3.push('G'); seq4.push('T');
        }
        for _ in 0..20 {
            seq1.push('A'); seq2.push('A'); seq3.push('A'); seq4.push('A');
        }

        let aln = make_alignment(&[
            ("T1", &seq1),
            ("T2", &seq2),
            ("T3", &seq3),
            ("T4", &seq4),
        ]);

        let result = lake_invariants(&aln).unwrap();
        assert_eq!(result.best_topology_index, 1, "T2 should be favored");
    }

    #[test]
    fn test_lake_clear_signal_t3() {
        // Create alignment with strong signal for T3: (1,4)|(2,3)
        // Use xyyx transversion patterns
        let mut seq1 = String::new();
        let mut seq2 = String::new();
        let mut seq3 = String::new();
        let mut seq4 = String::new();

        for _ in 0..50 {
            seq1.push('A'); seq2.push('C'); seq3.push('C'); seq4.push('A');
            seq1.push('G'); seq2.push('T'); seq3.push('T'); seq4.push('G');
        }
        for _ in 0..20 {
            seq1.push('A'); seq2.push('A'); seq3.push('A'); seq4.push('A');
        }

        let aln = make_alignment(&[
            ("T1", &seq1),
            ("T2", &seq2),
            ("T3", &seq3),
            ("T4", &seq4),
        ]);

        let result = lake_invariants(&aln).unwrap();
        assert_eq!(result.best_topology_index, 2, "T3 should be favored");
    }

    #[test]
    fn test_lake_equal_signal() {
        // Create alignment with equal support for all topologies
        let mut seq1 = String::new();
        let mut seq2 = String::new();
        let mut seq3 = String::new();
        let mut seq4 = String::new();

        // xxyy supporting T1
        for _ in 0..30 {
            seq1.push('A'); seq2.push('A'); seq3.push('C'); seq4.push('C');
        }
        // xyxy supporting T2
        for _ in 0..30 {
            seq1.push('A'); seq2.push('C'); seq3.push('A'); seq4.push('C');
        }
        // xyyx supporting T3
        for _ in 0..30 {
            seq1.push('A'); seq2.push('C'); seq3.push('C'); seq4.push('A');
        }

        let aln = make_alignment(&[
            ("T1", &seq1),
            ("T2", &seq2),
            ("T3", &seq3),
            ("T4", &seq4),
        ]);

        let result = lake_invariants(&aln).unwrap();
        // All should be approximately equal
        assert_eq!(result.topology_support[0], 30);
        assert_eq!(result.topology_support[1], 30);
        assert_eq!(result.topology_support[2], 30);
        // Chi-squared should be 0 for equal support
        assert!(result.chi_squared.abs() < 1e-10);
    }

    #[test]
    fn test_lake_no_informative_sites() {
        // All constant sites -> no informative patterns
        let aln = make_alignment(&[
            ("A", "AAAAAAAAAA"),
            ("B", "AAAAAAAAAA"),
            ("C", "AAAAAAAAAA"),
            ("D", "AAAAAAAAAA"),
        ]);

        let result = lake_invariants(&aln).unwrap();
        assert_eq!(result.topology_support, [0, 0, 0]);
        assert!(result.chi_squared.abs() < 1e-10);
    }

    #[test]
    fn test_lake_wrong_ntaxa() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "ACGT"),
        ]);
        assert!(lake_invariants(&aln).is_err());
    }

    #[test]
    fn test_lake_transitions_not_counted() {
        // Transitions (A<->G, C<->T) should NOT be counted as informative
        // Pattern AAGG: xxyy but A->G is a transition, not transversion
        let aln = make_alignment(&[
            ("A", "AAAA"),
            ("B", "AAAA"),
            ("C", "GGGG"),
            ("D", "GGGG"),
        ]);

        let result = lake_invariants(&aln).unwrap();
        // These are xxyy patterns but A->G is transition, not transversion
        assert_eq!(
            result.topology_support,
            [0, 0, 0],
            "transition-only patterns should not be counted"
        );
    }

    #[test]
    fn test_lake_chi_squared_nonzero() {
        // Unequal support -> chi-squared > 0
        let mut seq1 = String::new();
        let mut seq2 = String::new();
        let mut seq3 = String::new();
        let mut seq4 = String::new();

        // 60 sites supporting T1, 20 T2, 10 T3
        for _ in 0..60 {
            seq1.push('A'); seq2.push('A'); seq3.push('C'); seq4.push('C');
        }
        for _ in 0..20 {
            seq1.push('A'); seq2.push('C'); seq3.push('A'); seq4.push('C');
        }
        for _ in 0..10 {
            seq1.push('A'); seq2.push('C'); seq3.push('C'); seq4.push('A');
        }

        let aln = make_alignment(&[
            ("T1", &seq1),
            ("T2", &seq2),
            ("T3", &seq3),
            ("T4", &seq4),
        ]);

        let result = lake_invariants(&aln).unwrap();
        assert!(result.chi_squared > 0.0, "chi-squared should be positive for unequal support");
        assert_eq!(result.best_topology_index, 0);
    }

    // --- Cavender's invariants tests ---

    #[test]
    fn test_cavender_clear_signal_t1() {
        // Data consistent with T1: (1,2)|(3,4)
        // Strong xxyy pattern using transversions
        let mut seq1 = String::new();
        let mut seq2 = String::new();
        let mut seq3 = String::new();
        let mut seq4 = String::new();

        // RRYY patterns: taxa 1,2 share purine, taxa 3,4 share pyrimidine
        for _ in 0..80 {
            seq1.push('A'); seq2.push('A'); seq3.push('C'); seq4.push('C');
        }
        // YYRR patterns: complement
        for _ in 0..80 {
            seq1.push('C'); seq2.push('C'); seq3.push('A'); seq4.push('A');
        }
        // Some constant sites
        for _ in 0..40 {
            seq1.push('A'); seq2.push('A'); seq3.push('A'); seq4.push('A');
        }

        let aln = make_alignment(&[
            ("T1", &seq1),
            ("T2", &seq2),
            ("T3", &seq3),
            ("T4", &seq4),
        ]);

        let result = cavender_invariants(&aln).unwrap();
        assert_eq!(
            result.best_topology_index, 0,
            "T1 should be favored, invariants = {:?}",
            result.invariant_values
        );
    }

    #[test]
    fn test_cavender_clear_signal_t2() {
        // Data consistent with T2: (1,3)|(2,4)
        let mut seq1 = String::new();
        let mut seq2 = String::new();
        let mut seq3 = String::new();
        let mut seq4 = String::new();

        // RYRY patterns
        for _ in 0..80 {
            seq1.push('A'); seq2.push('C'); seq3.push('A'); seq4.push('C');
        }
        // YRYR patterns
        for _ in 0..80 {
            seq1.push('C'); seq2.push('A'); seq3.push('C'); seq4.push('A');
        }
        // Constant sites
        for _ in 0..40 {
            seq1.push('A'); seq2.push('A'); seq3.push('A'); seq4.push('A');
        }

        let aln = make_alignment(&[
            ("T1", &seq1),
            ("T2", &seq2),
            ("T3", &seq3),
            ("T4", &seq4),
        ]);

        let result = cavender_invariants(&aln).unwrap();
        assert_eq!(
            result.best_topology_index, 1,
            "T2 should be favored, invariants = {:?}",
            result.invariant_values
        );
    }

    #[test]
    fn test_cavender_clear_signal_t3() {
        // Data consistent with T3: (1,4)|(2,3)
        let mut seq1 = String::new();
        let mut seq2 = String::new();
        let mut seq3 = String::new();
        let mut seq4 = String::new();

        // RYYR patterns
        for _ in 0..80 {
            seq1.push('A'); seq2.push('C'); seq3.push('C'); seq4.push('A');
        }
        // YRRY patterns
        for _ in 0..80 {
            seq1.push('C'); seq2.push('A'); seq3.push('A'); seq4.push('C');
        }
        // Constant sites
        for _ in 0..40 {
            seq1.push('A'); seq2.push('A'); seq3.push('A'); seq4.push('A');
        }

        let aln = make_alignment(&[
            ("T1", &seq1),
            ("T2", &seq2),
            ("T3", &seq3),
            ("T4", &seq4),
        ]);

        let result = cavender_invariants(&aln).unwrap();
        assert_eq!(
            result.best_topology_index, 2,
            "T3 should be favored, invariants = {:?}",
            result.invariant_values
        );
    }

    #[test]
    fn test_cavender_wrong_ntaxa() {
        let aln = make_alignment(&[
            ("A", "ACGT"),
            ("B", "ACGT"),
            ("C", "ACGT"),
            ("D", "ACGT"),
            ("E", "ACGT"),
        ]);
        assert!(cavender_invariants(&aln).is_err());
    }

    #[test]
    fn test_cavender_all_constant() {
        let aln = make_alignment(&[
            ("A", "AAAAAAAAAA"),
            ("B", "AAAAAAAAAA"),
            ("C", "AAAAAAAAAA"),
            ("D", "AAAAAAAAAA"),
        ]);

        let result = cavender_invariants(&aln).unwrap();
        // All invariants should be zero for constant data
        for val in &result.invariant_values {
            assert!(
                val.abs() < 1e-10,
                "invariant should be zero for constant data, got {}",
                val
            );
        }
    }

    #[test]
    fn test_cavender_invariant_relationship() {
        // The three Cavender invariants satisfy:
        // inv1 - inv2 + inv3 = (n2-n3) - (n1-n3) + (n1-n2) = 0
        // i.e., inv1 + inv3 = inv2 always
        let aln = make_alignment(&[
            ("A", "ACGTACGTACGTACGTACGT"),
            ("B", "AACCGGTTAACCGGTTAACC"),
            ("C", "ACACACACACGTGTGTACAC"),
            ("D", "AACCAACCGGTTGGTTAACC"),
        ]);

        let result = cavender_invariants(&aln).unwrap();
        let check = result.invariant_values[0] - result.invariant_values[1]
            + result.invariant_values[2];
        assert!(
            check.abs() < 1e-10,
            "inv1 - inv2 + inv3 should be zero, got {} (values: {:?})",
            check,
            result.invariant_values
        );
    }

    #[test]
    fn test_lake_and_cavender_agree() {
        // Both methods should identify the same best topology when signal is strong
        let mut seq1 = String::new();
        let mut seq2 = String::new();
        let mut seq3 = String::new();
        let mut seq4 = String::new();

        for _ in 0..100 {
            seq1.push('A'); seq2.push('A'); seq3.push('C'); seq4.push('C');
            seq1.push('G'); seq2.push('G'); seq3.push('T'); seq4.push('T');
        }
        for _ in 0..50 {
            seq1.push('A'); seq2.push('A'); seq3.push('A'); seq4.push('A');
        }

        let aln = make_alignment(&[
            ("T1", &seq1),
            ("T2", &seq2),
            ("T3", &seq3),
            ("T4", &seq4),
        ]);

        let lake_result = lake_invariants(&aln).unwrap();
        let cav_result = cavender_invariants(&aln).unwrap();

        assert_eq!(
            lake_result.best_topology_index,
            cav_result.best_topology_index,
            "Lake and Cavender should agree on strong signal (Lake={}, Cavender={})",
            lake_result.best_topology_index,
            cav_result.best_topology_index
        );
    }
}

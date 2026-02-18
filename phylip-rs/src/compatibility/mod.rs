//! Compatibility analysis and clique finding for phylogenetic characters.
//!
//! This module implements the Le Quesne (1969) compatibility method for
//! binary characters and the PHYLIP dnacomp-style DNA compatibility analysis.
//!
//! Two binary characters are **compatible** if they can both be perfectly
//! explained on a single tree without homoplasy. The compatibility test:
//! form the 2x2 table of presence/absence patterns; characters are compatible
//! if and only if at most 3 of the 4 possible combinations are observed.
//!
//! The maximum clique of mutually compatible characters identifies the
//! largest set of characters that can all evolve on a single tree without
//! any homoplasy. Finding this clique is NP-hard in general, but the
//! Bron-Kerbosch algorithm with pivoting handles typical character counts
//! efficiently.
//!
//! # References
//!
//! - Le Quesne, W.J. (1969). A method of selection of characters in
//!   numerical taxonomy. *Systematic Zoology*, 18, 201-205.
//! - Felsenstein, J. (1979). Alternative methods of phylogenetic inference
//!   and their interrelationship. *Systematic Zoology*, 28, 49-62.
//! - Bron, C. & Kerbosch, J. (1973). Algorithm 457: finding all cliques
//!   of an undirected graph. *Communications of the ACM*, 16, 575-577.

pub mod clique;
pub mod dna_compat;

use crate::tree::types::Alignment;

/// Binary character data matrix.
///
/// Each taxon is scored for each character as either ancestral (false/0) or
/// derived (true/1). This representation is used for the Le Quesne
/// compatibility test and maximum clique analysis.
#[derive(Debug, Clone)]
pub struct BinaryMatrix {
    /// Taxon names.
    pub taxa: Vec<String>,
    /// characters[taxon][character] = true (derived/1) or false (ancestral/0).
    pub characters: Vec<Vec<bool>>,
    /// Number of characters.
    pub n_chars: usize,
}

impl BinaryMatrix {
    /// Create a new binary matrix from taxon names and character data.
    ///
    /// Returns an error if the dimensions are inconsistent.
    pub fn new(taxa: Vec<String>, characters: Vec<Vec<bool>>) -> Result<Self, String> {
        if taxa.is_empty() {
            return Err("no taxa provided".to_string());
        }
        if characters.is_empty() {
            return Err("no character data provided".to_string());
        }
        if taxa.len() != characters.len() {
            return Err(format!(
                "number of taxa ({}) does not match number of character rows ({})",
                taxa.len(),
                characters.len()
            ));
        }
        let n_chars = characters[0].len();
        if n_chars == 0 {
            return Err("no characters provided".to_string());
        }
        for (i, row) in characters.iter().enumerate() {
            if row.len() != n_chars {
                return Err(format!(
                    "taxon {} has {} characters, expected {}",
                    taxa[i],
                    row.len(),
                    n_chars
                ));
            }
        }
        Ok(Self {
            taxa,
            characters,
            n_chars,
        })
    }

    /// Number of taxa.
    pub fn n_taxa(&self) -> usize {
        self.taxa.len()
    }

    /// Create a binary matrix from an Alignment by treating each site
    /// as a binary character.
    ///
    /// Encoding: for each site, the most common unambiguous base is treated
    /// as ancestral (0) and any other unambiguous base is derived (1).
    /// Ambiguous bases and gaps are treated as ancestral (0).
    pub fn from_alignment(alignment: &Alignment) -> Self {
        let ntaxa = alignment.ntaxa();
        let nsites = alignment.nsites();
        let taxa: Vec<String> = (0..ntaxa)
            .map(|i| alignment.taxon_name(i).to_string())
            .collect();

        let mut characters = vec![vec![false; nsites]; ntaxa];

        for site in 0..nsites {
            // Count occurrences of each definite base at this site
            let mut counts = [0usize; 4]; // A, C, G, T
            for taxon in 0..ntaxa {
                let base = alignment.base(taxon, site);
                if let Some(idx) = base.index() {
                    counts[idx] += 1;
                }
            }

            // Find the most common base (ancestral state)
            let max_count = counts.iter().copied().max().unwrap_or(0);
            if max_count == 0 {
                // All ambiguous or gap; leave all as false
                continue;
            }
            // Pick the first base with the max count as ancestral
            let ancestral_idx = counts.iter().position(|&c| c == max_count).unwrap();

            // Encode: any definite base that differs from ancestral is derived (1)
            for taxon in 0..ntaxa {
                let base = alignment.base(taxon, site);
                if let Some(idx) = base.index() {
                    if idx != ancestral_idx {
                        characters[taxon][site] = true;
                    }
                }
                // Ambiguous/gap: leave as false (ancestral)
            }
        }

        Self {
            taxa,
            characters,
            n_chars: nsites,
        }
    }
}

/// Test whether two binary characters are compatible.
///
/// Two binary (0/1) characters are compatible if and only if at most 3
/// of the 4 possible taxon state combinations (0,0), (0,1), (1,0), (1,1)
/// are observed. If all 4 are present, the characters cannot both evolve
/// on any tree without homoplasy.
pub fn are_compatible(
    matrix: &BinaryMatrix,
    char_i: usize,
    char_j: usize,
) -> bool {
    let mut seen = [false; 4]; // (0,0), (0,1), (1,0), (1,1)

    for taxon in 0..matrix.n_taxa() {
        let a = matrix.characters[taxon][char_i];
        let b = matrix.characters[taxon][char_j];
        let idx = (a as usize) * 2 + (b as usize);
        seen[idx] = true;
    }

    let count = seen.iter().filter(|&&s| s).count();
    count <= 3
}

/// Build the compatibility graph: nodes = characters, edges = compatible pairs.
///
/// Returns an adjacency matrix where `graph[i][j]` is true if characters
/// i and j are compatible.
pub fn compatibility_graph(matrix: &BinaryMatrix) -> Vec<Vec<bool>> {
    let n = matrix.n_chars;
    let mut graph = vec![vec![false; n]; n];

    for i in 0..n {
        // Note: diagonal is false (vertices are not their own neighbors
        // in graph-theoretic sense, which is needed for clique finding).
        for j in (i + 1)..n {
            let compat = are_compatible(matrix, i, j);
            graph[i][j] = compat;
            graph[j][i] = compat;
        }
    }

    graph
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::types::{Alignment, Base, Sequence};

    fn make_binary_matrix(
        taxa: &[&str],
        data: &[Vec<bool>],
    ) -> BinaryMatrix {
        let taxa_vec: Vec<String> = taxa.iter().map(|s| s.to_string()).collect();
        BinaryMatrix::new(taxa_vec, data.to_vec()).unwrap()
    }

    // --- BinaryMatrix tests ---

    #[test]
    fn test_binary_matrix_new_valid() {
        let m = BinaryMatrix::new(
            vec!["A".into(), "B".into(), "C".into()],
            vec![
                vec![true, false],
                vec![false, true],
                vec![true, true],
            ],
        );
        assert!(m.is_ok());
        let m = m.unwrap();
        assert_eq!(m.n_taxa(), 3);
        assert_eq!(m.n_chars, 2);
    }

    #[test]
    fn test_binary_matrix_new_empty_taxa() {
        let m = BinaryMatrix::new(vec![], vec![]);
        assert!(m.is_err());
    }

    #[test]
    fn test_binary_matrix_new_mismatched_rows() {
        let m = BinaryMatrix::new(
            vec!["A".into(), "B".into()],
            vec![vec![true, false]],
        );
        assert!(m.is_err());
    }

    #[test]
    fn test_binary_matrix_new_unequal_chars() {
        let m = BinaryMatrix::new(
            vec!["A".into(), "B".into()],
            vec![vec![true, false], vec![true]],
        );
        assert!(m.is_err());
    }

    #[test]
    fn test_binary_matrix_from_alignment() {
        let aln = Alignment::new(vec![
            Sequence::new("A", vec![Base::A, Base::C]),
            Sequence::new("B", vec![Base::A, Base::T]),
            Sequence::new("C", vec![Base::T, Base::C]),
        ])
        .unwrap();

        let m = BinaryMatrix::from_alignment(&aln);
        assert_eq!(m.n_taxa(), 3);
        assert_eq!(m.n_chars, 2);
        // Site 0: A appears twice (ancestral), T is derived
        // A: false, B: false, C: true
        assert!(!m.characters[0][0]); // A has 'A' (ancestral)
        assert!(!m.characters[1][0]); // B has 'A' (ancestral)
        assert!(m.characters[2][0]); // C has 'T' (derived)
    }

    // --- Compatibility tests ---

    #[test]
    fn test_incompatible_all_four_patterns() {
        // Two characters with all 4 patterns (0,0), (0,1), (1,0), (1,1): incompatible
        let m = make_binary_matrix(
            &["A", "B", "C", "D"],
            &[
                vec![false, false], // (0,0)
                vec![false, true],  // (0,1)
                vec![true, false],  // (1,0)
                vec![true, true],   // (1,1)
            ],
        );
        assert!(!are_compatible(&m, 0, 1));
    }

    #[test]
    fn test_compatible_three_patterns() {
        // Only 3 of 4 patterns -> compatible
        let m = make_binary_matrix(
            &["A", "B", "C"],
            &[
                vec![false, false], // (0,0)
                vec![false, true],  // (0,1)
                vec![true, true],   // (1,1)
            ],
        );
        assert!(are_compatible(&m, 0, 1));
    }

    #[test]
    fn test_compatible_identical_characters() {
        // Identical characters: only (0,0) and (1,1) observed
        let m = make_binary_matrix(
            &["A", "B", "C"],
            &[
                vec![false, false],
                vec![true, true],
                vec![false, false],
            ],
        );
        assert!(are_compatible(&m, 0, 1));
    }

    #[test]
    fn test_compatible_constant_character() {
        // Constant character (all 0): compatible with everything
        let m = make_binary_matrix(
            &["A", "B", "C"],
            &[
                vec![false, true],
                vec![false, false],
                vec![false, true],
            ],
        );
        // Char 0 is constant (all false), char 1 varies.
        // Patterns: (0,1), (0,0), (0,1) -> only 2 patterns -> compatible
        assert!(are_compatible(&m, 0, 1));
    }

    #[test]
    fn test_compatible_constant_all_one() {
        // Constant character (all 1): compatible with everything
        let m = make_binary_matrix(
            &["A", "B", "C"],
            &[
                vec![true, true],
                vec![true, false],
                vec![true, true],
            ],
        );
        // Char 0 is constant (all true), char 1 varies.
        // Patterns: (1,1), (1,0), (1,1) -> 2 patterns -> compatible
        assert!(are_compatible(&m, 0, 1));
    }

    #[test]
    fn test_compatible_complementary_characters() {
        // One character is the complement of the other
        // Patterns: (0,1) and (1,0) -> only 2 patterns -> compatible
        let m = make_binary_matrix(
            &["A", "B", "C", "D"],
            &[
                vec![false, true],
                vec![false, true],
                vec![true, false],
                vec![true, false],
            ],
        );
        assert!(are_compatible(&m, 0, 1));
    }

    // --- Compatibility graph tests ---

    #[test]
    fn test_compatibility_graph_all_compatible() {
        // Three identical characters -> all pairwise compatible
        let m = make_binary_matrix(
            &["A", "B", "C"],
            &[
                vec![false, false, false],
                vec![true, true, true],
                vec![false, false, false],
            ],
        );
        let graph = compatibility_graph(&m);
        for i in 0..3 {
            for j in 0..3 {
                if i != j {
                    assert!(graph[i][j], "characters {} and {} should be compatible", i, j);
                } else {
                    // Diagonal is false (no self-loops in graph)
                    assert!(!graph[i][j]);
                }
            }
        }
    }

    #[test]
    fn test_compatibility_graph_incompatible_pair() {
        let m = make_binary_matrix(
            &["A", "B", "C", "D"],
            &[
                vec![false, false],
                vec![false, true],
                vec![true, false],
                vec![true, true],
            ],
        );
        let graph = compatibility_graph(&m);
        // Diagonal is false (graph-theoretic convention: no self-loops)
        assert!(!graph[0][0]);
        assert!(!graph[1][1]);
        // These two characters are incompatible (all 4 patterns)
        assert!(!graph[0][1]);
        assert!(!graph[1][0]);
    }

    #[test]
    fn test_compatibility_graph_mixed() {
        // 3 characters: 0 and 1 compatible, 0 and 2 compatible, 1 and 2 incompatible
        let m = make_binary_matrix(
            &["A", "B", "C", "D"],
            &[
                // char0, char1, char2
                vec![false, false, false],
                vec![true, true, false],
                vec![true, true, true],
                vec![false, false, true],
            ],
        );
        let graph = compatibility_graph(&m);
        // Char 0 and 1: patterns (0,0), (1,1), (1,1), (0,0) -> 2 patterns -> compatible
        assert!(graph[0][1]);
        // Char 0 and 2: patterns (0,0), (1,0), (1,1), (0,1) -> 4 patterns -> incompatible
        assert!(!graph[0][2]);
        // Char 1 and 2: patterns (0,0), (1,0), (1,1), (0,1) -> 4 patterns -> incompatible
        assert!(!graph[1][2]);
    }
}

//! Maximum clique finding via Bron-Kerbosch algorithm with pivoting.
//!
//! Given a compatibility graph (from [`super::compatibility_graph`]), finds
//! the largest clique of mutually compatible binary characters. This is the
//! set of characters that can all be simultaneously explained on a single
//! tree without any homoplasy.
//!
//! The maximum clique problem is NP-hard in general, but the Bron-Kerbosch
//! algorithm with pivot selection provides efficient pruning for the
//! typical sizes encountered in phylogenetic character analysis.
//!
//! # References
//!
//! - Bron, C. & Kerbosch, J. (1973). Algorithm 457: finding all cliques
//!   of an undirected graph. *Communications of the ACM*, 16, 575-577.
//! - Tomita, E., Tanaka, A. & Takahashi, H. (2006). The worst-case time
//!   complexity for generating all maximal cliques and computational
//!   experiments. *Theoretical Computer Science*, 363, 28-42.

use crate::tree::types::Tree;
use super::BinaryMatrix;

/// Result of a clique analysis.
#[derive(Debug, Clone)]
pub struct CliqueResult {
    /// Indices of characters in the largest clique.
    pub clique_characters: Vec<usize>,
    /// Size of the largest clique.
    pub clique_size: usize,
    /// The tree supported by the clique characters (built from compatible
    /// splits implied by the clique characters).
    pub tree: Option<Tree>,
    /// Total number of maximal cliques found during the search.
    pub total_cliques: usize,
}

/// Find the largest clique of mutually compatible binary characters.
///
/// Builds the compatibility graph from the binary matrix and then
/// finds the maximum clique using the Bron-Kerbosch algorithm with
/// pivoting.
///
/// The clique characters represent the largest set of characters that
/// can all evolve on a single tree without any homoplasy.
pub fn find_max_clique(matrix: &BinaryMatrix) -> CliqueResult {
    let adj = super::compatibility_graph(matrix);
    let n = matrix.n_chars;

    let mut r: Vec<usize> = Vec::new();
    let mut p: Vec<usize> = (0..n).collect();
    let mut x: Vec<usize> = Vec::new();
    let mut best: Vec<usize> = Vec::new();
    let mut total_cliques: usize = 0;

    bron_kerbosch(&adj, &mut r, &mut p, &mut x, &mut best, &mut total_cliques);

    best.sort();

    let clique_size = best.len();
    let tree = if clique_size >= 3 {
        Some(build_tree_from_clique(matrix, &best))
    } else {
        None
    };

    CliqueResult {
        clique_characters: best,
        clique_size,
        tree,
        total_cliques,
    }
}

/// Bron-Kerbosch algorithm with pivoting for maximum clique finding.
///
/// - `adj` — adjacency matrix (compatibility graph)
/// - `r` — current clique being built
/// - `p` — candidate vertices that can extend the clique
/// - `x` — excluded vertices (already processed)
/// - `best` — best (largest) clique found so far
/// - `total_cliques` — counter for the total number of maximal cliques found
///
/// The pivot is chosen as the vertex in P union X with the most neighbors
/// in P, to maximize pruning.
fn bron_kerbosch(
    adj: &[Vec<bool>],
    r: &mut Vec<usize>,
    p: &mut Vec<usize>,
    x: &mut Vec<usize>,
    best: &mut Vec<usize>,
    total_cliques: &mut usize,
) {
    if p.is_empty() && x.is_empty() {
        // r is a maximal clique
        *total_cliques += 1;
        if r.len() > best.len() {
            *best = r.clone();
        }
        return;
    }

    if p.is_empty() {
        return;
    }

    // Pivot selection: choose vertex u in P union X with most neighbors in P
    let pivot = choose_pivot(adj, p, x);

    // Candidates to explore: vertices in P that are NOT neighbors of the pivot
    let candidates: Vec<usize> = p
        .iter()
        .filter(|&&v| !adj[pivot][v])
        .copied()
        .collect();

    for v in candidates {
        // Add v to R
        r.push(v);

        // New P: intersection of P with neighbors of v
        let mut new_p: Vec<usize> = p
            .iter()
            .filter(|&&u| adj[v][u])
            .copied()
            .collect();

        // New X: intersection of X with neighbors of v
        let mut new_x: Vec<usize> = x
            .iter()
            .filter(|&&u| adj[v][u])
            .copied()
            .collect();

        bron_kerbosch(adj, r, &mut new_p, &mut new_x, best, total_cliques);

        // Remove v from R
        r.pop();

        // Move v from P to X
        p.retain(|&u| u != v);
        x.push(v);
    }
}

/// Choose the pivot vertex from P union X that has the most neighbors in P.
///
/// This maximizes pruning by reducing the number of recursive calls.
fn choose_pivot(
    adj: &[Vec<bool>],
    p: &[usize],
    x: &[usize],
) -> usize {
    let mut best_pivot = p[0];
    let mut best_count = 0;

    // Check vertices in P
    for &u in p {
        let count = p.iter().filter(|&&v| adj[u][v]).count();
        if count > best_count {
            best_count = count;
            best_pivot = u;
        }
    }

    // Check vertices in X
    for &u in x {
        let count = p.iter().filter(|&&v| adj[u][v]).count();
        if count > best_count {
            best_count = count;
            best_pivot = u;
        }
    }

    best_pivot
}

/// Build a tree from a set of clique characters.
///
/// Each binary character in the clique induces a bipartition (split) of
/// the taxa. Since all characters in the clique are mutually compatible,
/// these splits can be represented on a single tree.
///
/// We build the tree by iteratively refining partitions of the taxon set
/// based on the character state patterns.
fn build_tree_from_clique(matrix: &BinaryMatrix, clique: &[usize]) -> Tree {
    let ntaxa = matrix.n_taxa();

    // Collect the splits induced by each clique character.
    // Each character splits taxa into two groups: those with state 0
    // and those with state 1.
    let mut splits: Vec<Vec<bool>> = Vec::new();
    for &ch in clique {
        let split: Vec<bool> = (0..ntaxa)
            .map(|t| matrix.characters[t][ch])
            .collect();
        // Only add non-trivial splits (at least one 0 and one 1)
        let ones = split.iter().filter(|&&s| s).count();
        if ones > 0 && ones < ntaxa {
            splits.push(split);
        }
    }

    // Build tree by hierarchical clustering from compatible splits.
    // Start with each taxon in its own group and merge according to splits.
    build_tree_from_splits(&matrix.taxa, &splits)
}

/// Build a tree from a set of compatible binary splits using
/// a hierarchical refinement approach.
///
/// Each split divides the taxa into two groups. Compatible splits can
/// be nested or disjoint but never partially overlapping. We build the
/// tree by sorting splits by their size (number of taxa on the "1" side)
/// and inserting them one at a time.
fn build_tree_from_splits(taxa: &[String], splits: &[Vec<bool>]) -> Tree {
    let ntaxa = taxa.len();
    let mut tree = Tree::new();
    tree.rooted = false;

    // Create root
    let root = tree.add_node(None, None, None);
    tree.root = root;

    if splits.is_empty() {
        // Star tree: all taxa as children of root
        for name in taxa {
            tree.add_node(Some(name.clone()), None, Some(root));
        }
        return tree;
    }

    // Sort splits by number of "1" taxa (ascending) so we process
    // smaller groups first.
    let mut indexed_splits: Vec<(usize, &Vec<bool>)> =
        splits.iter().enumerate().collect();
    indexed_splits.sort_by_key(|(_, s)| s.iter().filter(|&&b| b).count());

    // group_of[taxon] = node_id of the internal node containing it
    let mut group_of = vec![root; ntaxa];

    for (_, split) in &indexed_splits {
        // Find the set of taxa on the "1" side
        let ones: Vec<usize> = (0..ntaxa).filter(|&t| split[t]).collect();
        if ones.is_empty() || ones.len() == ntaxa {
            continue;
        }

        // Find the parent node(s) for these taxa. All "1" taxa should
        // currently be assigned to the same node (since splits are compatible).
        let parent_node = group_of[ones[0]];

        // Check if all "1" taxa share the same current group
        let all_same = ones.iter().all(|&t| group_of[t] == parent_node);
        if !all_same {
            // Splits might not nest perfectly with our insertion order;
            // skip this split.
            continue;
        }

        // Check if the "1" taxa are a strict subset of the taxa at this node
        let taxa_at_parent: Vec<usize> = (0..ntaxa)
            .filter(|&t| group_of[t] == parent_node)
            .collect();
        if ones.len() == taxa_at_parent.len() {
            // The split doesn't subdivide this group further
            continue;
        }

        // Create a new internal node under parent_node for the "1" group
        let new_node = tree.add_node(None, None, Some(parent_node));

        // Update group assignments
        for &t in &ones {
            group_of[t] = new_node;
        }
    }

    // Add leaf nodes: each taxon becomes a child of its group node
    for (t, name) in taxa.iter().enumerate() {
        tree.add_node(Some(name.clone()), None, Some(group_of[t]));
    }

    tree
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compatibility::{are_compatible, compatibility_graph};

    fn make_binary_matrix(
        taxa: &[&str],
        data: &[Vec<bool>],
    ) -> BinaryMatrix {
        let taxa_vec: Vec<String> = taxa.iter().map(|s| s.to_string()).collect();
        BinaryMatrix::new(taxa_vec, data.to_vec()).unwrap()
    }

    // --- Clique tests ---

    #[test]
    fn test_all_compatible_clique_is_all() {
        // 3 identical characters -> all compatible -> clique = all
        let m = make_binary_matrix(
            &["A", "B", "C", "D"],
            &[
                vec![false, false, false],
                vec![true, true, true],
                vec![false, false, false],
                vec![true, true, true],
            ],
        );
        let result = find_max_clique(&m);
        assert_eq!(result.clique_size, 3);
        assert_eq!(result.clique_characters, vec![0, 1, 2]);
    }

    #[test]
    fn test_no_compatible_pairs_clique_size_1() {
        // Each pair of characters has all 4 patterns -> no two compatible
        // Need at least 4 taxa for all 4 patterns
        let m = make_binary_matrix(
            &["A", "B", "C", "D"],
            &[
                vec![false, false, false],
                vec![false, true, false],
                vec![true, false, true],
                vec![true, true, true],
            ],
        );
        // Check that chars 0 and 1 are incompatible
        assert!(!are_compatible(&m, 0, 1));
        // Chars 0 and 2: patterns (0,0), (0,0), (1,1), (1,1) -> 2 patterns -> compatible
        // Actually let me verify:
        // taxon A: (0,0), B: (0,0), C: (1,1), D: (1,1) -> only 2 patterns -> compatible
        // So chars 0 and 2 ARE compatible. Let me make truly all-incompatible.
        let m2 = make_binary_matrix(
            &["A", "B", "C", "D"],
            &[
                // Char 0: A=0, B=0, C=1, D=1
                // Char 1: A=0, B=1, C=0, D=1
                // Char 2: A=0, B=1, C=1, D=0
                vec![false, false, false],
                vec![false, true, true],
                vec![true, false, true],
                vec![true, true, false],
            ],
        );
        // 0 vs 1: (0,0),(0,1),(1,0),(1,1) -> 4 patterns -> incompatible
        assert!(!are_compatible(&m2, 0, 1));
        // 0 vs 2: (0,0),(0,1),(1,1),(1,0) -> 4 patterns -> incompatible
        assert!(!are_compatible(&m2, 0, 2));
        // 1 vs 2: (0,0),(1,1),(0,1),(1,0) -> 4 patterns -> incompatible
        assert!(!are_compatible(&m2, 1, 2));

        let result = find_max_clique(&m2);
        assert_eq!(result.clique_size, 1);
    }

    #[test]
    fn test_five_chars_known_clique() {
        // 5 characters with a known compatibility structure.
        // Characters 0,1,2 are pairwise compatible (form a clique of size 3).
        // Character 3 is incompatible with 0 but compatible with 1,2.
        // Character 4 is incompatible with 1 but compatible with 0,2.
        //
        // So the largest clique could be {0,1,2} or {1,2,3} or {0,2,4}.
        // All have size 3.
        let m = make_binary_matrix(
            &["A", "B", "C", "D", "E"],
            &[
                // c0  c1  c2  c3  c4
                vec![false, false, false, false, false], // A
                vec![true,  true,  true,  false, true],  // B
                vec![true,  true,  true,  true,  true],  // C
                vec![false, false, false, true,  false],  // D
                vec![false, true,  false, true,  false],  // E: makes c0-c3 incompatible, c1-c4 incompatible
            ],
        );

        let result = find_max_clique(&m);
        assert!(result.clique_size >= 3);

        // Verify that the clique is actually a clique (all pairs compatible)
        let graph = compatibility_graph(&m);
        for i in 0..result.clique_characters.len() {
            for j in (i + 1)..result.clique_characters.len() {
                let ci = result.clique_characters[i];
                let cj = result.clique_characters[j];
                assert!(
                    graph[ci][cj],
                    "clique chars {} and {} should be compatible",
                    ci, cj
                );
            }
        }

        // Verify maximality: no additional character can be added
        let non_clique: Vec<usize> = (0..m.n_chars)
            .filter(|c| !result.clique_characters.contains(c))
            .collect();
        for &c in &non_clique {
            let all_compat = result
                .clique_characters
                .iter()
                .all(|&cc| graph[c][cc]);
            // If c is compatible with all clique members, the clique
            // isn't maximal. This should NOT happen.
            assert!(
                !all_compat,
                "character {} is compatible with all clique members, clique is not maximal",
                c
            );
        }
    }

    #[test]
    fn test_clique_single_character() {
        let m = make_binary_matrix(
            &["A", "B", "C"],
            &[
                vec![false],
                vec![true],
                vec![false],
            ],
        );
        let result = find_max_clique(&m);
        assert_eq!(result.clique_size, 1);
        assert_eq!(result.clique_characters, vec![0]);
    }

    #[test]
    fn test_clique_two_compatible_chars() {
        let m = make_binary_matrix(
            &["A", "B", "C"],
            &[
                vec![false, false],
                vec![true, true],
                vec![false, false],
            ],
        );
        let result = find_max_clique(&m);
        assert_eq!(result.clique_size, 2);
    }

    #[test]
    fn test_clique_builds_tree_for_large_clique() {
        // 4 taxa, 3 compatible characters -> should build a tree
        let m = make_binary_matrix(
            &["A", "B", "C", "D"],
            &[
                // c0: A=0,B=0,C=1,D=1 (split {C,D})
                // c1: A=0,B=0,C=1,D=1 (same split)
                // c2: A=0,B=0,C=1,D=1 (same split)
                vec![false, false, false],
                vec![false, false, false],
                vec![true, true, true],
                vec![true, true, true],
            ],
        );
        let result = find_max_clique(&m);
        assert_eq!(result.clique_size, 3);
        assert!(result.tree.is_some());
        let tree = result.tree.unwrap();
        assert_eq!(tree.num_leaves(), 4);
    }

    #[test]
    fn test_clique_total_cliques_counted() {
        let m = make_binary_matrix(
            &["A", "B", "C", "D"],
            &[
                vec![false, false],
                vec![false, true],
                vec![true, false],
                vec![true, true],
            ],
        );
        // chars 0 and 1 are incompatible -> each forms its own maximal clique
        let result = find_max_clique(&m);
        assert!(result.total_cliques >= 2);
    }

    #[test]
    fn test_clique_with_nested_splits() {
        // Characters that induce nested splits:
        // c0: {C,D} vs {A,B}
        // c1: {D} vs {A,B,C}
        // These are compatible (nested: {D} subset of {C,D})
        let m = make_binary_matrix(
            &["A", "B", "C", "D"],
            &[
                vec![false, false],
                vec![false, false],
                vec![true, false],
                vec![true, true],
            ],
        );
        assert!(are_compatible(&m, 0, 1));
        let result = find_max_clique(&m);
        assert_eq!(result.clique_size, 2);
        assert_eq!(result.clique_characters, vec![0, 1]);
    }
}

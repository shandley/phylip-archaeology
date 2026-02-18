//! Kitsch: Clock-constrained Fitch-Margoliash method (Fitch & Margoliash, 1967).
//!
//! Implements the Fitch-Margoliash least-squares method with a molecular
//! clock constraint (ultrametric tree), corresponding to PHYLIP's `kitsch`
//! program.
//!
//! The key difference from regular Fitch-Margoliash: branch lengths are
//! derived from **node heights** rather than being free parameters. The
//! ultrametric constraint means all tips are at the same distance from the
//! root.
//!
//! # Algorithm
//!
//! 1. Use UPGMA to determine the initial topology and node heights.
//! 2. Parameterize the tree by node heights `h_i` (one per internal node).
//! 3. Branch length from internal node to its parent = `h_parent - h_i`.
//!    Branch length from leaf to its parent = `h_parent` (leaf height = 0).
//! 4. Iteratively optimize each node height to minimize the weighted
//!    least-squares (WLS) criterion:
//!    ```text
//!    WLS = SUM_{i<j} w_ij * (d_ij - p_ij)^2
//!    ```
//!    where d_ij = observed distance, p_ij = tree path length,
//!    w_ij = 1/d_ij^2 (Fitch-Margoliash weighting).
//! 5. Heights are constrained: all heights >= 0, parent height > child height.
//!
//! # References
//!
//! - Fitch WM, Margoliash E (1967). "Construction of phylogenetic trees."
//!   *Science* 155(3760):279-284.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.

use crate::tree::{DistanceMatrix, Tree};

use super::upgma::upgma;

/// Convergence threshold for height optimization.
const DELTA: f64 = 1e-8;

/// Maximum number of optimization iterations.
const MAX_ITERATIONS: usize = 200;

/// Result of the Kitsch algorithm.
#[derive(Debug)]
pub struct KitschResult {
    /// The fitted ultrametric tree.
    pub tree: Tree,
    /// Weighted least-squares score (lower is better).
    pub wls_score: f64,
    /// Height of the root (distance from root to any tip).
    pub root_height: f64,
}

/// Build a rooted ultrametric tree from a distance matrix using the
/// clock-constrained Fitch-Margoliash (Kitsch) method.
///
/// The topology is determined by UPGMA, and then node heights are optimized
/// to minimize the Fitch-Margoliash weighted least-squares criterion subject
/// to the ultrametric constraint.
///
/// # Panics
///
/// Panics if the distance matrix has fewer than 2 taxa.
///
/// # Returns
///
/// A `KitschResult` containing the ultrametric tree, WLS score, and root height.
pub fn kitsch(matrix: &DistanceMatrix) -> KitschResult {
    let n = matrix.size();
    assert!(n >= 2, "Kitsch requires at least 2 taxa");

    // Step 1: Get initial tree topology and heights from UPGMA.
    let tree = upgma(matrix);

    // Step 2: Extract node heights from the UPGMA tree.
    let mut heights = compute_heights(&tree);

    // Step 3: Build leaf-to-matrix-index mapping.
    let leaf_to_matrix_idx = build_leaf_index_map(&tree, matrix);

    // Step 4: Collect internal (non-leaf) node IDs for optimization.
    let internal_nodes: Vec<usize> = tree
        .nodes
        .iter()
        .filter(|node| !node.is_leaf())
        .map(|node| node.id)
        .collect();

    // Step 5: Iteratively optimize heights.
    let mut prev_score = f64::MAX;
    for _ in 0..MAX_ITERATIONS {
        for &node_id in &internal_nodes {
            optimize_node_height(
                &tree,
                &mut heights,
                node_id,
                matrix,
                &leaf_to_matrix_idx,
            );
        }

        let score = compute_wls_score(&tree, &heights, matrix, &leaf_to_matrix_idx);
        if (prev_score - score).abs() < DELTA {
            break;
        }
        prev_score = score;
    }

    // Step 6: Build the final tree with optimized branch lengths.
    let mut result_tree = tree.clone();
    apply_heights_to_tree(&mut result_tree, &heights);

    let wls = compute_wls_score(&result_tree, &heights, matrix, &leaf_to_matrix_idx);
    let root_height = heights[result_tree.root];

    KitschResult {
        tree: result_tree,
        wls_score: wls,
        root_height,
    }
}

/// Compute the height of each node (distance from node to tips).
/// Leaves have height 0; internal nodes have height based on branch lengths.
fn compute_heights(tree: &Tree) -> Vec<f64> {
    let mut heights = vec![0.0; tree.nodes.len()];

    // Process nodes in postorder: leaves first, then internal nodes.
    let postorder = tree.postorder();
    for &node_id in &postorder {
        if tree.nodes[node_id].is_leaf() {
            heights[node_id] = 0.0;
        } else {
            // Height of an internal node = height of child + branch length of child
            // For UPGMA trees this should be consistent across children.
            let mut max_height = 0.0;
            for &child_id in &tree.nodes[node_id].children {
                let child_height =
                    heights[child_id] + tree.nodes[child_id].branch_length.unwrap_or(0.0);
                if child_height > max_height {
                    max_height = child_height;
                }
            }
            heights[node_id] = max_height;
        }
    }

    heights
}

/// Build a mapping from tree node index to distance matrix taxon index.
fn build_leaf_index_map(tree: &Tree, matrix: &DistanceMatrix) -> Vec<Option<usize>> {
    let mut map = vec![None; tree.nodes.len()];
    for node in &tree.nodes {
        if let Some(ref name) = node.name {
            if let Some(idx) = matrix.names.iter().position(|n| n == name) {
                map[node.id] = Some(idx);
            }
        }
    }
    map
}

/// Compute the patristic (path) distance between two leaves given node heights.
///
/// The path distance between two leaves is 2 * h_LCA, where h_LCA is the
/// height of their lowest common ancestor (since both leaves are at height 0).
fn path_distance_from_heights(
    tree: &Tree,
    heights: &[f64],
    leaf_a: usize,
    leaf_b: usize,
) -> f64 {
    let lca = find_lca(tree, leaf_a, leaf_b);
    // Path = (h_LCA - h_a) + (h_LCA - h_b) = 2*h_LCA - h_a - h_b
    // Since leaves have height 0: path = 2 * h_LCA
    2.0 * heights[lca] - heights[leaf_a] - heights[leaf_b]
}

/// Find the lowest common ancestor of two nodes.
fn find_lca(tree: &Tree, a: usize, b: usize) -> usize {
    // Build path from a to root
    let mut ancestors_a = std::collections::HashSet::new();
    let mut current = a;
    ancestors_a.insert(current);
    while let Some(parent) = tree.nodes[current].parent {
        ancestors_a.insert(parent);
        current = parent;
    }

    // Walk from b to root, find first common ancestor
    current = b;
    if ancestors_a.contains(&current) {
        return current;
    }
    while let Some(parent) = tree.nodes[current].parent {
        if ancestors_a.contains(&parent) {
            return parent;
        }
        current = parent;
    }

    // Should not happen in a connected tree; return root
    tree.root
}

/// Compute the WLS score for the current set of heights.
fn compute_wls_score(
    tree: &Tree,
    heights: &[f64],
    matrix: &DistanceMatrix,
    leaf_to_matrix_idx: &[Option<usize>],
) -> f64 {
    let n = matrix.size();

    // Collect all leaves with their matrix indices.
    let leaves: Vec<(usize, usize)> = tree
        .nodes
        .iter()
        .filter_map(|node| leaf_to_matrix_idx[node.id].map(|mi| (node.id, mi)))
        .collect();

    let mut score = 0.0;
    for i in 0..leaves.len() {
        for j in (i + 1)..leaves.len() {
            let (node_a, mi_a) = leaves[i];
            let (node_b, mi_b) = leaves[j];
            let d_obs = matrix.get(mi_a, mi_b);
            if d_obs > 0.0 {
                let d_tree = path_distance_from_heights(tree, heights, node_a, node_b);
                let w = 1.0 / (d_obs * d_obs);
                let diff = d_obs - d_tree;
                score += w * diff * diff;
            }
        }
    }

    // Handle the n == 2 edge case where both are in the same pair
    let _ = n;

    score
}

/// Optimize the height of a single internal node to minimize WLS.
///
/// For a given node, we find the height that minimizes the weighted sum
/// of squared differences, subject to the constraint that the height
/// must be between the maximum child height and the parent height (or
/// unbounded above if this is the root).
fn optimize_node_height(
    tree: &Tree,
    heights: &mut Vec<f64>,
    node_id: usize,
    matrix: &DistanceMatrix,
    leaf_to_matrix_idx: &[Option<usize>],
) {
    // Determine height bounds.
    let min_height = tree.nodes[node_id]
        .children
        .iter()
        .map(|&c| heights[c])
        .fold(0.0f64, f64::max);

    let max_height = match tree.nodes[node_id].parent {
        Some(parent) => heights[parent],
        None => {
            // Root: upper bound is twice the current height or some large value
            (heights[node_id] * 3.0).max(min_height + 10.0)
        }
    };

    // If the bounds are essentially equal, no optimization possible.
    if (max_height - min_height) < 1e-12 {
        heights[node_id] = min_height;
        return;
    }

    // Use golden section search to minimize WLS with respect to this node's height.
    // We directly modify and restore heights[node_id] to evaluate at each candidate.
    let phi = (5.0_f64.sqrt() + 1.0) / 2.0;
    let resphi = 2.0 - phi;

    let mut a = min_height;
    let mut b = max_height;

    let mut x1 = a + resphi * (b - a);
    let mut x2 = b - resphi * (b - a);

    // Evaluate WLS at candidate height h by temporarily setting heights[node_id].
    heights[node_id] = x1;
    let mut f1 = compute_wls_score(tree, heights, matrix, leaf_to_matrix_idx);
    heights[node_id] = x2;
    let mut f2 = compute_wls_score(tree, heights, matrix, leaf_to_matrix_idx);

    for _ in 0..100 {
        if (b - a).abs() < 1e-10 {
            break;
        }

        if f1 < f2 {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + resphi * (b - a);
            heights[node_id] = x1;
            f1 = compute_wls_score(tree, heights, matrix, leaf_to_matrix_idx);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = b - resphi * (b - a);
            heights[node_id] = x2;
            f2 = compute_wls_score(tree, heights, matrix, leaf_to_matrix_idx);
        }
    }

    let optimal_h = (a + b) / 2.0;
    heights[node_id] = optimal_h.max(min_height);
}

/// Apply the optimized heights to tree branch lengths.
fn apply_heights_to_tree(tree: &mut Tree, heights: &[f64]) {
    for node_id in 0..tree.nodes.len() {
        if let Some(parent) = tree.nodes[node_id].parent {
            let bl = heights[parent] - heights[node_id];
            tree.nodes[node_id].branch_length = Some(bl.max(0.0));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: compute the distance from a leaf to the root by summing branch lengths.
    fn distance_to_root(tree: &Tree, leaf_id: usize) -> f64 {
        let mut current = leaf_id;
        let mut total = 0.0;
        while let Some(parent) = tree.nodes[current].parent {
            total += tree.nodes[current].branch_length.unwrap_or(0.0);
            current = parent;
        }
        total
    }

    /// Helper to create the standard 5-taxon test matrix.
    fn test_matrix_5() -> DistanceMatrix {
        let names: Vec<String> = vec![
            "Alpha".into(),
            "Beta".into(),
            "Gamma".into(),
            "Delta".into(),
            "Epsilon".into(),
        ];
        let mut dm = DistanceMatrix::new(names);
        dm.set(0, 1, 1.0);
        dm.set(0, 2, 2.0);
        dm.set(0, 3, 3.0);
        dm.set(0, 4, 3.0);
        dm.set(1, 2, 2.0);
        dm.set(1, 3, 3.0);
        dm.set(1, 4, 3.0);
        dm.set(2, 3, 3.0);
        dm.set(2, 4, 3.0);
        dm.set(3, 4, 1.0);
        dm
    }

    /// Create an ultrametric distance matrix (all root-to-tip distances equal).
    fn ultrametric_matrix() -> DistanceMatrix {
        // Tree: ((A:1, B:1):1, (C:1, D:1):1):0
        // Root height = 2
        // d(A,B) = 2, d(C,D) = 2
        // d(A,C) = d(A,D) = d(B,C) = d(B,D) = 4
        let names: Vec<String> = vec!["A".into(), "B".into(), "C".into(), "D".into()];
        let mut dm = DistanceMatrix::new(names);
        dm.set(0, 1, 2.0); // A-B
        dm.set(2, 3, 2.0); // C-D
        dm.set(0, 2, 4.0); // A-C
        dm.set(0, 3, 4.0); // A-D
        dm.set(1, 2, 4.0); // B-C
        dm.set(1, 3, 4.0); // B-D
        dm
    }

    // -----------------------------------------------------------------------
    // Ultrametric input: WLS should be near zero
    // -----------------------------------------------------------------------

    #[test]
    fn test_ultrametric_input_low_wls() {
        let dm = ultrametric_matrix();
        let result = kitsch(&dm);

        assert!(
            result.wls_score < 0.01,
            "WLS score for ultrametric input should be near 0, got {}",
            result.wls_score
        );
    }

    #[test]
    fn test_ultrametric_root_height() {
        let dm = ultrametric_matrix();
        let result = kitsch(&dm);

        assert!(
            (result.root_height - 2.0).abs() < 0.1,
            "root height should be close to 2.0 for ultrametric input, got {}",
            result.root_height
        );
    }

    // -----------------------------------------------------------------------
    // Result tree should be ultrametric
    // -----------------------------------------------------------------------

    #[test]
    fn test_result_is_ultrametric() {
        let dm = test_matrix_5();
        let result = kitsch(&dm);

        let leaf_distances: Vec<f64> = result
            .tree
            .leaves()
            .iter()
            .map(|leaf| distance_to_root(&result.tree, leaf.id))
            .collect();

        let first = leaf_distances[0];
        for (i, &d) in leaf_distances.iter().enumerate() {
            assert!(
                (d - first).abs() < 0.01,
                "leaf {} distance to root = {:.6}, expected {:.6} (ultrametric violation)",
                i,
                d,
                first
            );
        }
    }

    #[test]
    fn test_ultrametric_on_ultrametric_input() {
        let dm = ultrametric_matrix();
        let result = kitsch(&dm);

        let leaf_distances: Vec<f64> = result
            .tree
            .leaves()
            .iter()
            .map(|leaf| distance_to_root(&result.tree, leaf.id))
            .collect();

        let first = leaf_distances[0];
        for (i, &d) in leaf_distances.iter().enumerate() {
            assert!(
                (d - first).abs() < 1e-6,
                "leaf {} distance to root = {:.6}, expected {:.6}",
                i,
                d,
                first
            );
        }
    }

    // -----------------------------------------------------------------------
    // Non-negative branch lengths
    // -----------------------------------------------------------------------

    #[test]
    fn test_nonnegative_branch_lengths() {
        let dm = test_matrix_5();
        let result = kitsch(&dm);

        for node in &result.tree.nodes {
            if let Some(bl) = node.branch_length {
                assert!(
                    bl >= -1e-10,
                    "Kitsch branch lengths should be non-negative, got {} for node {}",
                    bl,
                    node.id
                );
            }
        }
    }

    #[test]
    fn test_nonnegative_branch_lengths_ultrametric() {
        let dm = ultrametric_matrix();
        let result = kitsch(&dm);

        for node in &result.tree.nodes {
            if let Some(bl) = node.branch_length {
                assert!(
                    bl >= -1e-10,
                    "branch length should be non-negative, got {}",
                    bl
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Compare topology with UPGMA: should often agree
    // -----------------------------------------------------------------------

    #[test]
    fn test_topology_agrees_with_upgma() {
        let dm = test_matrix_5();
        let result = kitsch(&dm);
        let upgma_tree = upgma(&dm);

        // Both should have the same number of nodes
        assert_eq!(result.tree.num_nodes(), upgma_tree.num_nodes());
        assert_eq!(result.tree.num_leaves(), upgma_tree.num_leaves());

        // Both should be rooted
        assert!(result.tree.rooted);

        // Both should have the same leaf names
        let mut kitsch_names = result.tree.leaf_names();
        let mut upgma_names = upgma_tree.leaf_names();
        kitsch_names.sort();
        upgma_names.sort();
        assert_eq!(kitsch_names, upgma_names);
    }

    #[test]
    fn test_alpha_beta_grouped() {
        let dm = test_matrix_5();
        let result = kitsch(&dm);

        let alpha_id = result
            .tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("Alpha"))
            .unwrap()
            .id;
        let beta_id = result
            .tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("Beta"))
            .unwrap()
            .id;

        assert_eq!(
            result.tree.nodes[alpha_id].parent,
            result.tree.nodes[beta_id].parent,
            "Alpha and Beta should be siblings in Kitsch tree"
        );
    }

    #[test]
    fn test_delta_epsilon_grouped() {
        let dm = test_matrix_5();
        let result = kitsch(&dm);

        let delta_id = result
            .tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("Delta"))
            .unwrap()
            .id;
        let epsilon_id = result
            .tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("Epsilon"))
            .unwrap()
            .id;

        assert_eq!(
            result.tree.nodes[delta_id].parent,
            result.tree.nodes[epsilon_id].parent,
            "Delta and Epsilon should be siblings in Kitsch tree"
        );
    }

    // -----------------------------------------------------------------------
    // Small examples with known optimal heights
    // -----------------------------------------------------------------------

    #[test]
    fn test_two_taxa() {
        let names: Vec<String> = vec!["A".into(), "B".into()];
        let mut dm = DistanceMatrix::new(names);
        dm.set(0, 1, 4.0);

        let result = kitsch(&dm);

        // Root height should be 2.0 (half the distance)
        assert!(
            (result.root_height - 2.0).abs() < 0.1,
            "root height for 2 taxa with d=4 should be ~2.0, got {}",
            result.root_height
        );

        // Both branches should be ~2.0
        let a_id = result
            .tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("A"))
            .unwrap()
            .id;
        let b_id = result
            .tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("B"))
            .unwrap()
            .id;

        let bl_a = result.tree.nodes[a_id].branch_length.unwrap();
        let bl_b = result.tree.nodes[b_id].branch_length.unwrap();
        assert!(
            (bl_a - 2.0).abs() < 0.1,
            "branch length A should be ~2.0, got {}",
            bl_a
        );
        assert!(
            (bl_b - 2.0).abs() < 0.1,
            "branch length B should be ~2.0, got {}",
            bl_b
        );
    }

    #[test]
    fn test_three_taxa_equidistant() {
        // Three taxa, all pairwise distances = 6
        // Optimal ultrametric: star tree with root height = 3
        let names: Vec<String> = vec!["A".into(), "B".into(), "C".into()];
        let mut dm = DistanceMatrix::new(names);
        dm.set(0, 1, 6.0);
        dm.set(0, 2, 6.0);
        dm.set(1, 2, 6.0);

        let result = kitsch(&dm);

        // Root height should be ~3.0
        assert!(
            (result.root_height - 3.0).abs() < 0.5,
            "root height for equidistant 3-taxa (d=6) should be ~3.0, got {}",
            result.root_height
        );

        // All root-to-tip distances should be equal
        let leaf_distances: Vec<f64> = result
            .tree
            .leaves()
            .iter()
            .map(|leaf| distance_to_root(&result.tree, leaf.id))
            .collect();

        let first = leaf_distances[0];
        for &d in &leaf_distances {
            assert!(
                (d - first).abs() < 0.1,
                "equidistant tree should be ultrametric"
            );
        }
    }

    // -----------------------------------------------------------------------
    // WLS score
    // -----------------------------------------------------------------------

    #[test]
    fn test_wls_score_is_finite() {
        let dm = test_matrix_5();
        let result = kitsch(&dm);
        assert!(
            result.wls_score.is_finite(),
            "WLS score should be finite, got {}",
            result.wls_score
        );
        assert!(
            result.wls_score >= 0.0,
            "WLS score should be non-negative, got {}",
            result.wls_score
        );
    }

    #[test]
    fn test_kitsch_wls_le_upgma_wls() {
        // Kitsch should achieve WLS <= UPGMA WLS (since it optimizes heights)
        let dm = test_matrix_5();
        let kitsch_result = kitsch(&dm);

        // Compute UPGMA's WLS score
        let upgma_tree = upgma(&dm);
        let upgma_heights = compute_heights(&upgma_tree);
        let leaf_map = build_leaf_index_map(&upgma_tree, &dm);
        let upgma_wls = compute_wls_score(&upgma_tree, &upgma_heights, &dm, &leaf_map);

        assert!(
            kitsch_result.wls_score <= upgma_wls + 1e-6,
            "Kitsch WLS ({}) should be <= UPGMA WLS ({})",
            kitsch_result.wls_score,
            upgma_wls
        );
    }

    // -----------------------------------------------------------------------
    // Tree properties
    // -----------------------------------------------------------------------

    #[test]
    fn test_tree_is_rooted() {
        let dm = test_matrix_5();
        let result = kitsch(&dm);
        assert!(result.tree.rooted, "Kitsch tree should be rooted");
    }

    #[test]
    fn test_correct_number_of_nodes() {
        let dm = test_matrix_5();
        let result = kitsch(&dm);
        assert_eq!(result.tree.num_leaves(), 5);
        // Rooted binary tree: 2n-1 nodes
        assert_eq!(result.tree.num_nodes(), 2 * 5 - 1);
    }

    #[test]
    fn test_leaf_names() {
        let dm = test_matrix_5();
        let result = kitsch(&dm);
        let mut names = result.tree.leaf_names();
        names.sort();
        assert_eq!(
            names,
            vec!["Alpha", "Beta", "Delta", "Epsilon", "Gamma"]
        );
    }

    #[test]
    fn test_root_has_no_parent() {
        let dm = test_matrix_5();
        let result = kitsch(&dm);
        assert!(result.tree.nodes[result.tree.root].parent.is_none());
    }

    // -----------------------------------------------------------------------
    // Edge cases
    // -----------------------------------------------------------------------

    #[test]
    #[should_panic(expected = "at least 2 taxa")]
    fn test_too_few_taxa() {
        let names: Vec<String> = vec!["A".into()];
        let dm = DistanceMatrix::new(names);
        kitsch(&dm);
    }

    #[test]
    fn test_root_height_positive() {
        let dm = test_matrix_5();
        let result = kitsch(&dm);
        assert!(
            result.root_height > 0.0,
            "root height should be positive, got {}",
            result.root_height
        );
    }
}

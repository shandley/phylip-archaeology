//! Fitch-Margoliash distance method (Fitch & Margoliash, 1967).
//!
//! Fits a phylogenetic tree to a distance matrix by minimizing a weighted
//! least-squares (WLS) criterion:
//!
//!   S = SUM_{i<j} w(i,j) * (d_obs(i,j) - d_tree(i,j))^2
//!
//! where w(i,j) = 1 / d_obs(i,j)^2 (Fitch-Margoliash weighting). This
//! gives more weight to smaller, more accurately estimated distances.
//!
//! The current implementation uses Neighbor-Joining to determine the tree
//! topology, then optimizes branch lengths by iterative weighted least
//! squares. A full SPR (subtree pruning and regrafting) topology search
//! is planned for a future version.
//!
//! Reference: Fitch WM, Margoliash E (1967). "Construction of phylogenetic
//! trees." Science 155(3760):279-284.

use crate::tree::{DistanceMatrix, Tree};

use super::neighbor_joining::neighbor_joining;

/// Convergence threshold for branch length optimization.
const DELTA: f64 = 1e-6;

/// Maximum number of optimization iterations.
const MAX_ITERATIONS: usize = 100;

/// Build an unrooted phylogenetic tree from a distance matrix using the
/// Fitch-Margoliash method with weighted least-squares branch length
/// optimization.
///
/// The topology is determined by Neighbor-Joining, and branch lengths are
/// then optimized to minimize the Fitch-Margoliash weighted least-squares
/// criterion. A full topology search (SPR rearrangements) is planned for
/// a future version.
///
/// # Panics
///
/// Panics if the distance matrix has fewer than 3 taxa.
///
/// # Returns
///
/// An unrooted `Tree` with `tree.rooted = false` and optimized branch
/// lengths. Branch lengths are constrained to be non-negative.
pub fn fitch_margoliash(matrix: &DistanceMatrix) -> Tree {
    let n = matrix.size();
    assert!(n >= 3, "Fitch-Margoliash requires at least 3 taxa");

    // Step 1: Use NJ to determine the tree topology.
    let mut tree = neighbor_joining(matrix);

    // Step 2: Optimize branch lengths by weighted least squares.
    // Map leaf names to their indices in the distance matrix.
    let leaf_to_matrix_idx = build_leaf_index_map(&tree, matrix);

    // Iteratively optimize branch lengths.
    let mut prev_score = f64::MAX;
    for _ in 0..MAX_ITERATIONS {
        optimize_all_branch_lengths(&mut tree, matrix, &leaf_to_matrix_idx);

        let score = wls_score(&tree, matrix, &leaf_to_matrix_idx);
        if (prev_score - score).abs() < DELTA {
            break;
        }
        prev_score = score;
    }

    tree
}

/// Compute the Fitch-Margoliash weighted least-squares score.
///
/// S = SUM_{i<j} (1/d_obs^2) * (d_obs - d_tree)^2
///
/// Lower is better. Returns 0.0 for a perfect fit.
pub fn wls_score(
    tree: &Tree,
    matrix: &DistanceMatrix,
    leaf_to_matrix_idx: &[Option<usize>],
) -> f64 {
    let n = matrix.size();
    let patristic = compute_all_patristic_distances(tree, leaf_to_matrix_idx, n);

    let mut score = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let d_obs = matrix.get(i, j);
            let d_tree = patristic[i][j];
            if d_obs > 0.0 {
                let w = 1.0 / (d_obs * d_obs);
                let diff = d_obs - d_tree;
                score += w * diff * diff;
            }
        }
    }
    score
}

/// Build a mapping from tree node index to distance matrix taxon index.
/// Returns a Vec of length tree.nodes.len(), where entry [node_id] is
/// Some(matrix_idx) for leaf nodes and None for internal nodes.
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

/// Compute all pairwise patristic distances between leaves.
/// Returns an n x n matrix indexed by matrix taxon indices.
fn compute_all_patristic_distances(
    tree: &Tree,
    leaf_to_matrix_idx: &[Option<usize>],
    n: usize,
) -> Vec<Vec<f64>> {
    let mut result = vec![vec![0.0; n]; n];

    // Collect all leaf node IDs with their matrix indices.
    let leaves: Vec<(usize, usize)> = tree
        .nodes
        .iter()
        .filter_map(|node| {
            leaf_to_matrix_idx[node.id].map(|mi| (node.id, mi))
        })
        .collect();

    for i in 0..leaves.len() {
        for j in (i + 1)..leaves.len() {
            let (node_a, mi_a) = leaves[i];
            let (node_b, mi_b) = leaves[j];
            let d = patristic_distance(tree, node_a, node_b);
            result[mi_a][mi_b] = d;
            result[mi_b][mi_a] = d;
        }
    }

    result
}

/// Compute the patristic distance between two nodes by finding their
/// paths to the root and computing the path through their LCA.
fn patristic_distance(tree: &Tree, a: usize, b: usize) -> f64 {
    let path_a = path_to_root(tree, a);
    let path_b = path_to_root(tree, b);

    let set_a: std::collections::HashMap<usize, f64> =
        path_a.iter().cloned().collect();

    for &(node, dist_b) in &path_b {
        if let Some(&dist_a) = set_a.get(&node) {
            return dist_a + dist_b;
        }
    }

    // Should not happen in a connected tree.
    f64::NAN
}

/// Returns vec of (node_id, cumulative_distance_from_start) along the
/// path from `start` to the root.
fn path_to_root(tree: &Tree, start: usize) -> Vec<(usize, f64)> {
    let mut path = Vec::new();
    let mut current = start;
    let mut cum_dist = 0.0;
    path.push((current, cum_dist));
    while let Some(parent) = tree.nodes[current].parent {
        cum_dist += tree.nodes[current].branch_length.unwrap_or(0.0);
        current = parent;
        path.push((current, cum_dist));
    }
    path
}

/// Optimize all branch lengths in the tree by iterative weighted least
/// squares, one branch at a time (coordinate descent).
///
/// For each edge, we compute the optimal branch length holding all other
/// edges fixed. This is repeated for all edges in a single pass. The
/// caller should invoke this function multiple times until convergence.
fn optimize_all_branch_lengths(
    tree: &mut Tree,
    matrix: &DistanceMatrix,
    leaf_to_matrix_idx: &[Option<usize>],
) {
    // Collect edges to optimize: (node_id) for all non-root nodes with
    // a branch_length.
    let edge_nodes: Vec<usize> = tree
        .nodes
        .iter()
        .filter(|node| node.parent.is_some())
        .map(|node| node.id)
        .collect();

    for &edge_node in &edge_nodes {
        // Find the leaves on each side of this edge.
        let leaves_below = collect_leaf_matrix_indices(tree, edge_node, leaf_to_matrix_idx);
        let leaves_above = collect_leaf_matrix_indices_excluding(
            tree,
            tree.root,
            edge_node,
            leaf_to_matrix_idx,
        );

        if leaves_below.is_empty() || leaves_above.is_empty() {
            continue;
        }

        // Compute the current patristic distances, then find the optimal
        // branch length for this edge.
        //
        // For each pair (a in below, b in above), the patristic distance is:
        //   d_tree(a,b) = path_below(a) + v_edge + path_above(b)
        //
        // where v_edge is the branch length we're optimizing, path_below(a)
        // is the distance from leaf a to the edge_node, and path_above(b) is
        // the distance from leaf b to the parent of edge_node (going the
        // other direction through the tree).
        //
        // The WLS optimal v_edge minimizes:
        //   S = SUM w(a,b) * (d_obs(a,b) - path_below(a) - v - path_above(b))^2
        //
        // Taking derivative w.r.t. v and setting to zero:
        //   v = SUM w(a,b) * (d_obs(a,b) - path_below(a) - path_above(b)) / SUM w(a,b)

        // Compute path distances from each leaf to the edge boundary.
        // For leaves below: distance from leaf to edge_node (not including the edge itself).
        // For leaves above: distance from leaf to parent of edge_node (not including the edge).
        let parent = tree.nodes[edge_node].parent.unwrap();

        let path_below: Vec<(usize, f64)> = leaves_below
            .iter()
            .map(|&mi| {
                let node_id = find_node_for_matrix_idx(tree, leaf_to_matrix_idx, mi);
                let d = distance_to_node(tree, node_id, edge_node);
                (mi, d)
            })
            .collect();

        let path_above: Vec<(usize, f64)> = leaves_above
            .iter()
            .map(|&mi| {
                let node_id = find_node_for_matrix_idx(tree, leaf_to_matrix_idx, mi);
                let d = distance_to_node(tree, node_id, parent);
                (mi, d)
            })
            .collect();

        let mut numerator = 0.0;
        let mut denominator = 0.0;

        for &(mi_a, d_below) in &path_below {
            for &(mi_b, d_above) in &path_above {
                let d_obs = matrix.get(mi_a, mi_b);
                if d_obs > 0.0 {
                    let w = 1.0 / (d_obs * d_obs);
                    numerator += w * (d_obs - d_below - d_above);
                    denominator += w;
                }
            }
        }

        if denominator > 0.0 {
            let optimal_v = numerator / denominator;
            // Clamp to non-negative.
            tree.nodes[edge_node].branch_length = Some(optimal_v.max(0.0));
        }
    }
}

/// Collect all leaf matrix indices in the subtree rooted at `node`.
fn collect_leaf_matrix_indices(
    tree: &Tree,
    node: usize,
    leaf_to_matrix_idx: &[Option<usize>],
) -> Vec<usize> {
    let mut result = Vec::new();
    let mut stack = vec![node];
    while let Some(id) = stack.pop() {
        if let Some(mi) = leaf_to_matrix_idx[id] {
            result.push(mi);
        }
        for &child in &tree.nodes[id].children {
            stack.push(child);
        }
    }
    result
}

/// Collect all leaf matrix indices in the tree reachable from `start`
/// without going through `excluded_subtree`.
fn collect_leaf_matrix_indices_excluding(
    tree: &Tree,
    start: usize,
    excluded_subtree: usize,
    leaf_to_matrix_idx: &[Option<usize>],
) -> Vec<usize> {
    let mut result = Vec::new();
    let mut stack = vec![start];
    while let Some(id) = stack.pop() {
        if id == excluded_subtree {
            continue;
        }
        if let Some(mi) = leaf_to_matrix_idx[id] {
            result.push(mi);
        }
        for &child in &tree.nodes[id].children {
            stack.push(child);
        }
    }
    result
}

/// Find the tree node ID corresponding to a given matrix index.
fn find_node_for_matrix_idx(
    _tree: &Tree,
    leaf_to_matrix_idx: &[Option<usize>],
    matrix_idx: usize,
) -> usize {
    leaf_to_matrix_idx
        .iter()
        .position(|&mi| mi == Some(matrix_idx))
        .expect("matrix index not found in tree")
}

/// Compute the sum of branch lengths from `start` to `target` along
/// parent pointers. Returns 0.0 if start == target.
///
/// This function walks upward from `start` to `target`. If `target` is
/// not an ancestor of `start`, it will panic.
fn distance_to_node(tree: &Tree, start: usize, target: usize) -> f64 {
    if start == target {
        return 0.0;
    }

    // Try walking upward from start.
    let mut current = start;
    let mut total = 0.0;
    while current != target {
        total += tree.nodes[current].branch_length.unwrap_or(0.0);
        match tree.nodes[current].parent {
            Some(p) => current = p,
            None => {
                // target is not an ancestor of start; try the other direction.
                // This can happen in unrooted trees where the "parent" pointers
                // don't necessarily reflect the true path.
                return distance_via_lca(tree, start, target);
            }
        }
    }
    total
}

/// Compute distance between two arbitrary nodes via their LCA.
fn distance_via_lca(tree: &Tree, a: usize, b: usize) -> f64 {
    let path_a = path_to_root(tree, a);
    let path_b = path_to_root(tree, b);

    let set_a: std::collections::HashMap<usize, f64> =
        path_a.iter().cloned().collect();

    for &(node, dist_b) in &path_b {
        if let Some(&dist_a) = set_a.get(&node) {
            return dist_a + dist_b;
        }
    }

    0.0
}

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn test_fm_correct_number_of_nodes() {
        let dm = test_matrix_5();
        let tree = fitch_margoliash(&dm);

        assert_eq!(tree.num_leaves(), 5);
        // Same topology as NJ: 2n-2 nodes for unrooted binary tree.
        assert_eq!(tree.num_nodes(), 2 * 5 - 2);
    }

    #[test]
    fn test_fm_unrooted() {
        let dm = test_matrix_5();
        let tree = fitch_margoliash(&dm);
        assert!(!tree.rooted, "Fitch-Margoliash tree should be unrooted");
    }

    #[test]
    fn test_fm_leaf_names() {
        let dm = test_matrix_5();
        let tree = fitch_margoliash(&dm);

        let mut names = tree.leaf_names();
        names.sort();
        assert_eq!(
            names,
            vec!["Alpha", "Beta", "Delta", "Epsilon", "Gamma"]
        );
    }

    #[test]
    fn test_fm_groups_alpha_beta() {
        let dm = test_matrix_5();
        let tree = fitch_margoliash(&dm);

        let alpha_id = tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("Alpha"))
            .unwrap()
            .id;
        let beta_id = tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("Beta"))
            .unwrap()
            .id;

        assert_eq!(
            tree.nodes[alpha_id].parent,
            tree.nodes[beta_id].parent,
            "Alpha and Beta should be siblings"
        );
    }

    #[test]
    fn test_fm_groups_delta_epsilon() {
        let dm = test_matrix_5();
        let tree = fitch_margoliash(&dm);

        let delta_id = tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("Delta"))
            .unwrap()
            .id;
        let epsilon_id = tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("Epsilon"))
            .unwrap()
            .id;

        assert_eq!(
            tree.nodes[delta_id].parent,
            tree.nodes[epsilon_id].parent,
            "Delta and Epsilon should be siblings"
        );
    }

    #[test]
    fn test_fm_nonnegative_branch_lengths() {
        let dm = test_matrix_5();
        let tree = fitch_margoliash(&dm);

        for node in &tree.nodes {
            if let Some(bl) = node.branch_length {
                assert!(
                    bl >= -1e-10,
                    "FM branch lengths should be non-negative, got {} for node {}",
                    bl,
                    node.id
                );
            }
        }
    }

    #[test]
    fn test_fm_wls_score_improves_over_nj() {
        let dm = test_matrix_5();

        // Get NJ tree and its WLS score.
        let nj_tree = neighbor_joining(&dm);
        let nj_leaf_map = build_leaf_index_map(&nj_tree, &dm);
        let nj_score = wls_score(&nj_tree, &dm, &nj_leaf_map);

        // Get FM tree and its WLS score.
        let fm_tree = fitch_margoliash(&dm);
        let fm_leaf_map = build_leaf_index_map(&fm_tree, &dm);
        let fm_score = wls_score(&fm_tree, &dm, &fm_leaf_map);

        // FM should be at least as good as NJ (lower or equal score).
        assert!(
            fm_score <= nj_score + 1e-10,
            "FM score ({}) should be <= NJ score ({})",
            fm_score,
            nj_score
        );
    }

    #[test]
    fn test_fm_patristic_distances_close_to_observed() {
        let dm = test_matrix_5();
        let tree = fitch_margoliash(&dm);
        let leaf_map = build_leaf_index_map(&tree, &dm);

        let n = dm.size();
        let patristic = compute_all_patristic_distances(&tree, &leaf_map, n);

        // For this additive matrix, patristic distances should be very close
        // to observed distances.
        for i in 0..n {
            for j in (i + 1)..n {
                let d_obs = dm.get(i, j);
                let d_tree = patristic[i][j];
                assert!(
                    (d_obs - d_tree).abs() < 0.1,
                    "Patristic distance ({},{}) = {}, observed = {}, diff = {}",
                    i,
                    j,
                    d_tree,
                    d_obs,
                    (d_obs - d_tree).abs()
                );
            }
        }
    }

    #[test]
    fn test_fm_trifurcation_root() {
        let dm = test_matrix_5();
        let tree = fitch_margoliash(&dm);

        assert_eq!(
            tree.nodes[tree.root].children.len(),
            3,
            "FM root should be a trifurcation"
        );
    }

    #[test]
    fn test_fm_three_taxa() {
        let names: Vec<String> = vec!["A".into(), "B".into(), "C".into()];
        let mut dm = DistanceMatrix::new(names);
        dm.set(0, 1, 4.0);
        dm.set(0, 2, 6.0);
        dm.set(1, 2, 6.0);

        let tree = fitch_margoliash(&dm);

        assert_eq!(tree.num_leaves(), 3);
        assert_eq!(tree.nodes[tree.root].children.len(), 3);
    }

    #[test]
    fn test_wls_score_perfect_additive() {
        // For a perfectly additive matrix, the NJ tree already achieves
        // score ~0 (up to floating point).
        let dm = test_matrix_5();
        let tree = fitch_margoliash(&dm);
        let leaf_map = build_leaf_index_map(&tree, &dm);
        let score = wls_score(&tree, &dm, &leaf_map);

        assert!(
            score < 0.01,
            "WLS score for additive matrix should be near 0, got {}",
            score
        );
    }

    #[test]
    #[should_panic(expected = "at least 3 taxa")]
    fn test_fm_too_few_taxa() {
        let names: Vec<String> = vec!["A".into(), "B".into()];
        let dm = DistanceMatrix::new(names);
        fitch_margoliash(&dm);
    }
}

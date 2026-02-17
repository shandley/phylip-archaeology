//! Neighbor-joining algorithm (Saitou & Nei, 1987).
//!
//! Constructs an unrooted phylogenetic tree from a distance matrix by
//! iteratively joining the pair of nodes that minimizes the Q-criterion:
//!
//!   Q(i,j) = (n-2) * d(i,j) - R(i) - R(j)
//!
//! where R(i) is the sum of distances from taxon i to all other active taxa.
//!
//! The algorithm is O(n^3) time and O(n^2) space, and produces an unrooted
//! tree that may contain negative branch lengths.
//!
//! Reference: Saitou N, Nei M (1987). "The neighbor-joining method: a new
//! method for reconstructing phylogenetic trees." Molecular Biology and
//! Evolution 4(4):406-425.

use crate::tree::{DistanceMatrix, Tree};

/// Build an unrooted phylogenetic tree from a distance matrix using the
/// Neighbor-Joining algorithm.
///
/// # Panics
///
/// Panics if the distance matrix has fewer than 3 taxa.
///
/// # Returns
///
/// An unrooted `Tree` with `tree.rooted = false`. The root node is the
/// central trifurcation connecting the final three clusters. NJ can produce
/// negative branch lengths; these are preserved as-is.
pub fn neighbor_joining(matrix: &DistanceMatrix) -> Tree {
    let n = matrix.size();
    assert!(n >= 3, "Neighbor-joining requires at least 3 taxa");

    // Work with a full distance matrix that we modify in place.
    // dist[i][j] gives distance between active clusters i and j.
    let mut dist = matrix.to_full_matrix();

    // Track which clusters are still active.
    let mut active: Vec<bool> = vec![true; n];

    // The tree node index corresponding to each cluster slot.
    // Initially each slot i maps to leaf node i in the tree.
    let mut cluster_node: Vec<usize> = Vec::with_capacity(n);

    // Accumulated branch length offsets (av[] in PHYLIP).
    let mut av: Vec<f64> = vec![0.0; n];

    // Build the tree.
    let mut tree = Tree::new();
    tree.rooted = false;

    // Create leaf nodes (indices 0..n-1).
    for i in 0..n {
        let id = tree.add_node(Some(matrix.names[i].clone()), None, None);
        cluster_node.push(id);
    }

    // Number of active clusters.
    let mut num_active = n;

    // Main NJ loop: iterate until 3 clusters remain.
    while num_active > 3 {
        // Step 1: Compute R[i] = sum of distances from i to all other active clusters.
        let mut r: Vec<f64> = vec![0.0; n];
        for i in 0..n {
            if !active[i] {
                continue;
            }
            let mut sum = 0.0;
            for j in 0..n {
                if j != i && active[j] {
                    sum += dist[i][j];
                }
            }
            r[i] = sum;
        }

        // Step 2: Find pair (mini, minj) minimizing Q(i,j) = (num_active-2)*d(i,j) - R[i] - R[j].
        let mut mini = 0;
        let mut minj = 0;
        let mut q_min = f64::MAX;
        let fotu2 = (num_active - 2) as f64;

        for i in 0..n {
            if !active[i] {
                continue;
            }
            for j in (i + 1)..n {
                if !active[j] {
                    continue;
                }
                let q = fotu2 * dist[i][j] - r[i] - r[j];
                if q < q_min {
                    q_min = q;
                    mini = i;
                    minj = j;
                }
            }
        }

        // Step 3: Compute branch lengths.
        let d_ij = dist[mini][minj];
        let bi = d_ij / 2.0 + (r[mini] - r[minj]) / (2.0 * fotu2);
        let bj = d_ij - bi;

        // Adjust for accumulated offsets.
        let branch_i = bi - av[mini];
        let branch_j = bj - av[minj];

        // Step 4: Create a new internal node joining mini and minj.
        let new_node = tree.add_node(None, None, None);

        // Attach children with branch lengths.
        let child_i = cluster_node[mini];
        let child_j = cluster_node[minj];
        tree.nodes[child_i].parent = Some(new_node);
        tree.nodes[child_i].branch_length = Some(branch_i);
        tree.nodes[child_j].parent = Some(new_node);
        tree.nodes[child_j].branch_length = Some(branch_j);
        tree.nodes[new_node].children.push(child_i);
        tree.nodes[new_node].children.push(child_j);

        // Step 5: Update distances using the PHYLIP formula.
        // d(new, k) = (d(i,k) + d(j,k)) / 2
        //
        // This is a simple average (not the textbook NJ formula which subtracts
        // d(i,j)). The accumulated branch lengths in av[] handle the adjustment
        // that the textbook formula bakes into the distance. This formulation is
        // equivalent and matches the original PHYLIP neighbor.c implementation.
        for k in 0..n {
            if !active[k] || k == mini || k == minj {
                continue;
            }
            let d_new = (dist[mini][k] + dist[minj][k]) / 2.0;
            dist[mini][k] = d_new;
            dist[k][mini] = d_new;
        }

        // Step 6: Update bookkeeping.
        // mini's slot now represents the new cluster.
        active[minj] = false;
        cluster_node[mini] = new_node;
        av[mini] = d_ij / 2.0;
        num_active -= 1;
    }

    // Final step: connect the last 3 clusters to a trifurcation.
    // Collect the remaining active indices.
    let remaining: Vec<usize> = (0..n).filter(|&i| active[i]).collect();
    assert_eq!(remaining.len(), 3);

    let a = remaining[0];
    let b = remaining[1];
    let c = remaining[2];

    // Compute branch lengths for the trifurcation using the triangle formula:
    //   bi = (d(a,b) + d(a,c) - d(b,c)) / 2 - av[a]
    //   bj = (d(a,b) + d(b,c) - d(a,c)) / 2 - av[b]
    //   bk = (d(a,c) + d(b,c) - d(a,b)) / 2 - av[c]
    let d_ab = dist[a][b];
    let d_ac = dist[a][c];
    let d_bc = dist[b][c];

    let branch_a = (d_ab + d_ac - d_bc) / 2.0 - av[a];
    let branch_b = (d_ab + d_bc - d_ac) / 2.0 - av[b];
    let branch_c = (d_ac + d_bc - d_ab) / 2.0 - av[c];

    // Create the central trifurcation node.
    let root = tree.add_node(None, None, None);
    tree.root = root;

    let node_a = cluster_node[a];
    let node_b = cluster_node[b];
    let node_c = cluster_node[c];

    tree.nodes[node_a].parent = Some(root);
    tree.nodes[node_a].branch_length = Some(branch_a);
    tree.nodes[node_b].parent = Some(root);
    tree.nodes[node_b].branch_length = Some(branch_b);
    tree.nodes[node_c].parent = Some(root);
    tree.nodes[node_c].branch_length = Some(branch_c);

    tree.nodes[root].children.push(node_a);
    tree.nodes[root].children.push(node_b);
    tree.nodes[root].children.push(node_c);

    tree
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
        // Alpha-Beta: 1, Alpha-Gamma: 2, Alpha-Delta: 3, Alpha-Epsilon: 3
        dm.set(0, 1, 1.0);
        dm.set(0, 2, 2.0);
        dm.set(0, 3, 3.0);
        dm.set(0, 4, 3.0);
        // Beta-Gamma: 2, Beta-Delta: 3, Beta-Epsilon: 3
        dm.set(1, 2, 2.0);
        dm.set(1, 3, 3.0);
        dm.set(1, 4, 3.0);
        // Gamma-Delta: 3, Gamma-Epsilon: 3
        dm.set(2, 3, 3.0);
        dm.set(2, 4, 3.0);
        // Delta-Epsilon: 1
        dm.set(3, 4, 1.0);
        dm
    }

    /// Find which leaf names are descendants of a given node.
    fn leaf_names_under(tree: &Tree, node_id: usize) -> Vec<String> {
        let mut result = Vec::new();
        let mut stack = vec![node_id];
        while let Some(id) = stack.pop() {
            let node = &tree.nodes[id];
            if node.is_leaf() {
                if let Some(ref name) = node.name {
                    result.push(name.clone());
                }
            }
            for &child in &node.children {
                stack.push(child);
            }
        }
        result.sort();
        result
    }

    #[test]
    fn test_nj_correct_number_of_nodes() {
        let dm = test_matrix_5();
        let tree = neighbor_joining(&dm);

        // 5 taxa => 5 leaves + internal nodes
        assert_eq!(tree.num_leaves(), 5);

        // For NJ with n taxa: n leaves + (n-2) internal nodes = 2n-2 total
        // (the trifurcation means one fewer internal node than a fully resolved
        // binary tree).
        assert_eq!(tree.num_nodes(), 2 * 5 - 2);
    }

    #[test]
    fn test_nj_leaf_names() {
        let dm = test_matrix_5();
        let tree = neighbor_joining(&dm);

        let mut names = tree.leaf_names();
        names.sort();
        assert_eq!(
            names,
            vec!["Alpha", "Beta", "Delta", "Epsilon", "Gamma"]
        );
    }

    #[test]
    fn test_nj_groups_alpha_beta() {
        let dm = test_matrix_5();
        let tree = neighbor_joining(&dm);

        // Alpha and Beta (distance 1) should be grouped together as siblings.
        // Find Alpha and Beta leaf nodes.
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

        // They should share the same parent.
        assert_eq!(
            tree.nodes[alpha_id].parent,
            tree.nodes[beta_id].parent,
            "Alpha and Beta should be siblings"
        );
    }

    #[test]
    fn test_nj_groups_delta_epsilon() {
        let dm = test_matrix_5();
        let tree = neighbor_joining(&dm);

        // Delta and Epsilon (distance 1) should be grouped together.
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
    fn test_nj_branch_lengths_sum_to_distances() {
        let dm = test_matrix_5();
        let tree = neighbor_joining(&dm);

        // For NJ, the patristic distance (sum of branch lengths on path)
        // between any pair of leaves should approximately equal the input
        // distance (exactly for an additive matrix).

        // This matrix is additive, so the distances should be exact.
        // Compute patristic distance between Alpha and Beta.
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

        let patristic = patristic_distance(&tree, alpha_id, beta_id);
        assert!(
            (patristic - 1.0).abs() < 1e-10,
            "Patristic distance Alpha-Beta should be 1.0, got {}",
            patristic
        );

        // Delta-Epsilon should be 1.0.
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

        let patristic = patristic_distance(&tree, delta_id, epsilon_id);
        assert!(
            (patristic - 1.0).abs() < 1e-10,
            "Patristic distance Delta-Epsilon should be 1.0, got {}",
            patristic
        );

        // Alpha-Gamma should be 2.0.
        let gamma_id = tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("Gamma"))
            .unwrap()
            .id;
        let patristic = patristic_distance(&tree, alpha_id, gamma_id);
        assert!(
            (patristic - 2.0).abs() < 1e-10,
            "Patristic distance Alpha-Gamma should be 2.0, got {}",
            patristic
        );

        // Alpha-Delta should be 3.0.
        let patristic = patristic_distance(&tree, alpha_id, delta_id);
        assert!(
            (patristic - 3.0).abs() < 1e-10,
            "Patristic distance Alpha-Delta should be 3.0, got {}",
            patristic
        );
    }

    #[test]
    fn test_nj_unrooted() {
        let dm = test_matrix_5();
        let tree = neighbor_joining(&dm);
        assert!(!tree.rooted, "NJ tree should be unrooted");
    }

    #[test]
    fn test_nj_trifurcation_root() {
        let dm = test_matrix_5();
        let tree = neighbor_joining(&dm);

        // The root should have exactly 3 children (trifurcation).
        assert_eq!(
            tree.nodes[tree.root].children.len(),
            3,
            "NJ root should be a trifurcation"
        );
    }

    #[test]
    fn test_nj_three_taxa() {
        // Minimal case: 3 taxa directly form a trifurcation.
        let names: Vec<String> = vec!["A".into(), "B".into(), "C".into()];
        let mut dm = DistanceMatrix::new(names);
        dm.set(0, 1, 4.0);
        dm.set(0, 2, 6.0);
        dm.set(1, 2, 6.0);

        let tree = neighbor_joining(&dm);

        assert_eq!(tree.num_leaves(), 3);
        assert_eq!(tree.nodes[tree.root].children.len(), 3);

        // Branch lengths from the triangle formula:
        // bA = (4+6-6)/2 = 2, bB = (4+6-6)/2 = 2, bC = (6+6-4)/2 = 4
        let a_id = tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("A"))
            .unwrap()
            .id;
        let b_id = tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("B"))
            .unwrap()
            .id;
        let c_id = tree
            .nodes
            .iter()
            .find(|n| n.name.as_deref() == Some("C"))
            .unwrap()
            .id;

        assert!(
            (tree.nodes[a_id].branch_length.unwrap() - 2.0).abs() < 1e-10
        );
        assert!(
            (tree.nodes[b_id].branch_length.unwrap() - 2.0).abs() < 1e-10
        );
        assert!(
            (tree.nodes[c_id].branch_length.unwrap() - 4.0).abs() < 1e-10
        );
    }

    #[test]
    #[should_panic(expected = "at least 3 taxa")]
    fn test_nj_too_few_taxa() {
        let names: Vec<String> = vec!["A".into(), "B".into()];
        let dm = DistanceMatrix::new(names);
        neighbor_joining(&dm);
    }

    /// Compute the patristic distance between two nodes by finding their
    /// paths to the root and computing the path through their LCA.
    fn patristic_distance(tree: &Tree, a: usize, b: usize) -> f64 {
        // Get path from a to root with cumulative distances.
        let path_a = path_to_root(tree, a);
        let path_b = path_to_root(tree, b);

        // Find the LCA (lowest common ancestor).
        let set_a: std::collections::HashMap<usize, f64> =
            path_a.iter().cloned().collect();

        for &(node, dist_b) in &path_b {
            if let Some(&dist_a) = set_a.get(&node) {
                return dist_a + dist_b;
            }
        }

        panic!("No common ancestor found");
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
}

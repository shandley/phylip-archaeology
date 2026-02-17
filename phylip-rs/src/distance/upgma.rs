//! UPGMA algorithm (Sokal & Michener, 1958).
//!
//! Unweighted Pair Group Method with Arithmetic Mean. Produces a rooted
//! ultrametric tree from a distance matrix by iteratively joining the pair
//! of clusters with the smallest distance.
//!
//! UPGMA assumes a molecular clock (equal rates across all lineages), which
//! means all tips are equidistant from the root. The resulting tree is
//! ultrametric: for any internal node, the distance from that node to all
//! of its descendant leaves is the same.
//!
//! The algorithm is O(n^3) time and O(n^2) space.
//!
//! Reference: Sokal RR, Michener CD (1958). "A statistical method for
//! evaluating systematic relationships." University of Kansas Science
//! Bulletin 38:1409-1438.

use crate::tree::{DistanceMatrix, Tree};

/// Build a rooted ultrametric tree from a distance matrix using the UPGMA
/// algorithm.
///
/// # Panics
///
/// Panics if the distance matrix has fewer than 2 taxa.
///
/// # Returns
///
/// A rooted `Tree` with `tree.rooted = true`. Branch lengths are always
/// non-negative (assuming a valid distance matrix), and the tree satisfies
/// the ultrametric property.
pub fn upgma(matrix: &DistanceMatrix) -> Tree {
    let n = matrix.size();
    assert!(n >= 2, "UPGMA requires at least 2 taxa");

    // Work with a full distance matrix that we modify in place.
    let mut dist = matrix.to_full_matrix();

    // Track which clusters are still active.
    let mut active: Vec<bool> = vec![true; n];

    // The tree node index corresponding to each cluster slot.
    let mut cluster_node: Vec<usize> = Vec::with_capacity(n);

    // Height (distance from tips) of each cluster.
    let mut height: Vec<f64> = vec![0.0; n];

    // Size (number of OTUs) of each cluster.
    let mut cluster_size: Vec<usize> = vec![1; n];

    // Build the tree.
    let mut tree = Tree::new();
    tree.rooted = true;

    // Create leaf nodes (indices 0..n-1).
    for i in 0..n {
        let id = tree.add_node(Some(matrix.names[i].clone()), None, None);
        cluster_node.push(id);
    }

    // Number of active clusters.
    let mut num_active = n;

    // Main UPGMA loop: iterate until 1 cluster remains.
    while num_active > 1 {
        // Step 1: Find pair (mini, minj) with minimum d(i,j).
        let mut mini = 0;
        let mut minj = 0;
        let mut d_min = f64::MAX;

        for i in 0..n {
            if !active[i] {
                continue;
            }
            for j in (i + 1)..n {
                if !active[j] {
                    continue;
                }
                if dist[i][j] < d_min {
                    d_min = dist[i][j];
                    mini = i;
                    minj = j;
                }
            }
        }

        // Step 2: Compute branch lengths.
        // The new node has height = d_min / 2.
        let new_height = d_min / 2.0;
        let branch_i = new_height - height[mini];
        let branch_j = new_height - height[minj];

        // Step 3: Create a new internal node joining mini and minj.
        let new_node = tree.add_node(None, None, None);

        let child_i = cluster_node[mini];
        let child_j = cluster_node[minj];
        tree.nodes[child_i].parent = Some(new_node);
        tree.nodes[child_i].branch_length = Some(branch_i);
        tree.nodes[child_j].parent = Some(new_node);
        tree.nodes[child_j].branch_length = Some(branch_j);
        tree.nodes[new_node].children.push(child_i);
        tree.nodes[new_node].children.push(child_j);

        // Step 4: Update distances using weighted average.
        // d(new, k) = (d(i,k) * n_i + d(j,k) * n_j) / (n_i + n_j)
        let ni = cluster_size[mini] as f64;
        let nj = cluster_size[minj] as f64;
        for k in 0..n {
            if !active[k] || k == mini || k == minj {
                continue;
            }
            let d_new = (dist[mini][k] * ni + dist[minj][k] * nj) / (ni + nj);
            dist[mini][k] = d_new;
            dist[k][mini] = d_new;
        }

        // Step 5: Update bookkeeping.
        active[minj] = false;
        cluster_node[mini] = new_node;
        height[mini] = new_height;
        cluster_size[mini] += cluster_size[minj];
        num_active -= 1;
    }

    // The last remaining cluster is the root.
    let root_idx = (0..n).find(|&i| active[i]).unwrap();
    tree.root = cluster_node[root_idx];

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
    fn test_upgma_correct_number_of_nodes() {
        let dm = test_matrix_5();
        let tree = upgma(&dm);

        // 5 leaves.
        assert_eq!(tree.num_leaves(), 5);

        // For UPGMA with n taxa: n leaves + (n-1) internal nodes = 2n-1 total.
        assert_eq!(tree.num_nodes(), 2 * 5 - 1);
    }

    #[test]
    fn test_upgma_rooted() {
        let dm = test_matrix_5();
        let tree = upgma(&dm);
        assert!(tree.rooted, "UPGMA tree should be rooted");
    }

    #[test]
    fn test_upgma_leaf_names() {
        let dm = test_matrix_5();
        let tree = upgma(&dm);

        let mut names = tree.leaf_names();
        names.sort();
        assert_eq!(
            names,
            vec!["Alpha", "Beta", "Delta", "Epsilon", "Gamma"]
        );
    }

    #[test]
    fn test_upgma_groups_alpha_beta() {
        let dm = test_matrix_5();
        let tree = upgma(&dm);

        // Alpha and Beta (distance 1) should be grouped first.
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
    fn test_upgma_groups_delta_epsilon() {
        let dm = test_matrix_5();
        let tree = upgma(&dm);

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
    fn test_upgma_ultrametric_branch_lengths() {
        let dm = test_matrix_5();
        let tree = upgma(&dm);

        // Alpha and Beta have distance 1.0, so their common ancestor
        // should be at height 0.5. Branch lengths should each be 0.5.
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

        assert!(
            (tree.nodes[alpha_id].branch_length.unwrap() - 0.5).abs() < 1e-10,
            "Alpha branch length should be 0.5, got {}",
            tree.nodes[alpha_id].branch_length.unwrap()
        );
        assert!(
            (tree.nodes[beta_id].branch_length.unwrap() - 0.5).abs() < 1e-10,
            "Beta branch length should be 0.5, got {}",
            tree.nodes[beta_id].branch_length.unwrap()
        );

        // Similarly for Delta-Epsilon.
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

        assert!(
            (tree.nodes[delta_id].branch_length.unwrap() - 0.5).abs() < 1e-10,
        );
        assert!(
            (tree.nodes[epsilon_id].branch_length.unwrap() - 0.5).abs() < 1e-10,
        );
    }

    #[test]
    fn test_upgma_ultrametric_property() {
        let dm = test_matrix_5();
        let tree = upgma(&dm);

        // Ultrametric property: the distance from the root to every leaf
        // should be the same.
        let root = tree.root;
        let leaf_distances: Vec<f64> = tree
            .leaves()
            .iter()
            .map(|leaf| distance_to_ancestor(&tree, leaf.id, root))
            .collect();

        let first = leaf_distances[0];
        for (i, &d) in leaf_distances.iter().enumerate() {
            assert!(
                (d - first).abs() < 1e-10,
                "Leaf {} has distance {} to root, expected {} (ultrametric violation)",
                i,
                d,
                first
            );
        }
    }

    #[test]
    fn test_upgma_nonnegative_branch_lengths() {
        let dm = test_matrix_5();
        let tree = upgma(&dm);

        for node in &tree.nodes {
            if let Some(bl) = node.branch_length {
                assert!(
                    bl >= -1e-10,
                    "UPGMA should produce non-negative branch lengths, got {} for node {}",
                    bl,
                    node.id
                );
            }
        }
    }

    #[test]
    fn test_upgma_two_taxa() {
        let names: Vec<String> = vec!["A".into(), "B".into()];
        let mut dm = DistanceMatrix::new(names);
        dm.set(0, 1, 4.0);

        let tree = upgma(&dm);

        assert_eq!(tree.num_leaves(), 2);
        assert_eq!(tree.num_nodes(), 3); // 2 leaves + 1 internal

        // Each branch should be 2.0 (half the distance).
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

        assert!((tree.nodes[a_id].branch_length.unwrap() - 2.0).abs() < 1e-10);
        assert!((tree.nodes[b_id].branch_length.unwrap() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_upgma_root_has_no_parent() {
        let dm = test_matrix_5();
        let tree = upgma(&dm);
        assert!(tree.nodes[tree.root].parent.is_none());
    }

    #[test]
    #[should_panic(expected = "at least 2 taxa")]
    fn test_upgma_too_few_taxa() {
        let names: Vec<String> = vec!["A".into()];
        let dm = DistanceMatrix::new(names);
        upgma(&dm);
    }

    /// Compute the sum of branch lengths from `node` up to `ancestor`.
    fn distance_to_ancestor(tree: &Tree, node: usize, ancestor: usize) -> f64 {
        let mut current = node;
        let mut total = 0.0;
        while current != ancestor {
            total += tree.nodes[current].branch_length.unwrap_or(0.0);
            current = tree.nodes[current]
                .parent
                .expect("reached root without finding ancestor");
        }
        total
    }
}

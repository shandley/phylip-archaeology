//! Felsenstein's phylogenetically independent contrasts (1985).
//!
//! Species data points are not statistically independent because species
//! share evolutionary history. Independent contrasts transform correlated
//! species data into a set of statistically independent values that can
//! be used in standard statistical analyses (regression, correlation, ANOVA).
//!
//! For n taxa, the algorithm produces exactly n-1 contrasts via postorder
//! traversal. At each internal node with children L and R:
//!
//! ```text
//! contrast = X_L - X_R
//! variance = v_L + v_R
//! standardized_contrast = contrast / sqrt(variance)
//!
//! X_parent = (X_L * v_R + X_R * v_L) / (v_L + v_R)
//! v_additional = (v_L * v_R) / (v_L + v_R)
//! ```
//!
//! The branch length from the parent to its own parent is extended by
//! `v_additional`, propagating the variance upward through the tree.
//!
//! Under Brownian motion, these standardized contrasts are independent
//! standard normal variates, enabling valid statistical tests on
//! phylogenetically structured data.
//!
//! # References
//!
//! - Felsenstein, J. (1985). Phylogenies and the comparative method.
//!   *American Naturalist*, 125, 1-15.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//!   Chapter 25: Quantitative characters, phylogenies, and morphometrics.

use crate::tree::Tree;

use super::ContinuousData;

/// Result of independent contrasts analysis.
#[derive(Debug, Clone)]
pub struct ContrastResult {
    /// Standardized contrasts for each character.
    /// `contrasts[character]` has `n_taxa - 1` values.
    pub contrasts: Vec<Vec<f64>>,
    /// Variance (sum of branch lengths) for each contrast.
    /// Length is `n_taxa - 1`.
    pub variances: Vec<f64>,
    /// Estimated root value for each character (weighted average at root).
    pub root_values: Vec<f64>,
    /// Correlation matrix between characters (if multiple characters).
    /// `correlations[i][j]` is the PIC correlation between characters i and j.
    pub correlations: Option<Vec<Vec<f64>>>,
}

/// Default branch length used when a node has no branch length specified.
const DEFAULT_BRANCH_LENGTH: f64 = 1.0;

/// Minimum branch length for numerical stability.
const MIN_BRANCH_LENGTH: f64 = 1e-10;

/// Compute phylogenetically independent contrasts.
///
/// Given a rooted phylogenetic tree and continuous character data for the
/// leaf taxa, computes Felsenstein's (1985) independent contrasts via
/// postorder traversal.
///
/// The tree must be binary (each internal node has exactly 2 children).
/// Leaf names must match taxon names in the data. Branch lengths are
/// required; a default of 1.0 is used if missing.
///
/// # Errors
///
/// Returns an error if:
/// - The tree has no leaves
/// - Leaf names do not match data taxa
/// - An internal node has other than 2 children (non-binary tree)
/// - Branch length sum is zero or negative at any contrast
pub fn independent_contrasts(
    tree: &Tree,
    data: &ContinuousData,
) -> Result<ContrastResult, String> {
    let n_leaves = tree.num_leaves();
    if n_leaves == 0 {
        return Err("tree has no leaves".to_string());
    }
    if n_leaves < 2 {
        return Err("need at least 2 taxa for contrasts".to_string());
    }

    let n_chars = data.n_chars;

    // Build a mapping from leaf name to taxon index in data.
    let name_to_idx: std::collections::HashMap<&str, usize> = data
        .taxa
        .iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i))
        .collect();

    let n_nodes = tree.num_nodes();

    // Per-node character values (propagated upward during postorder).
    // For leaves, these are the observed data values.
    // For internal nodes, these become the weighted averages.
    let mut node_values: Vec<Vec<f64>> = vec![vec![0.0; n_chars]; n_nodes];

    // Effective branch length for each node, modified during traversal.
    // This starts as the original branch length and gets augmented with
    // the additional variance term at each internal node.
    let mut effective_length: Vec<f64> = vec![0.0; n_nodes];

    // Initialize leaf node values and effective lengths.
    for node in &tree.nodes {
        if node.is_leaf() {
            let name = node.name.as_deref().ok_or_else(|| {
                format!("leaf node {} has no name", node.id)
            })?;
            let taxon_idx = name_to_idx.get(name).ok_or_else(|| {
                format!("leaf '{}' not found in data", name)
            })?;
            for c in 0..n_chars {
                node_values[node.id][c] = data.values[*taxon_idx][c];
            }
        }
        effective_length[node.id] = node
            .branch_length
            .unwrap_or(DEFAULT_BRANCH_LENGTH)
            .max(MIN_BRANCH_LENGTH);
    }

    // Collect contrasts during postorder traversal.
    let mut contrasts: Vec<Vec<f64>> = vec![Vec::new(); n_chars];
    let mut variances: Vec<f64> = Vec::new();

    let postorder = tree.postorder();

    for &node_id in &postorder {
        let children = &tree.nodes[node_id].children;
        if children.is_empty() {
            continue; // Leaf node: already initialized.
        }

        if children.len() != 2 {
            return Err(format!(
                "independent contrasts requires a binary tree, but node {} has {} children",
                node_id,
                children.len()
            ));
        }

        let left_id = children[0];
        let right_id = children[1];

        let v_l = effective_length[left_id];
        let v_r = effective_length[right_id];
        let v_sum = v_l + v_r;

        if v_sum <= 0.0 {
            return Err(format!(
                "zero or negative variance at node {}: v_L={}, v_R={}",
                node_id, v_l, v_r
            ));
        }

        let sqrt_v = v_sum.sqrt();
        variances.push(v_sum);

        for c in 0..n_chars {
            let x_l = node_values[left_id][c];
            let x_r = node_values[right_id][c];

            // Raw contrast (difference).
            let contrast = x_l - x_r;

            // Standardized contrast (divided by sqrt of variance).
            let standardized = contrast / sqrt_v;
            contrasts[c].push(standardized);

            // Propagate weighted average upward.
            node_values[node_id][c] = (x_l * v_r + x_r * v_l) / v_sum;
        }

        // Additional variance from this node, propagated to its parent.
        let v_additional = (v_l * v_r) / v_sum;

        // The effective branch length from this node to its parent is
        // the original branch length plus the additional variance.
        effective_length[node_id] += v_additional;
    }

    // Root values are the weighted averages at the root.
    let root_values: Vec<f64> = (0..n_chars)
        .map(|c| node_values[tree.root][c])
        .collect();

    // Compute correlation matrix if there are multiple characters.
    let correlations = if n_chars >= 2 {
        let mut corr_matrix = vec![vec![0.0; n_chars]; n_chars];
        for i in 0..n_chars {
            for j in 0..n_chars {
                corr_matrix[i][j] = pic_correlation_from_contrasts(&contrasts[i], &contrasts[j]);
            }
        }
        Some(corr_matrix)
    } else {
        None
    };

    Ok(ContrastResult {
        contrasts,
        variances,
        root_values,
        correlations,
    })
}

/// Compute the correlation between two characters using phylogenetically
/// independent contrasts.
///
/// Given a `ContrastResult`, extracts the standardized contrasts for the
/// two specified characters and computes their Pearson correlation:
///
/// ```text
/// r = sum(c1_i * c2_i) / sqrt(sum(c1_i^2) * sum(c2_i^2))
/// ```
///
/// Returns 0.0 if either character has zero variance in contrasts.
pub fn pic_correlation(contrasts: &ContrastResult, char1: usize, char2: usize) -> f64 {
    if char1 >= contrasts.contrasts.len() || char2 >= contrasts.contrasts.len() {
        return 0.0;
    }
    pic_correlation_from_contrasts(&contrasts.contrasts[char1], &contrasts.contrasts[char2])
}

/// Internal helper: compute Pearson correlation between two vectors of
/// standardized contrasts.
fn pic_correlation_from_contrasts(c1: &[f64], c2: &[f64]) -> f64 {
    assert_eq!(c1.len(), c2.len());
    let n = c1.len();
    if n == 0 {
        return 0.0;
    }

    let mut sum_prod = 0.0;
    let mut sum_sq1 = 0.0;
    let mut sum_sq2 = 0.0;

    for i in 0..n {
        sum_prod += c1[i] * c2[i];
        sum_sq1 += c1[i] * c1[i];
        sum_sq2 += c2[i] * c2[i];
    }

    let denom = (sum_sq1 * sum_sq2).sqrt();
    if denom < 1e-15 {
        return 0.0;
    }

    sum_prod / denom
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::comparative::ContinuousData;
    use crate::tree::newick::parse_newick;

    /// Helper to build ContinuousData from slices.
    fn make_data(taxa: &[&str], values: &[&[f64]]) -> ContinuousData {
        ContinuousData::new(
            taxa.iter().map(|s| s.to_string()).collect(),
            values.iter().map(|v| v.to_vec()).collect(),
        )
        .unwrap()
    }

    // ------------------------------------------------------------------
    // Test: Two taxa produce exactly one contrast
    // ------------------------------------------------------------------
    #[test]
    fn test_two_taxa_one_contrast() {
        // Tree: (A:1.0,B:1.0);
        let tree = parse_newick("(A:1.0,B:1.0);").unwrap();
        let data = make_data(&["A", "B"], &[&[5.0], &[3.0]]);

        let result = independent_contrasts(&tree, &data).unwrap();

        // One contrast for 2 taxa.
        assert_eq!(result.contrasts[0].len(), 1);
        assert_eq!(result.variances.len(), 1);

        // contrast = (5.0 - 3.0) / sqrt(1.0 + 1.0) = 2.0 / sqrt(2.0)
        let expected = 2.0 / (2.0_f64).sqrt();
        assert!(
            (result.contrasts[0][0] - expected).abs() < 1e-10,
            "expected contrast {}, got {}",
            expected,
            result.contrasts[0][0]
        );

        // Variance = 1.0 + 1.0 = 2.0
        assert!((result.variances[0] - 2.0).abs() < 1e-10);

        // Root value = weighted average = (5*1 + 3*1) / (1+1) = 4.0
        assert!((result.root_values[0] - 4.0).abs() < 1e-10);
    }

    // ------------------------------------------------------------------
    // Test: Star tree (all branch lengths equal)
    // ------------------------------------------------------------------
    #[test]
    fn test_star_tree_three_taxa() {
        // Star tree: ((A:1,B:1):0.0001,C:1); -- nearly star-like
        // A true star tree needs a trifurcation which violates binary requirement,
        // so we use a near-zero internal branch.
        let tree = parse_newick("((A:1.0,B:1.0):0.0001,C:1.0);").unwrap();
        let data = make_data(&["A", "B", "C"], &[&[2.0], &[4.0], &[6.0]]);

        let result = independent_contrasts(&tree, &data).unwrap();

        // Should produce 2 contrasts for 3 taxa.
        assert_eq!(result.contrasts[0].len(), 2);
        assert_eq!(result.variances.len(), 2);

        // First contrast (A vs B): (2 - 4) / sqrt(1 + 1) = -2/sqrt(2)
        let expected_1 = -2.0 / (2.0_f64).sqrt();
        assert!(
            (result.contrasts[0][0] - expected_1).abs() < 1e-6,
            "first contrast: expected {}, got {}",
            expected_1,
            result.contrasts[0][0]
        );
    }

    // ------------------------------------------------------------------
    // Test: Number of contrasts = n_taxa - 1
    // ------------------------------------------------------------------
    #[test]
    fn test_number_of_contrasts() {
        let tree = parse_newick("((A:1.0,B:1.0):0.5,(C:1.0,D:1.0):0.5);").unwrap();
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[1.0], &[2.0], &[3.0], &[4.0]],
        );

        let result = independent_contrasts(&tree, &data).unwrap();

        // 4 taxa -> 3 contrasts.
        assert_eq!(result.contrasts[0].len(), 3);
        assert_eq!(result.variances.len(), 3);
    }

    // ------------------------------------------------------------------
    // Test: 4 taxa on known tree, verify hand-calculated contrasts
    // ------------------------------------------------------------------
    #[test]
    fn test_four_taxa_hand_calculated() {
        // Tree: ((A:1,B:1):1,(C:2,D:2):1);
        //
        // Values: A=2, B=6, C=3, D=5
        //
        // First internal node (A,B):
        //   contrast = 2 - 6 = -4
        //   variance = 1 + 1 = 2
        //   standardized = -4/sqrt(2)
        //   weighted avg = (2*1 + 6*1) / 2 = 4
        //   v_additional = (1*1)/(1+1) = 0.5
        //   effective_length of AB-node = 1.0 + 0.5 = 1.5
        //
        // Second internal node (C,D):
        //   contrast = 3 - 5 = -2
        //   variance = 2 + 2 = 4
        //   standardized = -2/sqrt(4) = -1
        //   weighted avg = (3*2 + 5*2) / 4 = 4
        //   v_additional = (2*2)/(2+2) = 1
        //   effective_length of CD-node = 1.0 + 1.0 = 2.0
        //
        // Root node (AB vs CD):
        //   contrast = 4 - 4 = 0
        //   variance = 1.5 + 2.0 = 3.5
        //   standardized = 0/sqrt(3.5) = 0
        //   weighted avg = (4*2 + 4*1.5) / 3.5 = 4
        //
        let tree = parse_newick("((A:1.0,B:1.0):1.0,(C:2.0,D:2.0):1.0);").unwrap();
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[2.0], &[6.0], &[3.0], &[5.0]],
        );

        let result = independent_contrasts(&tree, &data).unwrap();

        let c = &result.contrasts[0];
        let v = &result.variances;

        // Contrast 1: A vs B
        let expected_c1 = -4.0 / (2.0_f64).sqrt();
        assert!(
            (c[0] - expected_c1).abs() < 1e-10,
            "contrast 1: expected {}, got {}",
            expected_c1,
            c[0]
        );
        assert!((v[0] - 2.0).abs() < 1e-10, "variance 1: expected 2.0, got {}", v[0]);

        // Contrast 2: C vs D
        let expected_c2 = -2.0 / (4.0_f64).sqrt();
        assert!(
            (c[1] - expected_c2).abs() < 1e-10,
            "contrast 2: expected {}, got {}",
            expected_c2,
            c[1]
        );
        assert!((v[1] - 4.0).abs() < 1e-10, "variance 2: expected 4.0, got {}", v[1]);

        // Contrast 3: AB vs CD
        let expected_c3 = 0.0 / (3.5_f64).sqrt();
        assert!(
            (c[2] - expected_c3).abs() < 1e-10,
            "contrast 3: expected {}, got {}",
            expected_c3,
            c[2]
        );
        assert!((v[2] - 3.5).abs() < 1e-10, "variance 3: expected 3.5, got {}", v[2]);

        // Root value should be 4.0.
        assert!(
            (result.root_values[0] - 4.0).abs() < 1e-10,
            "root value: expected 4.0, got {}",
            result.root_values[0]
        );
    }

    // ------------------------------------------------------------------
    // Test: Two perfectly correlated characters -> PIC correlation ~ 1.0
    // ------------------------------------------------------------------
    #[test]
    fn test_perfectly_correlated_characters() {
        let tree = parse_newick("((A:1.0,B:1.0):0.5,(C:1.0,D:1.0):0.5);").unwrap();
        // Character 2 = 2 * character 1 + constant
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[1.0, 12.0], &[2.0, 14.0], &[3.0, 16.0], &[4.0, 18.0]],
        );

        let result = independent_contrasts(&tree, &data).unwrap();
        let r = pic_correlation(&result, 0, 1);

        assert!(
            (r - 1.0).abs() < 1e-10,
            "perfectly correlated characters should have r=1.0, got {}",
            r
        );
    }

    // ------------------------------------------------------------------
    // Test: Two independent characters -> PIC correlation ~ 0.0
    // ------------------------------------------------------------------
    #[test]
    fn test_independent_characters() {
        // Use a symmetric tree where the characters are designed to have
        // zero correlation through the contrasts.
        let tree = parse_newick("((A:1.0,B:1.0):0.5,(C:1.0,D:1.0):0.5);").unwrap();

        // Char 1: A=1, B=3, C=2, D=4 -> contrasts proportional to [-2, -2, 0]
        // Char 2: A=1, B=-1, C=1, D=-1 -> contrasts proportional to [2, 2, 0]
        // These give correlation = sum(-2*2, -2*2, 0*0) / sqrt(sum(4+4+0) * sum(4+4+0))
        //   = -8 / 8 = -1.0  (perfectly anti-correlated)
        //
        // Let's instead pick values where contrasts are orthogonal.
        // Char 1: A=1, B=3, C=1, D=3  -> contrast1=-2, contrast2=-2, contrast3=0
        // Char 2: A=1, B=-1, C=-1, D=1 -> contrast1=2, contrast2=-2, contrast3=0
        // Correlation = sum(-2*2, -2*(-2), 0*0) / sqrt(8 * 8)
        //             = (-4 + 4) / 8 = 0.0
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[1.0, 1.0], &[3.0, -1.0], &[1.0, -1.0], &[3.0, 1.0]],
        );

        let result = independent_contrasts(&tree, &data).unwrap();
        let r = pic_correlation(&result, 0, 1);

        assert!(
            r.abs() < 1e-10,
            "independent characters should have r=0.0, got {}",
            r
        );
    }

    // ------------------------------------------------------------------
    // Test: Root value is weighted average
    // ------------------------------------------------------------------
    #[test]
    fn test_root_value_weighted_average() {
        // With unequal branch lengths, root value should reflect
        // the inverse-variance weighting.
        let tree = parse_newick("(A:1.0,B:3.0);").unwrap();
        let data = make_data(&["A", "B"], &[&[10.0], &[20.0]]);

        let result = independent_contrasts(&tree, &data).unwrap();

        // Root value = (10*3 + 20*1) / (1+3) = (30+20)/4 = 12.5
        let expected = 12.5;
        assert!(
            (result.root_values[0] - expected).abs() < 1e-10,
            "root value: expected {}, got {}",
            expected,
            result.root_values[0]
        );
    }

    // ------------------------------------------------------------------
    // Test: Self-correlation is 1.0
    // ------------------------------------------------------------------
    #[test]
    fn test_self_correlation() {
        let tree = parse_newick("((A:1.0,B:1.0):0.5,(C:1.0,D:1.0):0.5);").unwrap();
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[1.0, 5.0], &[2.0, 6.0], &[3.0, 7.0], &[4.0, 8.0]],
        );

        let result = independent_contrasts(&tree, &data).unwrap();

        // Self-correlation should be 1.0.
        let r00 = pic_correlation(&result, 0, 0);
        let r11 = pic_correlation(&result, 1, 1);
        assert!(
            (r00 - 1.0).abs() < 1e-10,
            "self-correlation should be 1.0, got {}",
            r00
        );
        assert!(
            (r11 - 1.0).abs() < 1e-10,
            "self-correlation should be 1.0, got {}",
            r11
        );
    }

    // ------------------------------------------------------------------
    // Test: Correlation matrix is symmetric
    // ------------------------------------------------------------------
    #[test]
    fn test_correlation_matrix_symmetric() {
        let tree = parse_newick("((A:1.0,B:2.0):0.5,(C:1.5,D:0.5):1.0);").unwrap();
        let data = make_data(
            &["A", "B", "C", "D"],
            &[
                &[1.0, 5.0, 9.0],
                &[2.0, 7.0, 8.0],
                &[3.0, 4.0, 7.0],
                &[4.0, 6.0, 10.0],
            ],
        );

        let result = independent_contrasts(&tree, &data).unwrap();
        let corr = result.correlations.as_ref().unwrap();

        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (corr[i][j] - corr[j][i]).abs() < 1e-10,
                    "correlation matrix not symmetric: corr[{}][{}]={}, corr[{}][{}]={}",
                    i, j, corr[i][j], j, i, corr[j][i]
                );
            }
        }
    }

    // ------------------------------------------------------------------
    // Test: Contrasts invariant to tree root position
    // When re-rooting preserves all pairwise distances between leaves,
    // the sum of squared standardized contrasts (which equals the
    // Brownian motion log-likelihood kernel) is invariant.
    // ------------------------------------------------------------------
    #[test]
    fn test_contrasts_rootposition_independence() {
        // Tree 1: ((A:1,B:1):1,(C:1,D:1):1);
        //   root -> AB_int:1, root -> CD_int:1
        //   AB_int -> A:1, AB_int -> B:1
        //   CD_int -> C:1, CD_int -> D:1
        //
        // Pairwise distances: A-B=2, C-D=2, A-C=A-D=B-C=B-D=4
        //
        // Re-rooted on the A branch (preserving all pairwise distances):
        // Tree 2: (A:0.5,((C:1,D:1):2,B:1):0.5);
        //   new_root -> A:0.5, new_root -> int1:0.5
        //   int1 -> int2:2, int1 -> B:1
        //   int2 -> C:1, int2 -> D:1
        //
        // Distances in tree2:
        //   A-B = 0.5 + 0.5 + 1 = 2 (ok)
        //   C-D = 1 + 1 = 2 (ok)
        //   A-C = 0.5 + 0.5 + 2 + 1 = 4 (ok)
        //   B-C = 1 + 2 + 1 = 4 (ok)

        let tree1 = parse_newick("((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);").unwrap();
        let tree2 = parse_newick("(A:0.5,((C:1.0,D:1.0):2.0,B:1.0):0.5);").unwrap();

        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[2.0], &[6.0], &[3.0], &[5.0]],
        );

        let result1 = independent_contrasts(&tree1, &data).unwrap();
        let result2 = independent_contrasts(&tree2, &data).unwrap();

        // Both produce 3 contrasts.
        assert_eq!(result1.contrasts[0].len(), 3);
        assert_eq!(result2.contrasts[0].len(), 3);

        // The C vs D contrast should be identical in both trees since
        // the subtree containing C and D is unchanged.
        // In tree1 postorder: (A,B) first, then (C,D), then root.
        // In tree2 postorder: (C,D) first, then (B, CD_int), then root.
        let cd_contrast_1 = result1.contrasts[0][1]; // C vs D in tree1
        let cd_contrast_2 = result2.contrasts[0][0]; // C vs D in tree2

        assert!(
            (cd_contrast_1 - cd_contrast_2).abs() < 1e-10,
            "C-D contrast should be the same: {} vs {}",
            cd_contrast_1,
            cd_contrast_2
        );

        // The sum of squared standardized contrasts is the key
        // root-invariance property: sum(c_i^2) is invariant to re-rooting
        // when all pairwise distances are preserved.
        let ss1: f64 = result1.contrasts[0].iter().map(|c| c * c).sum();
        let ss2: f64 = result2.contrasts[0].iter().map(|c| c * c).sum();

        assert!(
            (ss1 - ss2).abs() < 1e-8,
            "sum of squared contrasts should be root-invariant: {} vs {}",
            ss1,
            ss2
        );
    }

    // ------------------------------------------------------------------
    // Test: Variances are positive
    // ------------------------------------------------------------------
    #[test]
    fn test_variances_positive() {
        let tree = parse_newick("((A:0.5,B:1.5):0.3,(C:0.7,D:2.0):0.8);").unwrap();
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[1.0], &[2.0], &[3.0], &[4.0]],
        );

        let result = independent_contrasts(&tree, &data).unwrap();

        for (i, &v) in result.variances.iter().enumerate() {
            assert!(v > 0.0, "variance {} should be positive, got {}", i, v);
        }
    }

    // ------------------------------------------------------------------
    // Test: Error on non-binary tree
    // ------------------------------------------------------------------
    #[test]
    fn test_error_on_nonbinary_tree() {
        // Trifurcating tree: (A:1,B:1,C:1);
        let tree = parse_newick("(A:1.0,B:1.0,C:1.0);").unwrap();
        let data = make_data(&["A", "B", "C"], &[&[1.0], &[2.0], &[3.0]]);

        let result = independent_contrasts(&tree, &data);
        assert!(result.is_err(), "should error on non-binary tree");
    }

    // ------------------------------------------------------------------
    // Test: Error on mismatched taxa
    // ------------------------------------------------------------------
    #[test]
    fn test_error_on_mismatched_taxa() {
        let tree = parse_newick("(A:1.0,B:1.0);").unwrap();
        let data = make_data(&["X", "Y"], &[&[1.0], &[2.0]]);

        let result = independent_contrasts(&tree, &data);
        assert!(result.is_err(), "should error on mismatched taxa");
    }

    // ------------------------------------------------------------------
    // Test: Multiple characters produce correct number of contrast vectors
    // ------------------------------------------------------------------
    #[test]
    fn test_multiple_characters() {
        let tree = parse_newick("((A:1.0,B:1.0):0.5,(C:1.0,D:1.0):0.5);").unwrap();
        let data = make_data(
            &["A", "B", "C", "D"],
            &[
                &[1.0, 10.0, 100.0],
                &[2.0, 20.0, 200.0],
                &[3.0, 30.0, 300.0],
                &[4.0, 40.0, 400.0],
            ],
        );

        let result = independent_contrasts(&tree, &data).unwrap();

        // 3 characters, each with 3 contrasts (4 taxa - 1).
        assert_eq!(result.contrasts.len(), 3);
        for c in 0..3 {
            assert_eq!(result.contrasts[c].len(), 3);
        }
        assert_eq!(result.root_values.len(), 3);
    }

    // ------------------------------------------------------------------
    // Test: Five taxa tree
    // ------------------------------------------------------------------
    #[test]
    fn test_five_taxa() {
        let tree = parse_newick("((A:1.0,B:1.0):0.5,((C:1.0,D:1.0):0.5,E:1.5):0.5);").unwrap();
        let data = make_data(
            &["A", "B", "C", "D", "E"],
            &[&[1.0], &[2.0], &[3.0], &[4.0], &[5.0]],
        );

        let result = independent_contrasts(&tree, &data).unwrap();

        // 5 taxa -> 4 contrasts.
        assert_eq!(result.contrasts[0].len(), 4);
    }

    // ------------------------------------------------------------------
    // Test: pic_correlation with out-of-bounds index returns 0.0
    // ------------------------------------------------------------------
    #[test]
    fn test_pic_correlation_out_of_bounds() {
        let tree = parse_newick("(A:1.0,B:1.0);").unwrap();
        let data = make_data(&["A", "B"], &[&[1.0], &[2.0]]);
        let result = independent_contrasts(&tree, &data).unwrap();

        assert_eq!(pic_correlation(&result, 0, 5), 0.0);
        assert_eq!(pic_correlation(&result, 5, 0), 0.0);
    }

    // ------------------------------------------------------------------
    // Test: Negative correlation
    // ------------------------------------------------------------------
    #[test]
    fn test_negative_correlation() {
        let tree = parse_newick("((A:1.0,B:1.0):0.5,(C:1.0,D:1.0):0.5);").unwrap();
        // Character 2 = -2 * character 1 + constant
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[1.0, 8.0], &[2.0, 6.0], &[3.0, 4.0], &[4.0, 2.0]],
        );

        let result = independent_contrasts(&tree, &data).unwrap();
        let r = pic_correlation(&result, 0, 1);

        assert!(
            (r - (-1.0)).abs() < 1e-10,
            "anti-correlated characters should have r=-1.0, got {}",
            r
        );
    }

    // ------------------------------------------------------------------
    // Test: Unequal branch lengths affect root value
    // ------------------------------------------------------------------
    #[test]
    fn test_unequal_branches_root_value() {
        // If branch to A is very short and branch to B is very long,
        // root value should be closer to A's value.
        let tree = parse_newick("(A:0.01,B:100.0);").unwrap();
        let data = make_data(&["A", "B"], &[&[10.0], &[20.0]]);

        let result = independent_contrasts(&tree, &data).unwrap();

        // Root = (10 * 100 + 20 * 0.01) / (0.01 + 100)
        //      = (1000 + 0.2) / 100.01 ≈ 10.002
        assert!(
            (result.root_values[0] - 10.0).abs() < 0.01,
            "root should be close to A's value (10), got {}",
            result.root_values[0]
        );
    }
}

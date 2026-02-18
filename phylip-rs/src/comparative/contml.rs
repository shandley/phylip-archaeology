//! Maximum likelihood estimation under Brownian motion (Felsenstein, 1973, 1981).
//!
//! Given continuous character data and a phylogenetic tree, the Brownian motion
//! model assumes characters evolve as independent random walks along branches.
//! The data at the tips follows a multivariate normal distribution whose
//! covariance structure is determined by the tree topology and branch lengths.
//!
//! Rather than constructing and inverting the full n x n variance-covariance
//! matrix, we use the computationally efficient independent contrasts approach.
//! The log-likelihood for a single character is:
//!
//! ```text
//! ln L = -n/2 * ln(2 * pi) - (n-1)/2 * ln(sigma^2)
//!        - 1/2 * sum(ln(v_i)) - 1/(2 * sigma^2) * sum(c_i^2)
//! ```
//!
//! where `c_i` are the unstandardized contrasts and `v_i` are their variances
//! (sum of branch lengths at each contrast). The ML estimate of the
//! evolutionary rate is:
//!
//! ```text
//! sigma^2_ML = 1/(n-1) * sum(c_i^2 / v_i)
//!            = 1/(n-1) * sum(standardized_contrast_i^2)
//! ```
//!
//! For tree search, we implement stepwise addition followed by subtree
//! pruning and regrafting (SPR) rearrangements.
//!
//! # References
//!
//! - Felsenstein, J. (1973). Maximum-likelihood estimation of evolutionary
//!   trees from continuous characters. *American Journal of Human Genetics*,
//!   25, 471-492.
//! - Felsenstein, J. (1981). Continuous-character models. In *Inferring
//!   Phylogenies*. Sinauer Associates. Chapter 23.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.

use crate::tree::Tree;

use super::contrasts::independent_contrasts;
use super::ContinuousData;

/// Result of maximum likelihood analysis under Brownian motion.
#[derive(Debug, Clone)]
pub struct ContmlResult {
    /// The estimated ML tree.
    pub tree: Tree,
    /// Total log-likelihood of the data on the tree.
    pub lnl: f64,
    /// ML estimate of the evolutionary rate (sigma^2) for each character.
    pub sigma_squared: Vec<f64>,
    /// Estimated ancestral root value for each character.
    pub root_values: Vec<f64>,
}

/// The constant ln(2*pi).
const LN_2PI: f64 = 1.8378770664093453;

/// Minimum branch length for numerical stability.
const MIN_BRANCH_LENGTH: f64 = 1e-8;

/// Default branch length for initial tree construction.
const DEFAULT_INIT_BRANCH_LENGTH: f64 = 0.1;

/// Compute the Brownian motion log-likelihood for continuous character data
/// on a given tree.
///
/// Uses the independent contrasts approach, which is O(n) rather than the
/// O(n^3) direct matrix approach.
///
/// The log-likelihood is summed across all characters, each assumed to evolve
/// independently under Brownian motion with its own rate sigma^2.
///
/// # Errors
///
/// Returns an error if independent contrasts cannot be computed (e.g.,
/// non-binary tree, taxa mismatch).
pub fn brownian_log_likelihood(
    tree: &Tree,
    data: &ContinuousData,
) -> Result<f64, String> {
    let cr = independent_contrasts(tree, data)?;
    let n = data.n_taxa();
    let n_contrasts = n - 1; // number of contrasts

    if n_contrasts == 0 {
        return Err("need at least 2 taxa".to_string());
    }

    let mut total_lnl = 0.0;

    for char_idx in 0..data.n_chars {
        let contrasts = &cr.contrasts[char_idx];
        let variances = &cr.variances;

        // Compute sigma^2 ML estimate for this character.
        // sigma^2 = (1/(n-1)) * sum(standardized_contrast_i^2)
        // where standardized contrasts are already in cr.contrasts
        let mut sum_sq_std = 0.0;
        for &sc in contrasts {
            sum_sq_std += sc * sc;
        }
        let sigma2 = sum_sq_std / n_contrasts as f64;

        if sigma2 <= 0.0 {
            // All contrasts are zero (constant character). The likelihood
            // is technically infinite (or concentrated at a single point).
            // We assign a large but finite log-likelihood.
            // This avoids -inf and nan.
            total_lnl += 0.0;
            continue;
        }

        // Log-likelihood formula:
        // ln L = -(n-1)/2 * ln(2*pi) - (n-1)/2 * ln(sigma^2)
        //        - 1/2 * sum(ln(v_i)) - 1/(2*sigma^2) * sum(c_i^2)
        //
        // Note: c_i (unstandardized) = standardized_i * sqrt(v_i)
        // So sum(c_i^2) = sum(standardized_i^2 * v_i)
        //
        // And sum(c_i^2 / v_i) = sum(standardized_i^2) = (n-1) * sigma^2
        // So 1/(2*sigma^2) * sum(c_i^2) = 1/(2*sigma^2) * sum(std_i^2 * v_i)

        let nc = n_contrasts as f64;

        // Term 1: -(n-1)/2 * ln(2*pi)
        let term1 = -nc / 2.0 * LN_2PI;

        // Term 2: -(n-1)/2 * ln(sigma^2)
        let term2 = -nc / 2.0 * sigma2.ln();

        // Term 3: -1/2 * sum(ln(v_i))
        let mut sum_ln_v = 0.0;
        for &v in variances {
            sum_ln_v += v.ln();
        }
        let term3 = -0.5 * sum_ln_v;

        // Term 4: -1/(2*sigma^2) * sum(c_i^2)
        // c_i = standardized_i * sqrt(v_i), so c_i^2 = std_i^2 * v_i
        let mut sum_c2 = 0.0;
        for (i, &sc) in contrasts.iter().enumerate() {
            sum_c2 += sc * sc * variances[i];
        }
        let term4 = -sum_c2 / (2.0 * sigma2);

        total_lnl += term1 + term2 + term3 + term4;
    }

    Ok(total_lnl)
}

/// Estimate the evolutionary rate sigma^2 for each character using
/// maximum likelihood on the given tree.
///
/// sigma^2_ML = (1/(n-1)) * sum(standardized_contrast_i^2)
///
/// # Errors
///
/// Returns an error if contrasts cannot be computed.
pub fn estimate_rates(
    tree: &Tree,
    data: &ContinuousData,
) -> Result<Vec<f64>, String> {
    let cr = independent_contrasts(tree, data)?;
    let n_contrasts = data.n_taxa() - 1;

    if n_contrasts == 0 {
        return Err("need at least 2 taxa".to_string());
    }

    let mut rates = Vec::with_capacity(data.n_chars);

    for char_idx in 0..data.n_chars {
        let contrasts = &cr.contrasts[char_idx];
        let mut sum_sq = 0.0;
        for &sc in contrasts {
            sum_sq += sc * sc;
        }
        rates.push(sum_sq / n_contrasts as f64);
    }

    Ok(rates)
}

/// Search for the maximum likelihood tree under Brownian motion.
///
/// Uses stepwise addition to build an initial tree, then attempts
/// SPR rearrangements to improve the log-likelihood.
///
/// # Arguments
///
/// * `data` - Continuous character data for the taxa
/// * `seed` - Optional seed for deterministic taxon addition order
///
/// # Returns
///
/// A `ContmlResult` with the best tree found, its log-likelihood,
/// rate estimates, and root values.
pub fn contml_search(
    data: &ContinuousData,
    seed: Option<u64>,
) -> Result<ContmlResult, String> {
    let n = data.n_taxa();
    if n < 3 {
        return Err("need at least 3 taxa for tree search".to_string());
    }

    // Determine taxon addition order.
    let order = addition_order(n, seed);

    // Step 1: Build initial tree by stepwise addition.
    let mut tree = stepwise_addition(data, &order)?;

    // Step 2: Optimize branch lengths on the initial tree.
    optimize_all_branch_lengths(&mut tree, data)?;

    // Step 3: Try SPR rearrangements to improve the tree.
    spr_search(&mut tree, data)?;

    // Final optimization of branch lengths.
    optimize_all_branch_lengths(&mut tree, data)?;

    // Compute final results.
    let lnl = brownian_log_likelihood(&tree, data)?;
    let sigma_squared = estimate_rates(&tree, data)?;
    let cr = independent_contrasts(&tree, data)?;

    Ok(ContmlResult {
        tree,
        lnl,
        sigma_squared,
        root_values: cr.root_values,
    })
}

/// Generate the taxon addition order. If a seed is provided, use a
/// simple deterministic permutation; otherwise, use sequential order.
fn addition_order(n: usize, seed: Option<u64>) -> Vec<usize> {
    let mut order: Vec<usize> = (0..n).collect();
    if let Some(s) = seed {
        // Simple deterministic shuffle using a linear congruential generator.
        let mut rng = s;
        for i in (1..n).rev() {
            // LCG: next = (a * rng + c) mod m
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let j = (rng >> 33) as usize % (i + 1);
            order.swap(i, j);
        }
    }
    order
}

/// Build an initial binary tree by stepwise addition of taxa.
///
/// Starts with a 2-taxon tree, then adds each subsequent taxon
/// at the position that maximizes the log-likelihood.
fn stepwise_addition(
    data: &ContinuousData,
    order: &[usize],
) -> Result<Tree, String> {
    let n = order.len();
    if n < 2 {
        return Err("need at least 2 taxa".to_string());
    }

    // Start with a 2-taxon tree: root -> (taxon0, taxon1).
    let mut tree = Tree::new();
    tree.rooted = true;

    let root = tree.add_node(None, None, None);
    tree.root = root;

    let t0 = order[0];
    let t1 = order[1];
    tree.add_node(
        Some(data.taxa[t0].clone()),
        Some(DEFAULT_INIT_BRANCH_LENGTH),
        Some(root),
    );
    tree.add_node(
        Some(data.taxa[t1].clone()),
        Some(DEFAULT_INIT_BRANCH_LENGTH),
        Some(root),
    );

    // Add taxa one at a time.
    for &taxon_idx in &order[2..] {
        let taxon_name = data.taxa[taxon_idx].clone();
        add_taxon_best_position(&mut tree, data, &taxon_name)?;
    }

    Ok(tree)
}

/// Add a taxon to the tree at the position that maximizes the Brownian
/// motion log-likelihood. Tries inserting on each existing branch.
fn add_taxon_best_position(
    tree: &mut Tree,
    data: &ContinuousData,
    taxon_name: &str,
) -> Result<(), String> {
    // Collect all edges (parent, child) where we could insert.
    let edges: Vec<(usize, usize)> = tree
        .nodes
        .iter()
        .filter(|n| n.parent.is_some())
        .map(|n| (n.parent.unwrap(), n.id))
        .collect();

    if edges.is_empty() {
        return Err("no edges available for insertion".to_string());
    }

    let mut best_lnl = f64::NEG_INFINITY;
    let mut best_tree = tree.clone();

    for &(parent_id, child_id) in &edges {
        let mut candidate = tree.clone();

        // Insert the new taxon on the branch (parent_id -> child_id).
        insert_on_branch(
            &mut candidate,
            parent_id,
            child_id,
            taxon_name,
            DEFAULT_INIT_BRANCH_LENGTH,
        );

        // Evaluate likelihood.
        if let Ok(lnl) = brownian_log_likelihood(&candidate, data) {
            if lnl > best_lnl {
                best_lnl = lnl;
                best_tree = candidate;
            }
        }
    }

    *tree = best_tree;
    Ok(())
}

/// Insert a new taxon on the branch between parent_id and child_id.
///
/// Creates a new internal node that becomes the parent of both the
/// existing child and the new leaf.
fn insert_on_branch(
    tree: &mut Tree,
    parent_id: usize,
    child_id: usize,
    taxon_name: &str,
    branch_length: f64,
) {
    let child_bl = tree.nodes[child_id]
        .branch_length
        .unwrap_or(DEFAULT_INIT_BRANCH_LENGTH);
    let half_bl = child_bl / 2.0;

    // Remove child from parent's children list.
    tree.nodes[parent_id]
        .children
        .retain(|&c| c != child_id);

    // Create new internal node as child of parent.
    let new_internal = tree.add_node(None, Some(half_bl), Some(parent_id));

    // Reattach the original child to the new internal node.
    tree.nodes[child_id].parent = Some(new_internal);
    tree.nodes[child_id].branch_length = Some(half_bl.max(MIN_BRANCH_LENGTH));
    tree.nodes[new_internal].children.push(child_id);

    // Add new leaf as child of the new internal node.
    tree.add_node(
        Some(taxon_name.to_string()),
        Some(branch_length),
        Some(new_internal),
    );
}

/// Optimize all branch lengths in the tree to maximize the Brownian
/// motion log-likelihood.
///
/// Uses coordinate ascent: optimize one branch at a time, cycling
/// through all branches until convergence.
fn optimize_all_branch_lengths(
    tree: &mut Tree,
    data: &ContinuousData,
) -> Result<f64, String> {
    let max_rounds = 20;
    let tolerance = 1e-6;

    let mut prev_lnl = brownian_log_likelihood(tree, data)?;

    for _ in 0..max_rounds {
        // Get all non-root nodes with branch lengths.
        let node_ids: Vec<usize> = tree
            .nodes
            .iter()
            .filter(|n| n.parent.is_some())
            .map(|n| n.id)
            .collect();

        for node_id in node_ids {
            optimize_one_branch(tree, node_id, data)?;
        }

        let new_lnl = brownian_log_likelihood(tree, data)?;
        if (new_lnl - prev_lnl).abs() < tolerance {
            return Ok(new_lnl);
        }
        prev_lnl = new_lnl;
    }

    Ok(prev_lnl)
}

/// Optimize a single branch length using golden section search.
fn optimize_one_branch(
    tree: &mut Tree,
    node_id: usize,
    data: &ContinuousData,
) -> Result<f64, String> {
    let max_bl = 10.0;
    let tol = 1e-6;
    let phi = (5.0_f64.sqrt() - 1.0) / 2.0; // golden ratio conjugate

    let mut a = MIN_BRANCH_LENGTH;
    let mut b = max_bl;

    let mut c = b - phi * (b - a);
    let mut d = a + phi * (b - a);

    for _ in 0..50 {
        tree.nodes[node_id].branch_length = Some(c);
        let fc = brownian_log_likelihood(tree, data).unwrap_or(f64::NEG_INFINITY);

        tree.nodes[node_id].branch_length = Some(d);
        let fd = brownian_log_likelihood(tree, data).unwrap_or(f64::NEG_INFINITY);

        if fc > fd {
            b = d;
        } else {
            a = c;
        }

        if (b - a).abs() < tol {
            break;
        }

        c = b - phi * (b - a);
        d = a + phi * (b - a);
    }

    let best_bl = (a + b) / 2.0;
    tree.nodes[node_id].branch_length = Some(best_bl.max(MIN_BRANCH_LENGTH));

    Ok(best_bl)
}

/// Attempt SPR (Subtree Pruning and Regrafting) moves to improve the tree.
///
/// For each subtree, try pruning it and reattaching at every other branch.
/// Accept the move if it improves the log-likelihood.
fn spr_search(
    tree: &mut Tree,
    data: &ContinuousData,
) -> Result<bool, String> {
    let mut improved = true;
    let mut any_improved = false;

    while improved {
        improved = false;
        let current_lnl = brownian_log_likelihood(tree, data)?;

        // Collect candidate subtrees (non-root nodes that are not
        // children of the root in a 2-child root).
        let candidates: Vec<usize> = tree
            .nodes
            .iter()
            .filter(|n| {
                n.parent.is_some() && n.id != tree.root
            })
            .map(|n| n.id)
            .collect();

        for prune_id in &candidates {
            let prune_id = *prune_id;
            let parent_id = match tree.nodes[prune_id].parent {
                Some(p) => p,
                None => continue,
            };

            // Don't prune if it would leave the parent with < 1 child
            // (besides the pruned subtree).
            if tree.nodes[parent_id].children.len() < 2 {
                continue;
            }

            // Get the sibling(s) of the pruned subtree.
            let siblings: Vec<usize> = tree.nodes[parent_id]
                .children
                .iter()
                .copied()
                .filter(|&c| c != prune_id)
                .collect();

            if siblings.is_empty() {
                continue;
            }

            // Collect all possible regraft edges (excluding those in
            // the pruned subtree).
            let subtree_nodes = collect_subtree_nodes(tree, prune_id);
            let regraft_edges: Vec<(usize, usize)> = tree
                .nodes
                .iter()
                .filter(|n| {
                    n.parent.is_some()
                        && !subtree_nodes.contains(&n.id)
                        && n.id != prune_id
                        // Don't regraft to the same position
                        && n.parent.unwrap() != parent_id
                })
                .map(|n| (n.parent.unwrap(), n.id))
                .collect();

            for &(regraft_parent, regraft_child) in &regraft_edges {
                let mut candidate = tree.clone();

                // Prune: remove subtree from current position.
                if !prune_subtree(&mut candidate, prune_id) {
                    continue;
                }

                // Regraft: insert subtree at new position.
                let prune_bl = candidate.nodes[prune_id]
                    .branch_length
                    .unwrap_or(DEFAULT_INIT_BRANCH_LENGTH);
                insert_subtree_on_branch(
                    &mut candidate,
                    prune_id,
                    regraft_parent,
                    regraft_child,
                    prune_bl,
                );

                // Evaluate.
                if let Ok(new_lnl) = brownian_log_likelihood(&candidate, data) {
                    if new_lnl > current_lnl + 1e-8 {
                        *tree = candidate;
                        improved = true;
                        any_improved = true;
                        break; // restart with the improved tree
                    }
                }
            }

            if improved {
                break; // restart outer loop
            }
        }
    }

    Ok(any_improved)
}

/// Collect all node IDs in the subtree rooted at `node_id`.
fn collect_subtree_nodes(tree: &Tree, node_id: usize) -> Vec<usize> {
    let mut result = Vec::new();
    let mut stack = vec![node_id];
    while let Some(id) = stack.pop() {
        result.push(id);
        for &child in &tree.nodes[id].children {
            stack.push(child);
        }
    }
    result
}

/// Prune a subtree from its current position in the tree.
/// If the parent becomes a unary node, collapse it.
/// Returns false if the prune is not possible.
fn prune_subtree(tree: &mut Tree, subtree_id: usize) -> bool {
    let parent_id = match tree.nodes[subtree_id].parent {
        Some(p) => p,
        None => return false,
    };

    // Remove subtree from parent's children.
    tree.nodes[parent_id]
        .children
        .retain(|&c| c != subtree_id);
    tree.nodes[subtree_id].parent = None;

    // If parent now has only one child, collapse it (suppress the node).
    if tree.nodes[parent_id].children.len() == 1 {
        let remaining_child = tree.nodes[parent_id].children[0];

        if let Some(grandparent) = tree.nodes[parent_id].parent {
            // Sum branch lengths through the collapsed node.
            let parent_bl = tree.nodes[parent_id]
                .branch_length
                .unwrap_or(0.0);
            let child_bl = tree.nodes[remaining_child]
                .branch_length
                .unwrap_or(0.0);
            tree.nodes[remaining_child].branch_length = Some(parent_bl + child_bl);

            // Reattach remaining child to grandparent.
            tree.nodes[remaining_child].parent = Some(grandparent);
            // Replace parent_id with remaining_child in grandparent's children list.
            if let Some(pos) = tree.nodes[grandparent]
                .children
                .iter()
                .position(|&c| c == parent_id)
            {
                tree.nodes[grandparent].children[pos] = remaining_child;
            }
            // Mark parent as disconnected.
            tree.nodes[parent_id].children.clear();
            tree.nodes[parent_id].parent = None;
        } else {
            // Parent is the root. Make remaining child the new root.
            tree.root = remaining_child;
            tree.nodes[remaining_child].parent = None;
            tree.nodes[remaining_child].branch_length = None;
            tree.nodes[parent_id].children.clear();
        }
    }

    true
}

/// Insert a subtree onto a branch between regraft_parent and regraft_child.
fn insert_subtree_on_branch(
    tree: &mut Tree,
    subtree_id: usize,
    regraft_parent: usize,
    regraft_child: usize,
    subtree_bl: f64,
) {
    let child_bl = tree.nodes[regraft_child]
        .branch_length
        .unwrap_or(DEFAULT_INIT_BRANCH_LENGTH);
    let half_bl = child_bl / 2.0;

    // Remove regraft_child from regraft_parent's children.
    tree.nodes[regraft_parent]
        .children
        .retain(|&c| c != regraft_child);

    // Create a new internal node.
    let new_internal = tree.add_node(None, Some(half_bl), Some(regraft_parent));

    // Reattach regraft_child to new internal node.
    tree.nodes[regraft_child].parent = Some(new_internal);
    tree.nodes[regraft_child].branch_length = Some(half_bl.max(MIN_BRANCH_LENGTH));
    tree.nodes[new_internal].children.push(regraft_child);

    // Attach subtree to new internal node.
    tree.nodes[subtree_id].parent = Some(new_internal);
    tree.nodes[subtree_id].branch_length = Some(subtree_bl);
    tree.nodes[new_internal].children.push(subtree_id);
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
    // Test: Log-likelihood is finite and negative
    // ------------------------------------------------------------------
    #[test]
    fn test_brownian_lnl_is_negative() {
        let tree = parse_newick("((A:1.0,B:1.0):0.5,(C:1.0,D:1.0):0.5);").unwrap();
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[1.0], &[2.0], &[3.0], &[4.0]],
        );

        let lnl = brownian_log_likelihood(&tree, &data).unwrap();
        assert!(lnl.is_finite(), "log-likelihood should be finite, got {}", lnl);
        assert!(lnl < 0.0, "log-likelihood should be negative, got {}", lnl);
    }

    // ------------------------------------------------------------------
    // Test: Two taxa log-likelihood
    // ------------------------------------------------------------------
    #[test]
    fn test_two_taxa_lnl() {
        let tree = parse_newick("(A:1.0,B:1.0);").unwrap();
        let data = make_data(&["A", "B"], &[&[0.0], &[1.0]]);

        let lnl = brownian_log_likelihood(&tree, &data).unwrap();
        assert!(lnl.is_finite());
        assert!(lnl < 0.0);

        // Hand calculation for 1 contrast:
        // contrast = (0 - 1) = -1, variance = 2
        // standardized = -1/sqrt(2)
        // sigma^2 = (1/sqrt(2))^2 / 1 = 0.5
        // ln L = -1/2 * ln(2*pi) - 1/2 * ln(0.5) - 1/2 * ln(2) - 1/(2*0.5) * 1
        //      = -0.9189 - (-0.3466) - 0.3466 - 1.0
        //      = -0.9189 + 0.3466 - 0.3466 - 1.0
        //      = -1.9189
        let expected = -0.5 * LN_2PI - 0.5 * (0.5_f64).ln() - 0.5 * (2.0_f64).ln() - 1.0;
        assert!(
            (lnl - expected).abs() < 1e-8,
            "expected lnl {}, got {}",
            expected,
            lnl
        );
    }

    // ------------------------------------------------------------------
    // Test: Rate estimates are positive
    // ------------------------------------------------------------------
    #[test]
    fn test_rates_positive() {
        let tree = parse_newick("((A:1.0,B:1.0):0.5,(C:1.0,D:1.0):0.5);").unwrap();
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[1.0, 10.0], &[2.0, 20.0], &[3.0, 30.0], &[4.0, 40.0]],
        );

        let rates = estimate_rates(&tree, &data).unwrap();
        assert_eq!(rates.len(), 2);
        for (i, &r) in rates.iter().enumerate() {
            assert!(r > 0.0, "rate {} should be positive, got {}", i, r);
        }
    }

    // ------------------------------------------------------------------
    // Test: Known sigma^2 recovery
    // ------------------------------------------------------------------
    #[test]
    fn test_known_sigma_recovery() {
        // Construct data consistent with Brownian motion on a star-like tree.
        // With equal branch lengths = 1.0 and known sigma^2 = 1.0, the
        // expected variance of contrasts is sigma^2 * (v_L + v_R) / (v_L + v_R) = 1.0
        // (after standardization).
        //
        // For a simple symmetric tree with moderate data, the ML estimate
        // should be reasonably close to the true value.

        let tree = parse_newick("((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);").unwrap();

        // Generate data where differences are consistent with sigma^2 = 1.0
        // on branches of length 1.0 (so expected variance per branch = 1.0).
        // A simple test: if sigma^2 = 1 and branch length = 1,
        // then the standard deviation of change = 1.
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[0.0], &[1.0], &[0.5], &[1.5]],
        );

        let rates = estimate_rates(&tree, &data).unwrap();
        assert!(rates[0] > 0.0);
        // Rate should be a reasonable positive value (exact value depends on
        // the specific data, but should not be wildly off).
        assert!(rates[0] < 10.0, "rate estimate {} seems too large", rates[0]);
    }

    // ------------------------------------------------------------------
    // Test: Variance proportional to time
    // ------------------------------------------------------------------
    #[test]
    fn test_variance_proportional_to_time() {
        // If we double all branch lengths, the estimated sigma^2 should
        // be approximately halved (since sigma^2 is per unit branch length,
        // and doubling time with same data means less rate per unit time).

        let tree1 = parse_newick("((A:1.0,B:1.0):0.5,(C:1.0,D:1.0):0.5);").unwrap();
        let tree2 = parse_newick("((A:2.0,B:2.0):1.0,(C:2.0,D:2.0):1.0);").unwrap();

        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[1.0], &[3.0], &[2.0], &[4.0]],
        );

        let rate1 = estimate_rates(&tree1, &data).unwrap()[0];
        let rate2 = estimate_rates(&tree2, &data).unwrap()[0];

        // rate2 should be approximately rate1 / 2
        let ratio = rate1 / rate2;
        assert!(
            (ratio - 2.0).abs() < 0.5,
            "rate ratio should be ~2.0 when branches are doubled, got {}",
            ratio
        );
    }

    // ------------------------------------------------------------------
    // Test: Higher likelihood for correct tree than random tree
    // ------------------------------------------------------------------
    #[test]
    fn test_correct_tree_higher_likelihood() {
        // Data generated with strong signal for ((A,B),(C,D))
        let good_tree = parse_newick("((A:0.5,B:0.5):1.0,(C:0.5,D:0.5):1.0);").unwrap();
        let bad_tree = parse_newick("((A:0.5,C:0.5):1.0,(B:0.5,D:0.5):1.0);").unwrap();

        // A and B are similar, C and D are similar, but (A,B) vs (C,D) are different
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[1.0], &[1.1], &[10.0], &[10.1]],
        );

        let lnl_good = brownian_log_likelihood(&good_tree, &data).unwrap();
        let lnl_bad = brownian_log_likelihood(&bad_tree, &data).unwrap();

        assert!(
            lnl_good > lnl_bad,
            "correct tree ({}) should have higher lnl than incorrect tree ({})",
            lnl_good,
            lnl_bad
        );
    }

    // ------------------------------------------------------------------
    // Test: Log-likelihood with multiple characters
    // ------------------------------------------------------------------
    #[test]
    fn test_multiple_characters_lnl() {
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

        let lnl = brownian_log_likelihood(&tree, &data).unwrap();
        assert!(lnl.is_finite());
        assert!(lnl < 0.0);

        // The lnl for 3 characters should be the sum of individual character lnls.
        let mut sum_lnl = 0.0;
        for c in 0..3 {
            let char_data = make_data(
                &["A", "B", "C", "D"],
                &[
                    &[data.values[0][c]],
                    &[data.values[1][c]],
                    &[data.values[2][c]],
                    &[data.values[3][c]],
                ],
            );
            sum_lnl += brownian_log_likelihood(&tree, &char_data).unwrap();
        }

        assert!(
            (lnl - sum_lnl).abs() < 1e-8,
            "multi-char lnl ({}) should equal sum of single-char lnls ({})",
            lnl,
            sum_lnl
        );
    }

    // ------------------------------------------------------------------
    // Test: contml_search produces a valid result
    // ------------------------------------------------------------------
    #[test]
    fn test_contml_search_basic() {
        let data = make_data(
            &["A", "B", "C", "D", "E"],
            &[
                &[1.0, 2.0],
                &[1.1, 2.1],
                &[5.0, 6.0],
                &[5.1, 6.1],
                &[3.0, 4.0],
            ],
        );

        let result = contml_search(&data, Some(42)).unwrap();

        assert!(result.lnl.is_finite(), "lnl should be finite");
        assert!(result.lnl < 0.0, "lnl should be negative");
        assert_eq!(result.sigma_squared.len(), 2);
        assert_eq!(result.root_values.len(), 2);

        // The tree should have 5 leaves.
        assert_eq!(result.tree.num_leaves(), 5);

        // All rates should be positive.
        for &r in &result.sigma_squared {
            assert!(r > 0.0, "rate should be positive, got {}", r);
        }
    }

    // ------------------------------------------------------------------
    // Test: contml_search with clear structure recovers groupings
    // ------------------------------------------------------------------
    #[test]
    fn test_contml_search_finds_structure() {
        // Two clear clusters: {A,B} near 0, {C,D} near 10
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[0.0], &[0.1], &[10.0], &[10.1]],
        );

        let result = contml_search(&data, Some(123)).unwrap();

        // The found tree should have better likelihood than the wrong grouping.
        let wrong_tree = parse_newick("((A:0.5,C:0.5):1.0,(B:0.5,D:0.5):1.0);").unwrap();
        let wrong_lnl = brownian_log_likelihood(&wrong_tree, &data).unwrap();

        assert!(
            result.lnl >= wrong_lnl,
            "search result lnl ({}) should be >= wrong tree lnl ({})",
            result.lnl,
            wrong_lnl
        );
    }

    // ------------------------------------------------------------------
    // Test: contml_search deterministic with same seed
    // ------------------------------------------------------------------
    #[test]
    fn test_contml_search_deterministic() {
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[1.0], &[2.0], &[5.0], &[6.0]],
        );

        let result1 = contml_search(&data, Some(42)).unwrap();
        let result2 = contml_search(&data, Some(42)).unwrap();

        assert!(
            (result1.lnl - result2.lnl).abs() < 1e-8,
            "same seed should give same result: {} vs {}",
            result1.lnl,
            result2.lnl
        );
    }

    // ------------------------------------------------------------------
    // Test: Optimize branch lengths improves likelihood
    // ------------------------------------------------------------------
    #[test]
    fn test_optimize_improves_lnl() {
        let mut tree = parse_newick("((A:5.0,B:5.0):5.0,(C:5.0,D:5.0):5.0);").unwrap();
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[0.0], &[0.1], &[10.0], &[10.1]],
        );

        let lnl_before = brownian_log_likelihood(&tree, &data).unwrap();
        optimize_all_branch_lengths(&mut tree, &data).unwrap();
        let lnl_after = brownian_log_likelihood(&tree, &data).unwrap();

        assert!(
            lnl_after >= lnl_before - 1e-6,
            "optimization should not decrease lnl: {} -> {}",
            lnl_before,
            lnl_after
        );
    }

    // ------------------------------------------------------------------
    // Test: Error on too few taxa
    // ------------------------------------------------------------------
    #[test]
    fn test_contml_search_error_few_taxa() {
        let data = make_data(&["A", "B"], &[&[1.0], &[2.0]]);
        let result = contml_search(&data, None);
        assert!(result.is_err());
    }

    // ------------------------------------------------------------------
    // Test: Constant character has zero rate
    // ------------------------------------------------------------------
    #[test]
    fn test_constant_character_zero_rate() {
        let tree = parse_newick("((A:1.0,B:1.0):0.5,(C:1.0,D:1.0):0.5);").unwrap();
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[5.0], &[5.0], &[5.0], &[5.0]],
        );

        let rates = estimate_rates(&tree, &data).unwrap();
        assert!(
            rates[0].abs() < 1e-15,
            "constant character should have rate 0, got {}",
            rates[0]
        );
    }

    // ------------------------------------------------------------------
    // Test: Log-likelihood on 3-taxon tree
    // ------------------------------------------------------------------
    #[test]
    fn test_three_taxa_lnl() {
        let tree = parse_newick("((A:1.0,B:1.0):1.0,C:2.0);").unwrap();
        let data = make_data(&["A", "B", "C"], &[&[1.0], &[3.0], &[5.0]]);

        let lnl = brownian_log_likelihood(&tree, &data).unwrap();
        assert!(lnl.is_finite());
        assert!(lnl < 0.0);
    }

    // ------------------------------------------------------------------
    // Test: contml_search with 6 taxa
    // ------------------------------------------------------------------
    #[test]
    fn test_contml_search_six_taxa() {
        let data = make_data(
            &["A", "B", "C", "D", "E", "F"],
            &[
                &[1.0, 1.0],
                &[1.5, 1.5],
                &[5.0, 5.0],
                &[5.5, 5.5],
                &[9.0, 9.0],
                &[9.5, 9.5],
            ],
        );

        let result = contml_search(&data, Some(7)).unwrap();

        assert!(result.lnl.is_finite());
        assert_eq!(result.tree.num_leaves(), 6);
    }

    // ------------------------------------------------------------------
    // Test: Branch length optimization converges
    // ------------------------------------------------------------------
    #[test]
    fn test_branch_optimization_convergence() {
        let mut tree = parse_newick("((A:0.001,B:0.001):0.001,(C:0.001,D:0.001):0.001);").unwrap();
        let data = make_data(
            &["A", "B", "C", "D"],
            &[&[0.0], &[10.0], &[0.0], &[10.0]],
        );

        let lnl = optimize_all_branch_lengths(&mut tree, &data).unwrap();
        assert!(lnl.is_finite());

        // After optimization, branches connecting different taxa should be longer
        // than the initial 0.001.
        let leaves: Vec<usize> = tree
            .nodes
            .iter()
            .filter(|n| n.is_leaf())
            .map(|n| n.id)
            .collect();
        let mut has_longer = false;
        for &leaf_id in &leaves {
            if let Some(bl) = tree.nodes[leaf_id].branch_length {
                if bl > 0.01 {
                    has_longer = true;
                }
            }
        }
        assert!(has_longer, "some branches should have grown from 0.001");
    }
}

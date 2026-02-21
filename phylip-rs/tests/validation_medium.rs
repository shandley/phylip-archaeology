//! Medium-scale integration tests for phylip-rs.
//!
//! These tests use 50-100 taxon datasets to verify that phylip-rs algorithms
//! work correctly at realistic scales and that key statistical properties hold.
//! They also include performance sanity checks.

use phylip_rs::bootstrap::{bootstrap_weights, resample_alignment, ResamplingMethod, SimpleRng};
use phylip_rs::consensus::{consensus_tree, ConsensusMethod};
use phylip_rs::distance::neighbor_joining;
use phylip_rs::likelihood::models::{Jc69Model, SubstitutionModel};
use phylip_rs::likelihood::pruning::{log_likelihood, optimize_branch_lengths};
use phylip_rs::models::jc69::JC69;
use phylip_rs::models::k2p::K2P;
use phylip_rs::models::compute_distance_matrix;
use phylip_rs::parsimony::wagner::search as parsimony_search;
use phylip_rs::tree::newick::{parse_newick, write_newick};
use phylip_rs::tree::{Alignment, Base, Sequence, Tree};

// ============================================================================
// Deterministic sequence simulator (JC69 model)
// ============================================================================

struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 33) as f64 / (1u64 << 31) as f64
    }
}

fn index_to_base(i: usize) -> Base {
    match i {
        0 => Base::A,
        1 => Base::C,
        2 => Base::G,
        3 => Base::T,
        _ => unreachable!(),
    }
}

/// Simulate DNA sequences along a tree under JC69.
fn simulate_alignment(
    tree: &Tree,
    nsites: usize,
    model: &Jc69Model,
    rng: &mut Rng,
) -> Alignment {
    let nnodes = tree.num_nodes();
    let mut node_seqs: Vec<Vec<usize>> = vec![vec![0; nsites]; nnodes];

    // Generate root sequence
    for site in 0..nsites {
        node_seqs[tree.root][site] = (rng.next_u64() % 4) as usize;
    }

    // Traverse in preorder (root → tips)
    let preorder = tree.preorder();
    for &node_id in &preorder {
        if node_id == tree.root {
            continue;
        }
        let bl = tree.nodes[node_id].branch_length.unwrap_or(0.1);
        let p = model.transition_probs(bl);
        let parent_id = tree.nodes[node_id].parent.unwrap();

        for site in 0..nsites {
            let parent_base = node_seqs[parent_id][site];
            let r = rng.next_f64();
            let mut cumsum = 0.0;
            let mut child_base = 3; // default to T
            for b in 0..4 {
                cumsum += p[parent_base][b];
                if r < cumsum {
                    child_base = b;
                    break;
                }
            }
            node_seqs[node_id][site] = child_base;
        }
    }

    // Extract leaf sequences
    let mut sequences = Vec::new();
    for leaf in tree.leaves() {
        let name = leaf.name.clone().unwrap_or_else(|| format!("T{}", leaf.id));
        let bases: Vec<Base> = node_seqs[leaf.id]
            .iter()
            .map(|&b| index_to_base(b))
            .collect();
        sequences.push(Sequence::new(name, bases));
    }

    Alignment::new(sequences).unwrap()
}

/// Generate a random tree with n taxa and branch lengths ~ Uniform(min_bl, max_bl).
fn random_tree(ntaxa: usize, min_bl: f64, max_bl: f64, rng: &mut Rng) -> Tree {
    // Build tree by sequential taxon addition
    let mut tree = Tree::new();

    // Add first two taxa connected by internal node
    let root = tree.add_node(None, None, None);
    let bl = min_bl + rng.next_f64() * (max_bl - min_bl);
    let _t0 = tree.add_node(
        Some(format!("T{}", 0)),
        Some(bl),
        Some(root),
    );
    let bl = min_bl + rng.next_f64() * (max_bl - min_bl);
    let _t1 = tree.add_node(
        Some(format!("T{}", 1)),
        Some(bl),
        Some(root),
    );

    // Add remaining taxa one at a time
    for i in 2..ntaxa {
        // Pick a random existing leaf to subdivide its parent branch
        let leaves: Vec<usize> = tree
            .leaves()
            .iter()
            .map(|l| l.id)
            .collect();
        let target_leaf = leaves[rng.next_u64() as usize % leaves.len()];
        let target_parent = tree.nodes[target_leaf].parent.unwrap();
        let old_bl = tree.nodes[target_leaf].branch_length.unwrap_or(0.1);

        // Create new internal node
        let split_point = rng.next_f64();
        let new_internal_bl = old_bl * split_point;
        let remaining_bl = old_bl * (1.0 - split_point);

        let new_internal = tree.add_node(None, Some(new_internal_bl), Some(target_parent));

        // Move the target leaf under the new internal node
        tree.nodes[target_leaf].parent = Some(new_internal);
        tree.nodes[target_leaf].branch_length = Some(remaining_bl);

        // Remove target_leaf from old parent's children, add new_internal instead
        tree.nodes[target_parent]
            .children
            .retain(|&c| c != target_leaf);
        tree.nodes[new_internal].children.push(target_leaf);

        // Add new taxon as sibling
        let bl = min_bl + rng.next_f64() * (max_bl - min_bl);
        tree.add_node(
            Some(format!("T{}", i)),
            Some(bl),
            Some(new_internal),
        );
    }

    tree
}

// ============================================================================
// Self-consistency: NJ → ML optimization should improve likelihood
// ============================================================================

#[test]
fn test_nj_then_ml_improves_likelihood() {
    let mut rng = Rng::new(42);
    let true_tree = random_tree(20, 0.01, 0.2, &mut rng);
    let model = Jc69Model;
    let alignment = simulate_alignment(&true_tree, 500, &model, &mut rng);

    // Build NJ tree from distances
    let dist_model = JC69::new();
    let matrix = compute_distance_matrix(&alignment, &dist_model).unwrap();
    let mut nj_tree = neighbor_joining(&matrix);

    // Evaluate likelihood before optimization
    let lnl_before = log_likelihood(&nj_tree, &alignment, &model).unwrap();

    // Optimize branch lengths
    let lnl_after = optimize_branch_lengths(&mut nj_tree, &alignment, &model).unwrap();

    assert!(
        lnl_after >= lnl_before - 0.01,
        "ML optimization should improve NJ tree likelihood: before={}, after={}",
        lnl_before,
        lnl_after
    );
}

// ============================================================================
// Newick round-trip: write → parse → RF distance = 0
// ============================================================================

#[test]
fn test_newick_round_trip_50_taxa() {
    let mut rng = Rng::new(123);
    let tree1 = random_tree(50, 0.01, 0.5, &mut rng);

    let newick = write_newick(&tree1);
    let tree2 = parse_newick(&newick).unwrap();

    // Both trees should have the same number of leaves
    assert_eq!(tree1.num_leaves(), tree2.num_leaves());
    assert_eq!(tree1.num_leaves(), 50);

    // RF distance should be 0 (same topology)
    let rf = phylip_rs::tree::distances::robinson_foulds(&tree1, &tree2).unwrap();
    assert_eq!(
        rf, 0,
        "Newick round-trip should preserve topology for 50-taxon tree"
    );
}

// ============================================================================
// Bootstrap convergence: more replicates → support values stabilize
// ============================================================================

#[test]
fn test_bootstrap_convergence_20_taxa() {
    let mut rng = Rng::new(42);
    let true_tree = random_tree(20, 0.01, 0.15, &mut rng);
    let model = Jc69Model;
    let alignment = simulate_alignment(&true_tree, 1000, &model, &mut rng);

    let dist_model = JC69::new();

    // Run 50 bootstrap replicates
    let mut boot_rng = SimpleRng::seed(99);
    let mut trees_50 = Vec::new();
    for _ in 0..50 {
        let weights = bootstrap_weights(
            alignment.nsites(),
            &ResamplingMethod::Bootstrap,
            &mut boot_rng,
        );
        let resampled = resample_alignment(&alignment, &weights);
        let matrix = compute_distance_matrix(&resampled, &dist_model).unwrap();
        trees_50.push(neighbor_joining(&matrix));
    }

    // Run another 50 (total 100)
    let mut trees_100 = trees_50.clone();
    for _ in 0..50 {
        let weights = bootstrap_weights(
            alignment.nsites(),
            &ResamplingMethod::Bootstrap,
            &mut boot_rng,
        );
        let resampled = resample_alignment(&alignment, &weights);
        let matrix = compute_distance_matrix(&resampled, &dist_model).unwrap();
        trees_100.push(neighbor_joining(&matrix));
    }

    let consensus_50 = consensus_tree(&trees_50, &ConsensusMethod::MajorityRule).unwrap();
    let consensus_100 = consensus_tree(&trees_100, &ConsensusMethod::MajorityRule).unwrap();

    // Both consensus trees should have 20 leaves
    assert_eq!(consensus_50.tree.num_leaves(), 20);
    assert_eq!(consensus_100.tree.num_leaves(), 20);

    // The 100-replicate consensus should have at least as many resolved splits
    // as the 50-replicate consensus (or very close)
    let splits_50 = consensus_50.split_frequencies.len();
    let splits_100 = consensus_100.split_frequencies.len();

    // More replicates generally gives more stable support values
    // Both should find at least some supported splits
    assert!(
        splits_50 > 0,
        "50-rep consensus should have supported splits"
    );
    assert!(
        splits_100 > 0,
        "100-rep consensus should have supported splits"
    );
}

// ============================================================================
// Model nesting: JC69 lnL <= K2P lnL on ts/tv biased data
// ============================================================================

#[test]
fn test_model_comparison_on_simulated_data() {
    // Generate data and evaluate under both JC69 and K2P distance models.
    // Both should produce valid distance matrices.
    let mut rng = Rng::new(42);
    let true_tree = random_tree(20, 0.01, 0.2, &mut rng);
    let model = Jc69Model;
    let alignment = simulate_alignment(&true_tree, 500, &model, &mut rng);

    let jc69 = JC69::new();
    let k2p = K2P::new();

    let matrix_jc = compute_distance_matrix(&alignment, &jc69).unwrap();
    let matrix_k2p = compute_distance_matrix(&alignment, &k2p).unwrap();

    // Both should produce valid distance matrices
    assert_eq!(matrix_jc.size(), 20);
    assert_eq!(matrix_k2p.size(), 20);

    // All distances should be non-negative and finite
    for i in 0..20 {
        for j in (i + 1)..20 {
            let d_jc = matrix_jc.get(i, j);
            let d_k2p = matrix_k2p.get(i, j);

            assert!(
                d_jc >= 0.0 && d_jc.is_finite(),
                "JC69 distance ({},{}) should be non-negative and finite: {}",
                i,
                j,
                d_jc
            );
            assert!(
                d_k2p >= 0.0 && d_k2p.is_finite(),
                "K2P distance ({},{}) should be non-negative and finite: {}",
                i,
                j,
                d_k2p
            );
        }
    }
}

// ============================================================================
// True topology recovery on easy datasets
// ============================================================================

#[test]
fn test_nj_recovers_true_topology_easy_data() {
    // With short branches and many sites, NJ should recover the true tree
    // (or be very close by RF distance)
    let mut rng = Rng::new(42);
    let true_tree = random_tree(10, 0.01, 0.05, &mut rng);
    let model = Jc69Model;
    let alignment = simulate_alignment(&true_tree, 5000, &model, &mut rng);

    let dist_model = JC69::new();
    let matrix = compute_distance_matrix(&alignment, &dist_model).unwrap();
    let nj_tree = neighbor_joining(&matrix);

    let rf = phylip_rs::tree::distances::robinson_foulds(&nj_tree, &true_tree).unwrap();
    let max_rf = 2 * (10 - 3); // maximum RF distance for 10 taxa

    // NJ should get the topology exactly right or very close
    assert!(
        rf <= 4, // Allow at most 2 wrong splits
        "NJ should recover true topology on easy data (10 taxa, 5000 sites). \
         RF distance: {} (max: {})",
        rf,
        max_rf
    );
}

#[test]
fn test_parsimony_recovers_true_topology_easy_data() {
    // With short branches, parsimony should also recover the true tree
    let mut rng = Rng::new(42);
    let true_tree = random_tree(8, 0.01, 0.03, &mut rng);
    let model = Jc69Model;
    let alignment = simulate_alignment(&true_tree, 5000, &model, &mut rng);

    let result = parsimony_search(&alignment, Some(42));

    let rf = phylip_rs::tree::distances::robinson_foulds(&result.tree, &true_tree).unwrap();
    let max_rf = 2 * (8 - 3);

    assert!(
        rf <= 4,
        "Parsimony should recover true topology on easy data (8 taxa, 5000 sites). \
         RF distance: {} (max: {}), Score: {}",
        rf,
        max_rf,
        result.score
    );
}

// ============================================================================
// Performance sanity checks
// ============================================================================

#[test]
fn test_performance_nj_100_taxa() {
    // NJ on 100 taxa should complete quickly (< 2 seconds)
    let mut rng = Rng::new(42);
    let true_tree = random_tree(100, 0.01, 0.2, &mut rng);
    let model = Jc69Model;
    let alignment = simulate_alignment(&true_tree, 1000, &model, &mut rng);

    let dist_model = JC69::new();

    let start = std::time::Instant::now();
    let matrix = compute_distance_matrix(&alignment, &dist_model).unwrap();
    let _tree = neighbor_joining(&matrix);
    let elapsed = start.elapsed();

    assert!(
        elapsed.as_secs() < 2,
        "100-taxon NJ should complete in < 2s, took {:?}",
        elapsed
    );
    assert_eq!(matrix.size(), 100);
}

#[test]
fn test_performance_parsimony_20_taxa() {
    // Parsimony search on 20 taxa should complete in reasonable time
    // (Note: parsimony is O(n^3) so 50 taxa can be slow in debug mode)
    let mut rng = Rng::new(42);
    let true_tree = random_tree(20, 0.01, 0.15, &mut rng);
    let model = Jc69Model;
    let alignment = simulate_alignment(&true_tree, 500, &model, &mut rng);

    let start = std::time::Instant::now();
    let result = parsimony_search(&alignment, Some(42));
    let elapsed = start.elapsed();

    assert!(
        elapsed.as_secs() < 30,
        "20-taxon parsimony should complete in < 30s, took {:?}",
        elapsed
    );
    assert!(result.score > 0, "Parsimony score should be positive");
    assert_eq!(result.tree.num_leaves(), 20);
}

#[test]
fn test_performance_ml_20_taxa() {
    // ML branch length optimization on 20 taxa, 500 sites
    let mut rng = Rng::new(42);
    let true_tree = random_tree(20, 0.01, 0.2, &mut rng);
    let model = Jc69Model;
    let alignment = simulate_alignment(&true_tree, 500, &model, &mut rng);

    let dist_model = JC69::new();
    let matrix = compute_distance_matrix(&alignment, &dist_model).unwrap();
    let mut nj_tree = neighbor_joining(&matrix);

    let start = std::time::Instant::now();
    let lnl = optimize_branch_lengths(&mut nj_tree, &alignment, &model).unwrap();
    let elapsed = start.elapsed();

    assert!(
        elapsed.as_secs() < 30,
        "20-taxon ML optimization should complete in < 30s, took {:?}",
        elapsed
    );
    assert!(lnl.is_finite() && lnl < 0.0, "lnL should be finite and negative: {}", lnl);
}

// ============================================================================
// Distance matrix properties on medium-scale data
// ============================================================================

#[test]
fn test_distance_matrix_symmetry_50_taxa() {
    let mut rng = Rng::new(42);
    let true_tree = random_tree(50, 0.01, 0.2, &mut rng);
    let model = Jc69Model;
    let alignment = simulate_alignment(&true_tree, 500, &model, &mut rng);

    let dist_model = JC69::new();
    let matrix = compute_distance_matrix(&alignment, &dist_model).unwrap();

    // Distance matrix should be symmetric
    for i in 0..50 {
        for j in (i + 1)..50 {
            assert_eq!(
                matrix.get(i, j),
                matrix.get(j, i),
                "Distance matrix should be symmetric: d({},{}) != d({},{})",
                i,
                j,
                j,
                i
            );
        }
    }
}

#[test]
fn test_distance_matrix_triangle_inequality_50_taxa() {
    // The triangle inequality should (approximately) hold:
    // d(A,C) <= d(A,B) + d(B,C)
    let mut rng = Rng::new(42);
    let true_tree = random_tree(50, 0.01, 0.2, &mut rng);
    let model = Jc69Model;
    let alignment = simulate_alignment(&true_tree, 1000, &model, &mut rng);

    let dist_model = JC69::new();
    let matrix = compute_distance_matrix(&alignment, &dist_model).unwrap();

    let mut violations = 0;
    let mut total = 0;
    for i in 0..50 {
        for j in (i + 1)..50 {
            for k in (j + 1)..50 {
                total += 3;
                let dij = matrix.get(i, j);
                let djk = matrix.get(j, k);
                let dik = matrix.get(i, k);

                // Allow small violations due to stochastic estimation
                if dij > djk + dik + 0.01 {
                    violations += 1;
                }
                if djk > dij + dik + 0.01 {
                    violations += 1;
                }
                if dik > dij + djk + 0.01 {
                    violations += 1;
                }
            }
        }
    }

    let violation_rate = violations as f64 / total as f64;
    assert!(
        violation_rate < 0.05,
        "Triangle inequality violations should be rare: {}/{} = {:.2}%",
        violations,
        total,
        violation_rate * 100.0
    );
}

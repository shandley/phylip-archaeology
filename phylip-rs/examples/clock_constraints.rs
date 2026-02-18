// clock_constraints.rs
//
// Demonstrates how the Kitsch algorithm handles the molecular clock constraint
// and how it optimizes node heights to fit a distance matrix to an ultrametric
// tree. Compares UPGMA and Kitsch on two scenarios:
//
//   1. A "violated clock" case where one lineage (taxon E) evolves 3x faster
//      than the others, breaking the molecular clock assumption.
//   2. A "perfect clock" case where all lineages evolve at the same rate,
//      so UPGMA and Kitsch should agree.
//
// This illustrates one of the most important trade-offs in molecular
// phylogenetics: the molecular clock assumption simplifies tree inference
// (we estimate node heights instead of individual branch lengths), but
// when the assumption is violated, the resulting tree can be distorted.
//
// Kitsch optimizes node heights using weighted least squares (WLS) and
// generally achieves a better fit than UPGMA when the clock is violated.
//
// Usage:
//   cargo run --example clock_constraints

use phylip_rs::distance::kitsch::kitsch;
use phylip_rs::distance::upgma::upgma;
use phylip_rs::tree::newick::write_newick;
use phylip_rs::tree::DistanceMatrix;

// ---------------------------------------------------------------------------
// Helpers for analyzing ultrametric trees
// ---------------------------------------------------------------------------

/// Compute the distance from a leaf to the root by summing branch lengths.
fn distance_to_root(tree: &phylip_rs::tree::Tree, leaf_id: usize) -> f64 {
    let mut current = leaf_id;
    let mut total = 0.0;
    while let Some(parent) = tree.nodes[current].parent {
        total += tree.nodes[current].branch_length.unwrap_or(0.0);
        current = parent;
    }
    total
}

/// Compute the WLS (weighted least-squares) score for a tree against a
/// distance matrix, using Fitch-Margoliash weighting w_ij = 1/d_ij^2.
///
/// WLS = sum_{i<j} w_ij * (d_ij - p_ij)^2
///
/// where d_ij is the observed distance and p_ij is the tree path distance.
fn compute_wls(tree: &phylip_rs::tree::Tree, matrix: &DistanceMatrix) -> f64 {
    let n = matrix.size();

    // Build leaf name -> matrix index map
    let leaves: Vec<(usize, usize)> = tree
        .nodes
        .iter()
        .filter_map(|node| {
            if let Some(ref name) = node.name {
                matrix
                    .names
                    .iter()
                    .position(|n| n == name)
                    .map(|mi| (node.id, mi))
            } else {
                None
            }
        })
        .collect();

    let mut score = 0.0;
    for i in 0..leaves.len() {
        for j in (i + 1)..leaves.len() {
            let (node_a, mi_a) = leaves[i];
            let (node_b, mi_b) = leaves[j];
            let d_obs = matrix.get(mi_a, mi_b);
            if d_obs > 0.0 {
                let d_tree = path_distance(tree, node_a, node_b);
                let w = 1.0 / (d_obs * d_obs);
                let diff = d_obs - d_tree;
                score += w * diff * diff;
            }
        }
    }

    let _ = n;
    score
}

/// Compute the path distance between two nodes in a rooted tree.
fn path_distance(tree: &phylip_rs::tree::Tree, a: usize, b: usize) -> f64 {
    // Find LCA
    let mut ancestors_a = std::collections::HashSet::new();
    let mut current = a;
    ancestors_a.insert(current);
    while let Some(parent) = tree.nodes[current].parent {
        ancestors_a.insert(parent);
        current = parent;
    }

    let lca;
    current = b;
    if ancestors_a.contains(&current) {
        lca = current;
    } else {
        lca = 'find_lca: {
            while let Some(parent) = tree.nodes[current].parent {
                if ancestors_a.contains(&parent) {
                    break 'find_lca parent;
                }
                current = parent;
            }
            tree.root
        };
    }

    // Sum branch lengths from a to LCA + from b to LCA
    let dist_a = distance_between(tree, a, lca);
    let dist_b = distance_between(tree, b, lca);
    dist_a + dist_b
}

/// Sum branch lengths from `node` up to `ancestor`.
fn distance_between(tree: &phylip_rs::tree::Tree, node: usize, ancestor: usize) -> f64 {
    let mut current = node;
    let mut total = 0.0;
    while current != ancestor {
        total += tree.nodes[current].branch_length.unwrap_or(0.0);
        current = tree
            .nodes[current]
            .parent
            .expect("reached root without finding ancestor");
    }
    total
}

// ---------------------------------------------------------------------------
// Main demonstration
// ---------------------------------------------------------------------------

fn main() {
    println!();
    println!("CLOCK CONSTRAINTS: UPGMA vs. Kitsch Under Rate Heterogeneity");
    println!("==============================================================");
    println!();
    println!("The molecular clock hypothesis assumes that all lineages evolve");
    println!("at the same rate, which means all tips of the tree are equidistant");
    println!("from the root (the tree is 'ultrametric'). This is a powerful");
    println!("constraint: instead of estimating 2n-3 branch lengths for n taxa,");
    println!("we only need to estimate n-1 node heights.");
    println!();
    println!("Two methods enforce this constraint:");
    println!("  UPGMA  -- builds the tree by sequential clustering, assuming");
    println!("            that the smallest distance corresponds to the most");
    println!("            recent divergence. Heights are set directly from");
    println!("            average pairwise distances.");
    println!("  Kitsch -- uses the UPGMA topology but then OPTIMIZES node");
    println!("            heights to minimize a weighted least-squares (WLS)");
    println!("            criterion. This can improve the fit substantially.");
    println!();
    println!("When the clock holds, both methods should agree. When the clock");
    println!("is violated (some lineages evolve faster), Kitsch typically");
    println!("achieves a better fit because it optimizes heights rather than");
    println!("setting them from raw averages.");
    println!();

    // =====================================================================
    // Scenario 1: Violated molecular clock
    // =====================================================================

    println!("Scenario 1: Violated Molecular Clock (Rate Heterogeneity)");
    println!("==========================================================");
    println!();
    println!("True tree: ((A,B),(C,(D,E)))");
    println!("  - True divergence times: root = 1.0, AB split = 0.3, CDE split = 0.6, DE split = 0.2");
    println!("  - All lineages evolve at rate 1.0 EXCEPT taxon E at rate 3.0");
    println!("  - E's distances to everything are inflated by its fast rate");
    println!();

    // True tree structure:
    //   root (height 1.0) -> AB_ancestor (height 0.3), CDE_ancestor (height 0.6)
    //   AB_ancestor -> A (height 0), B (height 0)
    //   CDE_ancestor -> C (height 0), DE_ancestor (height 0.2)
    //   DE_ancestor -> D (height 0), E (height 0)
    //
    // True path distances (with equal rates):
    //   d(A,B) = 2 * 0.3 = 0.6       (time from A to AB_ancestor + back to B)
    //   Wait -- path distances use divergence times differently.
    //
    // For a clock tree, d(i,j) = 2 * height(LCA(i,j))
    //   d(A,B) = 2 * 0.3 = 0.6
    //   d(D,E) = 2 * 0.2 = 0.4
    //   d(C,D) = d(C,E) = 2 * 0.6 = 1.2 (C and D/E share CDE_ancestor)
    //   d(A,C) = d(A,D) = d(A,E) = d(B,C) = d(B,D) = d(B,E) = 2 * 1.0 = 2.0
    //
    // With rate heterogeneity (E at 3x):
    // E's branch from DE_ancestor is 0.2 time units * rate 3.0 = 0.6 (instead of 0.2)
    // E's branch from CDE_ancestor through DE_ancestor:
    //   path(E to CDE) = rate_E * 0.2 + rate_normal * 0.4 = 0.6 + 0.4 = 1.0
    //   Wait, the branch from DE_anc to CDE_anc has the "normal" rate.
    //
    // Actually, the rate applies to the terminal branch of E only.
    // E's terminal branch: time = 0.2, rate = 3.0, distance = 0.6
    // All other terminal branches: time = their height difference, rate = 1.0

    // Let's be explicit about branch lengths as distances (rate * time):
    // A: branch = 1.0 * 1.0 = 1.0 (from root height 1.0 to A, through AB_anc)
    //    Actually A's total path to root = A->AB_anc (0.3) + AB_anc->root (0.7) = 1.0
    // Let me think in terms of branch distances:
    //
    // Normal rate = 1.0:
    //   A -> AB_anc: 0.3 * 1.0 = 0.3
    //   B -> AB_anc: 0.3 * 1.0 = 0.3
    //   AB_anc -> root: 0.7 * 1.0 = 0.7
    //   C -> CDE_anc: 0.6 * 1.0 = 0.6
    //   D -> DE_anc: 0.2 * 1.0 = 0.2
    //   DE_anc -> CDE_anc: 0.4 * 1.0 = 0.4
    //   CDE_anc -> root: 0.4 * 1.0 = 0.4
    //
    // E has rate 3.0, so:
    //   E -> DE_anc: 0.2 * 3.0 = 0.6  (inflated!)
    //
    // Pairwise distances:
    //   d(A,B) = 0.3 + 0.3 = 0.6
    //   d(D,E) = 0.2 + 0.6 = 0.8  (not 0.4 anymore!)
    //   d(C,D) = 0.6 + 0.4 + 0.2 = 1.2
    //   d(C,E) = 0.6 + 0.4 + 0.6 = 1.6  (inflated!)
    //   d(A,C) = 0.3 + 0.7 + 0.4 + 0.6 = 2.0
    //   d(A,D) = 0.3 + 0.7 + 0.4 + 0.4 + 0.2 = 2.0
    //   d(A,E) = 0.3 + 0.7 + 0.4 + 0.4 + 0.6 = 2.4  (inflated!)
    //   d(B,C) = 0.3 + 0.7 + 0.4 + 0.6 = 2.0
    //   d(B,D) = 0.3 + 0.7 + 0.4 + 0.4 + 0.2 = 2.0
    //   d(B,E) = 0.3 + 0.7 + 0.4 + 0.4 + 0.6 = 2.4  (inflated!)

    let mut matrix1 = DistanceMatrix::new(vec![
        "A".into(),
        "B".into(),
        "C".into(),
        "D".into(),
        "E".into(),
    ]);
    matrix1.set(0, 1, 0.6); // A-B
    matrix1.set(0, 2, 2.0); // A-C
    matrix1.set(0, 3, 2.0); // A-D
    matrix1.set(0, 4, 2.4); // A-E (inflated)
    matrix1.set(1, 2, 2.0); // B-C
    matrix1.set(1, 3, 2.0); // B-D
    matrix1.set(1, 4, 2.4); // B-E (inflated)
    matrix1.set(2, 3, 1.2); // C-D
    matrix1.set(2, 4, 1.6); // C-E (inflated)
    matrix1.set(3, 4, 0.8); // D-E (inflated)

    println!("Distance matrix (E's distances are inflated by 3x rate):");
    print_distance_matrix(&matrix1);
    println!();

    // Run UPGMA
    let upgma_tree1 = upgma(&matrix1);
    let upgma_newick1 = write_newick(&upgma_tree1);
    let upgma_wls1 = compute_wls(&upgma_tree1, &matrix1);

    // Run Kitsch
    let kitsch_result1 = kitsch(&matrix1);
    let kitsch_newick1 = write_newick(&kitsch_result1.tree);

    println!("UPGMA result:");
    println!("  Tree:      {}", upgma_newick1);
    println!("  WLS score: {:.6}", upgma_wls1);
    println!();

    println!("Kitsch result:");
    println!("  Tree:       {}", kitsch_newick1);
    println!("  WLS score:  {:.6}", kitsch_result1.wls_score);
    println!("  Root height: {:.4}", kitsch_result1.root_height);
    println!();

    // Analyze ultrametric property
    println!("Ultrametric property (root-to-tip distances should be equal):");
    println!();
    println!("  UPGMA:");
    for leaf in upgma_tree1.leaves() {
        let d = distance_to_root(&upgma_tree1, leaf.id);
        println!(
            "    {} -> root: {:.4}",
            leaf.name.as_deref().unwrap_or("?"),
            d
        );
    }

    println!();
    println!("  Kitsch:");
    for leaf in kitsch_result1.tree.leaves() {
        let d = distance_to_root(&kitsch_result1.tree, leaf.id);
        println!(
            "    {} -> root: {:.4}",
            leaf.name.as_deref().unwrap_or("?"),
            d
        );
    }
    println!();

    // Compare WLS scores
    if kitsch_result1.wls_score < upgma_wls1 {
        println!(
            "  Kitsch WLS ({:.6}) < UPGMA WLS ({:.6})",
            kitsch_result1.wls_score, upgma_wls1
        );
        println!(
            "  Kitsch achieves a {:.1}% improvement in fit.",
            100.0 * (1.0 - kitsch_result1.wls_score / upgma_wls1)
        );
    } else {
        println!(
            "  UPGMA WLS ({:.6}) <= Kitsch WLS ({:.6})",
            upgma_wls1, kitsch_result1.wls_score
        );
        println!("  Both methods achieve similar fit.");
    }
    println!();

    // Show how the clock constraint handles the rate heterogeneity
    println!("How the clock constraint handles E's fast rate:");
    println!("-----------------------------------------------");
    println!();
    println!("  In the true tree, E has a 3x faster rate, making its branch");
    println!("  appear longer. But a clock-constrained tree MUST place all");
    println!("  tips at height 0 and assign each internal node a single height.");
    println!();
    println!("  This means the algorithm must find node heights that are a");
    println!("  COMPROMISE: E's distances pull the DE ancestor and root heights");
    println!("  upward, while the other taxa pull them downward. Kitsch finds");
    println!("  the best compromise by minimizing the weighted squared deviations.");
    println!();

    // =====================================================================
    // Scenario 2: Perfect molecular clock
    // =====================================================================

    println!("Scenario 2: Perfect Molecular Clock (Equal Rates)");
    println!("==================================================");
    println!();
    println!("Same tree topology: ((A,B),(C,(D,E)))");
    println!("  - All lineages evolve at rate 1.0 (clock holds perfectly)");
    println!("  - The distance matrix is perfectly ultrametric");
    println!();

    // With equal rates:
    //   d(A,B) = 0.6
    //   d(D,E) = 0.4
    //   d(C,D) = d(C,E) = 1.2
    //   d(A,C) = d(A,D) = d(A,E) = d(B,C) = d(B,D) = d(B,E) = 2.0

    let mut matrix2 = DistanceMatrix::new(vec![
        "A".into(),
        "B".into(),
        "C".into(),
        "D".into(),
        "E".into(),
    ]);
    matrix2.set(0, 1, 0.6); // A-B
    matrix2.set(0, 2, 2.0); // A-C
    matrix2.set(0, 3, 2.0); // A-D
    matrix2.set(0, 4, 2.0); // A-E (NOT inflated)
    matrix2.set(1, 2, 2.0); // B-C
    matrix2.set(1, 3, 2.0); // B-D
    matrix2.set(1, 4, 2.0); // B-E (NOT inflated)
    matrix2.set(2, 3, 1.2); // C-D
    matrix2.set(2, 4, 1.2); // C-E (NOT inflated)
    matrix2.set(3, 4, 0.4); // D-E (NOT inflated)

    println!("Distance matrix (perfect clock, all rates equal):");
    print_distance_matrix(&matrix2);
    println!();

    // Run UPGMA
    let upgma_tree2 = upgma(&matrix2);
    let upgma_newick2 = write_newick(&upgma_tree2);
    let upgma_wls2 = compute_wls(&upgma_tree2, &matrix2);

    // Run Kitsch
    let kitsch_result2 = kitsch(&matrix2);
    let kitsch_newick2 = write_newick(&kitsch_result2.tree);

    println!("UPGMA result:");
    println!("  Tree:      {}", upgma_newick2);
    println!("  WLS score: {:.6}", upgma_wls2);
    println!();

    println!("Kitsch result:");
    println!("  Tree:       {}", kitsch_newick2);
    println!("  WLS score:  {:.6}", kitsch_result2.wls_score);
    println!("  Root height: {:.4}", kitsch_result2.root_height);
    println!();

    // Both should have essentially zero WLS
    println!("  Both WLS scores should be near zero (perfect clock):");
    println!(
        "    UPGMA WLS:  {:.6}  {}",
        upgma_wls2,
        if upgma_wls2 < 0.001 {
            "(excellent)"
        } else {
            ""
        }
    );
    println!(
        "    Kitsch WLS: {:.6}  {}",
        kitsch_result2.wls_score,
        if kitsch_result2.wls_score < 0.001 {
            "(excellent)"
        } else {
            ""
        }
    );
    println!();

    // Verify ultrametric property for perfect clock
    println!("  Root-to-tip distances (should all be 1.0000):");
    println!();
    println!("  UPGMA:");
    for leaf in upgma_tree2.leaves() {
        let d = distance_to_root(&upgma_tree2, leaf.id);
        println!(
            "    {} -> root: {:.4}",
            leaf.name.as_deref().unwrap_or("?"),
            d
        );
    }
    println!();
    println!("  Kitsch:");
    for leaf in kitsch_result2.tree.leaves() {
        let d = distance_to_root(&kitsch_result2.tree, leaf.id);
        println!(
            "    {} -> root: {:.4}",
            leaf.name.as_deref().unwrap_or("?"),
            d
        );
    }
    println!();

    // Verify topology: both should recover ((A,B),(C,(D,E)))
    println!("  Topology verification:");

    // Check UPGMA topology
    let upgma_ab = check_siblings(&upgma_tree2, "A", "B");
    let upgma_de = check_siblings(&upgma_tree2, "D", "E");
    println!(
        "    UPGMA:  A-B siblings? {}  D-E siblings? {}",
        if upgma_ab { "YES" } else { "no" },
        if upgma_de { "YES" } else { "no" },
    );

    // Check Kitsch topology
    let kitsch_ab = check_siblings(&kitsch_result2.tree, "A", "B");
    let kitsch_de = check_siblings(&kitsch_result2.tree, "D", "E");
    println!(
        "    Kitsch: A-B siblings? {}  D-E siblings? {}",
        if kitsch_ab { "YES" } else { "no" },
        if kitsch_de { "YES" } else { "no" },
    );

    if upgma_ab && upgma_de && kitsch_ab && kitsch_de {
        println!("    --> Both methods recover the correct topology!");
    }
    println!();

    // =====================================================================
    // Summary comparison
    // =====================================================================

    println!("Summary: How Node Height Optimization Helps");
    println!("============================================");
    println!();
    println!("                   | Violated Clock | Perfect Clock");
    println!("  -----------------+----------------+--------------");
    println!(
        "  UPGMA WLS score  | {:>14.6} | {:>13.6}",
        upgma_wls1, upgma_wls2
    );
    println!(
        "  Kitsch WLS score | {:>14.6} | {:>13.6}",
        kitsch_result1.wls_score, kitsch_result2.wls_score
    );
    println!(
        "  Root height (K)  | {:>14.4} | {:>13.4}",
        kitsch_result1.root_height, kitsch_result2.root_height
    );
    println!();

    println!("Key insights:");
    println!();
    println!("  1. ULTRAMETRIC CONSTRAINT: Both UPGMA and Kitsch produce trees where");
    println!("     all tips are equidistant from the root. This is the molecular clock");
    println!("     constraint -- it assumes a constant rate of evolution.");
    println!();
    println!("  2. HEIGHT OPTIMIZATION: Kitsch improves on UPGMA by treating node");
    println!("     heights as free parameters and optimizing them to minimize the");
    println!("     weighted least-squares criterion. UPGMA simply averages distances.");
    println!();
    println!("  3. RATE HETEROGENEITY: When one lineage (E) evolves faster, the");
    println!("     clock is violated. Both methods must compromise, but Kitsch finds");
    println!("     a better compromise because it minimizes total error.");
    println!();
    println!("  4. PERFECT CLOCK: When the clock holds, the distance matrix is");
    println!("     already ultrametric, and both methods achieve near-zero WLS.");
    println!("     UPGMA's simple averaging is sufficient; Kitsch cannot improve.");
    println!();
    println!("  5. PHYLIP'S APPROACH: Felsenstein implemented both methods because");
    println!("     he understood that the molecular clock is sometimes useful (it");
    println!("     reduces parameters and enables divergence time estimation) but");
    println!("     also often violated. The kitsch program (named after the German");
    println!("     word for 'tacky art') reflects his characteristic dry humor.");
    println!();

    println!("References:");
    println!("  Sokal, R.R. & Michener, C.D. (1958). A statistical method for");
    println!("  evaluating systematic relationships. University of Kansas Science");
    println!("  Bulletin, 38:1409-1438.");
    println!();
    println!("  Fitch, W.M. & Margoliash, E. (1967). Construction of phylogenetic");
    println!("  trees. Science, 155(3760):279-284.");
    println!();
    println!("  Felsenstein, J. (2004). Inferring Phylogenies. Sinauer Associates.");
    println!("  Chapter 11: Distance matrix methods.");
    println!();
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Print a labeled distance matrix.
fn print_distance_matrix(matrix: &DistanceMatrix) {
    let n = matrix.size();

    print!("          ");
    for name in &matrix.names {
        print!("{:>8}", name);
    }
    println!();

    for i in 0..n {
        print!("  {:>6}  ", matrix.names[i]);
        for j in 0..n {
            if i == j {
                print!("       -");
            } else {
                print!("  {:.4}", matrix.get(i, j));
            }
        }
        println!();
    }
}

/// Check whether two named taxa are siblings (share the same parent).
fn check_siblings(tree: &phylip_rs::tree::Tree, name_a: &str, name_b: &str) -> bool {
    let node_a = tree
        .nodes
        .iter()
        .find(|n| n.name.as_deref() == Some(name_a));
    let node_b = tree
        .nodes
        .iter()
        .find(|n| n.name.as_deref() == Some(name_b));

    match (node_a, node_b) {
        (Some(a), Some(b)) => a.parent.is_some() && a.parent == b.parent,
        _ => false,
    }
}

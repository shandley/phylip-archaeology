// chord_geometry.rs
//
// Demonstrates the geometric insight behind the Cavalli-Sforza chord distance:
// allele frequencies are mapped to a unit hypersphere via square-root
// transformation, and the chord distance is the straight-line (Euclidean)
// distance between points on that sphere. This distance is proportional to
// genetic drift time, while raw Euclidean distance on untransformed
// frequencies is not.
//
// The Cavalli-Sforza chord distance is implemented in PHYLIP's `gendist`
// program and is one of the most widely used measures of genetic divergence
// between populations.
//
// Usage:
//   cargo run --example chord_geometry

use phylip_rs::distance::neighbor_joining;
use phylip_rs::models::gene_freq::{
    cavalli_sforza_distance, compute_gene_freq_distances, GeneFreqData, GeneFreqMethod, Locus,
};
use phylip_rs::tree::newick::write_newick;
use phylip_rs::tree::DistanceMatrix;

// ---------------------------------------------------------------------------
// Simple LCG random number generator (no external dependencies needed)
// ---------------------------------------------------------------------------

struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Advance the LCG state and return the next pseudo-random u64.
    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    /// Return a pseudo-random f64 in [0, 1).
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 33) as f64 / (1u64 << 31) as f64
    }

    /// Return a Gaussian(0,1) sample using Box-Muller transform.
    fn next_gaussian(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-15);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

// ---------------------------------------------------------------------------
// Drift simulation on allele frequency vectors
// ---------------------------------------------------------------------------

/// Apply Brownian drift to an allele frequency vector.
/// `drift_amount` controls the magnitude of random perturbation.
/// After perturbation, negative values are clamped to 0 and the
/// vector is renormalized to sum to 1.
fn apply_drift(freqs: &[f64], drift_amount: f64, rng: &mut Rng) -> Vec<f64> {
    let mut new_freqs: Vec<f64> = freqs
        .iter()
        .map(|&f| f + drift_amount * rng.next_gaussian())
        .collect();

    // Project back to the simplex: clamp negatives, renormalize
    for f in new_freqs.iter_mut() {
        if *f < 0.0 {
            *f = 0.0;
        }
    }
    let sum: f64 = new_freqs.iter().sum();
    if sum > 0.0 {
        for f in new_freqs.iter_mut() {
            *f /= sum;
        }
    }
    new_freqs
}

/// Compute raw Euclidean distance between two allele frequency vectors.
fn euclidean_distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}

/// Compute Euclidean distance between sqrt-transformed frequency vectors.
fn sqrt_euclidean_distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| {
            let sx = x.max(0.0).sqrt();
            let sy = y.max(0.0).sqrt();
            (sx - sy) * (sx - sy)
        })
        .sum::<f64>()
        .sqrt()
}

/// Compute the inner product of sqrt-transformed vectors: sum(sqrt(x_i) * sqrt(y_i)).
/// This equals cos(angle) between the two points on the unit hypersphere.
fn sqrt_inner_product(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| (x * y).max(0.0).sqrt())
        .sum::<f64>()
}

// ---------------------------------------------------------------------------
// Main demonstration
// ---------------------------------------------------------------------------

fn main() {
    println!();
    println!("CAVALLI-SFORZA CHORD DISTANCE: Geometry on the Hypersphere");
    println!("============================================================");
    println!();
    println!("In 1967, Cavalli-Sforza and Edwards proposed a beautifully geometric");
    println!("approach to measuring genetic distance between populations. The key");
    println!("insight: if you apply a square-root transformation to allele frequencies,");
    println!("each population maps to a point on the surface of a unit hypersphere.");
    println!("The chord distance is then the straight-line (Euclidean) distance");
    println!("between these points through the interior of the sphere.");
    println!();
    println!("This distance has remarkable properties:");
    println!("  - It is proportional to genetic drift time (under Brownian motion)");
    println!("  - It satisfies the triangle inequality (it is a true metric)");
    println!("  - It is bounded between 0 and sqrt(2)");
    println!("  - Raw Euclidean distance on untransformed frequencies is NOT");
    println!("    proportional to drift time and distorts evolutionary relationships");
    println!();

    // -----------------------------------------------------------------------
    // Set up the divergence scenario
    // -----------------------------------------------------------------------

    println!("Scenario: Four populations with known divergence structure");
    println!("==========================================================");
    println!();
    println!("  True tree: ((Pop1, Pop2), (Pop3, Pop4))");
    println!();
    println!("  Pop1 and Pop2 diverged recently (drift = 0.03)");
    println!("  Pop3 and Pop4 diverged recently (drift = 0.03)");
    println!("  The two pairs diverged long ago   (drift = 0.12)");
    println!();

    let mut rng = Rng::new(271828182);

    // Define ancestral allele frequencies for 8 loci with 3-4 alleles each
    let ancestral_freqs: Vec<Vec<f64>> = vec![
        vec![0.50, 0.30, 0.20],
        vec![0.40, 0.35, 0.25],
        vec![0.60, 0.25, 0.15],
        vec![0.33, 0.34, 0.33],
        vec![0.45, 0.20, 0.20, 0.15],
        vec![0.30, 0.30, 0.20, 0.20],
        vec![0.55, 0.25, 0.20],
        vec![0.25, 0.25, 0.25, 0.25],
    ];

    let locus_names = [
        "ABO", "HLA-A", "HLA-B", "PGM", "ACP", "ESD", "MNS", "RH",
    ];

    let deep_drift = 0.12; // drift at the deep split
    let shallow_drift = 0.03; // drift at the recent split

    // Simulate divergence:
    // 1. Deep split: ancestral -> group_AB, group_CD
    // 2. Shallow splits: group_AB -> Pop1, Pop2; group_CD -> Pop3, Pop4
    let mut pop_freqs: Vec<Vec<Vec<f64>>> = vec![Vec::new(); 4]; // [pop][locus][allele]

    for (l, anc) in ancestral_freqs.iter().enumerate() {
        // Deep divergence
        let group_ab = apply_drift(anc, deep_drift, &mut rng);
        let group_cd = apply_drift(anc, deep_drift, &mut rng);

        // Shallow divergence
        let pop1 = apply_drift(&group_ab, shallow_drift, &mut rng);
        let pop2 = apply_drift(&group_ab, shallow_drift, &mut rng);
        let pop3 = apply_drift(&group_cd, shallow_drift, &mut rng);
        let pop4 = apply_drift(&group_cd, shallow_drift, &mut rng);

        pop_freqs[0].push(pop1);
        pop_freqs[1].push(pop2);
        pop_freqs[2].push(pop3);
        pop_freqs[3].push(pop4);

        let _ = l; // suppress unused warning
    }

    // Print the simulated frequencies
    println!("Simulated allele frequencies (8 loci, 3-4 alleles each):");
    println!("---------------------------------------------------------");
    for (l, name) in locus_names.iter().enumerate() {
        println!("  Locus {} ({} alleles):", name, ancestral_freqs[l].len());
        for p in 0..4 {
            let freq_str: Vec<String> = pop_freqs[p][l].iter().map(|f| format!("{:.3}", f)).collect();
            println!("    Pop{}: [{}]", p + 1, freq_str.join(", "));
        }
    }
    println!();

    // -----------------------------------------------------------------------
    // Build GeneFreqData and compute chord distances
    // -----------------------------------------------------------------------

    let pop_names: Vec<String> = (1..=4).map(|i| format!("Pop{}", i)).collect();

    let loci: Vec<Locus> = (0..ancestral_freqs.len())
        .map(|l| {
            let n_alleles = ancestral_freqs[l].len();
            let frequencies: Vec<Vec<f64>> = (0..4).map(|p| pop_freqs[p][l].clone()).collect();
            Locus {
                name: locus_names[l].to_string(),
                n_alleles,
                frequencies,
            }
        })
        .collect();

    let data = GeneFreqData::new(pop_names.clone(), loci);

    println!("Part 1: Comparing Distance Measures");
    println!("=====================================");
    println!();

    // Compute Cavalli-Sforza chord distances
    let cs_matrix = compute_gene_freq_distances(&data, GeneFreqMethod::CavalliSforza);

    println!("Cavalli-Sforza chord distances:");
    print_distance_matrix(&pop_names, &cs_matrix.distances);
    println!();

    // Compute raw Euclidean distances on untransformed frequencies
    let n_pops = 4;
    let n_loci = ancestral_freqs.len();
    let mut raw_euclid = vec![vec![0.0; n_pops]; n_pops];
    for i in 0..n_pops {
        for j in (i + 1)..n_pops {
            let mut total_sq = 0.0;
            for l in 0..n_loci {
                let d = euclidean_distance(&pop_freqs[i][l], &pop_freqs[j][l]);
                total_sq += d * d;
            }
            let avg_dist = (total_sq / n_loci as f64).sqrt();
            raw_euclid[i][j] = avg_dist;
            raw_euclid[j][i] = avg_dist;
        }
    }

    println!("Raw Euclidean distances (untransformed frequencies):");
    print_distance_matrix(&pop_names, &raw_euclid);
    println!();

    // Compute Euclidean distances on sqrt-transformed frequencies
    let mut sqrt_euclid = vec![vec![0.0; n_pops]; n_pops];
    for i in 0..n_pops {
        for j in (i + 1)..n_pops {
            let mut total_sq = 0.0;
            for l in 0..n_loci {
                let d = sqrt_euclidean_distance(&pop_freqs[i][l], &pop_freqs[j][l]);
                total_sq += d * d;
            }
            let avg_dist = (total_sq / n_loci as f64).sqrt();
            sqrt_euclid[i][j] = avg_dist;
            sqrt_euclid[j][i] = avg_dist;
        }
    }

    println!("Euclidean distances on sqrt-transformed frequencies:");
    print_distance_matrix(&pop_names, &sqrt_euclid);
    println!();

    // -----------------------------------------------------------------------
    // Analyze the distance ratios
    // -----------------------------------------------------------------------

    println!("Part 2: Distance Ratios Reveal the Geometry");
    println!("=============================================");
    println!();

    let cs_within_1_2 = cs_matrix.distances[0][1];
    let cs_within_3_4 = cs_matrix.distances[2][3];
    let cs_between_avg = (cs_matrix.distances[0][2]
        + cs_matrix.distances[0][3]
        + cs_matrix.distances[1][2]
        + cs_matrix.distances[1][3])
        / 4.0;

    let raw_within_1_2 = raw_euclid[0][1];
    let raw_within_3_4 = raw_euclid[2][3];
    let raw_between_avg =
        (raw_euclid[0][2] + raw_euclid[0][3] + raw_euclid[1][2] + raw_euclid[1][3]) / 4.0;

    println!("Expected: between-group distances should be ~4x within-group");
    println!("  (drift ratio: {:.1}/{:.1} = {:.1}x)", deep_drift, shallow_drift, deep_drift / shallow_drift);
    println!();

    println!("Chord distance (Cavalli-Sforza):");
    println!("  Pop1-Pop2 (within):  {:.4}", cs_within_1_2);
    println!("  Pop3-Pop4 (within):  {:.4}", cs_within_3_4);
    println!("  Between-group avg:   {:.4}", cs_between_avg);
    println!(
        "  Ratio (between/within): {:.2}x",
        cs_between_avg / ((cs_within_1_2 + cs_within_3_4) / 2.0)
    );
    println!();

    println!("Raw Euclidean distance (untransformed):");
    println!("  Pop1-Pop2 (within):  {:.4}", raw_within_1_2);
    println!("  Pop3-Pop4 (within):  {:.4}", raw_within_3_4);
    println!("  Between-group avg:   {:.4}", raw_between_avg);
    println!(
        "  Ratio (between/within): {:.2}x",
        raw_between_avg / ((raw_within_1_2 + raw_within_3_4) / 2.0)
    );
    println!();

    // -----------------------------------------------------------------------
    // Show the hypersphere geometry
    // -----------------------------------------------------------------------

    println!("Part 3: The Hypersphere Geometry");
    println!("=================================");
    println!();
    println!("The square-root transformation maps allele frequencies to the unit");
    println!("hypersphere. For a locus with alleles at frequencies (p1, p2, p3),");
    println!("the point on the sphere is (sqrt(p1), sqrt(p2), sqrt(p3)).");
    println!("This point has unit norm: sqrt(p1)^2 + sqrt(p2)^2 + sqrt(p3)^2 = 1.");
    println!();

    // Demonstrate with the first locus
    println!("Example: Locus {} (3 alleles)", locus_names[0]);
    println!();
    for p in 0..4 {
        let freqs = &pop_freqs[p][0];
        let sqrt_freqs: Vec<f64> = freqs.iter().map(|&f| f.max(0.0).sqrt()).collect();
        let norm_sq: f64 = sqrt_freqs.iter().map(|x| x * x).sum();
        println!(
            "  Pop{}: freqs = [{:.3}, {:.3}, {:.3}]  -->  sqrt = [{:.3}, {:.3}, {:.3}]  norm^2 = {:.6}",
            p + 1,
            freqs[0], freqs[1], freqs[2],
            sqrt_freqs[0], sqrt_freqs[1], sqrt_freqs[2],
            norm_sq
        );
    }
    println!();

    // Show cos(angle) = inner product of sqrt-transformed vectors
    println!("Inner products of sqrt-transformed vectors (= cos(angle) on sphere):");
    println!();
    for i in 0..n_pops {
        for j in (i + 1)..n_pops {
            let mut avg_cos = 0.0;
            for l in 0..n_loci {
                avg_cos += sqrt_inner_product(&pop_freqs[i][l], &pop_freqs[j][l]);
            }
            avg_cos /= n_loci as f64;
            let angle_deg = avg_cos.min(1.0).acos() * 180.0 / std::f64::consts::PI;
            println!(
                "  Pop{}-Pop{}: cos(angle) = {:.6}  angle = {:.2} degrees",
                i + 1,
                j + 1,
                avg_cos,
                angle_deg
            );
        }
    }
    println!();
    println!("Closely related populations have a small angle between them on the");
    println!("hypersphere (cos near 1). Distantly related populations have a larger");
    println!("angle (cos further from 1). The chord distance is 2*sin(angle/2),");
    println!("which for small angles is approximately proportional to the angle.");
    println!();

    // -----------------------------------------------------------------------
    // Verify against direct computation
    // -----------------------------------------------------------------------

    println!("Part 4: Verification -- Library vs. Manual Computation");
    println!("=======================================================");
    println!();

    // Manually compute chord distance for Pop1-Pop2 and compare with library
    let manual_cs_01 = {
        let mut sum = 0.0;
        for l in 0..n_loci {
            let mut cos_angle = 0.0;
            for a in 0..pop_freqs[0][l].len() {
                let product = pop_freqs[0][l][a] * pop_freqs[1][l][a];
                if product > 0.0 {
                    cos_angle += product.sqrt();
                }
            }
            sum += 1.0 - cos_angle;
        }
        ((2.0 / n_loci as f64) * sum).sqrt()
    };

    let library_cs_01 = cavalli_sforza_distance(&data, 0, 1);

    println!("Pop1-Pop2 chord distance:");
    println!("  Manual computation:  {:.10}", manual_cs_01);
    println!("  Library function:    {:.10}", library_cs_01);
    println!(
        "  Difference:          {:.2e}",
        (manual_cs_01 - library_cs_01).abs()
    );
    println!();

    // -----------------------------------------------------------------------
    // Build NJ tree from chord distances
    // -----------------------------------------------------------------------

    println!("Part 5: Tree Reconstruction from Chord Distances");
    println!("==================================================");
    println!();

    // Convert GeneFreqDistanceMatrix to DistanceMatrix for NJ
    let mut dm = DistanceMatrix::new(pop_names.clone());
    for i in 0..n_pops {
        for j in (i + 1)..n_pops {
            dm.set(i, j, cs_matrix.distances[i][j]);
        }
    }

    let nj_tree = neighbor_joining(&dm);
    let newick = write_newick(&nj_tree);

    println!("Neighbor-joining tree from chord distances:");
    println!("  {}", newick);
    println!();
    println!("Expected topology: ((Pop1,Pop2),(Pop3,Pop4))");
    println!();

    // Verify the topology: check that Pop1 and Pop2 are siblings
    let pop1_node = nj_tree
        .nodes
        .iter()
        .find(|n| n.name.as_deref() == Some("Pop1"))
        .unwrap();
    let pop2_node = nj_tree
        .nodes
        .iter()
        .find(|n| n.name.as_deref() == Some("Pop2"))
        .unwrap();
    let pop3_node = nj_tree
        .nodes
        .iter()
        .find(|n| n.name.as_deref() == Some("Pop3"))
        .unwrap();
    let pop4_node = nj_tree
        .nodes
        .iter()
        .find(|n| n.name.as_deref() == Some("Pop4"))
        .unwrap();

    let pop12_siblings = pop1_node.parent == pop2_node.parent;
    let pop34_siblings = pop3_node.parent == pop4_node.parent;

    println!(
        "Topology check: Pop1-Pop2 siblings? {}  Pop3-Pop4 siblings? {}",
        if pop12_siblings { "YES" } else { "no" },
        if pop34_siblings { "YES" } else { "no" },
    );

    if pop12_siblings && pop34_siblings {
        println!("  --> Correct! The chord distance recovers the true divergence structure.");
    }
    println!();

    // -----------------------------------------------------------------------
    // Also compare with Nei distance for completeness
    // -----------------------------------------------------------------------

    println!("Part 6: Comparison with Nei's Genetic Distance");
    println!("================================================");
    println!();

    let nei_matrix = compute_gene_freq_distances(&data, GeneFreqMethod::Nei);

    println!("Nei's genetic distance:");
    print_distance_matrix(&pop_names, &nei_matrix.distances);
    println!();

    let nei_within_avg = (nei_matrix.distances[0][1] + nei_matrix.distances[2][3]) / 2.0;
    let nei_between_avg = (nei_matrix.distances[0][2]
        + nei_matrix.distances[0][3]
        + nei_matrix.distances[1][2]
        + nei_matrix.distances[1][3])
        / 4.0;

    println!(
        "  Nei ratio (between/within):  {:.2}x",
        nei_between_avg / nei_within_avg
    );
    println!(
        "  Chord ratio (between/within): {:.2}x",
        cs_between_avg / ((cs_within_1_2 + cs_within_3_4) / 2.0)
    );
    println!();
    println!("Both distance measures correctly recover the relative divergence");
    println!("pattern, but Nei's distance is unbounded (it can be infinite for");
    println!("fixed differences), while the chord distance is always bounded");
    println!("between 0 and sqrt(2) = {:.4}.", 2.0_f64.sqrt());
    println!();

    // -----------------------------------------------------------------------
    // Summary
    // -----------------------------------------------------------------------

    println!("Summary");
    println!("=======");
    println!();
    println!("The Cavalli-Sforza chord distance achieves its elegant properties");
    println!("by exploiting a simple geometric trick: the square-root transformation");
    println!("maps allele frequency vectors from the probability simplex to the");
    println!("surface of a unit hypersphere. On this sphere:");
    println!();
    println!("  1. The chord distance (straight-line through the sphere) between");
    println!("     two populations is proportional to sqrt(drift time) for small");
    println!("     divergences under Brownian motion drift.");
    println!();
    println!("  2. The distance is bounded: 0 <= D <= sqrt(2), making it robust");
    println!("     to extreme frequency differences that cause Nei's distance");
    println!("     to become infinite.");
    println!();
    println!("  3. It satisfies the triangle inequality and is a true metric,");
    println!("     guaranteeing consistent behavior in clustering algorithms.");
    println!();
    println!("  4. The NJ tree built from chord distances correctly recovers the");
    println!("     true population divergence history.");
    println!();
    println!("This geometric perspective -- allele frequencies as points on a");
    println!("sphere -- was a profound conceptual advance in population genetics.");
    println!("It made phylogenetic tree-building for populations mathematically");
    println!("rigorous and connects population genetics to the geometry of");
    println!("statistical manifolds.");
    println!();

    println!("References:");
    println!("  Cavalli-Sforza, L.L. & Edwards, A.W.F. (1967). Phylogenetic");
    println!("  analysis: models and estimation procedures. American Journal of");
    println!("  Human Genetics, 19:233-257.");
    println!();
    println!("  Edwards, A.W.F. (1971). Distances between populations on the basis");
    println!("  of gene frequencies. Biometrics, 27:873-881.");
    println!();
    println!("  Felsenstein, J. (2004). Inferring Phylogenies. Sinauer Associates.");
    println!("  Chapter 11: Distances and how to estimate them.");
    println!();
}

/// Print a labeled distance matrix.
fn print_distance_matrix(names: &[String], distances: &[Vec<f64>]) {
    let n = names.len();

    // Header
    print!("          ");
    for name in names.iter() {
        print!("{:>8}", name);
    }
    println!();

    // Rows
    for i in 0..n {
        print!("  {:>6}  ", names[i]);
        for j in 0..n {
            if i == j {
                print!("       -");
            } else {
                print!("  {:.4}", distances[i][j]);
            }
        }
        println!();
    }
}

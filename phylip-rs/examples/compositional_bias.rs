// compositional_bias.rs
//
// Demonstrates that the LogDet distance is robust to base composition bias
// while JC69 is not. When GC-content varies across taxa, JC69 can group
// taxa with similar composition (convergent bias) rather than by true
// evolutionary relationship, while LogDet correctly recovers the true tree.
//
// This is a fundamental result in molecular phylogenetics: the LogDet/
// Paralinear distance (Lake 1994, Lockhart et al. 1994) was specifically
// developed to handle non-stationary base composition across lineages.
//
// Usage:
//   cargo run --example compositional_bias

use phylip_rs::distance::neighbor_joining;
use phylip_rs::likelihood::models::{Jc69Model, SubstitutionModel};
use phylip_rs::models::jc69::JC69;
use phylip_rs::models::logdet::LogDet;
use phylip_rs::models::compute_distance_matrix;
use phylip_rs::tree::newick::{parse_newick, write_newick};
use phylip_rs::tree::{Alignment, Base, Sequence};

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
}

// ---------------------------------------------------------------------------
// Sequence simulator along a phylogenetic tree under JC69
// ---------------------------------------------------------------------------

fn index_to_base(i: usize) -> Base {
    match i {
        0 => Base::A,
        1 => Base::C,
        2 => Base::G,
        3 => Base::T,
        _ => unreachable!(),
    }
}

fn base_to_index(b: Base) -> usize {
    match b {
        Base::A => 0,
        Base::C => 1,
        Base::G => 2,
        Base::T => 3,
        _ => 0,
    }
}

/// Simulate DNA sequences along a tree using the JC69 model.
fn simulate_alignment(
    tree: &phylip_rs::tree::Tree,
    nsites: usize,
    model: &Jc69Model,
    rng: &mut Rng,
) -> Alignment {
    let nnodes = tree.num_nodes();
    let mut node_seqs: Vec<Vec<usize>> = vec![vec![0; nsites]; nnodes];

    // Generate root sequence
    for site in 0..nsites {
        let root_base = (rng.next_u64() >> 33) as usize % 4;
        node_seqs[tree.root][site] = root_base;
    }

    // Traverse the tree in preorder (parent before children)
    let preorder = tree.preorder();
    for &node_id in &preorder {
        if node_id == tree.root {
            continue;
        }

        let branch_len = tree.nodes[node_id].branch_length.unwrap_or(0.1);
        let trans_probs = model.transition_probs(branch_len);
        let parent_id = tree.nodes[node_id].parent.unwrap();

        for site in 0..nsites {
            let parent_base = node_seqs[parent_id][site];
            let row = &trans_probs[parent_base];
            let u = rng.next_f64();
            let mut cumulative = 0.0;
            let mut child_base = 3;
            for j in 0..4 {
                cumulative += row[j];
                if u < cumulative {
                    child_base = j;
                    break;
                }
            }
            node_seqs[node_id][site] = child_base;
        }
    }

    // Collect leaf sequences
    let mut sequences = Vec::new();
    for node in &tree.nodes {
        if node.is_leaf() {
            if let Some(ref name) = node.name {
                let bases: Vec<Base> = node_seqs[node.id]
                    .iter()
                    .map(|&i| index_to_base(i))
                    .collect();
                sequences.push(Sequence::new(name.clone(), bases));
            }
        }
    }

    sequences.sort_by(|a, b| a.name.cmp(&b.name));
    Alignment::new(sequences).expect("simulated alignment should be valid")
}

/// Introduce GC-bias into a sequence by randomly converting A->G and T->C.
/// `bias_rate` is the fraction of AT bases that get shifted to GC.
fn apply_gc_bias(alignment: &Alignment, taxon_name: &str, bias_rate: f64, rng: &mut Rng) -> Vec<Base> {
    let seq = alignment
        .sequences
        .iter()
        .find(|s| s.name == taxon_name)
        .expect("taxon not found in alignment");

    seq.bases
        .iter()
        .map(|&b| {
            match b {
                Base::A => {
                    if rng.next_f64() < bias_rate {
                        Base::G // A -> G (purine stays purine, but GC increases)
                    } else {
                        b
                    }
                }
                Base::T => {
                    if rng.next_f64() < bias_rate {
                        Base::C // T -> C (pyrimidine stays pyrimidine, GC increases)
                    } else {
                        b
                    }
                }
                _ => b,
            }
        })
        .collect()
}

/// Compute GC content of a base sequence.
fn gc_content(bases: &[Base]) -> f64 {
    let gc = bases
        .iter()
        .filter(|&&b| b == Base::G || b == Base::C)
        .count();
    gc as f64 / bases.len() as f64
}

/// Count base frequencies for a sequence.
fn base_frequencies(bases: &[Base]) -> [f64; 4] {
    let mut counts = [0usize; 4];
    for &b in bases {
        if b.is_definite() {
            counts[base_to_index(b)] += 1;
        }
    }
    let total = counts.iter().sum::<usize>() as f64;
    [
        counts[0] as f64 / total,
        counts[1] as f64 / total,
        counts[2] as f64 / total,
        counts[3] as f64 / total,
    ]
}

/// Check if a newick string groups the two specified taxa as sisters.
fn are_sisters_in_tree(newick: &str, taxon1: &str, taxon2: &str) -> bool {
    // Simple heuristic: check if the newick contains (taxon1...,taxon2...) or
    // (taxon2...,taxon1...) as a cherry (sister pair).
    // We look for patterns like (A:...,B:...) in the newick string.
    let newick_clean = newick.replace(' ', "");

    // Find all cherry pairs: substrings of the form (X:num,Y:num)
    // where X and Y are taxon names
    let taxa = [taxon1, taxon2];
    for window_start in 0..newick_clean.len() {
        if newick_clean.as_bytes()[window_start] == b'(' {
            // Find matching close paren at same nesting level
            let mut depth = 0;
            let mut end = window_start;
            for (i, ch) in newick_clean[window_start..].chars().enumerate() {
                match ch {
                    '(' => depth += 1,
                    ')' => {
                        depth -= 1;
                        if depth == 0 {
                            end = window_start + i;
                            break;
                        }
                    }
                    _ => {}
                }
            }

            let inner = &newick_clean[window_start + 1..end];
            // Check if this is a cherry (no nested parens)
            if !inner.contains('(') && !inner.contains(')') {
                // Split by comma
                let parts: Vec<&str> = inner.split(',').collect();
                if parts.len() == 2 {
                    let name1 = parts[0].split(':').next().unwrap_or("");
                    let name2 = parts[1].split(':').next().unwrap_or("");
                    if (name1 == taxa[0] && name2 == taxa[1])
                        || (name1 == taxa[1] && name2 == taxa[0])
                    {
                        return true;
                    }
                }
            }
        }
    }
    false
}

fn main() {
    println!();
    println!("COMPOSITIONAL BIAS: LogDet vs JC69 Under Non-Stationary Evolution");
    println!("===================================================================");
    println!();
    println!("When base composition varies across lineages (e.g., some organisms are");
    println!("GC-rich while others are AT-rich), distance methods that assume equal");
    println!("composition -- like JC69 -- can be misled. Taxa with similar GC-content");
    println!("appear artificially close, even when they are not closely related.");
    println!();
    println!("The LogDet/Paralinear distance (Lake 1994, Lockhart et al. 1994) was");
    println!("specifically designed to handle this situation. It uses the determinant");
    println!("of the divergence matrix to factor out composition effects, giving");
    println!("consistent distance estimates regardless of base frequency variation.");
    println!();

    // =========================================================================
    // Step 1: Define the true tree and simulate sequences
    // =========================================================================
    println!("=========================================================================");
    println!("  STEP 1: Simulate sequences on the TRUE tree under JC69");
    println!("=========================================================================");
    println!();

    let true_tree_newick = "((A:0.15,B:0.15):0.10,(C:0.15,D:0.15):0.10);";
    let true_tree = parse_newick(true_tree_newick).expect("true tree should parse");

    println!("  True tree: ((A,B),(C,D))");
    println!("  A and B are sisters (same clade)");
    println!("  C and D are sisters (same clade)");
    println!("  Branch lengths: 0.15 to tips, 0.10 internal");
    println!();

    let model = Jc69Model;
    let nsites = 2000;
    let mut rng = Rng::new(271828182845);

    let original_alignment = simulate_alignment(&true_tree, nsites, &model, &mut rng);

    println!("  Simulated {} sites on 4 taxa under JC69 (equal base frequencies).", nsites);
    println!();

    // Show original base frequencies
    println!("  Original base frequencies (before bias):");
    for seq in &original_alignment.sequences {
        let freq = base_frequencies(&seq.bases);
        println!(
            "    {}: A={:.3}  C={:.3}  G={:.3}  T={:.3}  (GC={:.1}%)",
            seq.name,
            freq[0],
            freq[1],
            freq[2],
            freq[3],
            gc_content(&seq.bases) * 100.0
        );
    }
    println!();

    // =========================================================================
    // Step 2: Introduce GC-bias in taxa A and C
    // =========================================================================
    println!("=========================================================================");
    println!("  STEP 2: Introduce GC-bias in taxa A and C (NOT sisters!)");
    println!("=========================================================================");
    println!();
    println!("  We post-hoc shift ~30% of A->G and T->C bases in taxa A and C.");
    println!("  This simulates a convergent shift to high GC-content in two");
    println!("  UNRELATED lineages. Taxa B and D retain their original composition.");
    println!();

    let bias_rate = 0.30;
    let mut rng2 = Rng::new(314159265358);
    let biased_a = apply_gc_bias(&original_alignment, "A", bias_rate, &mut rng2);
    let biased_c = apply_gc_bias(&original_alignment, "C", bias_rate, &mut rng2);

    // Get unmodified sequences for B and D
    let seq_b = original_alignment
        .sequences
        .iter()
        .find(|s| s.name == "B")
        .unwrap()
        .bases
        .clone();
    let seq_d = original_alignment
        .sequences
        .iter()
        .find(|s| s.name == "D")
        .unwrap()
        .bases
        .clone();

    let biased_alignment = Alignment::new(vec![
        Sequence::new("A", biased_a.clone()),
        Sequence::new("B", seq_b.clone()),
        Sequence::new("C", biased_c.clone()),
        Sequence::new("D", seq_d.clone()),
    ])
    .expect("biased alignment should be valid");

    println!("  Base frequencies AFTER compositional bias:");
    for seq in &biased_alignment.sequences {
        let freq = base_frequencies(&seq.bases);
        let gc = gc_content(&seq.bases);
        let marker = if seq.name == "A" || seq.name == "C" {
            " <-- GC-biased"
        } else {
            ""
        };
        println!(
            "    {}: A={:.3}  C={:.3}  G={:.3}  T={:.3}  (GC={:.1}%){}",
            seq.name, freq[0], freq[1], freq[2], freq[3], gc * 100.0, marker
        );
    }
    println!();
    println!("  Note: A and C now share similar high GC-content (~{:.0}%),",
             gc_content(&biased_a) * 100.0);
    println!("  even though they are NOT sister taxa in the true tree.");
    println!();

    // =========================================================================
    // Step 3: Compute distance matrices
    // =========================================================================
    println!("=========================================================================");
    println!("  STEP 3: Compute pairwise distances using JC69 and LogDet");
    println!("=========================================================================");
    println!();

    let jc69_model = JC69::new();
    let logdet_model = LogDet::new();

    let jc69_matrix = compute_distance_matrix(&biased_alignment, &jc69_model)
        .expect("JC69 distance computation should succeed");
    let logdet_matrix = compute_distance_matrix(&biased_alignment, &logdet_model)
        .expect("LogDet distance computation should succeed");

    let names = ["A", "B", "C", "D"];

    println!("  JC69 distance matrix:");
    print!("          ");
    for name in &names {
        print!("{:>8}", name);
    }
    println!();
    for i in 0..4 {
        print!("    {:>4}  ", names[i]);
        for j in 0..4 {
            print!("{:>8.4}", jc69_matrix.get(i, j));
        }
        println!();
    }
    println!();

    println!("  LogDet distance matrix:");
    print!("          ");
    for name in &names {
        print!("{:>8}", name);
    }
    println!();
    for i in 0..4 {
        print!("    {:>4}  ", names[i]);
        for j in 0..4 {
            print!("{:>8.4}", logdet_matrix.get(i, j));
        }
        println!();
    }
    println!();

    // Highlight the key distances
    let jc69_ab = jc69_matrix.get(0, 1);
    let jc69_ac = jc69_matrix.get(0, 2);
    let logdet_ab = logdet_matrix.get(0, 1);
    let logdet_ac = logdet_matrix.get(0, 2);

    println!("  KEY COMPARISON:");
    println!("    JC69:   d(A,B) = {:.4}   d(A,C) = {:.4}   {} A-B than A-C",
             jc69_ab, jc69_ac,
             if jc69_ab < jc69_ac { "CLOSER" } else { "FARTHER" });
    println!("    LogDet: d(A,B) = {:.4}   d(A,C) = {:.4}   {} A-B than A-C",
             logdet_ab, logdet_ac,
             if logdet_ab < logdet_ac { "CLOSER" } else { "FARTHER" });
    println!();

    if jc69_ac < jc69_ab {
        println!("  JC69 thinks A is CLOSER to C -- WRONG! (attracted by similar GC-content)");
    }
    if logdet_ab < logdet_ac {
        println!("  LogDet correctly finds A is closer to B -- RIGHT!");
    }
    println!();

    // =========================================================================
    // Step 4: Build NJ trees from each distance matrix
    // =========================================================================
    println!("=========================================================================");
    println!("  STEP 4: Build Neighbor-Joining trees from each distance matrix");
    println!("=========================================================================");
    println!();

    let jc69_tree = neighbor_joining(&jc69_matrix);
    let logdet_tree = neighbor_joining(&logdet_matrix);

    let jc69_newick = write_newick(&jc69_tree);
    let logdet_newick = write_newick(&logdet_tree);

    println!("  True tree:        ((A,B),(C,D))");
    println!("  NJ tree (JC69):   {}", jc69_newick);
    println!("  NJ tree (LogDet): {}", logdet_newick);
    println!();

    // Check which taxa are grouped together
    let jc69_ab_sisters = are_sisters_in_tree(&jc69_newick, "A", "B");
    let jc69_ac_sisters = are_sisters_in_tree(&jc69_newick, "A", "C");
    let logdet_ab_sisters = are_sisters_in_tree(&logdet_newick, "A", "B");

    println!("  Tree topology analysis:");
    if jc69_ac_sisters {
        println!("    JC69 groups A with C  --> WRONG (attracted by compositional similarity)");
    } else if jc69_ab_sisters {
        println!("    JC69 groups A with B  --> correct (but lucky -- bias was not strong enough)");
    } else {
        println!("    JC69 produces an unexpected grouping");
    }

    if logdet_ab_sisters {
        println!("    LogDet groups A with B --> CORRECT (robust to compositional bias)");
    } else {
        println!("    LogDet does not group A with B (bias may have overwhelmed the signal)");
    }
    println!();

    // =========================================================================
    // Step 5: Replicated experiment
    // =========================================================================
    println!("=========================================================================");
    println!("  STEP 5: Replicated experiment (50 replicates)");
    println!("=========================================================================");
    println!();
    println!("  For each replicate: simulate 2000 sites on ((A,B),(C,D)), then apply");
    println!("  GC-bias to A and C, build NJ trees from JC69 and LogDet distances.");
    println!();

    let num_replicates = 50;
    let base_seed: u64 = 161803398874;
    let bias_rates = [0.20, 0.30, 0.40];

    println!(
        "  {:>10} | {:>20} | {:>20}",
        "Bias rate", "JC69 correct (%)", "LogDet correct (%)"
    );
    println!(
        "  -----------|----------------------|----------------------"
    );

    for &br in &bias_rates {
        let mut jc69_correct = 0;
        let mut logdet_correct = 0;

        for rep in 0..num_replicates {
            let seed = base_seed
                .wrapping_add((br * 1000.0) as u64 * 999979)
                .wrapping_add(rep as u64 * 1000003);
            let mut sim_rng = Rng::new(seed);

            // Simulate alignment on true tree
            let sim_aln = simulate_alignment(&true_tree, nsites, &model, &mut sim_rng);

            // Apply GC-bias to A and C
            let mut bias_rng = Rng::new(seed.wrapping_add(42));
            let biased_a_rep = apply_gc_bias(&sim_aln, "A", br, &mut bias_rng);
            let biased_c_rep = apply_gc_bias(&sim_aln, "C", br, &mut bias_rng);
            let seq_b_rep = sim_aln.sequences.iter().find(|s| s.name == "B").unwrap().bases.clone();
            let seq_d_rep = sim_aln.sequences.iter().find(|s| s.name == "D").unwrap().bases.clone();

            let biased_aln = Alignment::new(vec![
                Sequence::new("A", biased_a_rep),
                Sequence::new("B", seq_b_rep),
                Sequence::new("C", biased_c_rep),
                Sequence::new("D", seq_d_rep),
            ])
            .expect("alignment should be valid");

            // JC69 tree
            if let Ok(jc69_dm) = compute_distance_matrix(&biased_aln, &jc69_model) {
                let jc69_t = neighbor_joining(&jc69_dm);
                let nwk = write_newick(&jc69_t);
                if are_sisters_in_tree(&nwk, "A", "B") {
                    jc69_correct += 1;
                }
            }

            // LogDet tree
            if let Ok(logdet_dm) = compute_distance_matrix(&biased_aln, &logdet_model) {
                let logdet_t = neighbor_joining(&logdet_dm);
                let nwk = write_newick(&logdet_t);
                if are_sisters_in_tree(&nwk, "A", "B") {
                    logdet_correct += 1;
                }
            }
        }

        let jc69_pct = 100.0 * jc69_correct as f64 / num_replicates as f64;
        let logdet_pct = 100.0 * logdet_correct as f64 / num_replicates as f64;
        println!(
            "  {:>9.0}% | {:>18.0}% | {:>18.0}%",
            br * 100.0, jc69_pct, logdet_pct
        );
    }
    println!();

    // =========================================================================
    // Interpretation
    // =========================================================================
    println!("=========================================================================");
    println!("  INTERPRETATION");
    println!("=========================================================================");
    println!();
    println!("  The JC69 model assumes that all four bases have equal frequency (0.25)");
    println!("  in all lineages. When this assumption is violated -- as happens when");
    println!("  some lineages become GC-rich -- JC69 underestimates the distance");
    println!("  between taxa that share similar composition (A and C in our experiment)");
    println!("  and overestimates the distance between taxa with different composition.");
    println!();
    println!("  This causes JC69 to group A with C (the two GC-rich taxa) instead of");
    println!("  A with B (the true sisters). This is compositional attraction, analogous");
    println!("  to long-branch attraction but caused by base frequency convergence.");
    println!();
    println!("  The LogDet distance avoids this problem by using the DETERMINANT of");
    println!("  the 4x4 divergence matrix F, where F[a][b] = proportion of sites with");
    println!("  base a in sequence 1 and base b in sequence 2. The key identity is:");
    println!();
    println!("    d = -(1/4) * ln(det(F)) + (1/8) * [ln(prod(f_i)) + ln(prod(f_j))]");
    println!();
    println!("  The marginal frequency terms (f_i, f_j) explicitly account for the");
    println!("  different base compositions of each sequence. This makes LogDet a");
    println!("  CONSISTENT estimator even under non-stationary evolution where base");
    println!("  frequencies drift independently in each lineage.");
    println!();
    println!("  This is why LogDet/Paralinear distances are recommended whenever");
    println!("  there is evidence of compositional heterogeneity across taxa.");
    println!();

    println!("=========================================================================");
    println!("  REFERENCES");
    println!("=========================================================================");
    println!();
    println!("  Lake, J.A. (1994). Reconstructing evolutionary trees from DNA and");
    println!("  protein sequences: Paralinear distances. PNAS, 91:1455-1459.");
    println!();
    println!("  Lockhart, P.J., Steel, M.A., Hendy, M.D. & Penny, D. (1994).");
    println!("  Recovering evolutionary trees under a more realistic model of");
    println!("  sequence evolution. Mol. Biol. Evol., 11:605-612.");
    println!();
    println!("  Steel, M.A. (1994). Recovering a tree from the leaf colourations");
    println!("  it generates under a Markov model. Applied Mathematics Letters,");
    println!("  7:19-23.");
    println!();
    println!("  Jukes, T.H. & Cantor, C.R. (1969). Evolution of protein molecules.");
    println!("  In Mammalian Protein Metabolism, Vol. III, pp. 21-132.");
    println!();
}

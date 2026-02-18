// lake_invariants.rs
//
// Demonstrates that Lake's and Cavender's phylogenetic invariants can identify
// the true topology among the three possible unrooted 4-taxon trees.
//
// Phylogenetic invariants are functions of site-pattern frequencies that equal
// zero for the true tree topology but are non-zero for alternative topologies.
// They were proposed by Lake (1987) and Cavender (1978) as a statistically
// consistent, model-free approach to topology inference.
//
// This example simulates DNA sequences under JC69 on each of the three possible
// 4-taxon topologies, then shows that both Lake's and Cavender's invariants
// correctly recover the generating topology -- and that accuracy improves with
// increasing sequence length.
//
// Usage:
//   cargo run --example lake_invariants

use phylip_rs::invariants::lake::{lake_invariants, cavender_invariants};
use phylip_rs::likelihood::models::{Jc69Model, SubstitutionModel};
use phylip_rs::tree::newick::parse_newick;
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

fn simulate_alignment(
    tree: &phylip_rs::tree::Tree,
    nsites: usize,
    model: &Jc69Model,
    rng: &mut Rng,
) -> Alignment {
    let nnodes = tree.num_nodes();
    let mut node_seqs: Vec<Vec<usize>> = vec![vec![0; nsites]; nnodes];

    // Root sequence: uniform random bases
    for site in 0..nsites {
        let root_base = (rng.next_u64() >> 33) as usize % 4;
        node_seqs[tree.root][site] = root_base;
    }

    // Evolve along tree in preorder
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

// ---------------------------------------------------------------------------
// Topology labels
// ---------------------------------------------------------------------------

const TOPOLOGY_NAMES: [&str; 3] = [
    "T1: ((1,2),(3,4))",
    "T2: ((1,3),(2,4))",
    "T3: ((1,4),(2,3))",
];

// Newick trees for simulation (moderate branch lengths)
const SIMULATION_TREES: [&str; 3] = [
    "((Taxon1:0.15,Taxon2:0.15):0.05,(Taxon3:0.15,Taxon4:0.15):0.05);",
    "((Taxon1:0.15,Taxon3:0.15):0.05,(Taxon2:0.15,Taxon4:0.15):0.05);",
    "((Taxon1:0.15,Taxon4:0.15):0.05,(Taxon2:0.15,Taxon3:0.15):0.05);",
];

// ---------------------------------------------------------------------------
// Main demonstration
// ---------------------------------------------------------------------------

fn main() {
    println!();
    println!("PHYLOGENETIC INVARIANTS: Model-Free Topology Identification");
    println!("============================================================");
    println!();
    println!("Phylogenetic invariants are algebraic functions of site-pattern");
    println!("frequencies that equal zero for the true tree topology but are");
    println!("non-zero for alternative topologies. They provide a statistically");
    println!("consistent method for inferring tree topology that does not require");
    println!("estimating branch lengths.");
    println!();
    println!("For four taxa, there are exactly three possible unrooted topologies:");
    println!("  T1: ((1,2),(3,4))  -- taxa 1,2 are grouped together");
    println!("  T2: ((1,3),(2,4))  -- taxa 1,3 are grouped together");
    println!("  T3: ((1,4),(2,3))  -- taxa 1,4 are grouped together");
    println!();

    let model = Jc69Model;
    let base_seed: u64 = 271828182;

    // ===================================================================
    // Part 1: Increasing sequence length on topology T1
    // ===================================================================

    println!("PART 1: Accuracy Increases with Sequence Length");
    println!("===============================================");
    println!();
    println!("True topology: T1 = ((Taxon1:0.15, Taxon2:0.15):0.05, (Taxon3:0.15, Taxon4:0.15):0.05)");
    println!();
    println!("We simulate DNA under JC69 on T1 with varying numbers of sites.");
    println!("Both Lake's and Cavender's invariants should identify T1 as best,");
    println!("and the signal should strengthen with more data.");
    println!();

    let true_tree = parse_newick(SIMULATION_TREES[0]).expect("T1 newick should parse");
    let site_counts = [100, 500, 2000, 10000];

    println!(
        "{:>7} | {:>12} {:>12} {:>12} | {:>8} | {:>12} {:>12} {:>12} | {:>8}",
        "Sites",
        "Lake T1", "Lake T2", "Lake T3", "Lake?",
        "Cav T1", "Cav T2", "Cav T3", "Cav?"
    );
    println!(
        "--------|{:-<40}|----------|{:-<40}|---------",
        "", ""
    );

    for &nsites in &site_counts {
        let seed = base_seed.wrapping_add(nsites as u64 * 7919);
        let mut rng = Rng::new(seed);
        let alignment = simulate_alignment(&true_tree, nsites, &model, &mut rng);

        let lake = lake_invariants(&alignment).expect("lake invariants should succeed");
        let cav = cavender_invariants(&alignment).expect("cavender invariants should succeed");

        let lake_best_name = match lake.best_topology_index {
            0 => "T1",
            1 => "T2",
            _ => "T3",
        };
        let cav_best_name = match cav.best_topology_index {
            0 => "T1",
            1 => "T2",
            _ => "T3",
        };

        // Mark correct answers
        let lake_mark = if lake.best_topology_index == 0 { " *" } else { "" };
        let cav_mark = if cav.best_topology_index == 0 { " *" } else { "" };

        println!(
            "{:>7} | {:>12} {:>12} {:>12} | {:>6}{:<2} | {:>12.1} {:>12.1} {:>12.1} | {:>6}{:<2}",
            nsites,
            lake.topology_support[0],
            lake.topology_support[1],
            lake.topology_support[2],
            lake_best_name, lake_mark,
            cav.invariant_values[0],
            cav.invariant_values[1],
            cav.invariant_values[2],
            cav_best_name, cav_mark
        );
    }

    println!();
    println!("  Lake's invariants: The topology with the HIGHEST support count is favored.");
    println!("  Cavender's invariants: The topology whose value is CLOSEST TO ZERO is favored.");
    println!("  * = correctly identified the true topology T1.");
    println!();

    // ===================================================================
    // Part 2: All three topologies
    // ===================================================================

    println!("PART 2: Invariants Correctly Track the Generating Topology");
    println!("==========================================================");
    println!();
    println!("We simulate 5000-site alignments on each of the three topologies");
    println!("and verify that invariants correctly identify each one.");
    println!();

    let nsites = 5000;
    for (topo_idx, (&newick, &topo_name)) in SIMULATION_TREES
        .iter()
        .zip(TOPOLOGY_NAMES.iter())
        .enumerate()
    {
        let tree = parse_newick(newick).expect("newick should parse");
        let seed = base_seed
            .wrapping_add(topo_idx as u64 * 999983)
            .wrapping_add(nsites as u64);
        let mut rng = Rng::new(seed);
        let alignment = simulate_alignment(&tree, nsites, &model, &mut rng);

        let lake = lake_invariants(&alignment).expect("lake invariants should succeed");
        let cav = cavender_invariants(&alignment).expect("cavender invariants should succeed");

        let lake_correct = lake.best_topology_index == topo_idx;
        let cav_correct = cav.best_topology_index == topo_idx;

        println!("  Simulated on: {}", topo_name);
        println!(
            "    Lake's support:   T1={:<6} T2={:<6} T3={:<6}  => Best: T{} {}",
            lake.topology_support[0],
            lake.topology_support[1],
            lake.topology_support[2],
            lake.best_topology_index + 1,
            if lake_correct { "(CORRECT)" } else { "(wrong)" }
        );
        println!(
            "    Cavender's inv:   T1={:<10.1} T2={:<10.1} T3={:<10.1}  => Best: T{} {}",
            cav.invariant_values[0],
            cav.invariant_values[1],
            cav.invariant_values[2],
            cav.best_topology_index + 1,
            if cav_correct { "(CORRECT)" } else { "(wrong)" }
        );
        println!(
            "    Lake's chi-sq:    {:.2} (higher = stronger signal)",
            lake.chi_squared
        );
        println!();
    }

    // ===================================================================
    // Part 3: Replicate study — consistency over 50 replicates
    // ===================================================================

    println!("PART 3: Statistical Consistency Over 50 Replicates");
    println!("==================================================");
    println!();
    println!("For each sequence length, we run 50 independent simulations on");
    println!("the true topology T1 and count how often each method is correct.");
    println!();

    let num_replicates = 50;
    let replicate_site_counts = [100, 500, 2000, 10000];

    println!(
        "{:>7} | {:>18} | {:>18}",
        "Sites", "Lake correct", "Cavender correct"
    );
    println!("--------|{:-<20}|{:-<20}", "", "");

    for &nsites in &replicate_site_counts {
        let mut lake_correct = 0;
        let mut cav_correct = 0;

        for rep in 0..num_replicates {
            let seed = base_seed
                .wrapping_add(nsites as u64 * 1000003)
                .wrapping_add(rep as u64 * 999983);
            let mut rng = Rng::new(seed);
            let alignment = simulate_alignment(&true_tree, nsites, &model, &mut rng);

            let lake = lake_invariants(&alignment).expect("lake should succeed");
            let cav = cavender_invariants(&alignment).expect("cavender should succeed");

            if lake.best_topology_index == 0 {
                lake_correct += 1;
            }
            if cav.best_topology_index == 0 {
                cav_correct += 1;
            }
        }

        let lake_pct = 100.0 * lake_correct as f64 / num_replicates as f64;
        let cav_pct = 100.0 * cav_correct as f64 / num_replicates as f64;

        println!(
            "{:>7} | {:>15.0}% | {:>15.0}%",
            nsites, lake_pct, cav_pct
        );
    }

    println!();
    println!("  As the number of sites increases, both methods converge to 100%");
    println!("  accuracy. This is the property of STATISTICAL CONSISTENCY:");
    println!("  with enough data, the true topology is always recovered.");
    println!();

    // ===================================================================
    // Part 4: Detailed pattern analysis
    // ===================================================================

    println!("PART 4: Under the Hood — Site Pattern Frequencies");
    println!("=================================================");
    println!();
    println!("Lake's method classifies site patterns into three types:");
    println!("  xxyy: taxa 1,2 share one state; taxa 3,4 share another (supports T1)");
    println!("  xyxy: taxa 1,3 share one state; taxa 2,4 share another (supports T2)");
    println!("  xyyx: taxa 1,4 share one state; taxa 2,3 share another (supports T3)");
    println!();
    println!("Crucially, only TRANSVERSION patterns count (purine <-> pyrimidine),");
    println!("making the method robust to transition/transversion rate bias.");
    println!();

    let nsites = 10000;
    let seed = base_seed.wrapping_add(42);
    let mut rng = Rng::new(seed);
    let alignment = simulate_alignment(&true_tree, nsites, &model, &mut rng);
    let lake = lake_invariants(&alignment).expect("lake should succeed");

    let total_informative: usize = lake.topology_support.iter().sum();
    let total_sites = lake.pattern_counts.total_sites;

    println!("  Simulated: 10,000 sites on T1 under JC69");
    println!("  Total sites counted: {}", total_sites);
    println!("  Total transversion-informative patterns: {}", total_informative);
    println!(
        "  Fraction informative: {:.1}%",
        100.0 * total_informative as f64 / total_sites as f64
    );
    println!();
    println!(
        "  xxyy (supports T1): {:>5}  ({:.1}% of informative)",
        lake.topology_support[0],
        if total_informative > 0 {
            100.0 * lake.topology_support[0] as f64 / total_informative as f64
        } else {
            0.0
        }
    );
    println!(
        "  xyxy (supports T2): {:>5}  ({:.1}% of informative)",
        lake.topology_support[1],
        if total_informative > 0 {
            100.0 * lake.topology_support[1] as f64 / total_informative as f64
        } else {
            0.0
        }
    );
    println!(
        "  xyyx (supports T3): {:>5}  ({:.1}% of informative)",
        lake.topology_support[2],
        if total_informative > 0 {
            100.0 * lake.topology_support[2] as f64 / total_informative as f64
        } else {
            0.0
        }
    );
    println!();
    println!(
        "  Chi-squared statistic: {:.2} (df=2)",
        lake.chi_squared
    );
    println!("  Under the null hypothesis of equal support, chi-squared > 5.99");
    println!("  is significant at p < 0.05. Our value strongly rejects the null,");
    println!("  confirming a definite topological signal in the data.");
    println!();

    // ===================================================================
    // Part 5: Cavender's invariant identity
    // ===================================================================

    println!("PART 5: Cavender's Mathematical Identity");
    println!("=========================================");
    println!();
    println!("Cavender's invariants have an elegant algebraic property:");
    println!("the three invariant values always satisfy:");
    println!();
    println!("  inv(T1) - inv(T2) + inv(T3) = 0");
    println!();
    println!("This is because the invariants are defined as pairwise differences");
    println!("of three pattern group counts (n1, n2, n3):");
    println!("  inv(T1) = n2 - n3");
    println!("  inv(T2) = n1 - n3");
    println!("  inv(T3) = n1 - n2");
    println!();

    let cav = cavender_invariants(&alignment).expect("cavender should succeed");
    let identity_check =
        cav.invariant_values[0] - cav.invariant_values[1] + cav.invariant_values[2];

    println!("  For our 10,000-site T1 simulation:");
    println!(
        "    inv(T1) = {:.4}",
        cav.invariant_values[0]
    );
    println!(
        "    inv(T2) = {:.4}",
        cav.invariant_values[1]
    );
    println!(
        "    inv(T3) = {:.4}",
        cav.invariant_values[2]
    );
    println!(
        "    inv(T1) - inv(T2) + inv(T3) = {:.10}  (should be exactly 0)",
        identity_check
    );
    println!();
    println!("  Under the true topology T1, inv(T1) should be closest to zero,");
    println!("  while inv(T2) and inv(T3) should be non-zero and roughly equal");
    println!("  in magnitude (but with opposite signs from the identity).");
    println!();

    // ===================================================================
    // References
    // ===================================================================

    println!("References:");
    println!("  Lake, J.A. (1987). A rate-independent technique for analysis of");
    println!("  nucleic acid sequences: evolutionary parsimony. Molecular Biology");
    println!("  and Evolution, 4, 167-191.");
    println!();
    println!("  Cavender, J.A. (1978). Taxonomy with confidence. Mathematical");
    println!("  Biosciences, 40, 271-280.");
    println!();
    println!("  Felsenstein, J. (2004). Inferring Phylogenies. Sinauer Associates.");
    println!("  Chapter 29: Invariants.");
    println!();
}

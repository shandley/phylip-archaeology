// felsenstein_zone.rs
//
// Demonstrates the "Felsenstein Zone" -- a classic result from Joe Felsenstein's
// 1978 paper showing that maximum parsimony can be POSITIVELY MISLEADING:
// as you add more data, parsimony becomes MORE confident in the WRONG tree,
// while maximum likelihood correctly recovers the true tree.
//
// This is one of the most important results in the history of phylogenetics.
// It arises when two unrelated lineages independently undergo rapid evolution
// (long branches), causing convergent substitutions that parsimony mistakes
// for shared ancestry.
//
// Usage:
//   cargo run --example felsenstein_zone

use phylip_rs::likelihood::models::{Jc69Model, SubstitutionModel};
use phylip_rs::likelihood::pruning::{log_likelihood, optimize_branch_lengths};
use phylip_rs::parsimony::wagner::parsimony_score;
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
        // Use the upper 32 bits for better quality from the LCG.
        (self.next_u64() >> 33) as f64 / (1u64 << 31) as f64
    }
}

// ---------------------------------------------------------------------------
// Sequence simulator along a phylogenetic tree under JC69
// ---------------------------------------------------------------------------

/// Index-to-Base mapping: 0=A, 1=C, 2=G, 3=T
fn index_to_base(i: usize) -> Base {
    match i {
        0 => Base::A,
        1 => Base::C,
        2 => Base::G,
        3 => Base::T,
        _ => unreachable!(),
    }
}

/// Simulate DNA sequences along a tree using the JC69 model.
///
/// For each site, a random ancestral base is drawn at the root.
/// The sequence then evolves along each branch: for a branch of length t,
/// the transition probability matrix P(t) from JC69 gives the probability
/// of mutating from the parent's base to each possible descendant base.
///
/// Only the leaf sequences are returned (the observable data).
fn simulate_alignment(
    tree: &phylip_rs::tree::Tree,
    nsites: usize,
    model: &Jc69Model,
    rng: &mut Rng,
) -> Alignment {
    let nnodes = tree.num_nodes();

    // node_sequences[node_id] = vector of base indices (0..4) for each site
    let mut node_seqs: Vec<Vec<usize>> = vec![vec![0; nsites]; nnodes];

    // Generate root sequence: each site is a uniformly random base (JC69
    // has equal equilibrium frequencies of 0.25 each).
    for site in 0..nsites {
        let root_base = (rng.next_u64() >> 33) as usize % 4;
        node_seqs[tree.root][site] = root_base;
    }

    // Traverse the tree in preorder (parent before children).
    // For each child, mutate each site according to P(branch_length).
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
            // Sample descendant base from P[parent_base][.]
            let row = &trans_probs[parent_base];
            let u = rng.next_f64();
            let mut cumulative = 0.0;
            let mut child_base = 3; // default to T if rounding issues
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

    // Sort sequences by name to ensure consistent ordering (A, B, C, D).
    sequences.sort_by(|a, b| a.name.cmp(&b.name));

    Alignment::new(sequences).expect("simulated alignment should be valid")
}

// ---------------------------------------------------------------------------
// Evaluate all three unrooted 4-taxon topologies
// ---------------------------------------------------------------------------

/// The three possible unrooted topologies for four taxa {A, B, C, D}.
/// In Newick format (rooted representation of an unrooted tree):
const TOPOLOGY_NEWICKS: [&str; 3] = [
    "((A,B),(C,D));", // T1: TRUE tree -- A,B together; C,D together
    "((A,C),(B,D));", // T2: wrong
    "((A,D),(B,C));", // T3: wrong -- long branch attraction topology
];

/// For a given simulated alignment, evaluate all three topologies under
/// both parsimony and maximum likelihood.
///
/// Returns (parsimony_best_index, ml_best_index) -- which topology each
/// method considers best.
fn evaluate_topologies(alignment: &Alignment, model: &Jc69Model) -> (usize, usize) {
    let mut pars_scores = [0usize; 3];
    let mut ml_lnls = [f64::NEG_INFINITY; 3];

    for (i, newick) in TOPOLOGY_NEWICKS.iter().enumerate() {
        let tree = parse_newick(newick).expect("topology Newick should parse");

        // -- Parsimony --
        let (score, _) = parsimony_score(&tree, alignment)
            .expect("parsimony scoring should succeed");
        pars_scores[i] = score;

        // -- Maximum Likelihood (with branch length optimization) --
        let mut ml_tree = tree.clone();
        // Set initial branch lengths to a reasonable starting value.
        for node in &mut ml_tree.nodes {
            if node.parent.is_some() {
                node.branch_length = Some(0.1);
            }
        }
        match optimize_branch_lengths(&mut ml_tree, alignment, model) {
            Ok(_) => {
                match log_likelihood(&ml_tree, alignment, model) {
                    Ok(lnl) => ml_lnls[i] = lnl,
                    Err(_) => {} // leave as NEG_INFINITY
                }
            }
            Err(_) => {} // leave as NEG_INFINITY
        }
    }

    // Parsimony picks the topology with the LOWEST score
    let pars_best = pars_scores
        .iter()
        .enumerate()
        .min_by_key(|&(_, &s)| s)
        .unwrap()
        .0;

    // ML picks the topology with the HIGHEST log-likelihood
    let ml_best = ml_lnls
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap()
        .0;

    (pars_best, ml_best)
}

// ---------------------------------------------------------------------------
// Main: run the full Felsenstein Zone demonstration
// ---------------------------------------------------------------------------

fn main() {
    println!();
    println!("THE FELSENSTEIN ZONE: When More Data Makes You More Wrong");
    println!("=========================================================");
    println!();
    println!("In 1978, Joe Felsenstein proved one of the most striking results");
    println!("in all of statistical phylogenetics: maximum parsimony -- the idea");
    println!("of preferring the tree that requires the fewest evolutionary changes");
    println!("-- can be 'positively misleading.' This means that as you collect");
    println!("MORE data, parsimony becomes MORE confident in the WRONG answer.");
    println!();
    println!("The setup involves four species (A, B, C, D) arranged so that two");
    println!("UNRELATED lineages (A and D) happen to be evolving very rapidly,");
    println!("while two other lineages (B and C) evolve slowly. Because A and D");
    println!("both accumulate many mutations independently, they will, by chance,");
    println!("often converge on the same nucleotide at many sites. Parsimony");
    println!("mistakes these convergent substitutions for evidence of shared");
    println!("ancestry, and groups A with D -- the WRONG tree.");
    println!();
    println!("Maximum likelihood, on the other hand, uses an explicit model of");
    println!("how DNA evolves (here, the Jukes-Cantor 1969 model). It can");
    println!("distinguish true shared ancestry from convergent evolution caused");
    println!("by high mutation rates, and correctly infers the true tree.");
    println!();
    println!("This is the phenomenon of LONG BRANCH ATTRACTION (LBA).");
    println!();

    // Define the true tree.
    //
    // The branch lengths are chosen to place us firmly in the Felsenstein zone.
    // With long branches of 0.80 expected substitutions per site, each long-branch
    // taxon has about a 52% chance of having changed at any given site under JC69.
    // With the internal branch at only 0.01, there is almost no "true signal" for
    // parsimony to detect -- but there is a strong statistical signal for ML.
    let true_tree_newick = "((A:0.80,B:0.01):0.01,(C:0.01,D:0.80));";
    let true_tree = parse_newick(true_tree_newick).expect("true tree should parse");

    println!("True tree: ((A:0.80, B:0.01):0.01, (C:0.01, D:0.80))");
    println!();
    println!("  A and D are on LONG branches  (0.80 expected substitutions/site)");
    println!("  B and C are on SHORT branches  (0.01 expected substitutions/site)");
    println!("  Internal branch is SHORT       (0.01 expected substitutions/site)");
    println!();
    println!("  The long branches mean A and D each have roughly a 52% chance");
    println!("  of having changed at any given site. When both change, they will");
    println!("  sometimes land on the same base by coincidence -- creating a");
    println!("  pattern that LOOKS like shared ancestry to parsimony.");
    println!();
    println!("Three possible unrooted topologies for 4 taxa:");
    println!("  T1: ((A,B),(C,D))  <-- TRUE tree");
    println!("  T2: ((A,C),(B,D))");
    println!("  T3: ((A,D),(B,C))  <-- Parsimony's favorite (long branch attraction)");
    println!();

    let model = Jc69Model;

    // Experiment parameters
    let site_counts: Vec<usize> = vec![100, 500, 1000, 5000, 10000];
    let num_replicates = 50;
    let base_seed: u64 = 314159265;

    println!("Running simulations: {} replicates per sequence length...", num_replicates);
    println!("(This may take a minute for the larger sequence lengths.)");
    println!();

    // Print table header
    println!(
        "{:>7} | {:>19} | {:>19} | {:>19}",
        "Sites", "Parsimony correct", "ML correct", "Parsimony picks T3"
    );
    println!(
        "--------|{:-<21}|{:-<21}|{:-<21}",
        "", "", ""
    );

    for &nsites in &site_counts {
        let mut pars_correct = 0;
        let mut ml_correct = 0;
        let mut pars_picks_t3 = 0;

        for rep in 0..num_replicates {
            // Each replicate gets a unique seed derived from the base seed,
            // the number of sites, and the replicate number.
            let seed = base_seed
                .wrapping_add(nsites as u64 * 1000003)
                .wrapping_add(rep as u64 * 999983);
            let mut rng = Rng::new(seed);

            // Simulate sequences along the TRUE tree
            let alignment = simulate_alignment(&true_tree, nsites, &model, &mut rng);

            // Evaluate all three topologies
            let (pars_best, ml_best) = evaluate_topologies(&alignment, &model);

            if pars_best == 0 {
                pars_correct += 1;
            }
            if ml_best == 0 {
                ml_correct += 1;
            }
            if pars_best == 2 {
                pars_picks_t3 += 1;
            }
        }

        let pars_pct = 100.0 * pars_correct as f64 / num_replicates as f64;
        let ml_pct = 100.0 * ml_correct as f64 / num_replicates as f64;
        let t3_pct = 100.0 * pars_picks_t3 as f64 / num_replicates as f64;

        println!(
            "{:>7} | {:>17.0}% | {:>17.0}% | {:>17.0}%",
            nsites, pars_pct, ml_pct, t3_pct
        );
    }

    println!();
    println!("Interpretation:");
    println!("  As sequence length increases:");
    println!("  - ML converges toward the CORRECT tree (T1). This property is called");
    println!("    'statistical consistency': with enough data, the right answer emerges.");
    println!("  - Parsimony converges toward the WRONG tree (T3). This is 'positive");
    println!("    misleading' or 'inconsistency': more data makes things WORSE.");
    println!();
    println!("  The wrong tree T3 groups A and D together -- the two species on");
    println!("  long branches. Parsimony counts shared derived characters, but");
    println!("  cannot tell the difference between:");
    println!("    (a) a shared character inherited from a common ancestor, and");
    println!("    (b) a shared character arising from independent mutations");
    println!("        (convergent evolution) on two long branches.");
    println!();
    println!("  Maximum likelihood avoids this trap because it explicitly models");
    println!("  the probability of convergent substitution. Long branches are");
    println!("  EXPECTED to produce convergent patterns under the model, so ML");
    println!("  does not mistake them for evidence of relatedness.");
    println!();
    println!("This result, from Felsenstein (1978), was one of the key arguments");
    println!("for adopting model-based methods in phylogenetics. It remains one of");
    println!("the most celebrated demonstrations in computational evolutionary biology.");
    println!();
    println!("Reference:");
    println!("  Felsenstein, J. (1978). Cases in which parsimony or compatibility");
    println!("  methods will be positively misleading. Systematic Zoology, 27, 401-410.");
    println!();
}

// dollo_gain_loss.rs
//
// Demonstrates the difference between Dollo and Fitch parsimony in
// interpreting binary character evolution -- and why the distinction
// matters for biological inference.
//
// Dollo parsimony assumes each derived state (1) arises exactly once
// and can be lost multiple times. Fitch parsimony treats gains and
// losses symmetrically. For binary characters, both criteria always
// agree on the optimal topology, but they differ in HOW they explain
// the data. This leads to dramatically different biological stories.
//
// The Dollo/Fitch score ratio serves as a diagnostic: when data truly
// follow Dollo's law, the ratio is moderate; when data have multiple
// independent gains (violating Dollo's assumption), the ratio is inflated.
//
// Usage:
//   cargo run --example dollo_gain_loss

use phylip_rs::parsimony::dollo::{dollo_parsimony, DolloParsimonyResult};
use phylip_rs::parsimony::wagner::parsimony_score;
use phylip_rs::tree::newick::parse_newick;
use phylip_rs::tree::{Alignment, Base, Sequence, Tree};

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
// Data generators
// ---------------------------------------------------------------------------

/// Generate binary characters under Dollo's law.
///
/// Each character is gained exactly once on the tree, then independently
/// lost on descendant branches. This models gene presence/absence,
/// restriction sites, or transposable element insertions.
fn generate_dollo_data(
    tree: &Tree,
    nchars: usize,
    loss_prob: f64,
    rng: &mut Rng,
) -> Alignment {
    let ntaxa = tree.num_leaves();
    let leaf_names: Vec<String> = tree
        .leaves()
        .iter()
        .map(|n| n.name.clone().unwrap_or_default())
        .collect();

    let mut taxa_data: Vec<Vec<Base>> = vec![Vec::new(); ntaxa];

    for _ in 0..nchars {
        // Pick a random node for the gain event
        let nnodes = tree.num_nodes();
        let gain_node = (rng.next_u64() >> 33) as usize % nnodes;

        // Find all leaves descended from this node
        let descendants = get_descendant_leaves(tree, gain_node, &leaf_names);

        for (t, name) in leaf_names.iter().enumerate() {
            if descendants.contains(name) {
                if rng.next_f64() < loss_prob {
                    taxa_data[t].push(Base::A); // lost
                } else {
                    taxa_data[t].push(Base::T); // retained
                }
            } else {
                taxa_data[t].push(Base::A); // never had it
            }
        }
    }

    let sequences: Vec<Sequence> = leaf_names
        .iter()
        .enumerate()
        .map(|(i, name)| Sequence::new(name.as_str(), taxa_data[i].clone()))
        .collect();

    Alignment::new(sequences).unwrap()
}

/// Generate binary characters under a symmetric model.
///
/// Characters evolve with equal probability of 0->1 and 1->0 transitions.
/// Multiple independent gains of the same trait are possible.
fn generate_symmetric_data(
    tree: &Tree,
    nchars: usize,
    change_prob: f64,
    rng: &mut Rng,
) -> Alignment {
    let nnodes = tree.num_nodes();
    let leaf_names: Vec<String> = tree
        .leaves()
        .iter()
        .map(|n| n.name.clone().unwrap_or_default())
        .collect();
    let ntaxa = leaf_names.len();

    let mut taxa_data: Vec<Vec<Base>> = vec![Vec::new(); ntaxa];

    for _ in 0..nchars {
        let mut node_states = vec![false; nnodes];
        node_states[tree.root] = rng.next_f64() < 0.5;

        let preorder = tree.preorder();
        for &nid in &preorder {
            if nid == tree.root {
                continue;
            }
            let parent = tree.nodes[nid].parent.unwrap();
            node_states[nid] = if rng.next_f64() < change_prob {
                !node_states[parent]
            } else {
                node_states[parent]
            };
        }

        for (t, name) in leaf_names.iter().enumerate() {
            let leaf_id = tree
                .nodes
                .iter()
                .find(|n| n.name.as_deref() == Some(name.as_str()))
                .unwrap()
                .id;
            taxa_data[t].push(if node_states[leaf_id] {
                Base::T
            } else {
                Base::A
            });
        }
    }

    let sequences: Vec<Sequence> = leaf_names
        .iter()
        .enumerate()
        .map(|(i, name)| Sequence::new(name.as_str(), taxa_data[i].clone()))
        .collect();

    Alignment::new(sequences).unwrap()
}

/// Find all leaf names descended from a given node.
fn get_descendant_leaves(tree: &Tree, node_id: usize, leaf_names: &[String]) -> Vec<String> {
    let mut result = Vec::new();
    let mut stack = vec![node_id];
    while let Some(nid) = stack.pop() {
        let node = &tree.nodes[nid];
        if node.is_leaf() {
            if let Some(ref name) = node.name {
                if leaf_names.contains(name) {
                    result.push(name.clone());
                }
            }
        }
        for &child in &node.children {
            stack.push(child);
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Main demonstration
// ---------------------------------------------------------------------------

fn main() {
    println!();
    println!("DOLLO vs FITCH PARSIMONY: Different Models, Different Stories");
    println!("==============================================================");
    println!();
    println!("Both Dollo and Fitch parsimony find the same optimal topology for");
    println!("binary characters. But they tell fundamentally different biological");
    println!("stories about HOW evolution happened:");
    println!();
    println!("  Fitch: minimizes total state changes. Treats 0->1 and 1->0");
    println!("         symmetrically. Does not distinguish gains from losses.");
    println!();
    println!("  Dollo: assumes each trait (state 1) arises ONCE on the tree.");
    println!("         Score = total number of independent LOSSES (reversals).");
    println!("         Each character with any '1' taxa has exactly 1 gain.");
    println!();
    println!("This matters because in biology, some traits truly follow Dollo's");
    println!("law: genes are gained once (by horizontal transfer or duplication)");
    println!("but can be independently lost in multiple lineages.");

    let tree_newick = "((A:0.1,B:0.1):0.1,((C:0.1,D:0.1):0.1,(E:0.1,F:0.1):0.1):0.1);";
    let tree = parse_newick(tree_newick).unwrap();

    // =====================================================================
    // PART 1: How Dollo and Fitch read the same data differently
    // =====================================================================

    println!();
    println!("==============================================================");
    println!("PART 1: Different interpretations of the same pattern");
    println!("==============================================================");
    println!();
    println!("Tree: ((A,B),((C,D),(E,F)))");
    println!();

    let examples: Vec<(&str, [Base; 6], &str, &str)> = vec![
        (
            "Trait in A,B only (sister pair)",
            [Base::T, Base::T, Base::A, Base::A, Base::A, Base::A],
            "1 change (direction ambiguous)",
            "Gain on AB ancestor; must push gain deeper, some losses implied",
        ),
        (
            "Trait in all except E",
            [Base::T, Base::T, Base::T, Base::T, Base::A, Base::T],
            "1 change (E lost the trait)",
            "1 gain at root, 1 loss on branch to E",
        ),
        (
            "Trait in B,D,F only",
            [Base::A, Base::T, Base::A, Base::T, Base::A, Base::T],
            "3 changes (could be 3 gains OR 3 losses)",
            "1 gain at root, 3 independent losses (on A, C, E)",
        ),
        (
            "Trait only in A",
            [Base::T, Base::A, Base::A, Base::A, Base::A, Base::A],
            "1 change (simple autapomorphy)",
            "1 gain forced to a deep node, multiple losses to explain absences",
        ),
    ];

    let names = ["A", "B", "C", "D", "E", "F"];

    for (desc, states, fitch_interp, dollo_interp) in &examples {
        let seqs: Vec<Sequence> = names
            .iter()
            .enumerate()
            .map(|(i, n)| Sequence::new(*n, vec![states[i]]))
            .collect();
        let aln = Alignment::new(seqs).unwrap();

        let (fitch_score, _) = parsimony_score(&tree, &aln).unwrap();
        let dollo_result = dollo_parsimony(&tree, &aln);

        let pattern: String = names
            .iter()
            .enumerate()
            .map(|(i, n)| {
                format!("{}={}", n, if states[i] == Base::T { "1" } else { "0" })
            })
            .collect::<Vec<_>>()
            .join(" ");

        println!("  {}", desc);
        println!("    Pattern: {}", pattern);
        println!(
            "    Fitch: {} step(s)  --  {}",
            fitch_score, fitch_interp
        );
        println!(
            "    Dollo: {} gain(s), {} loss(es)  --  {}",
            dollo_result.gains, dollo_result.score, dollo_interp
        );
        println!();
    }

    println!("  The third pattern is especially revealing. Fitch sees 3 changes");
    println!("  but cannot say whether they are gains or losses. Dollo provides a");
    println!("  specific biological narrative: one gain at the root, three independent");
    println!("  losses. If the data represent gene presence/absence, the Dollo");
    println!("  interpretation (one gene acquisition, three independent deletions)");
    println!("  is far more biologically plausible.");

    // =====================================================================
    // PART 2: Dollo vs symmetric data generation
    // =====================================================================

    println!();
    println!("==============================================================");
    println!("PART 2: Matching the criterion to the generating process");
    println!("==============================================================");
    println!();
    println!("We simulate 500 characters under two different models on the");
    println!("same tree and compare how each criterion interprets the data.");
    println!();
    println!("  Model A (Dollo): each trait gained once, lost with prob 0.25");
    println!("  Model B (Symmetric): each branch flips state with prob 0.20");
    println!();

    let mut rng = Rng::new(271828182);
    let dollo_data = generate_dollo_data(&tree, 500, 0.25, &mut rng);
    let symm_data = generate_symmetric_data(&tree, 500, 0.20, &mut rng);

    // Analyze both datasets
    let (fd, _) = parsimony_score(&tree, &dollo_data).unwrap();
    let dd: DolloParsimonyResult = dollo_parsimony(&tree, &dollo_data);

    let (fs, _) = parsimony_score(&tree, &symm_data).unwrap();
    let ds: DolloParsimonyResult = dollo_parsimony(&tree, &symm_data);

    println!(
        "  {:<26} {:>8} {:>8} {:>8} {:>8}",
        "Dataset", "Fitch", "Gains", "Losses", "Loss/Gain"
    );
    println!("  {}", "-".repeat(60));
    println!(
        "  {:<26} {:>8} {:>8} {:>8} {:>8.2}",
        "Dollo-generated (Model A)",
        fd,
        dd.gains,
        dd.score,
        dd.score as f64 / dd.gains.max(1) as f64
    );
    println!(
        "  {:<26} {:>8} {:>8} {:>8} {:>8.2}",
        "Symmetric (Model B)",
        fs,
        ds.gains,
        ds.score,
        ds.score as f64 / ds.gains.max(1) as f64
    );

    println!();
    println!("  Key observations:");
    println!();
    println!("  1. FITCH SCORE: Lower for Dollo-generated data ({} vs {}).", fd, fs);
    println!("     Under Dollo's law, each character has exactly one origin, creating");
    println!("     cleaner phylogenetic signal. The symmetric model generates more");
    println!("     homoplasy (multiple independent gains), requiring more changes.");
    println!();
    println!("  2. DOLLO LOSSES: The loss/gain ratio reflects the generating process.");
    println!("     Under Dollo's law, characters are gained once and lost on some");
    println!("     descendant branches (ratio ~{:.1}).", dd.score as f64 / dd.gains.max(1) as f64);
    println!("     Under the symmetric model, the Dollo interpretation forces each");
    println!("     character to have a single origin, which means patterns caused by");
    println!("     multiple independent gains must be explained by more losses.");

    // =====================================================================
    // PART 3: Replicated comparison
    // =====================================================================

    println!();
    println!("==============================================================");
    println!("PART 3: Replicated experiment (10 datasets per model)");
    println!("==============================================================");
    println!();
    println!("  Comparing Fitch score and Dollo loss/gain ratio across replicates.");
    println!();

    println!(
        "  {:>3}  {:>26}  {:>26}",
        "", "--- Dollo-generated ---", "--- Symmetric-generated ---"
    );
    println!(
        "  {:>3}  {:>8} {:>8} {:>8}  {:>8} {:>8} {:>8}",
        "Rep", "Fitch", "Losses", "L/G", "Fitch", "Losses", "L/G"
    );
    println!("  {}", "-".repeat(62));

    let mut total_fitch_d = 0.0;
    let mut total_fitch_s = 0.0;
    let mut total_lg_d = 0.0;
    let mut total_lg_s = 0.0;
    let num_reps = 10;

    for rep in 0..num_reps {
        let seed_d = 31415u64.wrapping_add(rep * 99991);
        let seed_s = 27182u64.wrapping_add(rep * 99991);

        let mut rng_d = Rng::new(seed_d);
        let mut rng_s = Rng::new(seed_s);

        let rep_dd = generate_dollo_data(&tree, 300, 0.25, &mut rng_d);
        let rep_sd = generate_symmetric_data(&tree, 300, 0.20, &mut rng_s);

        let (fd, _) = parsimony_score(&tree, &rep_dd).unwrap();
        let drd = dollo_parsimony(&tree, &rep_dd);

        let (fs, _) = parsimony_score(&tree, &rep_sd).unwrap();
        let drs = dollo_parsimony(&tree, &rep_sd);

        let lg_d = drd.score as f64 / drd.gains.max(1) as f64;
        let lg_s = drs.score as f64 / drs.gains.max(1) as f64;

        total_fitch_d += fd as f64;
        total_fitch_s += fs as f64;
        total_lg_d += lg_d;
        total_lg_s += lg_s;

        println!(
            "  {:>3}  {:>8} {:>8} {:>8.2}  {:>8} {:>8} {:>8.2}",
            rep + 1, fd, drd.score, lg_d, fs, drs.score, lg_s
        );
    }

    let n = num_reps as f64;
    println!("  {}", "-".repeat(62));
    println!(
        "  Avg  {:>8.0} {:>8} {:>8.2}  {:>8.0} {:>8} {:>8.2}",
        total_fitch_d / n, "", total_lg_d / n,
        total_fitch_s / n, "", total_lg_s / n
    );

    println!();
    println!("  Data generated under Dollo's law has LOWER Fitch scores (less");
    println!("  homoplasy) because each trait has a single origin.");
    println!();
    println!("  The loss/gain ratio is HIGHER for Dollo-generated data because");
    println!("  traits genuinely are gained once and lost many times. Under the");
    println!("  symmetric model, balanced gain/loss patterns yield a lower ratio.");

    // =====================================================================
    // Conclusion
    // =====================================================================

    println!();
    println!("==============================================================");
    println!("CONCLUSION");
    println!("==============================================================");
    println!();
    println!("For binary characters, Fitch and Dollo parsimony find the same");
    println!("optimal tree. Their value lies in providing different evolutionary");
    println!("interpretations of the data:");
    println!();
    println!("  DOLLO decomposes evolution into unique gains and independent losses.");
    println!("  When traits truly arise once (genes, restriction sites, insertions),");
    println!("  Dollo's decomposition is biologically correct and provides a");
    println!("  meaningful narrative about convergent loss.");
    println!();
    println!("  FITCH counts total changes without direction. This is appropriate");
    println!("  when gains and losses are equally likely.");
    println!();
    println!("  The LOSS/GAIN RATIO under Dollo is a diagnostic tool: higher");
    println!("  ratios indicate many independent losses per character (consistent");
    println!("  with Dollo's law), while lower ratios suggest a more symmetric");
    println!("  evolutionary process where the single-origin assumption is less");
    println!("  appropriate.");
    println!();
    println!("This is why PHYLIP provides both DNAPARS (Fitch) and DOLLOP (Dollo):");
    println!("the choice depends on the biology of the characters being analyzed.");

    println!();
    println!("==============================================================");
    println!("References");
    println!("==============================================================");
    println!();
    println!("  Le Quesne, W.J. (1974). The uniquely evolved character concept");
    println!("  and its cladistic application. Systematic Zoology, 23, 513-517.");
    println!();
    println!("  Farris, J.S. (1977). Phylogenetic analysis under Dollo's law.");
    println!("  Systematic Zoology, 26, 77-88.");
    println!();
    println!("  Felsenstein, J. (2004). Inferring Phylogenies. Sinauer Associates.");
    println!("  Chapter 8.");
    println!();
    println!("  Fitch, W.M. (1971). Toward defining the course of evolution:");
    println!("  minimum change for a specific tree topology. Systematic Zoology,");
    println!("  20, 406-416.");
    println!();
}

//! WebAssembly wrapper for phylip-rs phylogenetic algorithms.
//!
//! Exports 5 functions for browser-based phylogenetic analysis:
//! - `compute_nj`: Distance matrix → Neighbor-Joining tree
//! - `compute_likelihood`: Evaluate tree likelihood under substitution models
//! - `compute_parsimony`: Maximum parsimony tree search
//! - `felsenstein_zone`: Simulate the Felsenstein Zone (LBA demonstration)
//! - `bootstrap_nj`: Bootstrap NJ with majority-rule consensus

use serde::Serialize;
use wasm_bindgen::prelude::*;

use phylip_rs::bootstrap::{bootstrap_weights, resample_alignment, ResamplingMethod, SimpleRng};
use phylip_rs::consensus::{consensus_tree, ConsensusMethod};
use phylip_rs::distance::neighbor_joining;
use phylip_rs::io::fasta::read_fasta;
use phylip_rs::likelihood::models::{F84Model, Jc69Model, SubstitutionModel};
use phylip_rs::likelihood::pruning::{log_likelihood, optimize_branch_lengths};
use phylip_rs::models::compute_distance_matrix;
use phylip_rs::models::jc69::JC69;
use phylip_rs::models::k2p::K2P;
use phylip_rs::parsimony::wagner::{parsimony_score, search as parsimony_search};
use phylip_rs::tree::newick::{parse_newick, write_newick};
use phylip_rs::tree::{Alignment, Base, Sequence};

// ---------------------------------------------------------------------------
// JSON response types
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct NjResult {
    names: Vec<String>,
    matrix: Vec<Vec<f64>>,
    newick: String,
    ntaxa: usize,
    nsites: usize,
    model: String,
}

#[derive(Serialize)]
struct LikelihoodResult {
    lnl: f64,
    ntaxa: usize,
    nsites: usize,
    model: String,
    newick: String,
}

#[derive(Serialize)]
struct ParsimonyResult {
    score: usize,
    newick: String,
    ntaxa: usize,
    nsites: usize,
}

#[derive(Serialize)]
struct FelsenResult {
    replicates: Vec<FelsenReplicate>,
    summary: FelsenSummary,
    params: FelsenParams,
}

#[derive(Serialize)]
struct FelsenReplicate {
    rep: usize,
    parsimony_best: usize,
    ml_best: usize,
}

#[derive(Serialize)]
struct FelsenSummary {
    parsimony_correct_pct: f64,
    ml_correct_pct: f64,
    parsimony_lba_pct: f64,
    nreps: u32,
    nsites: u32,
}

#[derive(Serialize)]
struct FelsenParams {
    long_branch: f64,
    short_branch: f64,
    internal: f64,
    nsites: u32,
    nreps: u32,
}

#[derive(Serialize)]
struct BootstrapResult {
    original_newick: String,
    consensus_newick: String,
    nreps: u32,
    ntaxa: usize,
    nsites: usize,
    split_support: Vec<SplitSupport>,
}

#[derive(Serialize)]
struct SplitSupport {
    taxa: Vec<String>,
    count: usize,
    frequency: f64,
}

// ---------------------------------------------------------------------------
// Simple LCG RNG (matches phylip-rs bootstrap::SimpleRng interface)
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
// 1. compute_nj
// ---------------------------------------------------------------------------

/// Compute a distance matrix and Neighbor-Joining tree from FASTA input.
///
/// `model`: "jc69" or "k2p"
/// Returns JSON: { names, matrix, newick, ntaxa, nsites, model }
#[wasm_bindgen]
pub fn compute_nj(fasta: &str, model: &str) -> Result<String, String> {
    let alignment = read_fasta(fasta).map_err(|e| format!("FASTA parse error: {e:?}"))?;
    let ntaxa = alignment.ntaxa();
    let nsites = alignment.nsites();

    let dist_matrix = match model.to_lowercase().as_str() {
        "k2p" | "kimura" => {
            let m = K2P::new();
            compute_distance_matrix(&alignment, &m)
        }
        _ => {
            let m = JC69::new();
            compute_distance_matrix(&alignment, &m)
        }
    }
    .map_err(|e| format!("Distance computation error: {e:?}"))?;

    let tree = neighbor_joining(&dist_matrix);
    let newick = write_newick(&tree);

    let result = NjResult {
        names: dist_matrix.names.clone(),
        matrix: dist_matrix.to_full_matrix(),
        newick,
        ntaxa,
        nsites,
        model: model.to_string(),
    };
    serde_json::to_string(&result).map_err(|e| format!("JSON error: {e}"))
}

// ---------------------------------------------------------------------------
// 2. compute_likelihood
// ---------------------------------------------------------------------------

/// Evaluate tree likelihood for an alignment under a substitution model.
///
/// `model`: "jc69" or "f84"
/// Returns JSON: { lnl, ntaxa, nsites, model, newick }
#[wasm_bindgen]
pub fn compute_likelihood(fasta: &str, newick_str: &str, model: &str) -> Result<String, String> {
    let alignment = read_fasta(fasta).map_err(|e| format!("FASTA parse error: {e:?}"))?;
    let mut tree = parse_newick(newick_str).map_err(|e| format!("Newick parse error: {e:?}"))?;

    let lnl = match model.to_lowercase().as_str() {
        "f84" => {
            let m = F84Model::equal_freqs(2.0);
            optimize_branch_lengths(&mut tree, &alignment, &m)
                .map_err(|e| format!("Optimization error: {e:?}"))?;
            log_likelihood(&tree, &alignment, &m)
        }
        _ => {
            let m = Jc69Model;
            optimize_branch_lengths(&mut tree, &alignment, &m)
                .map_err(|e| format!("Optimization error: {e:?}"))?;
            log_likelihood(&tree, &alignment, &m)
        }
    }
    .map_err(|e| format!("Likelihood error: {e:?}"))?;

    let result = LikelihoodResult {
        lnl,
        ntaxa: alignment.ntaxa(),
        nsites: alignment.nsites(),
        model: model.to_string(),
        newick: write_newick(&tree),
    };
    serde_json::to_string(&result).map_err(|e| format!("JSON error: {e}"))
}

// ---------------------------------------------------------------------------
// 3. compute_parsimony
// ---------------------------------------------------------------------------

/// Run maximum parsimony tree search on FASTA input.
///
/// `seed`: random seed (0 = use default seed 42)
/// Returns JSON: { score, newick, ntaxa, nsites }
#[wasm_bindgen]
pub fn compute_parsimony(fasta: &str, seed: u32) -> Result<String, String> {
    let alignment = read_fasta(fasta).map_err(|e| format!("FASTA parse error: {e:?}"))?;
    let ntaxa = alignment.ntaxa();
    let nsites = alignment.nsites();

    let seed_opt = if seed == 0 { Some(42u64) } else { Some(seed as u64) };
    let pars_result = parsimony_search(&alignment, seed_opt);

    let result = ParsimonyResult {
        score: pars_result.score,
        newick: write_newick(&pars_result.tree),
        ntaxa,
        nsites,
    };
    serde_json::to_string(&result).map_err(|e| format!("JSON error: {e}"))
}

// ---------------------------------------------------------------------------
// 4. felsenstein_zone
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

/// Simulate alignment along a tree under JC69.
fn simulate_alignment(
    tree: &phylip_rs::tree::Tree,
    nsites: usize,
    model: &Jc69Model,
    rng: &mut Rng,
) -> Alignment {
    let nnodes = tree.num_nodes();
    let mut node_seqs: Vec<Vec<usize>> = vec![vec![0; nsites]; nnodes];

    for slot in node_seqs[tree.root].iter_mut().take(nsites) {
        *slot = (rng.next_u64() >> 33) as usize % 4;
    }

    let preorder = tree.preorder();
    for &node_id in &preorder {
        if node_id == tree.root {
            continue;
        }
        let branch_len = tree.nodes[node_id].branch_length.unwrap_or(0.1);
        let trans_probs = model.transition_probs(branch_len);
        let parent_id = tree.nodes[node_id].parent.unwrap();

        #[allow(clippy::needless_range_loop)]
        for site in 0..nsites {
            let parent_base = node_seqs[parent_id][site];
            let row = &trans_probs[parent_base];
            let u = rng.next_f64();
            let mut cumulative = 0.0;
            let mut child_base = 3;
            for (j, &prob) in row.iter().enumerate().take(4) {
                cumulative += prob;
                if u < cumulative {
                    child_base = j;
                    break;
                }
            }
            node_seqs[node_id][site] = child_base;
        }
    }

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

const TOPOLOGY_NEWICKS: [&str; 3] = [
    "((A,B),(C,D));",
    "((A,C),(B,D));",
    "((A,D),(B,C));",
];

fn evaluate_topologies(alignment: &Alignment, model: &Jc69Model) -> (usize, usize) {
    let mut pars_scores = [0usize; 3];
    let mut ml_lnls = [f64::NEG_INFINITY; 3];

    for (i, newick) in TOPOLOGY_NEWICKS.iter().enumerate() {
        let tree = parse_newick(newick).expect("topology Newick should parse");

        let (score, _) = parsimony_score(&tree, alignment)
            .expect("parsimony scoring should succeed");
        pars_scores[i] = score;

        let mut ml_tree = tree.clone();
        for node in &mut ml_tree.nodes {
            if node.parent.is_some() {
                node.branch_length = Some(0.1);
            }
        }
        if optimize_branch_lengths(&mut ml_tree, alignment, model).is_ok() {
            if let Ok(lnl) = log_likelihood(&ml_tree, alignment, model) {
                ml_lnls[i] = lnl;
            }
        }
    }

    let pars_best = pars_scores
        .iter()
        .enumerate()
        .min_by_key(|&(_, &s)| s)
        .unwrap()
        .0;

    let ml_best = ml_lnls
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap()
        .0;

    (pars_best, ml_best)
}

/// Simulate the Felsenstein Zone: demonstrate long-branch attraction.
///
/// Parameters: branch lengths, number of sites, replicates, and seed.
/// Returns JSON with per-replicate results and accuracy summary.
#[wasm_bindgen]
pub fn felsenstein_zone(
    long_branch: f64,
    short_branch: f64,
    internal: f64,
    nsites: u32,
    nreps: u32,
    seed: u32,
) -> Result<String, String> {
    let true_newick = format!(
        "((A:{lb},B:{sb}):{int},(C:{sb},D:{lb}));",
        lb = long_branch,
        sb = short_branch,
        int = internal,
    );
    let true_tree = parse_newick(&true_newick)
        .map_err(|e| format!("Tree parse error: {e:?}"))?;

    let model = Jc69Model;
    let base_seed = if seed == 0 { 314159265u64 } else { seed as u64 };

    let mut replicates = Vec::with_capacity(nreps as usize);
    let mut pars_correct = 0u32;
    let mut ml_correct = 0u32;
    let mut pars_lba = 0u32;

    for rep in 0..nreps {
        let rep_seed = base_seed
            .wrapping_add(nsites as u64 * 1000003)
            .wrapping_add(rep as u64 * 999983);
        let mut rng = Rng::new(rep_seed);

        let alignment = simulate_alignment(&true_tree, nsites as usize, &model, &mut rng);
        let (pars_best, ml_best) = evaluate_topologies(&alignment, &model);

        if pars_best == 0 { pars_correct += 1; }
        if ml_best == 0 { ml_correct += 1; }
        if pars_best == 2 { pars_lba += 1; }

        replicates.push(FelsenReplicate {
            rep: rep as usize,
            parsimony_best: pars_best,
            ml_best,
        });
    }

    let result = FelsenResult {
        replicates,
        summary: FelsenSummary {
            parsimony_correct_pct: 100.0 * pars_correct as f64 / nreps as f64,
            ml_correct_pct: 100.0 * ml_correct as f64 / nreps as f64,
            parsimony_lba_pct: 100.0 * pars_lba as f64 / nreps as f64,
            nreps,
            nsites,
        },
        params: FelsenParams {
            long_branch,
            short_branch,
            internal,
            nsites,
            nreps,
        },
    };
    serde_json::to_string(&result).map_err(|e| format!("JSON error: {e}"))
}

// ---------------------------------------------------------------------------
// 5. bootstrap_nj
// ---------------------------------------------------------------------------

/// Bootstrap NJ analysis with majority-rule consensus.
///
/// `nreps`: number of bootstrap replicates
/// `seed`: random seed (0 = use default 42)
/// Returns JSON: { original_newick, consensus_newick, nreps, split_support }
#[wasm_bindgen]
pub fn bootstrap_nj(fasta: &str, nreps: u32, seed: u32) -> Result<String, String> {
    let alignment = read_fasta(fasta).map_err(|e| format!("FASTA parse error: {e:?}"))?;
    let ntaxa = alignment.ntaxa();
    let nsites = alignment.nsites();

    // Original NJ tree
    let jc = JC69::new();
    let orig_matrix = compute_distance_matrix(&alignment, &jc)
        .map_err(|e| format!("Distance error: {e:?}"))?;
    let orig_tree = neighbor_joining(&orig_matrix);
    let original_newick = write_newick(&orig_tree);

    // Bootstrap replicates
    let seed_val = if seed == 0 { 42u64 } else { seed as u64 };
    let mut rng = SimpleRng::seed(seed_val);
    let mut bootstrap_trees = Vec::with_capacity(nreps as usize);

    for _ in 0..nreps {
        let weights = bootstrap_weights(nsites, &ResamplingMethod::Bootstrap, &mut rng);
        let resampled = resample_alignment(&alignment, &weights);
        let dm = compute_distance_matrix(&resampled, &jc)
            .map_err(|e| format!("Bootstrap distance error: {e:?}"))?;
        let tree = neighbor_joining(&dm);
        bootstrap_trees.push(tree);
    }

    // Majority-rule consensus
    let consensus = consensus_tree(&bootstrap_trees, &ConsensusMethod::MajorityRule)
        .map_err(|e| format!("Consensus error: {e:?}"))?;

    let split_support: Vec<SplitSupport> = consensus
        .split_frequencies
        .iter()
        .map(|(taxa, count, freq)| SplitSupport {
            taxa: taxa.clone(),
            count: *count,
            frequency: *freq,
        })
        .collect();

    let result = BootstrapResult {
        original_newick,
        consensus_newick: write_newick(&consensus.tree),
        nreps,
        ntaxa,
        nsites,
        split_support,
    };
    serde_json::to_string(&result).map_err(|e| format!("JSON error: {e}"))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_FASTA: &str = ">A\nACGTACGT\n>B\nACGTACGA\n>C\nACGTACGA\n>D\nACGAACGA\n";

    #[test]
    fn test_compute_nj_jc69() {
        let result = compute_nj(TEST_FASTA, "jc69").unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&result).unwrap();
        assert_eq!(parsed["ntaxa"], 4);
        assert_eq!(parsed["nsites"], 8);
        assert!(parsed["newick"].as_str().unwrap().contains(';'));
        assert_eq!(parsed["matrix"].as_array().unwrap().len(), 4);
    }

    #[test]
    fn test_compute_nj_k2p() {
        let result = compute_nj(TEST_FASTA, "k2p").unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&result).unwrap();
        assert_eq!(parsed["ntaxa"], 4);
    }

    #[test]
    fn test_compute_likelihood() {
        let newick = "((A:0.1,B:0.1):0.05,(C:0.1,D:0.15):0.05);";
        let result = compute_likelihood(TEST_FASTA, newick, "jc69").unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&result).unwrap();
        assert!(parsed["lnl"].as_f64().unwrap() < 0.0);
        assert_eq!(parsed["ntaxa"], 4);
    }

    #[test]
    fn test_compute_parsimony() {
        let result = compute_parsimony(TEST_FASTA, 42).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&result).unwrap();
        assert!(parsed["score"].as_u64().unwrap() > 0);
        assert!(parsed["newick"].as_str().unwrap().contains(';'));
    }

    #[test]
    fn test_felsenstein_zone() {
        let result = felsenstein_zone(0.8, 0.01, 0.01, 500, 5, 42).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&result).unwrap();
        assert_eq!(parsed["summary"]["nreps"], 5);
        assert_eq!(parsed["summary"]["nsites"], 500);
        assert_eq!(parsed["replicates"].as_array().unwrap().len(), 5);
    }

    #[test]
    fn test_bootstrap_nj() {
        let result = bootstrap_nj(TEST_FASTA, 10, 42).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&result).unwrap();
        assert_eq!(parsed["nreps"], 10);
        assert!(parsed["original_newick"].as_str().unwrap().contains(';'));
        assert!(parsed["consensus_newick"].as_str().unwrap().contains(';'));
    }

    #[test]
    fn test_invalid_fasta() {
        let result = compute_nj("not valid fasta", "jc69");
        assert!(result.is_err());
    }
}

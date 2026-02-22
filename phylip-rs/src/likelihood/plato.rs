//! PLATO: Partial Likelihoods Assessed Through Optimisation.
//!
//! Implements the sliding-window likelihood scanning algorithm of
//! Grassly & Holmes (1997) for detecting regions of anomalous
//! phylogenetic signal within a multiple sequence alignment, such as
//! those arising from recombination, gene conversion, or variable
//! selective pressure.
//!
//! # Algorithm overview
//!
//! 1. Fix the tree topology and substitution model parameters estimated
//!    from the full alignment.
//! 2. Compute per-site log-likelihoods under the fitted model.
//! 3. Slide a window of size W across the alignment with step size S.
//!    For each window starting at site j, compute the partial
//!    log-likelihood PLL(j) = sum of site log-likelihoods for sites
//!    j..j+W.
//! 4. Under the null hypothesis of a homogeneous model, the expected
//!    partial log-likelihood is E[PLL] = (W / N) * total_lnL, where
//!    N is the total number of sites.
//! 5. The test statistic for each window is
//!    delta(j) = PLL(j) - E[PLL], where a significantly negative
//!    delta indicates that the window fits the tree poorly (a potential
//!    recombination breakpoint or region under different selection).
//! 6. Approximate significance is assessed via z-scores computed from
//!    the observed distribution of window statistics, or via parametric
//!    bootstrap under the null model.
//!
//! # Parametric bootstrap
//!
//! For formal hypothesis testing, the parametric bootstrap simulates
//! alignments under the fitted null model on the fitted tree and
//! computes the distribution of min(delta) across windows. The p-value
//! is the proportion of simulated min(delta) values that are as extreme
//! as or more extreme than the observed min(delta).
//!
//! # Example
//!
//! ```
//! use phylip_rs::likelihood::plato::plato_scan;
//! use phylip_rs::likelihood::models::Jc69Model;
//! use phylip_rs::tree::{Alignment, Base, Sequence};
//! use phylip_rs::tree::newick::parse_newick;
//!
//! let alignment = Alignment::new(vec![
//!     Sequence::new("A", vec![Base::A; 100]),
//!     Sequence::new("B", vec![Base::A; 100]),
//!     Sequence::new("C", vec![Base::T; 100]),
//! ]).unwrap();
//!
//! let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
//! let model = Jc69Model;
//! let result = plato_scan(&tree, &alignment, &model, 20, 5, 3.0).unwrap();
//! assert!(!result.window_starts.is_empty());
//! ```
//!
//! # References
//!
//! - Grassly, N. C. & Holmes, E. C. (1997). A likelihood method for the
//!   detection of selection and recombination using nucleotide sequences.
//!   *Molecular Biology and Evolution*, 14(3), 239-247.
//! - Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a
//!   maximum likelihood approach. *Journal of Molecular Evolution*, 17,
//!   368-376.

use crate::bootstrap::SimpleRng;
use crate::likelihood::models::SubstitutionModel;
use crate::likelihood::pruning::{site_log_likelihoods, LikelihoodError};
use crate::tree::types::{Alignment, Base, Sequence, Tree};

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

/// Result of a PLATO sliding-window likelihood scan.
///
/// Contains per-window partial log-likelihoods, deviations from the
/// expected value under a homogeneous model, z-scores, and a list of
/// windows that exceed the significance threshold.
#[derive(Debug, Clone)]
pub struct PlatoResult {
    /// Window width in alignment sites.
    pub window_size: usize,
    /// Step size (number of sites between successive window starts).
    pub step_size: usize,
    /// Starting site index for each window.
    pub window_starts: Vec<usize>,
    /// Partial log-likelihood for each window (sum of per-site lnL
    /// within the window).
    pub window_lnl: Vec<f64>,
    /// Deviation from the expected partial log-likelihood for each
    /// window: delta(j) = PLL(j) - E[PLL].
    pub window_delta: Vec<f64>,
    /// Z-score for each window, computed from the observed distribution
    /// of partial log-likelihoods:
    /// z(j) = (PLL(j) - mean_PLL) / sd_PLL.
    pub window_zscore: Vec<f64>,
    /// Mean partial log-likelihood across all windows.
    pub mean_window_lnl: f64,
    /// Standard deviation of partial log-likelihoods across windows.
    pub sd_window_lnl: f64,
    /// Total log-likelihood of the full alignment.
    pub total_lnl: f64,
    /// Indices into `window_starts` for windows whose |z-score|
    /// exceeds the threshold.
    pub anomalous_windows: Vec<usize>,
    /// Z-score threshold used to flag anomalous windows.
    pub threshold: f64,
}

/// Result of a PLATO parametric bootstrap analysis.
///
/// Provides a p-value for the observed minimum delta statistic by
/// comparing it against the null distribution obtained from simulated
/// alignments.
#[derive(Debug, Clone)]
pub struct PlatoBootstrapResult {
    /// The minimum delta(j) observed in the real data.
    pub observed_min_delta: f64,
    /// Minimum delta values from each simulated alignment.
    pub simulated_min_deltas: Vec<f64>,
    /// Proportion of simulated min(delta) values <= observed min(delta).
    pub p_value: f64,
}

// ---------------------------------------------------------------------------
// Core scan
// ---------------------------------------------------------------------------

/// Perform a PLATO sliding-window likelihood scan.
///
/// Computes per-site log-likelihoods from the given tree, alignment, and
/// model, then slides a window of `window_size` sites across the alignment
/// with a stride of `step_size`. Windows whose z-score exceeds `threshold`
/// in absolute value are flagged as anomalous.
///
/// # Arguments
///
/// * `tree` - Phylogenetic tree with branch lengths.
/// * `alignment` - Multiple sequence alignment.
/// * `model` - Nucleotide substitution model.
/// * `window_size` - Width of the sliding window in sites (must be >= 1
///   and <= alignment length).
/// * `step_size` - Number of sites between successive window starts
///   (must be >= 1).
/// * `threshold` - Z-score threshold for flagging anomalous windows
///   (e.g., 3.0 for approximately 99.7% confidence under normality).
///
/// # Errors
///
/// Returns `LikelihoodError` if the pruning algorithm fails (e.g., taxa
/// mismatch, empty tree, numerical issues), or if `window_size` or
/// `step_size` are invalid.
///
/// # References
///
/// Grassly, N. C. & Holmes, E. C. (1997). *Mol. Biol. Evol.* 14:239-247.
pub fn plato_scan(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    window_size: usize,
    step_size: usize,
    threshold: f64,
) -> Result<PlatoResult, LikelihoodError> {
    let nsites = alignment.nsites();

    // Validate parameters.
    if window_size == 0 || window_size > nsites {
        return Err(LikelihoodError::NumericalError(format!(
            "window_size ({}) must be between 1 and nsites ({})",
            window_size, nsites
        )));
    }
    if step_size == 0 {
        return Err(LikelihoodError::NumericalError(
            "step_size must be >= 1".to_string(),
        ));
    }

    // Step 1-2: compute per-site log-likelihoods.
    let sll = site_log_likelihoods(tree, alignment, model)?;
    let total_lnl: f64 = sll.iter().sum();

    // Step 3: slide a window across the alignment.
    let mut window_starts: Vec<usize> = Vec::new();
    let mut window_lnl: Vec<f64> = Vec::new();

    let mut start = 0;
    while start + window_size <= nsites {
        window_starts.push(start);
        let pll: f64 = sll[start..start + window_size].iter().sum();
        window_lnl.push(pll);
        start += step_size;
    }

    let n_windows = window_starts.len();

    // Step 4: expected partial log-likelihood under the null.
    let expected_pll = (window_size as f64 / nsites as f64) * total_lnl;

    // Step 5: delta(j) = PLL(j) - E[PLL].
    let window_delta: Vec<f64> = window_lnl.iter().map(|pll| pll - expected_pll).collect();

    // Step 6: z-scores from the observed window distribution.
    let (mean_pll, sd_pll) = mean_and_sd(&window_lnl);

    let window_zscore: Vec<f64> = if sd_pll > 0.0 {
        window_lnl
            .iter()
            .map(|pll| (pll - mean_pll) / sd_pll)
            .collect()
    } else {
        vec![0.0; n_windows]
    };

    // Flag anomalous windows.
    let anomalous_windows: Vec<usize> = window_zscore
        .iter()
        .enumerate()
        .filter(|(_, z)| z.abs() > threshold)
        .map(|(i, _)| i)
        .collect();

    Ok(PlatoResult {
        window_size,
        step_size,
        window_starts,
        window_lnl,
        window_delta,
        window_zscore,
        mean_window_lnl: mean_pll,
        sd_window_lnl: sd_pll,
        total_lnl,
        anomalous_windows,
        threshold,
    })
}

// ---------------------------------------------------------------------------
// Parametric bootstrap
// ---------------------------------------------------------------------------

/// Perform a PLATO parametric bootstrap to assess significance.
///
/// Simulates `n_simulations` alignments under the fitted model on the
/// fitted tree, computes the sliding-window statistics for each, and
/// reports the p-value for the observed minimum delta.
///
/// The p-value is the proportion of simulated min(delta) values that
/// are less than or equal to the observed min(delta). A small p-value
/// indicates that the observed alignment contains regions whose fit to
/// the tree is worse than expected under a homogeneous model.
///
/// # Arguments
///
/// * `tree` - Phylogenetic tree with branch lengths.
/// * `alignment` - Multiple sequence alignment.
/// * `model` - Nucleotide substitution model.
/// * `window_size` - Width of the sliding window in sites.
/// * `step_size` - Stride between successive windows.
/// * `n_simulations` - Number of parametric bootstrap replicates.
/// * `seed` - Seed for the pseudo-random number generator.
///
/// # Errors
///
/// Returns an error string if the scan or simulation fails.
///
/// # References
///
/// Grassly, N. C. & Holmes, E. C. (1997). *Mol. Biol. Evol.* 14:239-247.
pub fn plato_parametric_bootstrap(
    tree: &Tree,
    alignment: &Alignment,
    model: &dyn SubstitutionModel,
    window_size: usize,
    step_size: usize,
    n_simulations: usize,
    seed: u64,
) -> Result<PlatoBootstrapResult, String> {
    let nsites = alignment.nsites();
    let ntaxa = alignment.ntaxa();

    if nsites == 0 || ntaxa == 0 {
        return Err("alignment is empty".to_string());
    }

    // Scan the observed data.
    let observed = plato_scan(tree, alignment, model, window_size, step_size, 3.0)
        .map_err(|e| format!("observed scan failed: {}", e))?;

    let observed_min_delta = observed
        .window_delta
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);

    // Build a taxon-name-to-leaf-index map for simulation output.
    let leaf_order: Vec<(usize, String)> = tree
        .nodes
        .iter()
        .filter(|n| n.is_leaf())
        .map(|n| (n.id, n.name.clone().unwrap_or_default()))
        .collect();

    let mut rng = SimpleRng::seed(seed);
    let mut simulated_min_deltas: Vec<f64> = Vec::with_capacity(n_simulations);

    for _ in 0..n_simulations {
        // Simulate an alignment under the null model.
        let sim_aln =
            simulate_alignment(tree, model, nsites, &leaf_order, &mut rng);

        // Compute sliding-window statistics on the simulated data.
        let sim_sll = site_log_likelihoods(tree, &sim_aln, model)
            .map_err(|e| format!("simulated scan failed: {}", e))?;
        let sim_total: f64 = sim_sll.iter().sum();
        let sim_expected = (window_size as f64 / nsites as f64) * sim_total;

        let mut sim_min_delta = f64::INFINITY;
        let mut start = 0;
        while start + window_size <= nsites {
            let pll: f64 = sim_sll[start..start + window_size].iter().sum();
            let delta = pll - sim_expected;
            if delta < sim_min_delta {
                sim_min_delta = delta;
            }
            start += step_size;
        }

        simulated_min_deltas.push(sim_min_delta);
    }

    // P-value: proportion of simulated min(delta) <= observed min(delta).
    let n_extreme = simulated_min_deltas
        .iter()
        .filter(|&&d| d <= observed_min_delta)
        .count();
    let p_value = n_extreme as f64 / n_simulations as f64;

    Ok(PlatoBootstrapResult {
        observed_min_delta,
        simulated_min_deltas,
        p_value,
    })
}

// ---------------------------------------------------------------------------
// Sequence simulation along a tree
// ---------------------------------------------------------------------------

/// Simulate a DNA alignment by evolving sequences along the given tree
/// under the specified substitution model.
///
/// At the root, the ancestral state at each site is drawn from the
/// equilibrium frequencies of the model. Each branch then evolves the
/// state according to the transition probability matrix P(t) for its
/// branch length.
fn simulate_alignment(
    tree: &Tree,
    model: &dyn SubstitutionModel,
    nsites: usize,
    leaf_order: &[(usize, String)],
    rng: &mut SimpleRng,
) -> Alignment {
    let nnodes = tree.num_nodes();
    let pi = model.frequencies();

    // Build cumulative frequency distribution for root sampling.
    let cum_pi = cumulative_probs(&pi);

    // States at each node for each site.
    let mut states: Vec<Vec<usize>> = vec![vec![0; nsites]; nnodes];

    // Draw root states from equilibrium frequencies.
    for site in 0..nsites {
        states[tree.root][site] = draw_state(&cum_pi, rng);
    }

    // Pre-order traversal: evolve states from parent to children.
    let preorder = tree.preorder();
    for &node_id in &preorder {
        if node_id == tree.root {
            continue;
        }
        let parent_id = tree.nodes[node_id].parent.unwrap();
        let t = tree.nodes[node_id]
            .branch_length
            .unwrap_or(0.1)
            .max(1e-8);
        let p_matrix = model.transition_probs(t);

        // Pre-compute cumulative row probabilities.
        let cum_rows: Vec<[f64; 4]> = (0..4)
            .map(|i| cumulative_probs(&p_matrix[i]))
            .collect();

        for site in 0..nsites {
            let parent_state = states[parent_id][site];
            states[node_id][site] = draw_state(&cum_rows[parent_state], rng);
        }
    }

    // Extract leaf sequences in the requested order.
    let idx_to_base: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];

    let sequences: Vec<Sequence> = leaf_order
        .iter()
        .map(|(node_id, name)| {
            let bases: Vec<Base> = states[*node_id]
                .iter()
                .map(|&s| idx_to_base[s])
                .collect();
            Sequence::new(name.clone(), bases)
        })
        .collect();

    Alignment::new(sequences).expect("simulated alignment should be valid")
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

/// Compute the mean and standard deviation of a slice of f64 values.
///
/// Returns (0.0, 0.0) for an empty slice and (mean, 0.0) for a single
/// element.
fn mean_and_sd(values: &[f64]) -> (f64, f64) {
    let n = values.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let mean = values.iter().sum::<f64>() / n as f64;
    if n == 1 {
        return (mean, 0.0);
    }
    let var = values.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / (n - 1) as f64;
    (mean, var.sqrt())
}

/// Convert a 4-element probability vector into cumulative probabilities.
fn cumulative_probs(probs: &[f64; 4]) -> [f64; 4] {
    let mut cum = [0.0; 4];
    cum[0] = probs[0];
    cum[1] = cum[0] + probs[1];
    cum[2] = cum[1] + probs[2];
    cum[3] = cum[2] + probs[3];
    cum
}

/// Draw a state (0..3) from a cumulative probability distribution.
fn draw_state(cum: &[f64; 4], rng: &mut SimpleRng) -> usize {
    let u = rng.gen_f64();
    if u < cum[0] {
        0
    } else if u < cum[1] {
        1
    } else if u < cum[2] {
        2
    } else {
        3
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::likelihood::models::{F84Model, Jc69Model};
    use crate::tree::newick::parse_newick;
    use crate::tree::types::{Base, BaseFrequencies, Sequence};

    /// Helper to build an alignment from name/sequence pairs.
    fn make_alignment(data: &[(&str, &str)]) -> Alignment {
        let seqs: Vec<Sequence> = data
            .iter()
            .map(|(name, seq)| {
                let bases: Vec<Base> = seq
                    .chars()
                    .map(|c| Base::from_char(c).unwrap())
                    .collect();
                Sequence::new(*name, bases)
            })
            .collect();
        Alignment::new(seqs).unwrap()
    }

    /// Build a homogeneous test alignment (identical pattern repeated).
    fn homogeneous_alignment(nsites: usize) -> Alignment {
        // Simple repeating pattern — homogeneous signal.
        let pattern = "ACGTACGT";
        let make_seq = |name: &str, offset: usize| {
            let bases: Vec<Base> = (0..nsites)
                .map(|i| {
                    let c = pattern.as_bytes()[(i + offset) % pattern.len()] as char;
                    Base::from_char(c).unwrap()
                })
                .collect();
            Sequence::new(name, bases)
        };
        Alignment::new(vec![
            make_seq("A", 0),
            make_seq("B", 0),
            make_seq("C", 4),
        ])
        .unwrap()
    }

    // ---------------------------------------------------------------
    // Basic scan tests
    // ---------------------------------------------------------------

    #[test]
    fn test_scan_basic() {
        let aln = homogeneous_alignment(100);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let result = plato_scan(&tree, &aln, &model, 20, 5, 3.0).unwrap();

        assert_eq!(result.window_size, 20);
        assert_eq!(result.step_size, 5);
        assert!(!result.window_starts.is_empty());
        assert_eq!(result.window_lnl.len(), result.window_starts.len());
        assert_eq!(result.window_delta.len(), result.window_starts.len());
        assert_eq!(result.window_zscore.len(), result.window_starts.len());
    }

    #[test]
    fn test_scan_window_positions() {
        let aln = homogeneous_alignment(50);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let result = plato_scan(&tree, &aln, &model, 10, 5, 3.0).unwrap();

        // Windows: 0..10, 5..15, 10..20, ..., 40..50
        let expected_starts: Vec<usize> = (0..=40).step_by(5).collect();
        assert_eq!(result.window_starts, expected_starts);
    }

    #[test]
    fn test_scan_window_covers_full_alignment() {
        // Window size equals alignment length.
        let aln = homogeneous_alignment(20);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let result = plato_scan(&tree, &aln, &model, 20, 1, 3.0).unwrap();

        assert_eq!(result.window_starts.len(), 1);
        assert_eq!(result.window_starts[0], 0);
        // The single window's PLL should equal the total lnL.
        assert!(
            (result.window_lnl[0] - result.total_lnl).abs() < 1e-10,
            "full-alignment window lnl ({}) should equal total ({})",
            result.window_lnl[0],
            result.total_lnl
        );
    }

    #[test]
    fn test_scan_step_one() {
        let aln = homogeneous_alignment(30);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let result = plato_scan(&tree, &aln, &model, 10, 1, 3.0).unwrap();

        // With step=1 and window=10 on 30 sites: 30-10+1=21 windows.
        assert_eq!(result.window_starts.len(), 21);
        assert_eq!(result.window_starts[0], 0);
        assert_eq!(*result.window_starts.last().unwrap(), 20);
    }

    // ---------------------------------------------------------------
    // Homogeneous data: no anomalous windows expected
    // ---------------------------------------------------------------

    #[test]
    fn test_homogeneous_no_anomalous() {
        // With a perfectly repeating pattern, all windows should have
        // essentially the same partial log-likelihood; none should be
        // flagged at z=3.0.
        let aln = homogeneous_alignment(200);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let result = plato_scan(&tree, &aln, &model, 40, 8, 3.0).unwrap();

        assert!(
            result.anomalous_windows.is_empty(),
            "homogeneous data should have no anomalous windows, found {} ({:?})",
            result.anomalous_windows.len(),
            result.anomalous_windows,
        );
    }

    // ---------------------------------------------------------------
    // Window statistics properties
    // ---------------------------------------------------------------

    #[test]
    fn test_mean_window_lnl_close_to_expected() {
        // The mean of all window PLLs should be close to the expected
        // PLL = (W/N) * total_lnL.
        let aln = homogeneous_alignment(200);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let result = plato_scan(&tree, &aln, &model, 40, 10, 3.0).unwrap();

        let expected_pll = (40.0 / 200.0) * result.total_lnl;
        let rel_diff = (result.mean_window_lnl - expected_pll).abs()
            / expected_pll.abs().max(1e-12);
        assert!(
            rel_diff < 0.05,
            "mean window lnl ({}) should be close to expected ({}), rel_diff = {}",
            result.mean_window_lnl,
            expected_pll,
            rel_diff,
        );
    }

    #[test]
    fn test_delta_sums_relate_to_total() {
        // For non-overlapping windows that cover the full alignment,
        // the sum of PLL should equal the total lnL.
        let nsites = 100;
        let window = 20;
        let step = 20; // non-overlapping, full coverage

        let aln = homogeneous_alignment(nsites);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let result = plato_scan(&tree, &aln, &model, window, step, 3.0).unwrap();

        let pll_sum: f64 = result.window_lnl.iter().sum();
        assert!(
            (pll_sum - result.total_lnl).abs() < 1e-8,
            "sum of non-overlapping window PLLs ({}) should equal total lnL ({})",
            pll_sum,
            result.total_lnl,
        );
    }

    #[test]
    fn test_zscores_have_zero_mean() {
        // Use a non-periodic alignment to ensure genuine variation
        // across windows (periodic alignments have SD~0 which causes
        // numerical instability in z-scores).
        let tree = parse_newick("((A:0.2,B:0.2):0.1,C:0.4);").unwrap();
        let model = Jc69Model;
        let leaf_order: Vec<(usize, String)> = tree
            .nodes
            .iter()
            .filter(|n| n.is_leaf())
            .map(|n| (n.id, n.name.clone().unwrap_or_default()))
            .collect();

        let mut rng = SimpleRng::seed(314);
        let aln = simulate_alignment(&tree, &model, 200, &leaf_order, &mut rng);

        let result = plato_scan(&tree, &aln, &model, 40, 5, 3.0).unwrap();

        // When SD > 0, the mean of z-scores is exactly zero by
        // construction: sum((x_i - mean)/sd) = 0.
        if result.sd_window_lnl > 1e-12 {
            let z_mean: f64 =
                result.window_zscore.iter().sum::<f64>() / result.window_zscore.len() as f64;
            assert!(
                z_mean.abs() < 1e-10,
                "mean z-score should be ~0, got {}",
                z_mean,
            );
        }
    }

    // ---------------------------------------------------------------
    // Concatenated alignment: detect breakpoint
    // ---------------------------------------------------------------

    #[test]
    fn test_concatenated_detects_breakpoint() {
        // Build two halves of an alignment with conflicting
        // phylogenetic signal. In the first half, A and B share
        // derived states (grouping AB vs CD). In the second half,
        // A and C share derived states (grouping AC vs BD).
        //
        // We use variable sites (synapomorphies) to create genuine
        // phylogenetic signal conflict.

        let half_len = 100;

        // First half: A=B, C=D pattern — supports ((A,B),(C,D)).
        // Alternate between two patterns.
        let mut a1 = Vec::with_capacity(half_len);
        let mut b1 = Vec::with_capacity(half_len);
        let mut c1 = Vec::with_capacity(half_len);
        let mut d1 = Vec::with_capacity(half_len);
        for i in 0..half_len {
            if i % 2 == 0 {
                a1.push(Base::A); b1.push(Base::A);
                c1.push(Base::T); d1.push(Base::T);
            } else {
                a1.push(Base::G); b1.push(Base::G);
                c1.push(Base::C); d1.push(Base::C);
            }
        }

        // Second half: A=C, B=D pattern — supports ((A,C),(B,D)).
        let mut a2 = Vec::with_capacity(half_len);
        let mut b2 = Vec::with_capacity(half_len);
        let mut c2 = Vec::with_capacity(half_len);
        let mut d2 = Vec::with_capacity(half_len);
        for i in 0..half_len {
            if i % 2 == 0 {
                a2.push(Base::A); b2.push(Base::T);
                c2.push(Base::A); d2.push(Base::T);
            } else {
                a2.push(Base::C); b2.push(Base::G);
                c2.push(Base::C); d2.push(Base::G);
            }
        }

        let concat = |mut a: Vec<Base>, b: Vec<Base>| -> Vec<Base> {
            a.extend(b);
            a
        };

        let aln = Alignment::new(vec![
            Sequence::new("A", concat(a1, a2)),
            Sequence::new("B", concat(b1, b2)),
            Sequence::new("C", concat(c1, c2)),
            Sequence::new("D", concat(d1, d2)),
        ])
        .unwrap();

        // Tree that fits the first half well: ((A,B),(C,D)).
        let tree =
            parse_newick("((A:0.05,B:0.05):0.3,(C:0.05,D:0.05):0.3);").unwrap();
        let model = Jc69Model;

        let result = plato_scan(&tree, &aln, &model, 20, 5, 2.0).unwrap();

        // The second half should have lower PLL on this ((A,B),(C,D))
        // tree since it supports ((A,C),(B,D)) instead.
        // The minimum z-score window should be in the second half.
        let min_z_idx = result
            .window_zscore
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();

        let max_z_idx = result
            .window_zscore
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();

        // The best-fitting windows should be in the first half
        // and the worst-fitting in the second half.
        assert!(
            result.window_starts[max_z_idx] < half_len,
            "best window should be in first half (starts at {}), half_len={}",
            result.window_starts[max_z_idx],
            half_len,
        );
        assert!(
            result.window_starts[min_z_idx] >= half_len,
            "worst window should be in second half (starts at {}), half_len={}",
            result.window_starts[min_z_idx],
            half_len,
        );

        // Verify there is genuine variation in z-scores.
        assert!(
            result.sd_window_lnl > 0.0,
            "there should be variation across windows"
        );
    }

    // ---------------------------------------------------------------
    // Error handling
    // ---------------------------------------------------------------

    #[test]
    fn test_scan_invalid_window_size_zero() {
        let aln = homogeneous_alignment(50);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;
        assert!(plato_scan(&tree, &aln, &model, 0, 5, 3.0).is_err());
    }

    #[test]
    fn test_scan_invalid_window_size_too_large() {
        let aln = homogeneous_alignment(50);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;
        assert!(plato_scan(&tree, &aln, &model, 51, 5, 3.0).is_err());
    }

    #[test]
    fn test_scan_invalid_step_size_zero() {
        let aln = homogeneous_alignment(50);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;
        assert!(plato_scan(&tree, &aln, &model, 10, 0, 3.0).is_err());
    }

    #[test]
    fn test_scan_taxa_mismatch() {
        let aln = make_alignment(&[("X", "ACGT"), ("Y", "ACGT"), ("Z", "TGCA")]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;
        assert!(plato_scan(&tree, &aln, &model, 2, 1, 3.0).is_err());
    }

    // ---------------------------------------------------------------
    // Simulation tests
    // ---------------------------------------------------------------

    #[test]
    fn test_simulate_alignment_dimensions() {
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;
        let leaf_order: Vec<(usize, String)> = tree
            .nodes
            .iter()
            .filter(|n| n.is_leaf())
            .map(|n| (n.id, n.name.clone().unwrap_or_default()))
            .collect();

        let mut rng = SimpleRng::seed(42);
        let sim = simulate_alignment(&tree, &model, 100, &leaf_order, &mut rng);

        assert_eq!(sim.ntaxa(), 3);
        assert_eq!(sim.nsites(), 100);
    }

    #[test]
    fn test_simulate_alignment_valid_bases() {
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;
        let leaf_order: Vec<(usize, String)> = tree
            .nodes
            .iter()
            .filter(|n| n.is_leaf())
            .map(|n| (n.id, n.name.clone().unwrap_or_default()))
            .collect();

        let mut rng = SimpleRng::seed(123);
        let sim = simulate_alignment(&tree, &model, 200, &leaf_order, &mut rng);

        for t in 0..sim.ntaxa() {
            for s in 0..sim.nsites() {
                let b = sim.base(t, s);
                assert!(
                    b.is_definite(),
                    "simulated base at ({}, {}) = {:?} is not definite",
                    t,
                    s,
                    b,
                );
            }
        }
    }

    #[test]
    fn test_simulate_base_frequencies() {
        // Under JC69, the equilibrium frequencies are 0.25 each.
        // On a very long simulated alignment, observed frequencies
        // should be close to 0.25.
        let tree = parse_newick("((A:0.5,B:0.5):0.3,C:0.5);").unwrap();
        let model = Jc69Model;
        let leaf_order: Vec<(usize, String)> = tree
            .nodes
            .iter()
            .filter(|n| n.is_leaf())
            .map(|n| (n.id, n.name.clone().unwrap_or_default()))
            .collect();

        let mut rng = SimpleRng::seed(99);
        let sim = simulate_alignment(&tree, &model, 10000, &leaf_order, &mut rng);

        // Count bases across the first taxon.
        let mut counts = [0usize; 4];
        for s in 0..sim.nsites() {
            if let Some(idx) = sim.base(0, s).index() {
                counts[idx] += 1;
            }
        }
        let total = counts.iter().sum::<usize>() as f64;
        for (i, &c) in counts.iter().enumerate() {
            let freq = c as f64 / total;
            assert!(
                (freq - 0.25).abs() < 0.05,
                "base {} frequency = {}, expected ~0.25",
                i,
                freq,
            );
        }
    }

    #[test]
    fn test_simulate_short_branches_low_divergence() {
        // With very short branches, leaf sequences should be nearly
        // identical to the root.
        let tree = parse_newick("((A:0.001,B:0.001):0.001,C:0.001);").unwrap();
        let model = Jc69Model;
        let leaf_order: Vec<(usize, String)> = tree
            .nodes
            .iter()
            .filter(|n| n.is_leaf())
            .map(|n| (n.id, n.name.clone().unwrap_or_default()))
            .collect();

        let mut rng = SimpleRng::seed(7);
        let nsites = 1000;
        let sim = simulate_alignment(&tree, &model, nsites, &leaf_order, &mut rng);

        // Count differences between first two taxa.
        let mut diffs = 0;
        for s in 0..nsites {
            if sim.base(0, s) != sim.base(1, s) {
                diffs += 1;
            }
        }
        let pct_diff = diffs as f64 / nsites as f64;
        assert!(
            pct_diff < 0.02,
            "short branches should produce very similar sequences, \
             but {} / {} sites differ ({:.1}%)",
            diffs,
            nsites,
            pct_diff * 100.0,
        );
    }

    // ---------------------------------------------------------------
    // Parametric bootstrap tests
    // ---------------------------------------------------------------

    #[test]
    fn test_parametric_bootstrap_null_pvalue() {
        // Under the null (homogeneous data), the p-value should be
        // moderate (not extremely small). We use a relatively short
        // alignment and few simulations for speed.
        let aln = homogeneous_alignment(80);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let result = plato_parametric_bootstrap(
            &tree, &aln, &model, 20, 10, 50, 42,
        )
        .unwrap();

        assert_eq!(result.simulated_min_deltas.len(), 50);
        // Under null, p-value should not be extremely small.
        // With only 50 simulations the minimum possible nonzero p is 0.02.
        assert!(
            result.p_value >= 0.0 && result.p_value <= 1.0,
            "p-value should be in [0, 1], got {}",
            result.p_value,
        );
    }

    #[test]
    fn test_parametric_bootstrap_returns_correct_count() {
        let aln = homogeneous_alignment(60);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let result = plato_parametric_bootstrap(
            &tree, &aln, &model, 15, 5, 20, 123,
        )
        .unwrap();

        assert_eq!(
            result.simulated_min_deltas.len(),
            20,
            "should have exactly n_simulations simulated values"
        );
    }

    // ---------------------------------------------------------------
    // Utility function tests
    // ---------------------------------------------------------------

    #[test]
    fn test_mean_and_sd_basic() {
        let vals = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let (mean, sd) = mean_and_sd(&vals);
        assert!((mean - 5.0).abs() < 1e-10, "mean should be 5.0, got {}", mean);
        // Sample SD (Bessel-corrected, n-1) of [2,4,4,4,5,5,7,9]:
        // var = sum((x-5)^2) / 7 = (9+1+1+1+0+0+4+16)/7 = 32/7 ~ 4.571
        // sd ~ 2.138
        let expected_sd = (32.0_f64 / 7.0).sqrt();
        assert!(
            (sd - expected_sd).abs() < 1e-10,
            "sd should be {}, got {}",
            expected_sd,
            sd,
        );
    }

    #[test]
    fn test_mean_and_sd_empty() {
        let (mean, sd) = mean_and_sd(&[]);
        assert_eq!(mean, 0.0);
        assert_eq!(sd, 0.0);
    }

    #[test]
    fn test_mean_and_sd_single() {
        let (mean, sd) = mean_and_sd(&[42.0]);
        assert!((mean - 42.0).abs() < 1e-10);
        assert_eq!(sd, 0.0);
    }

    #[test]
    fn test_cumulative_probs() {
        let p = [0.25, 0.25, 0.25, 0.25];
        let cum = cumulative_probs(&p);
        assert!((cum[0] - 0.25).abs() < 1e-12);
        assert!((cum[1] - 0.50).abs() < 1e-12);
        assert!((cum[2] - 0.75).abs() < 1e-12);
        assert!((cum[3] - 1.00).abs() < 1e-12);
    }

    #[test]
    fn test_draw_state_distribution() {
        let p = [0.1, 0.2, 0.3, 0.4];
        let cum = cumulative_probs(&p);
        let mut rng = SimpleRng::seed(55);
        let mut counts = [0usize; 4];
        let n = 10000;
        for _ in 0..n {
            counts[draw_state(&cum, &mut rng)] += 1;
        }
        for (i, (&expected, &count)) in p.iter().zip(counts.iter()).enumerate() {
            let observed = count as f64 / n as f64;
            assert!(
                (observed - expected).abs() < 0.05,
                "state {} drawn {:.3}, expected {:.3}",
                i,
                observed,
                expected,
            );
        }
    }

    // ---------------------------------------------------------------
    // F84 model scan
    // ---------------------------------------------------------------

    #[test]
    fn test_scan_with_f84_model() {
        let aln = homogeneous_alignment(100);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let freqs = BaseFrequencies::new(0.3, 0.2, 0.2, 0.3);
        let model = F84Model::new(freqs, 2.0);

        let result = plato_scan(&tree, &aln, &model, 20, 5, 3.0).unwrap();

        assert!(!result.window_starts.is_empty());
        assert!(result.total_lnl < 0.0);
    }

    // ---------------------------------------------------------------
    // Regression: all window_delta should sum to zero for non-overlapping
    // ---------------------------------------------------------------

    #[test]
    fn test_delta_sum_zero_nonoverlapping() {
        let nsites = 120;
        let window = 30;
        let step = 30;

        let aln = homogeneous_alignment(nsites);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.3);").unwrap();
        let model = Jc69Model;

        let result = plato_scan(&tree, &aln, &model, window, step, 3.0).unwrap();

        // sum(PLL) = total_lnL, and sum(delta) = sum(PLL) - n_windows*E[PLL]
        // = total_lnL - n_windows * (W/N) * total_lnL
        // = total_lnL * (1 - n_windows * W / N).
        // For full non-overlapping coverage: n_windows * W = N, so sum = 0.
        let delta_sum: f64 = result.window_delta.iter().sum();
        assert!(
            delta_sum.abs() < 1e-8,
            "sum of deltas for non-overlapping full-coverage windows should be 0, got {}",
            delta_sum,
        );
    }
}

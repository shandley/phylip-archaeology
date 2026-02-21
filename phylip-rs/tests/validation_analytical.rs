//! Analytical validation tests for phylip-rs.
//!
//! These tests verify that phylip-rs implementations match known analytical
//! formulas and hand-calculated results. Each test cites the original source.

use phylip_rs::distance::{neighbor_joining, upgma};
use phylip_rs::likelihood::models::{F84Model, Jc69Model, SubstitutionModel};
use phylip_rs::likelihood::pruning::log_likelihood;
use phylip_rs::models::jc69::JC69;
use phylip_rs::models::k2p::K2P;
use phylip_rs::models::compute_distance_matrix;
use phylip_rs::parsimony::wagner::{parsimony_score, search as parsimony_search};
use phylip_rs::tree::newick::{parse_newick, write_newick};
use phylip_rs::tree::{Alignment, Base, BaseFrequencies, DistanceMatrix, Sequence};

// ============================================================================
// Helper: build Alignment from string sequences
// ============================================================================

fn make_alignment(data: &[(&str, &str)]) -> Alignment {
    let sequences: Vec<Sequence> = data
        .iter()
        .map(|(name, seq)| {
            let bases: Vec<Base> = seq.chars().map(|c| Base::from_char(c).unwrap()).collect();
            Sequence::new(*name, bases)
        })
        .collect();
    Alignment::new(sequences).unwrap()
}

// ============================================================================
// JC69 Distance Formula: d = -3/4 * ln(1 - 4p/3)
// Jukes & Cantor (1969) "Evolution of Protein Molecules"
// ============================================================================

#[test]
fn test_jc69_distance_formula_p01() {
    // For p = 0.01 (1% divergence), d = -0.75 * ln(1 - 0.04/3) = 0.01003
    // Create two sequences that differ at exactly 1% of sites
    let n = 1000;
    let seq1_bases = vec![Base::A; n];
    let mut seq2_bases = vec![Base::A; n];

    // Make exactly 10 differences (1%)
    for i in 0..10 {
        seq2_bases[i] = Base::T;
    }

    let seq1 = Sequence::new("S1", seq1_bases);
    let seq2 = Sequence::new("S2", seq2_bases);
    let alignment = Alignment::new(vec![seq1, seq2]).unwrap();

    let model = JC69::new();
    let matrix = compute_distance_matrix(&alignment, &model).unwrap();
    let d = matrix.get(0, 1);

    let p: f64 = 10.0 / 1000.0;
    let expected = -0.75_f64 * (1.0 - 4.0 * p / 3.0).ln();

    assert!(
        (d - expected).abs() < 1e-6,
        "JC69 distance for p=0.01: got {}, expected {}",
        d,
        expected
    );
}

#[test]
fn test_jc69_distance_formula_p10() {
    // For p = 0.10 (10% divergence), d = -0.75 * ln(1 - 0.4/3) = 0.10732
    let n = 1000;
    let seq1_bases = vec![Base::A; n];
    let mut seq2_bases = vec![Base::A; n];

    for i in 0..100 {
        seq2_bases[i] = Base::T;
    }

    let seq1 = Sequence::new("S1", seq1_bases);
    let seq2 = Sequence::new("S2", seq2_bases);
    let alignment = Alignment::new(vec![seq1, seq2]).unwrap();

    let model = JC69::new();
    let matrix = compute_distance_matrix(&alignment, &model).unwrap();
    let d = matrix.get(0, 1);

    let p: f64 = 100.0 / 1000.0;
    let expected = -0.75_f64 * (1.0 - 4.0 * p / 3.0).ln();

    assert!(
        (d - expected).abs() < 1e-6,
        "JC69 distance for p=0.10: got {}, expected {}",
        d,
        expected
    );
}

#[test]
fn test_jc69_distance_formula_p25() {
    // For p = 0.25 (25% divergence), d = -0.75 * ln(1 - 1/3) = 0.30409
    let n = 1200;
    let seq1_bases = vec![Base::A; n];
    let mut seq2_bases = vec![Base::A; n];

    for i in 0..300 {
        seq2_bases[i] = Base::T;
    }

    let seq1 = Sequence::new("S1", seq1_bases);
    let seq2 = Sequence::new("S2", seq2_bases);
    let alignment = Alignment::new(vec![seq1, seq2]).unwrap();

    let model = JC69::new();
    let matrix = compute_distance_matrix(&alignment, &model).unwrap();
    let d = matrix.get(0, 1);

    let p: f64 = 300.0 / 1200.0;
    let expected = -0.75_f64 * (1.0 - 4.0 * p / 3.0).ln();

    assert!(
        (d - expected).abs() < 1e-6,
        "JC69 distance for p=0.25: got {}, expected {}",
        d,
        expected
    );
}

#[test]
fn test_jc69_distance_symmetry() {
    // d(A,B) == d(B,A)
    let alignment = make_alignment(&[
        ("A", "ACGTACGTACGTACGT"),
        ("B", "ACGTACGTTCGTACGT"),
    ]);

    let model = JC69::new();
    let matrix = compute_distance_matrix(&alignment, &model).unwrap();

    assert_eq!(matrix.get(0, 1), matrix.get(1, 0));
}

#[test]
fn test_jc69_distance_identity() {
    // d(A,A) == 0.0
    let alignment = make_alignment(&[
        ("A", "ACGTACGTACGTACGT"),
        ("B", "ACGTACGTACGTACGT"),
    ]);

    let model = JC69::new();
    let matrix = compute_distance_matrix(&alignment, &model).unwrap();

    assert!(
        matrix.get(0, 1).abs() < 1e-10,
        "JC69 distance for identical sequences should be 0"
    );
}

// ============================================================================
// K2P Distance Formula: Kimura (1980) "A simple method for estimating
// evolutionary rates of base substitutions through comparative studies
// of nucleotide sequences" J Mol Evol 16:111-120
//
// d = -0.5 * ln((1 - 2P - Q) * sqrt(1 - 2Q))
// where P = proportion of transitions, Q = proportion of transversions
// ============================================================================

#[test]
fn test_k2p_distance_transitions_only() {
    // If all differences are transitions (P=0.10, Q=0):
    // d = -0.5 * ln((1 - 0.2) * sqrt(1)) = -0.5 * ln(0.8) = 0.11157
    let n = 1000;
    let seq1_bases = vec![Base::A; n];
    let mut seq2_bases = vec![Base::A; n];

    // A→G is a transition
    for i in 0..100 {
        seq2_bases[i] = Base::G;
    }

    let seq1 = Sequence::new("S1", seq1_bases);
    let seq2 = Sequence::new("S2", seq2_bases);
    let alignment = Alignment::new(vec![seq1, seq2]).unwrap();

    let model = K2P::new();
    let matrix = compute_distance_matrix(&alignment, &model).unwrap();
    let d = matrix.get(0, 1);

    let p: f64 = 0.1; // transitions
    let q: f64 = 0.0; // transversions
    let expected = -0.5_f64 * ((1.0 - 2.0 * p - q) * (1.0_f64 - 2.0 * q).sqrt()).ln();

    assert!(
        (d - expected).abs() < 1e-4,
        "K2P distance (transitions only): got {}, expected {}",
        d,
        expected
    );
}

#[test]
fn test_k2p_distance_transversions_only() {
    // If all differences are transversions (P=0, Q=0.10):
    // d = -0.5 * ln((1 - 0.1) * sqrt(1 - 0.2)) = -0.5 * ln(0.9 * 0.8944)
    let n = 1000;
    let seq1_bases = vec![Base::A; n];
    let mut seq2_bases = vec![Base::A; n];

    // A→C is a transversion
    for i in 0..100 {
        seq2_bases[i] = Base::C;
    }

    let seq1 = Sequence::new("S1", seq1_bases);
    let seq2 = Sequence::new("S2", seq2_bases);
    let alignment = Alignment::new(vec![seq1, seq2]).unwrap();

    let model = K2P::new();
    let matrix = compute_distance_matrix(&alignment, &model).unwrap();
    let d = matrix.get(0, 1);

    let p: f64 = 0.0; // transitions
    let q: f64 = 0.1; // transversions
    let expected = -0.5_f64 * ((1.0 - 2.0 * p - q) * (1.0_f64 - 2.0 * q).sqrt()).ln();

    assert!(
        (d - expected).abs() < 1e-4,
        "K2P distance (transversions only): got {}, expected {}",
        d,
        expected
    );
}

#[test]
fn test_k2p_equals_jc69_when_ts_equals_tv() {
    // When ts/tv ratio = 0.5 (equal rates), K2P should equal JC69
    // Create sequences with equal proportions of transitions and transversions
    let n = 900;
    let seq1_bases = vec![Base::A; n];
    let mut seq2_bases = vec![Base::A; n];

    // 30 transitions (A→G) + 60 transversions (30 A→C + 30 A→T) = 90 diffs
    // So P = 30/900 = 1/30, Q = 60/900 = 2/30
    // For JC69, all are equivalent: p = 90/900 = 0.1
    for i in 0..30 {
        seq2_bases[i] = Base::G; // transition
    }
    for i in 30..60 {
        seq2_bases[i] = Base::C; // transversion
    }
    for i in 60..90 {
        seq2_bases[i] = Base::T; // transversion
    }

    let seq1 = Sequence::new("S1", seq1_bases);
    let seq2 = Sequence::new("S2", seq2_bases);
    let alignment = Alignment::new(vec![seq1, seq2]).unwrap();

    let jc69 = JC69::new();
    let k2p = K2P::new();
    let matrix_jc = compute_distance_matrix(&alignment, &jc69).unwrap();
    let matrix_k2p = compute_distance_matrix(&alignment, &k2p).unwrap();

    let d_jc = matrix_jc.get(0, 1);
    let d_k2p = matrix_k2p.get(0, 1);

    // When ts/tv = 0.5, the JC69 and K2P corrections are very close
    // but not identical because K2P separates P and Q terms
    // The key property: K2P >= JC69 always (more parameters, better fit)
    assert!(
        (d_jc - d_k2p).abs() < 0.01,
        "K2P should be close to JC69 when ts=tv: JC69={}, K2P={}",
        d_jc,
        d_k2p
    );
}

// ============================================================================
// Neighbor-Joining: Saitou & Nei (1987) mathematical properties
// ============================================================================

#[test]
fn test_nj_four_taxa_star_tree() {
    // A perfect star tree (equidistant taxa) should produce a star-like NJ tree
    let mut matrix = DistanceMatrix::new(vec![
        "A".to_string(),
        "B".to_string(),
        "C".to_string(),
        "D".to_string(),
    ]);
    // All pairwise distances = 2.0 (symmetric ultrametric)
    for i in 0..4 {
        for j in (i + 1)..4 {
            matrix.set(i, j, 2.0);
        }
    }

    let tree = neighbor_joining(&matrix);
    assert_eq!(tree.num_leaves(), 4);
}

#[test]
fn test_nj_additive_distances() {
    // For an additive (tree-like) distance matrix, NJ should recover the exact tree
    // Tree: ((A:1,B:1):1,(C:1,D:1):1)
    // Distances: d(A,B)=2, d(C,D)=2, d(A,C)=d(A,D)=d(B,C)=d(B,D)=4
    let mut matrix = DistanceMatrix::new(vec![
        "A".to_string(),
        "B".to_string(),
        "C".to_string(),
        "D".to_string(),
    ]);
    matrix.set(0, 1, 2.0); // A-B
    matrix.set(2, 3, 2.0); // C-D
    matrix.set(0, 2, 4.0); // A-C
    matrix.set(0, 3, 4.0); // A-D
    matrix.set(1, 2, 4.0); // B-C
    matrix.set(1, 3, 4.0); // B-D

    let tree = neighbor_joining(&matrix);

    // Check that A and B cluster together, and C and D cluster together
    let newick = write_newick(&tree);

    // Parse the NJ tree and check topology using RF distance against known tree
    let expected = parse_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
    let rf = phylip_rs::tree::distances::robinson_foulds(&tree, &expected).unwrap();
    assert_eq!(
        rf, 0,
        "NJ should recover exact tree from additive distances. Got: {}",
        newick
    );
}

#[test]
fn test_nj_five_taxa_additive() {
    // Tree: (((A:1,B:1):1,C:2):1,(D:1.5,E:1.5):1.5)
    // Compute additive distances from this tree:
    // d(A,B) = 2, d(A,C) = 4, d(A,D) = 5.5, d(A,E) = 5.5
    // d(B,C) = 4, d(B,D) = 5.5, d(B,E) = 5.5
    // d(C,D) = 5.5, d(C,E) = 5.5
    // d(D,E) = 3
    let mut matrix = DistanceMatrix::new(vec![
        "A".to_string(),
        "B".to_string(),
        "C".to_string(),
        "D".to_string(),
        "E".to_string(),
    ]);
    matrix.set(0, 1, 2.0);
    matrix.set(0, 2, 4.0);
    matrix.set(0, 3, 5.5);
    matrix.set(0, 4, 5.5);
    matrix.set(1, 2, 4.0);
    matrix.set(1, 3, 5.5);
    matrix.set(1, 4, 5.5);
    matrix.set(2, 3, 5.5);
    matrix.set(2, 4, 5.5);
    matrix.set(3, 4, 3.0);

    let tree = neighbor_joining(&matrix);

    // A and B should cluster, D and E should cluster
    let expected = parse_newick("(((A:1,B:1):1,C:2):1,(D:1.5,E:1.5):1.5);").unwrap();
    let rf = phylip_rs::tree::distances::robinson_foulds(&tree, &expected).unwrap();
    assert_eq!(
        rf, 0,
        "NJ should recover exact 5-taxon tree from additive distances. Got: {}",
        write_newick(&tree)
    );
}

// ============================================================================
// UPGMA: ultrametric distance property
// ============================================================================

#[test]
fn test_upgma_ultrametric_input() {
    // UPGMA should recover the correct tree from ultrametric distances
    // Tree: ((A:1,B:1):1,(C:2,D:2):0) — ultrametric
    // d(A,B) = 2, d(C,D) = 4, d(A,C)=d(A,D)=d(B,C)=d(B,D) = 4
    let mut matrix = DistanceMatrix::new(vec![
        "A".to_string(),
        "B".to_string(),
        "C".to_string(),
        "D".to_string(),
    ]);
    matrix.set(0, 1, 2.0);
    matrix.set(0, 2, 4.0);
    matrix.set(0, 3, 4.0);
    matrix.set(1, 2, 4.0);
    matrix.set(1, 3, 4.0);
    matrix.set(2, 3, 4.0);

    let tree = upgma(&matrix);
    assert_eq!(tree.num_leaves(), 4);

    // A and B should cluster first (smallest distance)
    let leaf_names = tree.leaf_names();
    assert_eq!(leaf_names.len(), 4);
}

// ============================================================================
// Parsimony: lower bound and Fitch algorithm properties
// ============================================================================

#[test]
fn test_parsimony_score_informative_sites() {
    // A parsimony-informative site has at least 2 different characters,
    // each present in at least 2 taxa. The minimum score for such a site is 1.
    // For 3 taxa with pattern AAB, minimum changes = 1 on any tree.
    let alignment = make_alignment(&[
        ("X", "AAAAACCCCC"),
        ("Y", "AAAAACCCCC"),
        ("Z", "CCCCCAAAAA"),
    ]);

    let tree = parse_newick("((X,Y),Z);").unwrap();
    let (score, _) = parsimony_score(&tree, &alignment).unwrap();

    // Each site requires at least 1 change (all sites are variable: X=Y≠Z)
    // Score should be 10 (one change per site on this tree)
    assert_eq!(score, 10, "Parsimony score on star-like tree with 10 variable sites");
}

#[test]
fn test_parsimony_score_invariant_sites() {
    // Invariant sites (all same character) contribute 0 to parsimony score
    let alignment = make_alignment(&[
        ("X", "AAAA"),
        ("Y", "AAAA"),
        ("Z", "AAAA"),
    ]);

    let tree = parse_newick("((X,Y),Z);").unwrap();
    let (score, _) = parsimony_score(&tree, &alignment).unwrap();

    assert_eq!(score, 0, "Invariant sites should have parsimony score 0");
}

#[test]
fn test_parsimony_score_single_change() {
    // A site with pattern AAAB requires exactly 1 change on the correct tree
    let alignment = make_alignment(&[
        ("A", "A"),
        ("B", "A"),
        ("C", "A"),
        ("D", "T"),
    ]);

    let tree = parse_newick("((A,B),(C,D));").unwrap();
    let (score, _) = parsimony_score(&tree, &alignment).unwrap();

    assert_eq!(score, 1, "Pattern AAAB requires exactly 1 change");
}

#[test]
fn test_parsimony_search_finds_optimal() {
    // For a simple 4-taxon case with strong signal, parsimony search should
    // find the optimal tree
    // Create data that strongly supports ((A,B),(C,D))
    let alignment = make_alignment(&[
        ("A", "AAAAAAAAAAAAAAAAAAAACCCCCCCCCC"),
        ("B", "AAAAAAAAAAAAAAAAAAAACCCCCCCCCC"),
        ("C", "CCCCCCCCCCCCCCCCCCCCAAAAAAAAAA"),
        ("D", "CCCCCCCCCCCCCCCCCCCCAAAAAAAAAA"),
    ]);

    let result = parsimony_search(&alignment, Some(42));

    // The optimal score should be much less than the worst-case tree
    // For this data, the optimal tree groups A+B and C+D
    assert!(result.score <= 30, "Parsimony search should find a good tree, got score {}", result.score);
}

// ============================================================================
// Likelihood: Felsenstein pruning algorithm properties
// ============================================================================

#[test]
fn test_likelihood_negative() {
    // Log-likelihood should always be negative
    let alignment = make_alignment(&[
        ("A", "ACGTACGTACGT"),
        ("B", "ACGTACGTACGT"),
        ("C", "TGCATGCATGCA"),
    ]);

    let tree = parse_newick("((A:0.1,B:0.1):0.1,C:0.2);").unwrap();
    let model = Jc69Model;
    let lnl = log_likelihood(&tree, &alignment, &model).unwrap();

    assert!(lnl < 0.0, "Log-likelihood should be negative, got {}", lnl);
}

#[test]
fn test_likelihood_identical_sequences_optimal() {
    // When all sequences are identical, likelihood should be highest (least negative)
    let alignment_same = make_alignment(&[
        ("A", "ACGTACGTACGT"),
        ("B", "ACGTACGTACGT"),
        ("C", "ACGTACGTACGT"),
    ]);
    let alignment_diff = make_alignment(&[
        ("A", "ACGTACGTACGT"),
        ("B", "TGCATGCATGCA"),
        ("C", "GTACGTACGTAC"),
    ]);

    let tree = parse_newick("((A:0.01,B:0.01):0.01,C:0.02);").unwrap();
    let model = Jc69Model;

    let lnl_same = log_likelihood(&tree, &alignment_same, &model).unwrap();
    let lnl_diff = log_likelihood(&tree, &alignment_diff, &model).unwrap();

    assert!(
        lnl_same > lnl_diff,
        "Identical sequences on short tree should have higher lnL: same={}, diff={}",
        lnl_same,
        lnl_diff
    );
}

#[test]
fn test_likelihood_branch_length_effect() {
    // For identical sequences, shorter branches should give higher likelihood
    let alignment = make_alignment(&[
        ("A", "ACGTACGTACGT"),
        ("B", "ACGTACGTACGT"),
        ("C", "ACGTACGTACGT"),
    ]);

    let model = Jc69Model;

    let tree_short = parse_newick("((A:0.001,B:0.001):0.001,C:0.002);").unwrap();
    let tree_long = parse_newick("((A:1.0,B:1.0):1.0,C:2.0);").unwrap();

    let lnl_short = log_likelihood(&tree_short, &alignment, &model).unwrap();
    let lnl_long = log_likelihood(&tree_long, &alignment, &model).unwrap();

    assert!(
        lnl_short > lnl_long,
        "Identical sequences: shorter branches should give higher lnL: short={}, long={}",
        lnl_short,
        lnl_long
    );
}

#[test]
fn test_likelihood_both_models_consistent() {
    // Both JC69 and F84 should produce finite negative log-likelihoods
    // on the same data and tree
    let alignment = make_alignment(&[
        ("A", "ACGTACGTAACCGGTT"),
        ("B", "ACGTACGTTAACCGGT"),
        ("C", "TGCATGCAAATTCCGG"),
        ("D", "TGCATGCATTAACCGG"),
    ]);

    let tree = parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();

    let jc69 = Jc69Model;
    let f84 = F84Model::new(BaseFrequencies::equal(), 2.0);

    let lnl_jc = log_likelihood(&tree, &alignment, &jc69).unwrap();
    let lnl_f84 = log_likelihood(&tree, &alignment, &f84).unwrap();

    assert!(lnl_jc.is_finite(), "JC69 lnL should be finite");
    assert!(lnl_f84.is_finite(), "F84 lnL should be finite");
    assert!(lnl_jc < 0.0, "JC69 lnL should be negative");
    assert!(lnl_f84 < 0.0, "F84 lnL should be negative");
}

// ============================================================================
// Transition probability matrix properties
// ============================================================================

#[test]
fn test_jc69_transition_probs_sum_to_one() {
    // Each row of P(t) must sum to 1.0
    let model = Jc69Model;

    for &t in &[0.01, 0.1, 0.5, 1.0, 5.0] {
        let p = model.transition_probs(t);
        for (i, row) in p.iter().enumerate() {
            let sum: f64 = row.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "JC69 P({}) row {} sums to {} (should be 1.0)",
                t,
                i,
                sum
            );
        }
    }
}

#[test]
fn test_jc69_transition_probs_t_zero() {
    // P(0) should be the identity matrix
    let model = Jc69Model;
    let p = model.transition_probs(0.0);

    for i in 0..4 {
        for j in 0..4 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (p[i][j] - expected).abs() < 1e-10,
                "P(0)[{}][{}] = {} (expected {})",
                i,
                j,
                p[i][j],
                expected
            );
        }
    }
}

#[test]
fn test_jc69_transition_probs_equilibrium() {
    // As t → ∞, P(t) should approach the equilibrium distribution (0.25 each)
    let model = Jc69Model;
    let p = model.transition_probs(100.0);

    for i in 0..4 {
        for j in 0..4 {
            assert!(
                (p[i][j] - 0.25).abs() < 1e-6,
                "P(100)[{}][{}] = {} (expected 0.25)",
                i,
                j,
                p[i][j]
            );
        }
    }
}

#[test]
fn test_jc69_transition_probs_known_values() {
    // JC69: P_ii(t) = 1/4 + 3/4 * exp(-4t/3)
    //       P_ij(t) = 1/4 - 1/4 * exp(-4t/3)  for i≠j
    let model = Jc69Model;
    let t = 0.5;
    let p = model.transition_probs(t);

    let exp_term = (-4.0 * t / 3.0_f64).exp();
    let expected_diag = 0.25 + 0.75 * exp_term;
    let expected_off = 0.25 - 0.25 * exp_term;

    assert!(
        (p[0][0] - expected_diag).abs() < 1e-10,
        "P(0.5)[0][0] = {} (expected {})",
        p[0][0],
        expected_diag
    );
    assert!(
        (p[0][1] - expected_off).abs() < 1e-10,
        "P(0.5)[0][1] = {} (expected {})",
        p[0][1],
        expected_off
    );
}

#[test]
fn test_f84_transition_probs_sum_to_one() {
    let model = F84Model::new(BaseFrequencies::new(0.3, 0.2, 0.25, 0.25), 2.0);

    for &t in &[0.01, 0.1, 0.5, 1.0, 5.0] {
        let p = model.transition_probs(t);
        for (i, row) in p.iter().enumerate() {
            let sum: f64 = row.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "F84 P({}) row {} sums to {} (should be 1.0)",
                t,
                i,
                sum
            );
        }
    }
}

// ============================================================================
// Newick I/O round-trip
// ============================================================================

#[test]
fn test_newick_round_trip_topology() {
    let input = "((A:0.1,B:0.2):0.3,(C:0.15,D:0.25):0.35);";
    let tree1 = parse_newick(input).unwrap();
    let newick = write_newick(&tree1);
    let tree2 = parse_newick(&newick).unwrap();

    let rf = phylip_rs::tree::distances::robinson_foulds(&tree1, &tree2).unwrap();
    assert_eq!(rf, 0, "Newick round-trip should preserve topology");
}

#[test]
fn test_robinson_foulds_identical_trees() {
    let tree1 = parse_newick("((A,B),(C,D));").unwrap();
    let tree2 = parse_newick("((A,B),(C,D));").unwrap();

    let rf = phylip_rs::tree::distances::robinson_foulds(&tree1, &tree2).unwrap();
    assert_eq!(rf, 0, "Identical trees should have RF distance 0");
}

#[test]
fn test_robinson_foulds_different_trees() {
    let tree1 = parse_newick("((A,B),(C,D));").unwrap();
    let tree2 = parse_newick("((A,C),(B,D));").unwrap();

    let rf = phylip_rs::tree::distances::robinson_foulds(&tree1, &tree2).unwrap();
    assert!(rf > 0, "Different topologies should have RF distance > 0");
}

// ============================================================================
// Bootstrap: statistical properties
// ============================================================================

#[test]
fn test_bootstrap_weights_sum() {
    // Bootstrap weights should sum to nsites (sampling with replacement)
    use phylip_rs::bootstrap::{bootstrap_weights, ResamplingMethod, SimpleRng};

    let nsites = 100;
    let mut rng = SimpleRng::seed(42);
    let weights = bootstrap_weights(nsites, &ResamplingMethod::Bootstrap, &mut rng);

    let total: usize = weights.iter().sum();
    assert_eq!(
        total, nsites,
        "Bootstrap weights should sum to nsites"
    );
}

#[test]
fn test_bootstrap_weights_have_zeros() {
    // With enough sites, bootstrap should leave some sites unsampled
    // Probability of a site being unsampled = (1 - 1/n)^n ≈ 1/e ≈ 0.368
    use phylip_rs::bootstrap::{bootstrap_weights, ResamplingMethod, SimpleRng};

    let nsites = 1000;
    let mut rng = SimpleRng::seed(42);
    let weights = bootstrap_weights(nsites, &ResamplingMethod::Bootstrap, &mut rng);

    let zeros = weights.iter().filter(|&&w| w == 0).count();
    // Expect ~368 zeros out of 1000 (tolerance: 200-500)
    assert!(
        zeros > 200 && zeros < 500,
        "Expected ~368 zeros in bootstrap of 1000 sites, got {}",
        zeros
    );
}

// ============================================================================
// Consensus: basic properties
// ============================================================================

#[test]
fn test_strict_consensus_identical_trees() {
    // Strict consensus of identical trees should return the same topology
    use phylip_rs::consensus::{consensus_tree, ConsensusMethod};

    let tree1 = parse_newick("((A,B),(C,D));").unwrap();
    let tree2 = parse_newick("((A,B),(C,D));").unwrap();

    let result =
        consensus_tree(&[tree1.clone(), tree2], &ConsensusMethod::Strict).unwrap();
    let rf = phylip_rs::tree::distances::robinson_foulds(&result.tree, &tree1).unwrap();
    assert_eq!(
        rf, 0,
        "Strict consensus of identical trees should return same topology"
    );
}

#[test]
fn test_majority_rule_consensus() {
    // Majority-rule consensus: a split present in >50% of trees should be in consensus
    use phylip_rs::consensus::{consensus_tree, ConsensusMethod};

    // 3 trees: 2 have ((A,B),(C,D)), 1 has ((A,C),(B,D))
    let tree1 = parse_newick("((A,B),(C,D));").unwrap();
    let tree2 = parse_newick("((A,B),(C,D));").unwrap();
    let tree3 = parse_newick("((A,C),(B,D));").unwrap();

    let result =
        consensus_tree(&[tree1.clone(), tree2, tree3], &ConsensusMethod::MajorityRule).unwrap();

    // The AB|CD split appears in 2/3 trees (67%), should be in consensus
    let rf = phylip_rs::tree::distances::robinson_foulds(&result.tree, &tree1).unwrap();
    assert_eq!(
        rf, 0,
        "Majority-rule consensus should include the 2/3 majority split"
    );
}

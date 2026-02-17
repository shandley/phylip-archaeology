//! Protein (amino acid) substitution models for phylogenetic likelihood calculations.
//!
//! This module provides protein sequence types, substitution models, FASTA parsing,
//! and a pruning algorithm for computing the likelihood of protein alignments on
//! phylogenetic trees.
//!
//! # Amino acid types
//!
//! The [`AminoAcid`] enum represents the 20 standard amino acids plus `Gap` and
//! `Unknown`. The standard order used throughout this module is:
//!
//! A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V
//!
//! # Substitution models
//!
//! Two protein substitution models are provided:
//!
//! - [`PoissonModel`] — Equal exchangeabilities with configurable equilibrium
//!   frequencies. The simplest protein model.
//! - [`WagFreqModel`] — Poisson exchangeabilities with the published WAG
//!   equilibrium frequencies (Whelan & Goldman, 2001).
//!
//! Both implement the [`ProteinModel`] trait, which is analogous to the DNA
//! [`SubstitutionModel`](crate::likelihood::models::SubstitutionModel) trait
//! but operates on 20 states instead of 4.
//!
//! # Example
//!
//! ```
//! use phylip_rs::models::protein::{
//!     AminoAcid, ProteinSequence, ProteinAlignment,
//!     ProteinModel, PoissonModel, WagFreqModel,
//!     read_protein_fasta, protein_log_likelihood,
//! };
//! use phylip_rs::tree::newick::parse_newick;
//!
//! // Parse a protein FASTA alignment
//! let fasta = ">Seq1\nACDEFG\n>Seq2\nACDEYG\n>Seq3\nMCDEFG\n";
//! let alignment = read_protein_fasta(fasta).unwrap();
//! assert_eq!(alignment.n_sites, 6);
//!
//! // Evaluate likelihood under the Poisson model
//! let tree = parse_newick("((Seq1:0.1,Seq2:0.1):0.05,Seq3:0.2);").unwrap();
//! let model = PoissonModel::equal_frequencies();
//! let lnl = protein_log_likelihood(&tree, &alignment, &model).unwrap();
//! assert!(lnl < 0.0);
//! ```
//!
//! # References
//!
//! - Jones, D.T., Taylor, W.R. & Thornton, J.M. (1992). The rapid generation
//!   of mutation data matrices from protein sequences. *CABIOS*, 8, 275-282.
//! - Whelan, S. & Goldman, N. (2001). A general empirical model of protein
//!   evolution derived from multiple protein families using a maximum-likelihood
//!   approach. *Molecular Biology and Evolution*, 18, 691-699.
//! - Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum
//!   likelihood approach. *Journal of Molecular Evolution*, 17, 368-376.

use std::fmt;

use crate::tree::types::Tree;

/// Number of standard amino acid states.
const N_AA: usize = 20;

/// The 20 standard amino acids in canonical order, plus Gap and Unknown.
///
/// The canonical order is: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V.
/// This ordering matches the convention used by the JTT and WAG models.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AminoAcid {
    Ala, // A - Alanine
    Arg, // R - Arginine
    Asn, // N - Asparagine
    Asp, // D - Aspartic acid
    Cys, // C - Cysteine
    Gln, // Q - Glutamine
    Glu, // E - Glutamic acid
    Gly, // G - Glycine
    His, // H - Histidine
    Ile, // I - Isoleucine
    Leu, // L - Leucine
    Lys, // K - Lysine
    Met, // M - Methionine
    Phe, // F - Phenylalanine
    Pro, // P - Proline
    Ser, // S - Serine
    Thr, // T - Threonine
    Trp, // W - Tryptophan
    Tyr, // Y - Tyrosine
    Val, // V - Valine
    Gap,
    Unknown,
}

/// All 20 standard amino acids in canonical order.
const ALL_AMINO_ACIDS: [AminoAcid; N_AA] = [
    AminoAcid::Ala,
    AminoAcid::Arg,
    AminoAcid::Asn,
    AminoAcid::Asp,
    AminoAcid::Cys,
    AminoAcid::Gln,
    AminoAcid::Glu,
    AminoAcid::Gly,
    AminoAcid::His,
    AminoAcid::Ile,
    AminoAcid::Leu,
    AminoAcid::Lys,
    AminoAcid::Met,
    AminoAcid::Phe,
    AminoAcid::Pro,
    AminoAcid::Ser,
    AminoAcid::Thr,
    AminoAcid::Trp,
    AminoAcid::Tyr,
    AminoAcid::Val,
];

impl AminoAcid {
    /// Parse a single character into an AminoAcid.
    ///
    /// Recognizes all standard one-letter amino acid codes (case-insensitive),
    /// gap characters (`-` and `.`), and unknown/ambiguous characters
    /// (`X`, `?`, `*`, `B`, `Z`, `J`).
    pub fn from_char(c: char) -> AminoAcid {
        match c.to_ascii_uppercase() {
            'A' => AminoAcid::Ala,
            'R' => AminoAcid::Arg,
            'N' => AminoAcid::Asn,
            'D' => AminoAcid::Asp,
            'C' => AminoAcid::Cys,
            'Q' => AminoAcid::Gln,
            'E' => AminoAcid::Glu,
            'G' => AminoAcid::Gly,
            'H' => AminoAcid::His,
            'I' => AminoAcid::Ile,
            'L' => AminoAcid::Leu,
            'K' => AminoAcid::Lys,
            'M' => AminoAcid::Met,
            'F' => AminoAcid::Phe,
            'P' => AminoAcid::Pro,
            'S' => AminoAcid::Ser,
            'T' => AminoAcid::Thr,
            'W' => AminoAcid::Trp,
            'Y' => AminoAcid::Tyr,
            'V' => AminoAcid::Val,
            '-' | '.' => AminoAcid::Gap,
            // Ambiguity codes and stop codons treated as unknown
            'X' | '?' | '*' | 'B' | 'Z' | 'J' | 'U' | 'O' => AminoAcid::Unknown,
            _ => AminoAcid::Unknown,
        }
    }

    /// Convert to the standard one-letter code.
    pub fn to_char(self) -> char {
        match self {
            AminoAcid::Ala => 'A',
            AminoAcid::Arg => 'R',
            AminoAcid::Asn => 'N',
            AminoAcid::Asp => 'D',
            AminoAcid::Cys => 'C',
            AminoAcid::Gln => 'Q',
            AminoAcid::Glu => 'E',
            AminoAcid::Gly => 'G',
            AminoAcid::His => 'H',
            AminoAcid::Ile => 'I',
            AminoAcid::Leu => 'L',
            AminoAcid::Lys => 'K',
            AminoAcid::Met => 'M',
            AminoAcid::Phe => 'F',
            AminoAcid::Pro => 'P',
            AminoAcid::Ser => 'S',
            AminoAcid::Thr => 'T',
            AminoAcid::Trp => 'W',
            AminoAcid::Tyr => 'Y',
            AminoAcid::Val => 'V',
            AminoAcid::Gap => '-',
            AminoAcid::Unknown => 'X',
        }
    }

    /// Returns the index (0-19) for standard amino acids, or None for Gap/Unknown.
    pub fn index(self) -> Option<usize> {
        match self {
            AminoAcid::Ala => Some(0),
            AminoAcid::Arg => Some(1),
            AminoAcid::Asn => Some(2),
            AminoAcid::Asp => Some(3),
            AminoAcid::Cys => Some(4),
            AminoAcid::Gln => Some(5),
            AminoAcid::Glu => Some(6),
            AminoAcid::Gly => Some(7),
            AminoAcid::His => Some(8),
            AminoAcid::Ile => Some(9),
            AminoAcid::Leu => Some(10),
            AminoAcid::Lys => Some(11),
            AminoAcid::Met => Some(12),
            AminoAcid::Phe => Some(13),
            AminoAcid::Pro => Some(14),
            AminoAcid::Ser => Some(15),
            AminoAcid::Thr => Some(16),
            AminoAcid::Trp => Some(17),
            AminoAcid::Tyr => Some(18),
            AminoAcid::Val => Some(19),
            AminoAcid::Gap | AminoAcid::Unknown => None,
        }
    }

    /// Returns true if this is one of the 20 standard amino acids.
    pub fn is_standard(self) -> bool {
        self.index().is_some()
    }

    /// Returns a 20-element conditional likelihood vector.
    ///
    /// For standard amino acids, the vector has 1.0 at the amino acid's index
    /// and 0.0 elsewhere. For Gap and Unknown, all entries are 1.0 (fully
    /// ambiguous).
    pub fn likelihood_vector(self) -> [f64; N_AA] {
        match self.index() {
            Some(idx) => {
                let mut v = [0.0; N_AA];
                v[idx] = 1.0;
                v
            }
            None => [1.0; N_AA], // Gap or Unknown: all states possible
        }
    }
}

impl fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

/// A named protein sequence.
#[derive(Debug, Clone)]
pub struct ProteinSequence {
    pub name: String,
    pub residues: Vec<AminoAcid>,
}

impl ProteinSequence {
    pub fn new(name: impl Into<String>, residues: Vec<AminoAcid>) -> Self {
        Self {
            name: name.into(),
            residues,
        }
    }

    pub fn len(&self) -> usize {
        self.residues.len()
    }

    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }
}

/// A multiple protein sequence alignment.
#[derive(Debug, Clone)]
pub struct ProteinAlignment {
    pub sequences: Vec<ProteinSequence>,
    pub n_sites: usize,
}

impl ProteinAlignment {
    /// Create a new protein alignment from a vector of sequences.
    ///
    /// All sequences must have the same length. Returns an error if the
    /// sequences are empty or have unequal lengths.
    pub fn new(sequences: Vec<ProteinSequence>) -> Result<Self, String> {
        if sequences.is_empty() {
            return Err("alignment contains no sequences".to_string());
        }
        let n_sites = sequences[0].len();
        if sequences.iter().any(|s| s.len() != n_sites) {
            return Err("sequences have different lengths".to_string());
        }
        Ok(Self {
            sequences,
            n_sites,
        })
    }

    /// Number of taxa (sequences).
    pub fn ntaxa(&self) -> usize {
        self.sequences.len()
    }

    /// Get the name of taxon i.
    pub fn taxon_name(&self, i: usize) -> &str {
        &self.sequences[i].name
    }

    /// Get the amino acid at (taxon, site).
    pub fn residue(&self, taxon: usize, site: usize) -> AminoAcid {
        self.sequences[taxon].residues[site]
    }
}

// ---------------------------------------------------------------------------
// Protein substitution model trait
// ---------------------------------------------------------------------------

/// A protein substitution model for the pruning algorithm.
///
/// Provides transition probability matrices P(t) and equilibrium amino acid
/// frequencies needed for computing tree likelihoods with 20 amino acid states.
///
/// States are indexed in canonical order:
/// A(0), R(1), N(2), D(3), C(4), Q(5), E(6), G(7), H(8), I(9),
/// L(10), K(11), M(12), F(13), P(14), S(15), T(16), W(17), Y(18), V(19).
pub trait ProteinModel {
    /// Compute the 20x20 transition probability matrix P(t).
    ///
    /// `P[i][j]` = Pr(state j at descendant | state i at ancestor),
    /// given branch length `t` in expected substitutions per site.
    fn transition_probs(&self, t: f64) -> [[f64; N_AA]; N_AA];

    /// Equilibrium amino acid frequencies, in canonical order.
    fn frequencies(&self) -> [f64; N_AA];

    /// Model name.
    fn name(&self) -> &str;
}

// ---------------------------------------------------------------------------
// Poisson model
// ---------------------------------------------------------------------------

/// Poisson protein substitution model.
///
/// All amino acid exchange rates are equal. The rate matrix Q is defined as:
///
/// Q[i][j] = pi[j]   for i != j
/// Q[i][i] = -(1 - pi[i])
///
/// where pi[j] are the equilibrium frequencies. This model has a simple
/// eigenvalue structure: eigenvalue 0 (once) and eigenvalue -1 (19 times),
/// giving the transition probabilities:
///
/// P(i->j; t) = pi[j] + (delta(i,j) - pi[j]) * exp(-t)
#[derive(Debug, Clone)]
pub struct PoissonModel {
    frequencies: [f64; N_AA],
}

impl PoissonModel {
    /// Create a Poisson model with custom equilibrium frequencies.
    ///
    /// The frequencies should sum to 1.0.
    pub fn new(frequencies: [f64; N_AA]) -> Self {
        Self { frequencies }
    }

    /// Create a Poisson model with equal frequencies (1/20 each).
    pub fn equal_frequencies() -> Self {
        Self {
            frequencies: [1.0 / N_AA as f64; N_AA],
        }
    }
}

impl ProteinModel for PoissonModel {
    fn transition_probs(&self, t: f64) -> [[f64; N_AA]; N_AA] {
        let mut p = [[0.0f64; N_AA]; N_AA];
        let exp_neg_t = (-t).exp();

        for i in 0..N_AA {
            for j in 0..N_AA {
                if i == j {
                    p[i][j] = self.frequencies[j]
                        + (1.0 - self.frequencies[j]) * exp_neg_t;
                } else {
                    p[i][j] = self.frequencies[j] * (1.0 - exp_neg_t);
                }
            }
        }

        p
    }

    fn frequencies(&self) -> [f64; N_AA] {
        self.frequencies
    }

    fn name(&self) -> &str {
        "Poisson"
    }
}

// ---------------------------------------------------------------------------
// WAG equilibrium frequencies
// ---------------------------------------------------------------------------

/// Published WAG equilibrium frequencies (Whelan & Goldman, 2001).
///
/// Order: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V
const WAG_FREQUENCIES: [f64; N_AA] = [
    0.0866279, // A
    0.0439720, // R
    0.0390894, // N
    0.0570451, // D
    0.0193078, // C
    0.0367281, // Q
    0.0580589, // E
    0.0832518, // G
    0.0244313, // H
    0.0484660, // I
    0.0862090, // L
    0.0620286, // K
    0.0195027, // M
    0.0384319, // F
    0.0457631, // P
    0.0695179, // S
    0.0610127, // T
    0.0143859, // W
    0.0352742, // Y
    0.0708956, // V
];

/// WAG-frequencies protein model.
///
/// Uses Poisson (equal) exchangeabilities with the published WAG equilibrium
/// frequencies. This captures the amino acid composition bias observed in
/// real protein evolution without requiring the full WAG rate matrix.
///
/// Reference: Whelan, S. & Goldman, N. (2001). A general empirical model of
/// protein evolution derived from multiple protein families using a
/// maximum-likelihood approach. *Mol. Biol. Evol.*, 18, 691-699.
#[derive(Debug, Clone)]
pub struct WagFreqModel {
    inner: PoissonModel,
}

impl WagFreqModel {
    /// Create a new WAG-frequencies model.
    ///
    /// The published WAG frequencies are normalized to sum to exactly 1.0
    /// to ensure numerical correctness of transition probability matrices.
    pub fn new() -> Self {
        // Normalize the published frequencies to sum exactly to 1.0
        let sum: f64 = WAG_FREQUENCIES.iter().sum();
        let mut freqs = WAG_FREQUENCIES;
        for f in &mut freqs {
            *f /= sum;
        }
        Self {
            inner: PoissonModel::new(freqs),
        }
    }
}

impl Default for WagFreqModel {
    fn default() -> Self {
        Self::new()
    }
}

impl ProteinModel for WagFreqModel {
    fn transition_probs(&self, t: f64) -> [[f64; N_AA]; N_AA] {
        self.inner.transition_probs(t)
    }

    fn frequencies(&self) -> [f64; N_AA] {
        self.inner.frequencies()
    }

    fn name(&self) -> &str {
        "WAG-freq"
    }
}

// ---------------------------------------------------------------------------
// Protein FASTA parser
// ---------------------------------------------------------------------------

/// Parse a FASTA format string into a `ProteinAlignment`.
///
/// This function follows the same conventions as the DNA FASTA parser:
/// - Headers start with `>`
/// - Sequence name is the first whitespace-delimited token after `>`
/// - Sequence data may span multiple lines
/// - Whitespace and digits within sequence lines are ignored
/// - All sequences must have the same length
///
/// Any character that is not a standard amino acid code is treated as
/// `Unknown` (equivalent to fully ambiguous).
///
/// # Examples
///
/// ```
/// use phylip_rs::models::protein::read_protein_fasta;
///
/// let input = ">Seq1\nACDEFGHIKLMN\n>Seq2\nACDEFGHIKLMN\n";
/// let alignment = read_protein_fasta(input).unwrap();
/// assert_eq!(alignment.ntaxa(), 2);
/// assert_eq!(alignment.n_sites, 12);
/// ```
pub fn read_protein_fasta(input: &str) -> Result<ProteinAlignment, String> {
    let input = input.trim();
    if input.is_empty() {
        return Err("empty input".to_string());
    }

    let lines: Vec<&str> = input.lines().collect();

    // The first non-empty, non-comment line should start with '>'
    let first_content = lines.iter().find(|l| {
        let trimmed = l.trim();
        !trimmed.is_empty() && !trimmed.starts_with('#')
    });
    match first_content {
        Some(line) if !line.trim().starts_with('>') => {
            return Err("FASTA input does not start with '>' header".to_string());
        }
        None => return Err("empty input".to_string()),
        _ => {}
    }

    let mut sequences: Vec<ProteinSequence> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_residues: Vec<AminoAcid> = Vec::new();

    for line in &lines {
        let trimmed = line.trim();

        // Skip empty lines and comments
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed.starts_with('>') {
            // Save previous sequence if any
            if let Some(name) = current_name.take() {
                if current_residues.is_empty() {
                    return Err(format!("sequence '{}' has no data", name));
                }
                sequences.push(ProteinSequence::new(name, current_residues));
                current_residues = Vec::new();
            }

            // Parse header line
            let header = &trimmed[1..];
            let name = header
                .trim()
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            if name.is_empty() {
                return Err("empty sequence name".to_string());
            }
            current_name = Some(name);
        } else {
            // Sequence data line
            if current_name.is_none() {
                return Err("FASTA input does not start with '>' header".to_string());
            }

            for c in trimmed.chars() {
                if c.is_whitespace() || c.is_ascii_digit() {
                    continue;
                }
                current_residues.push(AminoAcid::from_char(c));
            }
        }
    }

    // Save the last sequence
    if let Some(name) = current_name.take() {
        if current_residues.is_empty() {
            return Err(format!("sequence '{}' has no data", name));
        }
        sequences.push(ProteinSequence::new(name, current_residues));
    }

    if sequences.is_empty() {
        return Err("empty input".to_string());
    }

    ProteinAlignment::new(sequences)
}

// ---------------------------------------------------------------------------
// Protein pruning algorithm
// ---------------------------------------------------------------------------

/// Minimum branch length for numerical stability.
const MIN_BRANCH_LENGTH: f64 = 1e-8;
/// Default branch length when none is specified.
const DEFAULT_BRANCH_LENGTH: f64 = 0.1;
/// Underflow threshold: rescale when max conditional likelihood drops below this.
const UNDERFLOW_THRESHOLD: f64 = 1e-100;

/// Conditional likelihood vector at a node for one site (20 amino acid states).
type CondLike = [f64; N_AA];

/// Per-node data for the protein pruning algorithm.
#[derive(Debug, Clone)]
struct NodeLikelihoods {
    /// Conditional likelihood vectors, one per site.
    cond: Vec<CondLike>,
    /// Log-scaling factor per site (for underflow prevention).
    log_scale: Vec<f64>,
}

impl NodeLikelihoods {
    fn new(nsites: usize) -> Self {
        Self {
            cond: vec![[0.0; N_AA]; nsites],
            log_scale: vec![0.0; nsites],
        }
    }
}

/// Get the effective branch length for a node, using default if missing.
fn effective_branch_length(tree: &Tree, node_id: usize) -> f64 {
    tree.nodes[node_id]
        .branch_length
        .unwrap_or(DEFAULT_BRANCH_LENGTH)
        .max(MIN_BRANCH_LENGTH)
}

/// Compute the total log-likelihood of a protein alignment on a phylogenetic tree.
///
/// This implements Felsenstein's pruning algorithm (1981) for 20 amino acid
/// states. The algorithm traverses the tree in postorder, computing conditional
/// likelihoods at each internal node, and combines them at the root using
/// the model's equilibrium frequencies.
///
/// The tree topology (from [`Tree`]) is reused from the DNA implementation;
/// only the data type and state space differ.
///
/// # Arguments
///
/// * `tree` - The phylogenetic tree with branch lengths
/// * `alignment` - The protein multiple sequence alignment
/// * `model` - The protein substitution model
///
/// # Returns
///
/// The total log-likelihood (sum over all sites), or an error string
/// if the tree leaves do not match the alignment taxa.
///
/// # Example
///
/// ```
/// use phylip_rs::models::protein::*;
/// use phylip_rs::tree::newick::parse_newick;
///
/// let alignment = ProteinAlignment::new(vec![
///     ProteinSequence::new("A", vec![AminoAcid::Ala, AminoAcid::Arg]),
///     ProteinSequence::new("B", vec![AminoAcid::Ala, AminoAcid::Asn]),
///     ProteinSequence::new("C", vec![AminoAcid::Cys, AminoAcid::Arg]),
/// ]).unwrap();
///
/// let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
/// let model = PoissonModel::equal_frequencies();
/// let lnl = protein_log_likelihood(&tree, &alignment, &model).unwrap();
/// assert!(lnl < 0.0);
/// ```
pub fn protein_log_likelihood(
    tree: &Tree,
    alignment: &ProteinAlignment,
    model: &dyn ProteinModel,
) -> Result<f64, String> {
    let nsites = alignment.n_sites;
    let nnodes = tree.num_nodes();

    if tree.num_leaves() == 0 {
        return Err("tree has no leaves".to_string());
    }

    // Build name -> alignment index map
    let name_to_idx: std::collections::HashMap<&str, usize> = alignment
        .sequences
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.as_str(), i))
        .collect();

    // Initialize node data
    let mut data = vec![NodeLikelihoods::new(nsites); nnodes];

    // Set leaf conditional likelihoods from observed data
    for node in &tree.nodes {
        if node.is_leaf() {
            if let Some(name) = &node.name {
                let taxon_idx = name_to_idx
                    .get(name.as_str())
                    .ok_or_else(|| format!("leaf '{}' not found in alignment", name))?;
                for site in 0..nsites {
                    data[node.id].cond[site] =
                        alignment.residue(*taxon_idx, site).likelihood_vector();
                }
            }
        }
    }

    // Postorder traversal: pruning algorithm
    let postorder = tree.postorder();
    for &node_id in &postorder {
        let children = &tree.nodes[node_id].children;
        if children.is_empty() {
            continue; // leaf: already initialized
        }

        // For each child, precompute the transition probability matrix
        let child_data: Vec<(usize, [[f64; N_AA]; N_AA])> = children
            .iter()
            .map(|&child_id| {
                let t = effective_branch_length(tree, child_id);
                let p = model.transition_probs(t);
                (child_id, p)
            })
            .collect();

        for site in 0..nsites {
            let mut node_cond = [1.0f64; N_AA];
            let mut node_log_scale = 0.0f64;

            for &(child_id, ref p_matrix) in &child_data {
                let child_cond = &data[child_id].cond[site];
                node_log_scale += data[child_id].log_scale[site];

                let mut partial = [0.0f64; N_AA];
                for b in 0..N_AA {
                    let mut sum = 0.0;
                    for j in 0..N_AA {
                        sum += p_matrix[b][j] * child_cond[j];
                    }
                    partial[b] = sum;
                }

                // Multiply into the running product
                for b in 0..N_AA {
                    node_cond[b] *= partial[b];
                }
            }

            // Underflow prevention: rescale if values are too small
            let max_val = node_cond.iter().copied().fold(0.0f64, f64::max);

            if max_val > 0.0 && max_val < UNDERFLOW_THRESHOLD {
                for b in 0..N_AA {
                    node_cond[b] /= max_val;
                }
                node_log_scale += max_val.ln();
            }

            data[node_id].cond[site] = node_cond;
            data[node_id].log_scale[site] = node_log_scale;
        }
    }

    // Compute total log-likelihood at the root
    let root = tree.root;
    let pi = model.frequencies();
    let mut total_lnl = 0.0;

    for site in 0..nsites {
        let root_cond = &data[root].cond[site];
        let log_scale = data[root].log_scale[site];

        let mut site_like = 0.0;
        for b in 0..N_AA {
            site_like += pi[b] * root_cond[b];
        }

        if site_like <= 0.0 {
            return Err(format!("site {} has zero likelihood", site));
        }

        total_lnl += site_like.ln() + log_scale;
    }

    Ok(total_lnl)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::newick::parse_newick;

    // -----------------------------------------------------------------------
    // AminoAcid tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_amino_acid_from_char_standard() {
        let chars = "ARNDCQEGHILKMFPSTWYV";
        let expected = ALL_AMINO_ACIDS;
        for (c, &aa) in chars.chars().zip(expected.iter()) {
            assert_eq!(AminoAcid::from_char(c), aa, "failed for char '{}'", c);
        }
    }

    #[test]
    fn test_amino_acid_from_char_lowercase() {
        let chars = "arndcqeghilkmfpstwyv";
        let expected = ALL_AMINO_ACIDS;
        for (c, &aa) in chars.chars().zip(expected.iter()) {
            assert_eq!(AminoAcid::from_char(c), aa, "failed for lowercase '{}'", c);
        }
    }

    #[test]
    fn test_amino_acid_from_char_gap() {
        assert_eq!(AminoAcid::from_char('-'), AminoAcid::Gap);
        assert_eq!(AminoAcid::from_char('.'), AminoAcid::Gap);
    }

    #[test]
    fn test_amino_acid_from_char_unknown() {
        assert_eq!(AminoAcid::from_char('X'), AminoAcid::Unknown);
        assert_eq!(AminoAcid::from_char('?'), AminoAcid::Unknown);
        assert_eq!(AminoAcid::from_char('*'), AminoAcid::Unknown);
        assert_eq!(AminoAcid::from_char('B'), AminoAcid::Unknown);
        assert_eq!(AminoAcid::from_char('Z'), AminoAcid::Unknown);
    }

    #[test]
    fn test_amino_acid_to_char_roundtrip() {
        for &aa in &ALL_AMINO_ACIDS {
            let c = aa.to_char();
            let aa2 = AminoAcid::from_char(c);
            assert_eq!(aa, aa2, "roundtrip failed for {:?}", aa);
        }
    }

    #[test]
    fn test_amino_acid_index() {
        for (i, &aa) in ALL_AMINO_ACIDS.iter().enumerate() {
            assert_eq!(aa.index(), Some(i), "index failed for {:?}", aa);
        }
        assert_eq!(AminoAcid::Gap.index(), None);
        assert_eq!(AminoAcid::Unknown.index(), None);
    }

    #[test]
    fn test_amino_acid_is_standard() {
        for &aa in &ALL_AMINO_ACIDS {
            assert!(aa.is_standard(), "{:?} should be standard", aa);
        }
        assert!(!AminoAcid::Gap.is_standard());
        assert!(!AminoAcid::Unknown.is_standard());
    }

    #[test]
    fn test_amino_acid_likelihood_vector_standard() {
        for (i, &aa) in ALL_AMINO_ACIDS.iter().enumerate() {
            let v = aa.likelihood_vector();
            for j in 0..N_AA {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert_eq!(
                    v[j], expected,
                    "likelihood_vector[{}] for {:?}: got {}, expected {}",
                    j, aa, v[j], expected
                );
            }
        }
    }

    #[test]
    fn test_amino_acid_likelihood_vector_gap_unknown() {
        let gap_v = AminoAcid::Gap.likelihood_vector();
        let unk_v = AminoAcid::Unknown.likelihood_vector();
        for j in 0..N_AA {
            assert_eq!(gap_v[j], 1.0);
            assert_eq!(unk_v[j], 1.0);
        }
    }

    #[test]
    fn test_amino_acid_display() {
        assert_eq!(format!("{}", AminoAcid::Ala), "A");
        assert_eq!(format!("{}", AminoAcid::Trp), "W");
        assert_eq!(format!("{}", AminoAcid::Gap), "-");
        assert_eq!(format!("{}", AminoAcid::Unknown), "X");
    }

    // -----------------------------------------------------------------------
    // ProteinSequence and ProteinAlignment tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_protein_sequence_basic() {
        let seq = ProteinSequence::new(
            "test",
            vec![AminoAcid::Ala, AminoAcid::Arg, AminoAcid::Asn],
        );
        assert_eq!(seq.name, "test");
        assert_eq!(seq.len(), 3);
        assert!(!seq.is_empty());
    }

    #[test]
    fn test_protein_alignment_new() {
        let s1 = ProteinSequence::new("A", vec![AminoAcid::Ala, AminoAcid::Arg]);
        let s2 = ProteinSequence::new("B", vec![AminoAcid::Cys, AminoAcid::Asp]);
        let aln = ProteinAlignment::new(vec![s1, s2]).unwrap();
        assert_eq!(aln.ntaxa(), 2);
        assert_eq!(aln.n_sites, 2);
        assert_eq!(aln.taxon_name(0), "A");
        assert_eq!(aln.residue(0, 0), AminoAcid::Ala);
        assert_eq!(aln.residue(1, 1), AminoAcid::Asp);
    }

    #[test]
    fn test_protein_alignment_empty() {
        let result = ProteinAlignment::new(vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn test_protein_alignment_unequal() {
        let s1 = ProteinSequence::new("A", vec![AminoAcid::Ala]);
        let s2 = ProteinSequence::new("B", vec![AminoAcid::Cys, AminoAcid::Asp]);
        let result = ProteinAlignment::new(vec![s1, s2]);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // Poisson model tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_poisson_identity_at_zero() {
        let model = PoissonModel::equal_frequencies();
        let p = model.transition_probs(0.0);
        for i in 0..N_AA {
            for j in 0..N_AA {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (p[i][j] - expected).abs() < 1e-12,
                    "P[{}][{}] = {}, expected {}",
                    i, j, p[i][j], expected
                );
            }
        }
    }

    #[test]
    fn test_poisson_equilibrium_at_large_t() {
        let model = PoissonModel::equal_frequencies();
        let p = model.transition_probs(100.0);
        let pi = model.frequencies();
        for i in 0..N_AA {
            for j in 0..N_AA {
                assert!(
                    (p[i][j] - pi[j]).abs() < 1e-10,
                    "P[{}][{}] = {}, expected pi[{}] = {}",
                    i, j, p[i][j], j, pi[j]
                );
            }
        }
    }

    #[test]
    fn test_poisson_rows_sum_to_one() {
        let model = PoissonModel::equal_frequencies();
        for &t in &[0.0, 0.01, 0.1, 0.5, 1.0, 5.0, 100.0] {
            let p = model.transition_probs(t);
            for i in 0..N_AA {
                let row_sum: f64 = p[i].iter().sum();
                assert!(
                    (row_sum - 1.0).abs() < 1e-12,
                    "Row {} sums to {} at t={}",
                    i, row_sum, t
                );
            }
        }
    }

    #[test]
    fn test_poisson_all_non_negative() {
        let model = PoissonModel::equal_frequencies();
        for &t in &[0.0, 0.001, 0.01, 0.1, 0.5, 1.0, 5.0] {
            let p = model.transition_probs(t);
            for i in 0..N_AA {
                for j in 0..N_AA {
                    assert!(
                        p[i][j] >= 0.0,
                        "P[{}][{}] = {} at t={} is negative",
                        i, j, p[i][j], t
                    );
                }
            }
        }
    }

    #[test]
    fn test_poisson_detailed_balance() {
        let model = PoissonModel::equal_frequencies();
        let p = model.transition_probs(0.3);
        let pi = model.frequencies();
        for i in 0..N_AA {
            for j in 0..N_AA {
                let lhs = pi[i] * p[i][j];
                let rhs = pi[j] * p[j][i];
                assert!(
                    (lhs - rhs).abs() < 1e-12,
                    "Detailed balance violated: pi[{}]*P[{}][{}] = {} != pi[{}]*P[{}][{}] = {}",
                    i, i, j, lhs, j, j, i, rhs
                );
            }
        }
    }

    #[test]
    fn test_poisson_name() {
        let model = PoissonModel::equal_frequencies();
        assert_eq!(model.name(), "Poisson");
    }

    // -----------------------------------------------------------------------
    // WAG-freq model tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_wag_frequencies_sum_to_one() {
        let model = WagFreqModel::new();
        let pi = model.frequencies();
        let sum: f64 = pi.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-6,
            "WAG frequencies sum to {}, expected 1.0",
            sum
        );
    }

    #[test]
    fn test_wag_frequencies_all_positive() {
        let model = WagFreqModel::new();
        let pi = model.frequencies();
        for (i, &f) in pi.iter().enumerate() {
            assert!(
                f > 0.0,
                "WAG frequency[{}] = {} should be positive",
                i, f
            );
        }
    }

    #[test]
    fn test_wag_identity_at_zero() {
        let model = WagFreqModel::new();
        let p = model.transition_probs(0.0);
        for i in 0..N_AA {
            for j in 0..N_AA {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (p[i][j] - expected).abs() < 1e-12,
                    "WAG P[{}][{}] = {}, expected {}",
                    i, j, p[i][j], expected
                );
            }
        }
    }

    #[test]
    fn test_wag_equilibrium_at_large_t() {
        let model = WagFreqModel::new();
        let p = model.transition_probs(100.0);
        let pi = model.frequencies();
        for i in 0..N_AA {
            for j in 0..N_AA {
                assert!(
                    (p[i][j] - pi[j]).abs() < 1e-10,
                    "WAG P[{}][{}] = {}, expected pi[{}] = {}",
                    i, j, p[i][j], j, pi[j]
                );
            }
        }
    }

    #[test]
    fn test_wag_rows_sum_to_one() {
        let model = WagFreqModel::new();
        for &t in &[0.0, 0.01, 0.1, 0.5, 1.0, 5.0] {
            let p = model.transition_probs(t);
            for i in 0..N_AA {
                let row_sum: f64 = p[i].iter().sum();
                assert!(
                    (row_sum - 1.0).abs() < 1e-12,
                    "WAG row {} sums to {} at t={}",
                    i, row_sum, t
                );
            }
        }
    }

    #[test]
    fn test_wag_all_non_negative() {
        let model = WagFreqModel::new();
        for &t in &[0.0, 0.001, 0.01, 0.1, 0.5, 1.0, 5.0] {
            let p = model.transition_probs(t);
            for i in 0..N_AA {
                for j in 0..N_AA {
                    assert!(
                        p[i][j] >= 0.0,
                        "WAG P[{}][{}] = {} at t={} is negative",
                        i, j, p[i][j], t
                    );
                }
            }
        }
    }

    #[test]
    fn test_wag_detailed_balance() {
        let model = WagFreqModel::new();
        let p = model.transition_probs(0.5);
        let pi = model.frequencies();
        for i in 0..N_AA {
            for j in 0..N_AA {
                let lhs = pi[i] * p[i][j];
                let rhs = pi[j] * p[j][i];
                assert!(
                    (lhs - rhs).abs() < 1e-12,
                    "WAG detailed balance violated at ({},{}): {} != {}",
                    i, j, lhs, rhs
                );
            }
        }
    }

    #[test]
    fn test_wag_name() {
        let model = WagFreqModel::new();
        assert_eq!(model.name(), "WAG-freq");
    }

    #[test]
    fn test_wag_default() {
        let m1 = WagFreqModel::new();
        let m2 = WagFreqModel::default();
        assert_eq!(m1.frequencies(), m2.frequencies());
    }

    // -----------------------------------------------------------------------
    // Poisson model with unequal frequencies
    // -----------------------------------------------------------------------

    #[test]
    fn test_poisson_unequal_freqs_rows_sum_to_one() {
        // Use the normalized WAG frequencies via WagFreqModel
        let model = WagFreqModel::new();
        for &t in &[0.0, 0.01, 0.1, 0.5, 1.0, 5.0] {
            let p = model.transition_probs(t);
            for i in 0..N_AA {
                let row_sum: f64 = p[i].iter().sum();
                assert!(
                    (row_sum - 1.0).abs() < 1e-12,
                    "Unequal freq row {} sums to {} at t={}",
                    i, row_sum, t
                );
            }
        }
    }

    #[test]
    fn test_poisson_unequal_freqs_equilibrium() {
        // Use the normalized WAG frequencies via WagFreqModel
        let model = WagFreqModel::new();
        let p = model.transition_probs(100.0);
        let pi = model.frequencies();
        for i in 0..N_AA {
            for j in 0..N_AA {
                assert!(
                    (p[i][j] - pi[j]).abs() < 1e-10,
                    "Unequal freq P[{}][{}] = {}, expected {}",
                    i, j, p[i][j], pi[j]
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Protein FASTA parser tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_protein_fasta_simple() {
        let input = ">Seq1\nACDEFG\n>Seq2\nACDEFG\n";
        let aln = read_protein_fasta(input).unwrap();
        assert_eq!(aln.ntaxa(), 2);
        assert_eq!(aln.n_sites, 6);
        assert_eq!(aln.taxon_name(0), "Seq1");
        assert_eq!(aln.residue(0, 0), AminoAcid::Ala);
        assert_eq!(aln.residue(0, 4), AminoAcid::Phe);
    }

    #[test]
    fn test_read_protein_fasta_multiline() {
        let input = ">Seq1\nACDE\nFGHI\n>Seq2\nKLMN\nPQRS\n";
        let aln = read_protein_fasta(input).unwrap();
        assert_eq!(aln.ntaxa(), 2);
        assert_eq!(aln.n_sites, 8);
    }

    #[test]
    fn test_read_protein_fasta_with_gaps() {
        let input = ">Seq1\nACD-FG\n>Seq2\nACDEFG\n";
        let aln = read_protein_fasta(input).unwrap();
        assert_eq!(aln.residue(0, 3), AminoAcid::Gap);
    }

    #[test]
    fn test_read_protein_fasta_empty() {
        assert!(read_protein_fasta("").is_err());
    }

    #[test]
    fn test_read_protein_fasta_no_header() {
        assert!(read_protein_fasta("ACDEFG").is_err());
    }

    #[test]
    fn test_read_protein_fasta_unequal_lengths() {
        let input = ">Seq1\nACDE\n>Seq2\nACDEFG\n";
        assert!(read_protein_fasta(input).is_err());
    }

    #[test]
    fn test_read_protein_fasta_header_with_description() {
        let input = ">Seq1 some protein\nACDE\n>Seq2 another protein\nACDE\n";
        let aln = read_protein_fasta(input).unwrap();
        assert_eq!(aln.taxon_name(0), "Seq1");
        assert_eq!(aln.taxon_name(1), "Seq2");
    }

    #[test]
    fn test_read_protein_fasta_all_amino_acids() {
        let input = ">Seq1\nARNDCQEGHILKMFPSTWYV\n>Seq2\nARNDCQEGHILKMFPSTWYV\n";
        let aln = read_protein_fasta(input).unwrap();
        assert_eq!(aln.n_sites, 20);
        for (i, &aa) in ALL_AMINO_ACIDS.iter().enumerate() {
            assert_eq!(aln.residue(0, i), aa, "failed at site {}", i);
        }
    }

    #[test]
    fn test_read_protein_fasta_lowercase() {
        let input = ">Seq1\nacdefg\n>Seq2\nacdefg\n";
        let aln = read_protein_fasta(input).unwrap();
        assert_eq!(aln.residue(0, 0), AminoAcid::Ala);
        assert_eq!(aln.residue(0, 1), AminoAcid::Cys);
    }

    // -----------------------------------------------------------------------
    // Protein likelihood tests
    // -----------------------------------------------------------------------

    fn make_protein_alignment(data: &[(&str, &str)]) -> ProteinAlignment {
        let seqs: Vec<ProteinSequence> = data
            .iter()
            .map(|(name, seq)| {
                let residues: Vec<AminoAcid> =
                    seq.chars().map(AminoAcid::from_char).collect();
                ProteinSequence::new(*name, residues)
            })
            .collect();
        ProteinAlignment::new(seqs).unwrap()
    }

    #[test]
    fn test_protein_likelihood_is_negative() {
        let aln = make_protein_alignment(&[
            ("A", "ACDE"),
            ("B", "ACDE"),
            ("C", "FGHI"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = PoissonModel::equal_frequencies();
        let lnl = protein_log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0, "log-likelihood should be negative, got {}", lnl);
    }

    #[test]
    fn test_protein_likelihood_identical_sequences() {
        let aln = make_protein_alignment(&[
            ("A", "ACDEFGHIKLMN"),
            ("B", "ACDEFGHIKLMN"),
            ("C", "ACDEFGHIKLMN"),
        ]);
        let tree = parse_newick("((A:0.001,B:0.001):0.001,C:0.001);").unwrap();
        let model = PoissonModel::equal_frequencies();
        let lnl = protein_log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
        // Per-site should be close to ln(1/20) = -2.996 for equal frequencies
        let per_site = lnl / 12.0;
        assert!(
            per_site > -4.0,
            "per-site lnl should be > -4.0 for identical seqs, got {}",
            per_site
        );
    }

    #[test]
    fn test_protein_likelihood_divergent_sequences() {
        let aln = make_protein_alignment(&[
            ("A", "AAAAAAAA"),
            ("B", "CCCCCCCC"),
            ("C", "DDDDDDDD"),
        ]);
        let tree = parse_newick("((A:0.01,B:0.01):0.01,C:0.01);").unwrap();
        let model = PoissonModel::equal_frequencies();
        let lnl = protein_log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < -30.0, "divergent sequences on short tree should have very low likelihood, got {}", lnl);
    }

    #[test]
    fn test_protein_likelihood_taxa_mismatch() {
        let aln = make_protein_alignment(&[
            ("X", "ACDE"),
            ("Y", "ACDE"),
            ("Z", "FGHI"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = PoissonModel::equal_frequencies();
        let result = protein_log_likelihood(&tree, &aln, &model);
        assert!(result.is_err());
    }

    #[test]
    fn test_protein_likelihood_with_wag() {
        let aln = make_protein_alignment(&[
            ("A", "ACDEFGHIKLMN"),
            ("B", "ACDEFGHIKLMN"),
            ("C", "FGHIKLMNPQRS"),
        ]);
        let tree = parse_newick("((A:0.2,B:0.2):0.1,C:0.3);").unwrap();
        let model = WagFreqModel::new();
        let lnl = protein_log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    #[test]
    fn test_protein_likelihood_four_taxa() {
        let aln = make_protein_alignment(&[
            ("A", "AACCDD"),
            ("B", "AACCDD"),
            ("C", "DDAACC"),
            ("D", "DDAACC"),
        ]);
        let tree =
            parse_newick("((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);").unwrap();
        let model = PoissonModel::equal_frequencies();
        let lnl = protein_log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    #[test]
    fn test_protein_likelihood_better_topology() {
        // Data has clear ((A,B),(C,D)) signal
        let aln = make_protein_alignment(&[
            ("A", "AAAAAACCCCCC"),
            ("B", "AAAAAACCCCCC"),
            ("C", "CCCCCCAAAAAA"),
            ("D", "CCCCCCAAAAAA"),
        ]);

        let good_tree =
            parse_newick("((A:0.05,B:0.05):0.2,(C:0.05,D:0.05):0.2);").unwrap();
        let bad_tree =
            parse_newick("((A:0.05,C:0.05):0.2,(B:0.05,D:0.05):0.2);").unwrap();

        let model = PoissonModel::equal_frequencies();
        let lnl_good = protein_log_likelihood(&good_tree, &aln, &model).unwrap();
        let lnl_bad = protein_log_likelihood(&bad_tree, &aln, &model).unwrap();

        assert!(
            lnl_good > lnl_bad,
            "correct topology ({}) should have higher likelihood than wrong ({})",
            lnl_good, lnl_bad
        );
    }

    #[test]
    fn test_protein_likelihood_with_gaps() {
        // Gaps should be treated as fully ambiguous (all states possible)
        let aln = make_protein_alignment(&[
            ("A", "A-DE"),
            ("B", "ACDE"),
            ("C", "FGHI"),
        ]);
        let tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.2);").unwrap();
        let model = PoissonModel::equal_frequencies();
        let lnl = protein_log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    #[test]
    fn test_protein_likelihood_tree_without_branch_lengths() {
        let aln = make_protein_alignment(&[
            ("A", "ACDE"),
            ("B", "ACDE"),
            ("C", "FGHI"),
        ]);
        // Tree without branch lengths: should use defaults
        let tree = parse_newick("((A,B),C);").unwrap();
        let model = PoissonModel::equal_frequencies();
        let lnl = protein_log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
    }

    #[test]
    fn test_protein_likelihood_poisson_vs_wag_differ() {
        // The Poisson (equal freq) and WAG-freq models should generally give
        // different likelihoods on the same data
        let aln = make_protein_alignment(&[
            ("A", "ACDEFGHIKLMN"),
            ("B", "ACDEFGHIRLMN"),
            ("C", "FGHIKLMNPQRS"),
        ]);
        let tree = parse_newick("((A:0.2,B:0.2):0.1,C:0.3);").unwrap();

        let poisson = PoissonModel::equal_frequencies();
        let wag = WagFreqModel::new();

        let lnl_poisson =
            protein_log_likelihood(&tree, &aln, &poisson).unwrap();
        let lnl_wag = protein_log_likelihood(&tree, &aln, &wag).unwrap();

        assert!(
            (lnl_poisson - lnl_wag).abs() > 0.01,
            "Poisson ({}) and WAG-freq ({}) should give different likelihoods",
            lnl_poisson, lnl_wag
        );
    }

    #[test]
    fn test_protein_likelihood_single_site_constant() {
        // Single constant site: all taxa have the same amino acid
        let aln = make_protein_alignment(&[("A", "A"), ("B", "A"), ("C", "A")]);
        let tree = parse_newick("((A:0.1,B:0.1):0.1,C:0.1);").unwrap();
        let model = PoissonModel::equal_frequencies();
        let lnl = protein_log_likelihood(&tree, &aln, &model).unwrap();
        assert!(lnl < 0.0);
        // Should be close to ln(1/20) = -2.996 but higher due to finite branch
        // lengths (most probability mass stays on the observed state)
        assert!(
            lnl > -5.0,
            "single constant site should have reasonable lnl, got {}",
            lnl
        );
    }

    #[test]
    fn test_protein_likelihood_empty_tree() {
        let aln = make_protein_alignment(&[("A", "ACDE")]);
        let tree = Tree::new(); // empty tree
        let model = PoissonModel::equal_frequencies();
        let result = protein_log_likelihood(&tree, &aln, &model);
        assert!(result.is_err());
    }
}

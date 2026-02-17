//! Core phylogenetic tree types.
//!
//! Provides representations for rooted and unrooted phylogenetic trees,
//! distance matrices, and DNA sequence alignments.

use std::fmt;

/// A nucleotide base, including IUPAC ambiguity codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Base {
    A,
    C,
    G,
    T,
    /// Two-fold ambiguity codes
    M, // A or C
    R, // A or G
    W, // A or T
    S, // C or G
    Y, // C or T
    K, // G or T
    /// Three-fold ambiguity codes
    V, // A, C, G
    H, // A, C, T
    D, // A, G, T
    B, // C, G, T
    /// Fully ambiguous / unknown / gap
    N,
    Gap,
}

impl Base {
    /// Parse a single character into a Base.
    pub fn from_char(c: char) -> Option<Base> {
        match c.to_ascii_uppercase() {
            'A' => Some(Base::A),
            'C' => Some(Base::C),
            'G' => Some(Base::G),
            'T' | 'U' => Some(Base::T),
            'M' => Some(Base::M),
            'R' => Some(Base::R),
            'W' => Some(Base::W),
            'S' => Some(Base::S),
            'Y' => Some(Base::Y),
            'K' => Some(Base::K),
            'V' => Some(Base::V),
            'H' => Some(Base::H),
            'D' => Some(Base::D),
            'B' => Some(Base::B),
            'N' | 'X' | '?' | 'O' => Some(Base::N),
            '-' => Some(Base::Gap),
            _ => None,
        }
    }

    /// Convert to a character.
    pub fn to_char(self) -> char {
        match self {
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
            Base::M => 'M',
            Base::R => 'R',
            Base::W => 'W',
            Base::S => 'S',
            Base::Y => 'Y',
            Base::K => 'K',
            Base::V => 'V',
            Base::H => 'H',
            Base::D => 'D',
            Base::B => 'B',
            Base::N => 'N',
            Base::Gap => '-',
        }
    }

    /// Returns the set of possible bases as a 4-element boolean array [A, C, G, T].
    pub fn base_set(self) -> [bool; 4] {
        match self {
            Base::A => [true, false, false, false],
            Base::C => [false, true, false, false],
            Base::G => [false, false, true, false],
            Base::T => [false, false, false, true],
            Base::M => [true, true, false, false],
            Base::R => [true, false, true, false],
            Base::W => [true, false, false, true],
            Base::S => [false, true, true, false],
            Base::Y => [false, true, false, true],
            Base::K => [false, false, true, true],
            Base::V => [true, true, true, false],
            Base::H => [true, true, false, true],
            Base::D => [true, false, true, true],
            Base::B => [false, true, true, true],
            Base::N | Base::Gap => [true, true, true, true],
        }
    }

    /// Returns the conditional likelihood vector [P(A), P(C), P(G), P(T)].
    /// For unambiguous bases, this is a 0/1 indicator.
    /// For ambiguous bases, multiple positions are 1.0.
    pub fn likelihood_vector(self) -> [f64; 4] {
        let bs = self.base_set();
        [
            if bs[0] { 1.0 } else { 0.0 },
            if bs[1] { 1.0 } else { 0.0 },
            if bs[2] { 1.0 } else { 0.0 },
            if bs[3] { 1.0 } else { 0.0 },
        ]
    }

    /// Returns true if this is an unambiguous base (A, C, G, or T).
    pub fn is_definite(self) -> bool {
        matches!(self, Base::A | Base::C | Base::G | Base::T)
    }

    /// Returns the index (0=A, 1=C, 2=G, 3=T) for definite bases.
    pub fn index(self) -> Option<usize> {
        match self {
            Base::A => Some(0),
            Base::C => Some(1),
            Base::G => Some(2),
            Base::T => Some(3),
            _ => None,
        }
    }
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

/// A named DNA sequence.
#[derive(Debug, Clone)]
pub struct Sequence {
    pub name: String,
    pub bases: Vec<Base>,
}

impl Sequence {
    pub fn new(name: impl Into<String>, bases: Vec<Base>) -> Self {
        Self {
            name: name.into(),
            bases,
        }
    }

    pub fn len(&self) -> usize {
        self.bases.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }
}

/// A multiple sequence alignment: a collection of sequences of equal length.
#[derive(Debug, Clone)]
pub struct Alignment {
    pub sequences: Vec<Sequence>,
}

impl Alignment {
    pub fn new(sequences: Vec<Sequence>) -> Result<Self, AlignmentError> {
        if sequences.is_empty() {
            return Err(AlignmentError::Empty);
        }
        let len = sequences[0].len();
        if sequences.iter().any(|s| s.len() != len) {
            return Err(AlignmentError::UnequalLengths);
        }
        Ok(Self { sequences })
    }

    /// Number of taxa (sequences).
    pub fn ntaxa(&self) -> usize {
        self.sequences.len()
    }

    /// Length of the alignment (number of sites).
    pub fn nsites(&self) -> usize {
        self.sequences.first().map_or(0, |s| s.len())
    }

    /// Get the name of taxon i.
    pub fn taxon_name(&self, i: usize) -> &str {
        &self.sequences[i].name
    }

    /// Get base at (taxon, site).
    pub fn base(&self, taxon: usize, site: usize) -> Base {
        self.sequences[taxon].bases[site]
    }
}

/// Errors in alignment construction.
#[derive(Debug, Clone)]
pub enum AlignmentError {
    Empty,
    UnequalLengths,
}

impl fmt::Display for AlignmentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AlignmentError::Empty => write!(f, "alignment contains no sequences"),
            AlignmentError::UnequalLengths => {
                write!(f, "sequences have different lengths")
            }
        }
    }
}

impl std::error::Error for AlignmentError {}

/// A symmetric distance matrix with taxon names.
#[derive(Debug, Clone)]
pub struct DistanceMatrix {
    /// Taxon names, length n.
    pub names: Vec<String>,
    /// Distances stored as a flat upper-triangular vector.
    /// Entry (i, j) where i < j is at index i*n - i*(i+1)/2 + j - i - 1.
    data: Vec<f64>,
    n: usize,
}

impl DistanceMatrix {
    /// Create a new zero-filled distance matrix for n taxa.
    pub fn new(names: Vec<String>) -> Self {
        let n = names.len();
        let size = n * (n - 1) / 2;
        Self {
            names,
            data: vec![0.0; size],
            n,
        }
    }

    /// Number of taxa.
    pub fn size(&self) -> usize {
        self.n
    }

    fn tri_index(&self, i: usize, j: usize) -> usize {
        debug_assert!(i < j && j < self.n);
        i * self.n - i * (i + 1) / 2 + j - i - 1
    }

    /// Get distance between taxa i and j.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        if i == j {
            0.0
        } else if i < j {
            self.data[self.tri_index(i, j)]
        } else {
            self.data[self.tri_index(j, i)]
        }
    }

    /// Set distance between taxa i and j.
    pub fn set(&mut self, i: usize, j: usize, value: f64) {
        if i == j {
            return;
        }
        let idx = if i < j {
            self.tri_index(i, j)
        } else {
            self.tri_index(j, i)
        };
        self.data[idx] = value;
    }

    /// Returns the full n x n matrix as a Vec<Vec<f64>>.
    pub fn to_full_matrix(&self) -> Vec<Vec<f64>> {
        let mut mat = vec![vec![0.0; self.n]; self.n];
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                let d = self.get(i, j);
                mat[i][j] = d;
                mat[j][i] = d;
            }
        }
        mat
    }
}

/// A node in a phylogenetic tree.
#[derive(Debug, Clone)]
pub struct TreeNode {
    /// Node identifier.
    pub id: usize,
    /// Taxon name (leaf nodes only).
    pub name: Option<String>,
    /// Branch length to parent.
    pub branch_length: Option<f64>,
    /// Children (empty for leaf nodes).
    pub children: Vec<usize>,
    /// Parent node index (None for root).
    pub parent: Option<usize>,
}

impl TreeNode {
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }
}

/// A phylogenetic tree stored as a vector of nodes.
///
/// Nodes are stored in a flat vector. The root is at index `root`.
/// Leaf nodes have `name` set and empty `children`.
/// Internal nodes have non-empty `children`.
#[derive(Debug, Clone)]
pub struct Tree {
    pub nodes: Vec<TreeNode>,
    pub root: usize,
    pub rooted: bool,
}

impl Tree {
    /// Create a new empty tree.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            root: 0,
            rooted: true,
        }
    }

    /// Add a node to the tree and return its index.
    pub fn add_node(
        &mut self,
        name: Option<String>,
        branch_length: Option<f64>,
        parent: Option<usize>,
    ) -> usize {
        let id = self.nodes.len();
        self.nodes.push(TreeNode {
            id,
            name,
            branch_length,
            children: Vec::new(),
            parent,
        });
        if let Some(p) = parent {
            self.nodes[p].children.push(id);
        }
        id
    }

    /// Number of nodes in the tree.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Returns all leaf nodes.
    pub fn leaves(&self) -> Vec<&TreeNode> {
        self.nodes.iter().filter(|n| n.is_leaf()).collect()
    }

    /// Number of leaf nodes (taxa).
    pub fn num_leaves(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_leaf()).count()
    }

    /// Get the names of all leaf nodes in order.
    pub fn leaf_names(&self) -> Vec<&str> {
        self.leaves()
            .iter()
            .filter_map(|n| n.name.as_deref())
            .collect()
    }

    /// Post-order traversal of node indices.
    pub fn postorder(&self) -> Vec<usize> {
        let mut result = Vec::new();
        let mut stack = vec![(self.root, false)];
        while let Some((node_id, visited)) = stack.pop() {
            if visited || self.nodes[node_id].children.is_empty() {
                result.push(node_id);
            } else {
                stack.push((node_id, true));
                for &child in self.nodes[node_id].children.iter().rev() {
                    stack.push((child, false));
                }
            }
        }
        result
    }

    /// Pre-order traversal of node indices.
    pub fn preorder(&self) -> Vec<usize> {
        let mut result = Vec::new();
        let mut stack = vec![self.root];
        while let Some(node_id) = stack.pop() {
            result.push(node_id);
            for &child in self.nodes[node_id].children.iter().rev() {
                stack.push(child);
            }
        }
        result
    }
}

impl Default for Tree {
    fn default() -> Self {
        Self::new()
    }
}

/// Base frequency parameters used by substitution models.
#[derive(Debug, Clone, Copy)]
pub struct BaseFrequencies {
    pub freq_a: f64,
    pub freq_c: f64,
    pub freq_g: f64,
    pub freq_t: f64,
}

impl BaseFrequencies {
    pub fn new(a: f64, c: f64, g: f64, t: f64) -> Self {
        Self {
            freq_a: a,
            freq_c: c,
            freq_g: g,
            freq_t: t,
        }
    }

    /// Equal base frequencies (0.25 each).
    pub fn equal() -> Self {
        Self::new(0.25, 0.25, 0.25, 0.25)
    }

    /// Estimate base frequencies from an alignment.
    pub fn from_alignment(alignment: &Alignment) -> Self {
        let mut counts = [0.0_f64; 4];
        let mut total = 0.0;
        for seq in &alignment.sequences {
            for &base in &seq.bases {
                if let Some(idx) = base.index() {
                    counts[idx] += 1.0;
                    total += 1.0;
                }
            }
        }
        if total == 0.0 {
            return Self::equal();
        }
        Self::new(
            counts[0] / total,
            counts[1] / total,
            counts[2] / total,
            counts[3] / total,
        )
    }

    /// Purine frequency (A + G).
    pub fn freq_r(&self) -> f64 {
        self.freq_a + self.freq_g
    }

    /// Pyrimidine frequency (C + T).
    pub fn freq_y(&self) -> f64 {
        self.freq_c + self.freq_t
    }
}

impl Default for BaseFrequencies {
    fn default() -> Self {
        Self::equal()
    }
}

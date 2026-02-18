//! Protein parsimony using genetic code step costs (protpars).
//!
//! Implements PHYLIP's `protpars` approach: the cost of changing between
//! amino acids is the minimum number of nucleotide substitutions required
//! by the standard genetic code.
//!
//! The 20x20 step matrix is computed from first principles by enumerating
//! all codon pairs for each amino acid pair and taking the minimum Hamming
//! distance. The Sankoff (1975) algorithm is then used to score trees under
//! this non-uniform cost matrix.
//!
//! # Example
//!
//! ```
//! use phylip_rs::parsimony::protein_parsimony::{
//!     protein_parsimony_score, protein_parsimony_search,
//!     genetic_code_step_matrix, AA_ORDER,
//! };
//! use phylip_rs::models::protein::{
//!     AminoAcid, ProteinSequence, ProteinAlignment,
//! };
//! use phylip_rs::tree::newick::parse_newick;
//!
//! let alignment = ProteinAlignment::new(vec![
//!     ProteinSequence::new("A", vec![AminoAcid::Ala, AminoAcid::Ala]),
//!     ProteinSequence::new("B", vec![AminoAcid::Ala, AminoAcid::Ala]),
//!     ProteinSequence::new("C", vec![AminoAcid::Ala, AminoAcid::Val]),
//! ]).unwrap();
//!
//! let tree = parse_newick("((A,B),C);").unwrap();
//! let score = protein_parsimony_score(&tree, &alignment);
//! assert_eq!(score, 1); // Ala->Val requires 1 nucleotide change
//! ```
//!
//! # References
//!
//! - Sankoff, D. (1975). Minimal mutation trees of sequences.
//!   *SIAM Journal on Applied Mathematics*, 28, 35-42.
//! - Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
//!   Chapter 8: The parsimony criterion.

use crate::models::protein::{AminoAcid, ProteinAlignment};
use crate::tree::types::{Tree, TreeNode};

/// The 20 standard amino acids in alphabetical order by one-letter code.
pub const AA_ORDER: [char; 20] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
];

/// Number of amino acid states.
const N_AA: usize = 20;

/// Map from amino acid one-letter code to index in `AA_ORDER` (0..19).
fn aa_char_to_index(c: char) -> Option<usize> {
    AA_ORDER.iter().position(|&x| x == c)
}

/// Map from `AminoAcid` enum to index in `AA_ORDER`.
fn amino_acid_to_aa_order_index(aa: AminoAcid) -> Option<usize> {
    aa_char_to_index(aa.to_char())
}

/// Standard genetic code: maps each amino acid one-letter code to its
/// RNA codons (using U instead of T).
///
/// Returns a list of (amino_acid_char, codons) pairs.
fn genetic_code_table() -> Vec<(char, Vec<[u8; 3]>)> {
    // RNA nucleotide encoding: U=0, C=1, A=2, G=3
    // Codons as [pos1, pos2, pos3] where each is 0..3
    //
    // The standard genetic code:
    // UUU=F UUC=F UUA=L UUG=L
    // UCU=S UCC=S UCA=S UCG=S
    // UAU=Y UAC=Y UAA=* UAG=*
    // UGU=C UGC=C UGA=* UGG=W
    // CUU=L CUC=L CUA=L CUG=L
    // CCU=P CCC=P CCA=P CCG=P
    // CAU=H CAC=H CAA=Q CAG=Q
    // CGU=R CGC=R CGA=R CGG=R
    // AUU=I AUC=I AUA=I AUG=M
    // ACU=T ACC=T ACA=T ACG=T
    // AAU=N AAC=N AAA=K AAG=K
    // AGU=S AGC=S AGA=R AGG=R
    // GUU=V GUC=V GUA=V GUG=V
    // GCU=A GCC=A GCA=A GCG=A
    // GAU=D GAC=D GAA=E GAG=E
    // GGU=G GGC=G GGA=G GGG=G

    let u = 0u8;
    let c = 1u8;
    let a = 2u8;
    let g = 3u8;

    vec![
        ('F', vec![[u, u, u], [u, u, c]]),
        (
            'L',
            vec![
                [u, u, a],
                [u, u, g],
                [c, u, u],
                [c, u, c],
                [c, u, a],
                [c, u, g],
            ],
        ),
        ('I', vec![[a, u, u], [a, u, c], [a, u, a]]),
        ('M', vec![[a, u, g]]),
        ('V', vec![[g, u, u], [g, u, c], [g, u, a], [g, u, g]]),
        (
            'S',
            vec![
                [u, c, u],
                [u, c, c],
                [u, c, a],
                [u, c, g],
                [a, g, u],
                [a, g, c],
            ],
        ),
        ('P', vec![[c, c, u], [c, c, c], [c, c, a], [c, c, g]]),
        ('T', vec![[a, c, u], [a, c, c], [a, c, a], [a, c, g]]),
        ('A', vec![[g, c, u], [g, c, c], [g, c, a], [g, c, g]]),
        ('Y', vec![[u, a, u], [u, a, c]]),
        ('H', vec![[c, a, u], [c, a, c]]),
        ('Q', vec![[c, a, a], [c, a, g]]),
        ('N', vec![[a, a, u], [a, a, c]]),
        ('K', vec![[a, a, a], [a, a, g]]),
        ('D', vec![[g, a, u], [g, a, c]]),
        ('E', vec![[g, a, a], [g, a, g]]),
        ('C', vec![[u, g, u], [u, g, c]]),
        ('W', vec![[u, g, g]]),
        (
            'R',
            vec![
                [c, g, u],
                [c, g, c],
                [c, g, a],
                [c, g, g],
                [a, g, a],
                [a, g, g],
            ],
        ),
        ('G', vec![[g, g, u], [g, g, c], [g, g, a], [g, g, g]]),
    ]
}

/// Compute the Hamming distance between two codons (number of
/// positions that differ).
fn codon_hamming(c1: &[u8; 3], c2: &[u8; 3]) -> u8 {
    let mut d = 0u8;
    for i in 0..3 {
        if c1[i] != c2[i] {
            d += 1;
        }
    }
    d
}

/// Build the 20x20 genetic code step matrix.
///
/// For each pair of amino acids (i, j), the cost is the minimum number
/// of nucleotide changes needed to convert any codon encoding amino acid i
/// to any codon encoding amino acid j.
///
/// The amino acids are ordered alphabetically by one-letter code as in
/// `AA_ORDER`: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y.
pub fn genetic_code_step_matrix() -> [[u8; 20]; 20] {
    let code = genetic_code_table();

    // Build a map from AA char to list of codons
    let mut codons_for_aa: [Vec<[u8; 3]>; 20] = Default::default();
    for (aa_char, codon_list) in &code {
        if let Some(idx) = aa_char_to_index(*aa_char) {
            codons_for_aa[idx] = codon_list.clone();
        }
    }

    let mut matrix = [[0u8; 20]; 20];

    for i in 0..N_AA {
        for j in 0..N_AA {
            if i == j {
                matrix[i][j] = 0;
            } else {
                let mut min_dist = u8::MAX;
                for c1 in &codons_for_aa[i] {
                    for c2 in &codons_for_aa[j] {
                        let d = codon_hamming(c1, c2);
                        if d < min_dist {
                            min_dist = d;
                        }
                    }
                }
                matrix[i][j] = min_dist;
            }
        }
    }

    matrix
}

/// Protein parsimony result.
#[derive(Debug, Clone)]
pub struct ProteinParsimonyResult {
    /// The best tree found.
    pub tree: Tree,
    /// The parsimony score (total steps across all sites).
    pub score: usize,
}

/// Score a tree using protein parsimony with genetic code costs.
///
/// Uses the Sankoff algorithm at each site:
/// - At each leaf, the cost for the observed amino acid is 0;
///   the cost for all other amino acids is infinity.
/// - At each internal node, for each possible amino acid state s:
///   cost[s] = sum over children c of min_over_states_r(cost_c[r] + step[s][r])
/// - The tree score at a site is the minimum cost at the root.
///
/// Gap and Unknown amino acids are treated as fully ambiguous (cost 0
/// for all states at the leaf).
pub fn protein_parsimony_score(tree: &Tree, alignment: &ProteinAlignment) -> usize {
    let step_matrix = genetic_code_step_matrix();
    protein_parsimony_score_with_matrix(tree, alignment, &step_matrix)
}

/// Score a tree using protein parsimony with a given step matrix.
///
/// This is the internal implementation that accepts a precomputed step matrix
/// to avoid recomputing it for every call during search.
fn protein_parsimony_score_with_matrix(
    tree: &Tree,
    alignment: &ProteinAlignment,
    step_matrix: &[[u8; 20]; 20],
) -> usize {
    let nsites = alignment.n_sites;
    let nnodes = tree.num_nodes();

    if nnodes == 0 {
        return 0;
    }

    // Build name -> alignment index map
    let name_to_idx: std::collections::HashMap<&str, usize> = alignment
        .sequences
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.as_str(), i))
        .collect();

    // Sankoff costs: cost[node][site][state] = minimum cost for subtree
    // rooted at this node, assuming this node has the given amino acid state.
    // We use u32 to avoid overflow (max cost per site is ~3 * nsites).
    const INF: u32 = u32::MAX / 2;

    // Allocate per-node cost arrays
    let mut costs: Vec<Vec<[u32; N_AA]>> = vec![vec![[INF; N_AA]; nsites]; nnodes];

    // Initialize leaves
    for node in &tree.nodes {
        if node.is_leaf() {
            if let Some(name) = &node.name {
                if let Some(&taxon_idx) = name_to_idx.get(name.as_str()) {
                    for site in 0..nsites {
                        let aa = alignment.residue(taxon_idx, site);
                        if let Some(idx) = amino_acid_to_aa_order_index(aa) {
                            // Known amino acid: cost 0 for this state, INF for others
                            costs[node.id][site] = [INF; N_AA];
                            costs[node.id][site][idx] = 0;
                        } else {
                            // Gap or Unknown: all states cost 0 (fully ambiguous)
                            costs[node.id][site] = [0; N_AA];
                        }
                    }
                }
            }
        }
    }

    // Postorder traversal: Sankoff algorithm
    let postorder = tree.postorder();
    for &node_id in &postorder {
        let children = &tree.nodes[node_id].children;
        if children.is_empty() {
            continue; // leaf: already initialized
        }

        for site in 0..nsites {
            // For each possible parent state s:
            // cost[s] = sum over children c of min_over_r(cost_c[r] + step[s][r])
            let mut parent_cost = [0u32; N_AA];

            for &child_id in children {
                let child_costs = &costs[child_id][site];

                for s in 0..N_AA {
                    let mut min_transition_cost = INF;
                    for r in 0..N_AA {
                        let tc = child_costs[r].saturating_add(step_matrix[s][r] as u32);
                        if tc < min_transition_cost {
                            min_transition_cost = tc;
                        }
                    }
                    parent_cost[s] = parent_cost[s].saturating_add(min_transition_cost);
                }
            }

            costs[node_id][site] = parent_cost;
        }
    }

    // Total score: sum of minimum costs at root over all sites
    let root = tree.root;
    let mut total_score: usize = 0;
    for site in 0..nsites {
        let min_cost = costs[root][site].iter().copied().min().unwrap_or(0);
        total_score += min_cost as usize;
    }

    total_score
}

/// Search for the most parsimonious tree for protein data.
///
/// Uses stepwise addition followed by SPR rearrangement, scoring
/// trees with the genetic code step matrix via the Sankoff algorithm.
///
/// # Arguments
///
/// * `alignment` - the protein alignment to analyze
/// * `seed` - optional random seed for taxon addition order
pub fn protein_parsimony_search(
    alignment: &ProteinAlignment,
    seed: Option<u64>,
) -> ProteinParsimonyResult {
    let ntaxa = alignment.ntaxa();
    assert!(ntaxa >= 3, "need at least 3 taxa for tree search");

    let step_matrix = genetic_code_step_matrix();

    // Taxon addition order
    let mut order: Vec<usize> = (0..ntaxa).collect();
    if let Some(s) = seed {
        let mut rng = SimpleRng::new(s);
        rng.shuffle(&mut order);
    }

    // Phase 1: Stepwise addition -- start with a 3-taxon star
    let mut tree = Tree::new();
    tree.rooted = false;
    let root = tree.add_node(None, None, None);
    tree.root = root;
    for &taxon_idx in &order[..3] {
        let name = alignment.taxon_name(taxon_idx).to_string();
        tree.add_node(Some(name), None, Some(root));
    }

    // Add remaining taxa one by one
    for i in 3..ntaxa {
        let taxon_name = alignment.taxon_name(order[i]);
        tree = pp_best_insertion(&tree, taxon_name, alignment, &step_matrix);
    }

    let mut score = protein_parsimony_score_with_matrix(&tree, alignment, &step_matrix);

    // Phase 2: SPR rearrangement
    loop {
        let (new_tree, new_score, improved) =
            pp_spr_round(&tree, alignment, &step_matrix, score);
        if !improved {
            break;
        }
        tree = new_tree;
        score = new_score;
    }

    ProteinParsimonyResult { tree, score }
}

// ============================================================
// Internal: tree manipulation (mirroring multistate.rs patterns)
// ============================================================

/// Simple LCG random number generator.
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1),
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    fn gen_range(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    fn shuffle<T>(&mut self, slice: &mut [T]) {
        for i in (1..slice.len()).rev() {
            let j = self.gen_range(i + 1);
            slice.swap(i, j);
        }
    }
}

/// Try inserting a new taxon at every possible edge; return the tree
/// with the best parsimony score.
fn pp_best_insertion(
    tree: &Tree,
    taxon_name: &str,
    alignment: &ProteinAlignment,
    step_matrix: &[[u8; 20]; 20],
) -> Tree {
    let mut best_tree = None;
    let mut best_score = usize::MAX;

    let edges: Vec<(usize, usize)> = tree
        .nodes
        .iter()
        .filter_map(|n| n.parent.map(|p| (p, n.id)))
        .collect();

    for (parent_id, child_id) in edges {
        let candidate = pp_insert_on_edge(tree, parent_id, child_id, taxon_name);
        let score = protein_parsimony_score_with_matrix(&candidate, alignment, step_matrix);
        if score < best_score {
            best_score = score;
            best_tree = Some(candidate);
        }
    }

    best_tree.unwrap_or_else(|| {
        let mut t = tree.clone();
        let root = t.root;
        t.add_node(Some(taxon_name.to_string()), None, Some(root));
        t
    })
}

/// Insert a new taxon on the edge between parent_id and child_id.
fn pp_insert_on_edge(
    tree: &Tree,
    parent_id: usize,
    child_id: usize,
    taxon_name: &str,
) -> Tree {
    let mut new_tree = Tree::new();
    new_tree.rooted = tree.rooted;

    let mut id_map = vec![0usize; tree.nodes.len()];

    for node in &tree.nodes {
        let new_id = new_tree.nodes.len();
        id_map[node.id] = new_id;
        new_tree.nodes.push(TreeNode {
            id: new_id,
            name: node.name.clone(),
            branch_length: node.branch_length,
            children: Vec::new(),
            parent: None,
        });
    }
    new_tree.root = id_map[tree.root];

    let new_internal_id = new_tree.nodes.len();
    new_tree.nodes.push(TreeNode {
        id: new_internal_id,
        name: None,
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    let new_leaf_id = new_tree.nodes.len();
    new_tree.nodes.push(TreeNode {
        id: new_leaf_id,
        name: Some(taxon_name.to_string()),
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    for node in &tree.nodes {
        let mapped_id = id_map[node.id];
        for &old_child in &node.children {
            if node.id == parent_id && old_child == child_id {
                new_tree.nodes[mapped_id].children.push(new_internal_id);
                new_tree.nodes[new_internal_id].parent = Some(mapped_id);

                let mapped_child = id_map[child_id];
                new_tree.nodes[new_internal_id].children.push(mapped_child);
                new_tree.nodes[mapped_child].parent = Some(new_internal_id);

                new_tree.nodes[new_internal_id].children.push(new_leaf_id);
                new_tree.nodes[new_leaf_id].parent = Some(new_internal_id);
            } else {
                let mapped_child = id_map[old_child];
                new_tree.nodes[mapped_id].children.push(mapped_child);
                new_tree.nodes[mapped_child].parent = Some(mapped_id);
            }
        }
    }

    new_tree
}

// --- SPR rearrangement for protein parsimony ---

/// Result of pruning a subtree.
struct PpPrunedTree {
    main_tree: Tree,
    subtree: Tree,
}

/// Collect all node IDs in the subtree rooted at `node_id`.
fn pp_collect_subtree_ids(tree: &Tree, node_id: usize) -> Vec<usize> {
    let mut ids = vec![node_id];
    let mut stack = vec![node_id];
    while let Some(nid) = stack.pop() {
        for &child in &tree.nodes[nid].children {
            ids.push(child);
            stack.push(child);
        }
    }
    ids
}

/// Prune a subtree rooted at `prune_id`.
fn pp_prune_subtree(tree: &Tree, prune_id: usize) -> Option<PpPrunedTree> {
    let parent_id = tree.nodes[prune_id].parent?;

    let subtree_ids: std::collections::HashSet<usize> =
        pp_collect_subtree_ids(tree, prune_id).into_iter().collect();

    // Build the subtree
    let mut subtree = Tree::new();
    subtree.rooted = false;
    let mut subtree_id_map = vec![usize::MAX; tree.nodes.len()];

    for &sid in &subtree_ids {
        let new_id = subtree.nodes.len();
        subtree_id_map[sid] = new_id;
        subtree.nodes.push(TreeNode {
            id: new_id,
            name: tree.nodes[sid].name.clone(),
            branch_length: tree.nodes[sid].branch_length,
            children: Vec::new(),
            parent: None,
        });
    }
    subtree.root = subtree_id_map[prune_id];

    for &sid in &subtree_ids {
        for &child in &tree.nodes[sid].children {
            if subtree_ids.contains(&child) {
                let mp = subtree_id_map[sid];
                let mc = subtree_id_map[child];
                subtree.nodes[mp].children.push(mc);
                subtree.nodes[mc].parent = Some(mp);
            }
        }
    }

    // Build the main tree (without subtree nodes)
    let mut main_tree = Tree::new();
    main_tree.rooted = tree.rooted;
    let mut main_id_map = vec![usize::MAX; tree.nodes.len()];

    let parent_remaining_children: Vec<usize> = tree.nodes[parent_id]
        .children
        .iter()
        .filter(|&&c| !subtree_ids.contains(&c))
        .copied()
        .collect();

    let collapse_parent =
        parent_remaining_children.len() == 1 && parent_id != tree.root;

    for node in &tree.nodes {
        if subtree_ids.contains(&node.id) {
            continue;
        }
        if collapse_parent && node.id == parent_id {
            continue;
        }
        let new_id = main_tree.nodes.len();
        main_id_map[node.id] = new_id;
        main_tree.nodes.push(TreeNode {
            id: new_id,
            name: node.name.clone(),
            branch_length: node.branch_length,
            children: Vec::new(),
            parent: None,
        });
    }

    main_tree.root = main_id_map[tree.root];

    for node in &tree.nodes {
        if subtree_ids.contains(&node.id) {
            continue;
        }
        if collapse_parent && node.id == parent_id {
            continue;
        }

        let mapped_id = main_id_map[node.id];

        for &child in &node.children {
            if subtree_ids.contains(&child) {
                continue;
            }

            if collapse_parent && child == parent_id {
                let surviving_child = parent_remaining_children[0];
                let mapped_child = main_id_map[surviving_child];
                main_tree.nodes[mapped_id].children.push(mapped_child);
                main_tree.nodes[mapped_child].parent = Some(mapped_id);
            } else {
                let mapped_child = main_id_map[child];
                main_tree.nodes[mapped_id].children.push(mapped_child);
                main_tree.nodes[mapped_child].parent = Some(mapped_id);
            }
        }
    }

    Some(PpPrunedTree {
        main_tree,
        subtree,
    })
}

/// Regraft a pruned subtree onto an edge in the main tree.
fn pp_regraft_subtree(
    pruned: &PpPrunedTree,
    graft_parent: usize,
    graft_child: usize,
) -> Tree {
    let main = &pruned.main_tree;
    let sub = &pruned.subtree;

    let mut new_tree = Tree::new();
    new_tree.rooted = main.rooted;

    let mut main_map = vec![0usize; main.nodes.len()];
    for node in &main.nodes {
        let new_id = new_tree.nodes.len();
        main_map[node.id] = new_id;
        new_tree.nodes.push(TreeNode {
            id: new_id,
            name: node.name.clone(),
            branch_length: node.branch_length,
            children: Vec::new(),
            parent: None,
        });
    }
    new_tree.root = main_map[main.root];

    let mut sub_map = vec![0usize; sub.nodes.len()];
    for node in &sub.nodes {
        let new_id = new_tree.nodes.len();
        sub_map[node.id] = new_id;
        new_tree.nodes.push(TreeNode {
            id: new_id,
            name: node.name.clone(),
            branch_length: node.branch_length,
            children: Vec::new(),
            parent: None,
        });
    }

    let graft_node_id = new_tree.nodes.len();
    new_tree.nodes.push(TreeNode {
        id: graft_node_id,
        name: None,
        branch_length: None,
        children: Vec::new(),
        parent: None,
    });

    for node in &main.nodes {
        let mapped_id = main_map[node.id];
        for &child in &node.children {
            if node.id == graft_parent && child == graft_child {
                new_tree.nodes[mapped_id].children.push(graft_node_id);
                new_tree.nodes[graft_node_id].parent = Some(mapped_id);

                let mapped_child = main_map[child];
                new_tree.nodes[graft_node_id].children.push(mapped_child);
                new_tree.nodes[mapped_child].parent = Some(graft_node_id);

                let subtree_root = sub_map[sub.root];
                new_tree.nodes[graft_node_id].children.push(subtree_root);
                new_tree.nodes[subtree_root].parent = Some(graft_node_id);
            } else {
                let mapped_child = main_map[child];
                new_tree.nodes[mapped_id].children.push(mapped_child);
                new_tree.nodes[mapped_child].parent = Some(mapped_id);
            }
        }
    }

    for node in &sub.nodes {
        let mapped_id = sub_map[node.id];
        for &child in &node.children {
            let mapped_child = sub_map[child];
            new_tree.nodes[mapped_id].children.push(mapped_child);
            new_tree.nodes[mapped_child].parent = Some(mapped_id);
        }
    }

    new_tree
}

/// One round of SPR rearrangement for protein parsimony.
fn pp_spr_round(
    tree: &Tree,
    alignment: &ProteinAlignment,
    step_matrix: &[[u8; 20]; 20],
    current_score: usize,
) -> (Tree, usize, bool) {
    let mut best_tree = tree.clone();
    let mut best_score = current_score;
    let mut improved = false;

    let nodes_to_try: Vec<usize> = tree
        .nodes
        .iter()
        .filter(|n| n.parent.is_some() && n.parent != Some(tree.root))
        .map(|n| n.id)
        .collect();

    for &prune_node in &nodes_to_try {
        if let Some(pruned) = pp_prune_subtree(tree, prune_node) {
            let edges: Vec<(usize, usize)> = pruned
                .main_tree
                .nodes
                .iter()
                .filter_map(|n| n.parent.map(|p| (p, n.id)))
                .collect();

            for (graft_parent, graft_child) in edges {
                let candidate = pp_regraft_subtree(&pruned, graft_parent, graft_child);
                let score =
                    protein_parsimony_score_with_matrix(&candidate, alignment, step_matrix);
                if score < best_score {
                    best_score = score;
                    best_tree = candidate;
                    improved = true;
                }
            }
        }
    }

    (best_tree, best_score, improved)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::protein::ProteinSequence;
    use crate::tree::newick::parse_newick;

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

    // --- Step matrix tests ---

    #[test]
    fn test_step_matrix_diagonal_zero() {
        let m = genetic_code_step_matrix();
        for i in 0..20 {
            assert_eq!(m[i][i], 0, "diagonal should be zero at ({}, {})", i, i);
        }
    }

    #[test]
    fn test_step_matrix_symmetric() {
        let m = genetic_code_step_matrix();
        for i in 0..20 {
            for j in 0..20 {
                assert_eq!(
                    m[i][j], m[j][i],
                    "step matrix should be symmetric at ({}, {})",
                    i, j
                );
            }
        }
    }

    #[test]
    fn test_step_matrix_all_positive_off_diagonal() {
        let m = genetic_code_step_matrix();
        for i in 0..20 {
            for j in 0..20 {
                if i != j {
                    assert!(
                        m[i][j] >= 1,
                        "off-diagonal should be >= 1 at ({}, {}), got {}",
                        i, j, m[i][j]
                    );
                }
            }
        }
    }

    #[test]
    fn test_step_matrix_max_is_3() {
        let m = genetic_code_step_matrix();
        let max = m.iter().flat_map(|row| row.iter()).copied().max().unwrap();
        assert_eq!(
            max, 3,
            "maximum step cost should be 3 (single codon = 3 positions)"
        );
    }

    #[test]
    fn test_step_matrix_f_to_l_is_1() {
        // F (Phe): UUU, UUC
        // L (Leu): UUA, UUG, CUU, CUC, CUA, CUG
        // UUU -> UUA = 1 change (3rd position)
        let m = genetic_code_step_matrix();
        let f_idx = aa_char_to_index('F').unwrap();
        let l_idx = aa_char_to_index('L').unwrap();
        assert_eq!(m[f_idx][l_idx], 1, "F -> L should require 1 nucleotide change");
    }

    #[test]
    fn test_step_matrix_a_to_v_is_1() {
        // A (Ala): GCU, GCC, GCA, GCG
        // V (Val): GUU, GUC, GUA, GUG
        // GCU -> GUU = 1 change (2nd position C->U)
        let m = genetic_code_step_matrix();
        let a_idx = aa_char_to_index('A').unwrap();
        let v_idx = aa_char_to_index('V').unwrap();
        assert_eq!(m[a_idx][v_idx], 1, "A -> V should require 1 nucleotide change");
    }

    #[test]
    fn test_step_matrix_d_to_e_is_1() {
        // D (Asp): GAU, GAC
        // E (Glu): GAA, GAG
        // GAU -> GAA = 1 change (3rd position)
        let m = genetic_code_step_matrix();
        let d_idx = aa_char_to_index('D').unwrap();
        let e_idx = aa_char_to_index('E').unwrap();
        assert_eq!(m[d_idx][e_idx], 1, "D -> E should require 1 nucleotide change");
    }

    #[test]
    fn test_step_matrix_known_2_change() {
        // F (Phe): UUU, UUC
        // W (Trp): UGG
        // UUU -> UGG = 2 changes (pos 2: U->G, pos 3: U->G)
        // UUC -> UGG = 2 changes (pos 2: U->G, pos 3: C->G)
        let m = genetic_code_step_matrix();
        let f_idx = aa_char_to_index('F').unwrap();
        let w_idx = aa_char_to_index('W').unwrap();
        assert_eq!(m[f_idx][w_idx], 2, "F -> W should require 2 nucleotide changes");
    }

    #[test]
    fn test_step_matrix_m_to_w_is_2() {
        // M (Met): AUG
        // W (Trp): UGG
        // AUG -> UGG = 2 changes (pos 1: A->U, pos 2: U->G)
        let m = genetic_code_step_matrix();
        let m_idx = aa_char_to_index('M').unwrap();
        let w_idx = aa_char_to_index('W').unwrap();
        assert_eq!(m[m_idx][w_idx], 2, "M -> W should require 2 nucleotide changes");
    }

    #[test]
    fn test_step_matrix_values_between_1_and_3() {
        let m = genetic_code_step_matrix();
        for i in 0..20 {
            for j in 0..20 {
                if i != j {
                    assert!(
                        m[i][j] >= 1 && m[i][j] <= 3,
                        "off-diagonal step cost should be in [1,3], got {} at ({}, {})",
                        m[i][j], i, j
                    );
                }
            }
        }
    }

    // --- Scoring tests ---

    #[test]
    fn test_identical_sequences_zero_score() {
        let aln = make_protein_alignment(&[
            ("A", "ACDE"),
            ("B", "ACDE"),
            ("C", "ACDE"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let score = protein_parsimony_score(&tree, &aln);
        assert_eq!(score, 0);
    }

    #[test]
    fn test_single_1_step_change() {
        // A (Ala, GCx) -> V (Val, GUx): 1 nucleotide change
        let aln = make_protein_alignment(&[
            ("A", "A"),
            ("B", "A"),
            ("C", "V"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let score = protein_parsimony_score(&tree, &aln);
        assert_eq!(score, 1, "A->V should cost 1 step");
    }

    #[test]
    fn test_single_2_step_change() {
        // F (Phe, UUU/UUC) -> W (Trp, UGG): 2 nucleotide changes
        let aln = make_protein_alignment(&[
            ("A", "F"),
            ("B", "F"),
            ("C", "W"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let score = protein_parsimony_score(&tree, &aln);
        assert_eq!(score, 2, "F->W should cost 2 steps");
    }

    #[test]
    fn test_single_3_step_change() {
        // For amino acids with minimum codon distance of 3, the Sankoff
        // algorithm on tree ((A,B),C) may find a cheaper path through an
        // intermediate state k: cost = min_k(m[i][k] + m[j][k]).
        // Verify the parsimony score matches this expected minimum.
        let m = genetic_code_step_matrix();
        let mut found = false;
        for i in 0..20 {
            for j in 0..20 {
                if m[i][j] == 3 {
                    // Compute expected Sankoff cost: min over all
                    // intermediate states k of m[i][k] + m[j][k].
                    let expected: u8 = (0..20)
                        .map(|k| m[i][k].saturating_add(m[j][k]))
                        .min()
                        .unwrap();
                    let aln = make_protein_alignment(&[
                        ("A", &AA_ORDER[i].to_string()),
                        ("B", &AA_ORDER[i].to_string()),
                        ("C", &AA_ORDER[j].to_string()),
                    ]);
                    let tree = parse_newick("((A,B),C);").unwrap();
                    let score = protein_parsimony_score(&tree, &aln);
                    assert_eq!(
                        score, expected as usize,
                        "Sankoff score for {} -> {} should be {} (direct cost {})",
                        AA_ORDER[i], AA_ORDER[j], expected, m[i][j]
                    );
                    // Score should be at most the direct cost
                    assert!(score <= m[i][j] as usize);
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }
        assert!(found, "should find at least one pair requiring 3 steps");
    }

    #[test]
    fn test_topology_matters_protein() {
        // A=A, B=A, C=V, D=V (A->V costs 1)
        // ((A,B),(C,D)): 1 step
        // ((A,C),(B,D)): 2 steps (both need to explain an A/V split)
        let aln = make_protein_alignment(&[
            ("A", "A"),
            ("B", "A"),
            ("C", "V"),
            ("D", "V"),
        ]);

        let t1 = parse_newick("((A,B),(C,D));").unwrap();
        let t2 = parse_newick("((A,C),(B,D));").unwrap();

        let s1 = protein_parsimony_score(&t1, &aln);
        let s2 = protein_parsimony_score(&t2, &aln);

        assert_eq!(s1, 1, "good topology should have score 1");
        assert!(s2 > s1, "bad topology should have worse score");
    }

    #[test]
    fn test_gap_handling_protein_parsimony() {
        // Gaps should be fully ambiguous (cost 0 for any state)
        let aln = make_protein_alignment(&[
            ("A", "A"),
            ("B", "A"),
            ("C", "-"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let score = protein_parsimony_score(&tree, &aln);
        assert_eq!(score, 0, "gap should be compatible with any amino acid");
    }

    #[test]
    fn test_multisite_protein_parsimony() {
        // Multiple sites; each site is scored independently
        let aln = make_protein_alignment(&[
            ("A", "AV"),
            ("B", "AV"),
            ("C", "VA"),
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let score = protein_parsimony_score(&tree, &aln);
        // Each site: A/A/V or V/V/A -> 1 step each = 2 total
        assert_eq!(score, 2);
    }

    #[test]
    fn test_protein_parsimony_vs_naive_count() {
        // The genetic code step matrix gives lower scores than simply
        // counting all amino acid differences as cost 1.
        let aln = make_protein_alignment(&[
            ("A", "A"),
            ("B", "A"),
            ("C", "V"), // A->V costs 1 nucleotide change
        ]);
        let tree = parse_newick("((A,B),C);").unwrap();
        let prot_score = protein_parsimony_score(&tree, &aln);

        // Naive approach (unordered, all changes cost 1): also 1 step
        // But for amino acids that are 2 or 3 apart, the protein score
        // will be higher than the naive score.
        assert_eq!(prot_score, 1, "for 1-step changes, scores should agree");

        // Now test a 2-step change
        let aln2 = make_protein_alignment(&[
            ("A", "F"),
            ("B", "F"),
            ("C", "W"), // F->W costs 2 nucleotide changes
        ]);
        let prot_score2 = protein_parsimony_score(&tree, &aln2);
        assert_eq!(
            prot_score2, 2,
            "protein parsimony should count nucleotide steps, not amino acid changes"
        );
    }

    // --- Search tests ---

    #[test]
    fn test_search_three_taxa() {
        let aln = make_protein_alignment(&[
            ("A", "ACDE"),
            ("B", "ACDE"),
            ("C", "FGHI"),
        ]);
        let result = protein_parsimony_search(&aln, Some(42));
        assert_eq!(result.tree.num_leaves(), 3);
        assert!(result.score > 0);
    }

    #[test]
    fn test_search_four_taxa() {
        let aln = make_protein_alignment(&[
            ("A", "AAAA"),
            ("B", "AAAA"),
            ("C", "VVVV"),
            ("D", "VVVV"),
        ]);
        let result = protein_parsimony_search(&aln, Some(42));
        assert_eq!(result.tree.num_leaves(), 4);
        // Best grouping: ((A,B),(C,D)) -> 1 step per site = 4
        assert_eq!(result.score, 4);
    }

    #[test]
    fn test_search_identical_sequences() {
        let aln = make_protein_alignment(&[
            ("A", "ACDE"),
            ("B", "ACDE"),
            ("C", "ACDE"),
            ("D", "ACDE"),
        ]);
        let result = protein_parsimony_search(&aln, None);
        assert_eq!(result.score, 0);
    }

    #[test]
    fn test_search_deterministic_with_seed() {
        let aln = make_protein_alignment(&[
            ("A", "ACDEFGHI"),
            ("B", "ACDEYGHI"),
            ("C", "FGHIACDE"),
            ("D", "FGHIVCDE"),
        ]);
        let r1 = protein_parsimony_search(&aln, Some(42));
        let r2 = protein_parsimony_search(&aln, Some(42));
        assert_eq!(r1.score, r2.score);
    }

    #[test]
    fn test_search_finds_good_tree() {
        // Clear signal: AB share states, CD share states
        let aln = make_protein_alignment(&[
            ("A", "AAAAAVVVV"),
            ("B", "AAAAAVVVV"),
            ("C", "VVVVVAAAA"),
            ("D", "VVVVVAAAA"),
        ]);
        let result = protein_parsimony_search(&aln, Some(7));
        // Optimal grouping gives 1 step per site = 9
        // (A->V is 1 nucleotide change, so each discordant site costs 1)
        assert_eq!(result.score, 9);
    }

    #[test]
    fn test_search_five_taxa() {
        let aln = make_protein_alignment(&[
            ("A", "AACCDD"),
            ("B", "AACCDD"),
            ("C", "DDAACC"),
            ("D", "DDAACC"),
            ("E", "CCDDAA"),
        ]);
        let result = protein_parsimony_search(&aln, Some(99));
        assert_eq!(result.tree.num_leaves(), 5);
        assert!(result.score > 0);
    }

    // --- Codon table completeness ---

    #[test]
    fn test_genetic_code_covers_all_61_codons() {
        let code = genetic_code_table();
        let total_codons: usize = code.iter().map(|(_, v)| v.len()).sum();
        // 64 total codons - 3 stop codons = 61 sense codons
        assert_eq!(total_codons, 61, "should have exactly 61 sense codons");
    }

    #[test]
    fn test_genetic_code_covers_all_20_amino_acids() {
        let code = genetic_code_table();
        let mut aa_chars: Vec<char> = code.iter().map(|(c, _)| *c).collect();
        aa_chars.sort();
        aa_chars.dedup();
        assert_eq!(
            aa_chars.len(),
            20,
            "should have exactly 20 amino acids"
        );
        // Check each one is in AA_ORDER
        for c in &aa_chars {
            assert!(
                AA_ORDER.contains(c),
                "amino acid '{}' not in AA_ORDER",
                c
            );
        }
    }
}

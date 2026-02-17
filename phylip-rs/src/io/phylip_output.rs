//! PHYLIP-style output formatting for analysis results.
//!
//! Produces human-readable reports resembling PHYLIP's classic `outfile`
//! output, including tree topology descriptions, branch length tables,
//! likelihood/parsimony scores, and model parameters.
//!
//! The `outtree` file (Newick format) is already handled by
//! [`crate::tree::newick::write_newick`].

use crate::likelihood::search::MLSearchResult;
use crate::parsimony::wagner::ParsimonyResult;
use crate::tree::types::{DistanceMatrix, Tree};

// ============================================================
// Constants
// ============================================================

/// Width of the separator line used in reports.
const SEPARATOR_WIDTH: usize = 70;

// ============================================================
// Public formatting functions
// ============================================================

/// Format a maximum likelihood analysis report.
///
/// Produces output resembling PHYLIP's DNAML `outfile`, including:
/// - Header with model name
/// - Tree topology description
/// - Branch length table
/// - Log-likelihood score
/// - Per-site log-likelihood summary statistics
///
/// # Arguments
/// * `result` - The ML search result containing the tree, log-likelihood,
///   and per-site log-likelihoods.
/// * `model_name` - The name of the substitution model used (e.g., "JC69", "F84").
pub fn format_ml_report(result: &MLSearchResult, model_name: &str) -> String {
    let mut out = String::new();

    // Header
    out.push_str(&format_header("Maximum Likelihood Analysis"));
    out.push('\n');
    out.push_str(&format!("Substitution model: {}\n", model_name));
    out.push_str(&format!(
        "Number of taxa:     {}\n",
        result.tree.num_leaves()
    ));
    out.push('\n');

    // Separator
    out.push_str(&separator());
    out.push('\n');

    // Tree topology description
    out.push_str("Tree topology:\n\n");
    out.push_str(&format_tree_description(&result.tree));
    out.push('\n');

    // Separator
    out.push_str(&separator());
    out.push('\n');

    // Branch length table
    out.push_str("Branch lengths:\n\n");
    out.push_str(&format_branch_table(&result.tree));
    out.push('\n');

    // Separator
    out.push_str(&separator());
    out.push('\n');

    // Log-likelihood
    out.push_str(&format!("Log-likelihood: {:.6}\n", result.lnl));
    out.push('\n');

    // Per-site log-likelihood summary
    if !result.site_lnl.is_empty() {
        out.push_str(&format_site_lnl_summary(&result.site_lnl));
    }

    out.push_str(&separator());
    out.push('\n');

    out
}

/// Format a parsimony analysis report.
///
/// Produces output resembling PHYLIP's DNAPARS `outfile`, including:
/// - Header
/// - Tree topology description
/// - Parsimony score
/// - Per-site step count summary
///
/// # Arguments
/// * `result` - The parsimony result containing the tree, score,
///   and per-site step counts.
pub fn format_parsimony_report(result: &ParsimonyResult) -> String {
    let mut out = String::new();

    // Header
    out.push_str(&format_header("Maximum Parsimony Analysis"));
    out.push('\n');
    out.push_str(&format!(
        "Number of taxa:     {}\n",
        result.tree.num_leaves()
    ));
    out.push('\n');

    // Separator
    out.push_str(&separator());
    out.push('\n');

    // Tree topology description
    out.push_str("Tree topology:\n\n");
    out.push_str(&format_tree_description(&result.tree));
    out.push('\n');

    // Separator
    out.push_str(&separator());
    out.push('\n');

    // Parsimony score
    out.push_str(&format!("Parsimony score (total steps): {}\n", result.score));
    out.push('\n');

    // Per-site step summary
    if !result.site_steps.is_empty() {
        out.push_str(&format_site_steps_summary(&result.site_steps));
    }

    out.push_str(&separator());
    out.push('\n');

    out
}

/// Format a distance analysis report.
///
/// Produces output resembling PHYLIP's NEIGHBOR/FITCH `outfile`, including:
/// - Header with method name
/// - Tree topology description
/// - Branch length table
/// - Input distance matrix
///
/// # Arguments
/// * `tree` - The inferred tree.
/// * `matrix` - The input distance matrix.
/// * `method` - The distance method name (e.g., "Neighbor-Joining", "UPGMA").
pub fn format_distance_report(tree: &Tree, matrix: &DistanceMatrix, method: &str) -> String {
    let mut out = String::new();

    // Header
    out.push_str(&format_header("Distance Matrix Analysis"));
    out.push('\n');
    out.push_str(&format!("Method:             {}\n", method));
    out.push_str(&format!("Number of taxa:     {}\n", matrix.size()));
    out.push('\n');

    // Separator
    out.push_str(&separator());
    out.push('\n');

    // Tree topology description
    out.push_str("Tree topology:\n\n");
    out.push_str(&format_tree_description(tree));
    out.push('\n');

    // Separator
    out.push_str(&separator());
    out.push('\n');

    // Branch length table
    out.push_str("Branch lengths:\n\n");
    out.push_str(&format_branch_table(tree));
    out.push('\n');

    // Separator
    out.push_str(&separator());
    out.push('\n');

    // Distance matrix
    out.push_str("Input distance matrix:\n\n");
    out.push_str(&format_distance_matrix(matrix));

    out.push_str(&separator());
    out.push('\n');

    out
}

/// Format a table of branch lengths for all edges in the tree.
///
/// Each row lists the branch endpoint (node name or internal node number),
/// its parent, and the branch length. The root node (which has no parent)
/// is omitted.
///
/// # Example output
/// ```text
/// From             To               Length
/// ----             --               ------
/// (internal 0)     A                0.100000
/// (internal 0)     B                0.200000
/// (root)           (internal 0)     0.050000
/// (root)           C                0.300000
/// ```
pub fn format_branch_table(tree: &Tree) -> String {
    let mut out = String::new();

    // Table header
    out.push_str(&format!(
        "{:<20} {:<20} {}\n",
        "From", "To", "Length"
    ));
    out.push_str(&format!(
        "{:<20} {:<20} {}\n",
        "----", "--", "------"
    ));

    // Collect branches in preorder for a logical display order
    let preorder = tree.preorder();
    for &node_id in &preorder {
        let node = &tree.nodes[node_id];
        if let Some(parent_id) = node.parent {
            let parent_label = node_label(tree, parent_id);
            let child_label = node_label(tree, node_id);
            let bl_str = match node.branch_length {
                Some(bl) => format!("{:.6}", bl),
                None => "n/a".to_string(),
            };
            out.push_str(&format!(
                "{:<20} {:<20} {}\n",
                parent_label, child_label, bl_str
            ));
        }
    }

    out
}

/// Format a human-readable description of the tree topology.
///
/// Produces a textual description that identifies clades and their
/// relationships, listing each internal node and its descendant taxa.
///
/// # Example output
/// ```text
/// Between        And            Length
/// -------        ---            ------
///    1           A              0.100000
///    1           B              0.100000
///    root        1              0.050000
///    root        C              0.200000
///
/// Tree ((A:0.1,B:0.1):0.05,C:0.2) has 3 taxa and 4 nodes.
///
/// Clade 1: { A, B }
/// ```
pub fn format_tree_description(tree: &Tree) -> String {
    let mut out = String::new();

    if tree.nodes.is_empty() {
        out.push_str("(empty tree)\n");
        return out;
    }

    // Numbered description table
    // Assign numbers to internal nodes for readability
    let internal_numbers = assign_internal_numbers(tree);

    out.push_str(&format!(
        "{:<14} {:<14} {}\n",
        "Between", "And", "Length"
    ));
    out.push_str(&format!(
        "{:<14} {:<14} {}\n",
        "-------", "---", "------"
    ));

    let preorder = tree.preorder();
    for &node_id in &preorder {
        let node = &tree.nodes[node_id];
        if let Some(parent_id) = node.parent {
            let parent_str = numbered_label(tree, parent_id, &internal_numbers);
            let child_str = numbered_label(tree, node_id, &internal_numbers);
            let bl_str = match node.branch_length {
                Some(bl) => format!("{:.6}", bl),
                None => "n/a".to_string(),
            };
            out.push_str(&format!(
                "{:<14} {:<14} {}\n",
                parent_str, child_str, bl_str
            ));
        }
    }

    // Summary line
    out.push_str(&format!(
        "\nTree has {} taxa and {} nodes.\n",
        tree.num_leaves(),
        tree.num_nodes()
    ));

    // List clades (internal nodes with their descendant leaves)
    let clades = collect_clades(tree);
    if !clades.is_empty() {
        out.push('\n');
        for (num, members) in &clades {
            out.push_str(&format!("Clade {}: {{ {} }}\n", num, members.join(", ")));
        }
    }

    out
}

// ============================================================
// Internal helpers
// ============================================================

/// Format a section header with a title line and decoration.
fn format_header(title: &str) -> String {
    let mut out = String::new();
    out.push_str(&separator());
    out.push('\n');
    // Center the title
    let padding = if SEPARATOR_WIDTH > title.len() {
        (SEPARATOR_WIDTH - title.len()) / 2
    } else {
        0
    };
    for _ in 0..padding {
        out.push(' ');
    }
    out.push_str(title);
    out.push('\n');
    out.push_str(&separator());
    out.push('\n');
    out
}

/// Produce a separator line of dashes.
fn separator() -> String {
    "-".repeat(SEPARATOR_WIDTH)
}

/// Get a display label for a node: its name if it has one, or a positional
/// label like "(internal 3)" or "(root)".
fn node_label(tree: &Tree, node_id: usize) -> String {
    let node = &tree.nodes[node_id];
    if let Some(ref name) = node.name {
        name.clone()
    } else if node_id == tree.root {
        "(root)".to_string()
    } else {
        format!("(internal {})", node_id)
    }
}

/// Assign sequential numbers to internal nodes for the topology description.
/// Returns a mapping from node_id to number. Leaves get no number.
fn assign_internal_numbers(tree: &Tree) -> Vec<Option<usize>> {
    let mut numbers = vec![None; tree.nodes.len()];
    let mut counter = 1usize;
    let preorder = tree.preorder();
    for &node_id in &preorder {
        if !tree.nodes[node_id].is_leaf() {
            if node_id == tree.root {
                // Root gets a special label
                numbers[node_id] = Some(0); // sentinel for "root"
            } else {
                numbers[node_id] = Some(counter);
                counter += 1;
            }
        }
    }
    numbers
}

/// Get a numbered label for the topology description.
fn numbered_label(tree: &Tree, node_id: usize, numbers: &[Option<usize>]) -> String {
    let node = &tree.nodes[node_id];
    if let Some(ref name) = node.name {
        name.clone()
    } else if node_id == tree.root {
        "root".to_string()
    } else {
        match numbers[node_id] {
            Some(n) => n.to_string(),
            None => format!("?{}", node_id),
        }
    }
}

/// Collect clades: for each internal non-root node, gather its descendant
/// leaf names. Returns a list of (clade_number, leaf_names).
fn collect_clades(tree: &Tree) -> Vec<(usize, Vec<String>)> {
    let internal_numbers = assign_internal_numbers(tree);
    let mut clades = Vec::new();

    let preorder = tree.preorder();
    for &node_id in &preorder {
        if tree.nodes[node_id].is_leaf() || node_id == tree.root {
            continue;
        }
        let leaves = descendant_leaves(tree, node_id);
        if let Some(Some(num)) = internal_numbers.get(node_id) {
            clades.push((*num, leaves));
        }
    }

    clades
}

/// Recursively collect the names of all leaf descendants of a node.
fn descendant_leaves(tree: &Tree, node_id: usize) -> Vec<String> {
    let node = &tree.nodes[node_id];
    if node.is_leaf() {
        return vec![node
            .name
            .clone()
            .unwrap_or_else(|| format!("unnamed_{}", node_id))];
    }
    let mut leaves = Vec::new();
    for &child in &node.children {
        leaves.extend(descendant_leaves(tree, child));
    }
    leaves
}

/// Format a summary of per-site log-likelihoods.
fn format_site_lnl_summary(site_lnl: &[f64]) -> String {
    let mut out = String::new();

    let nsites = site_lnl.len();
    let min = site_lnl
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let max = site_lnl
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let mean = site_lnl.iter().sum::<f64>() / nsites as f64;

    out.push_str("Per-site log-likelihood summary:\n");
    out.push_str(&format!("  Number of sites:  {}\n", nsites));
    out.push_str(&format!("  Mean:             {:.6}\n", mean));
    out.push_str(&format!("  Minimum:          {:.6}\n", min));
    out.push_str(&format!("  Maximum:          {:.6}\n", max));
    out.push('\n');

    out
}

/// Format a summary of per-site parsimony steps.
fn format_site_steps_summary(site_steps: &[usize]) -> String {
    let mut out = String::new();

    let nsites = site_steps.len();
    let min = site_steps.iter().copied().min().unwrap_or(0);
    let max = site_steps.iter().copied().max().unwrap_or(0);
    let total: usize = site_steps.iter().sum();
    let mean = if nsites > 0 {
        total as f64 / nsites as f64
    } else {
        0.0
    };

    // Count sites with zero steps (constant sites)
    let constant_sites = site_steps.iter().filter(|&&s| s == 0).count();
    let informative_sites = nsites - constant_sites;

    out.push_str("Per-site step summary:\n");
    out.push_str(&format!("  Number of sites:      {}\n", nsites));
    out.push_str(&format!("  Constant sites:       {}\n", constant_sites));
    out.push_str(&format!("  Informative sites:    {}\n", informative_sites));
    out.push_str(&format!("  Mean steps per site:  {:.4}\n", mean));
    out.push_str(&format!("  Minimum steps:        {}\n", min));
    out.push_str(&format!("  Maximum steps:        {}\n", max));
    out.push('\n');

    out
}

/// Format a distance matrix for inclusion in the report.
fn format_distance_matrix(matrix: &DistanceMatrix) -> String {
    let n = matrix.size();
    let mut out = String::new();

    // Header line: taxon count
    out.push_str(&format!("    {}\n", n));

    for i in 0..n {
        // Name padded to 10 characters
        let name = &matrix.names[i];
        let padded_name = if name.len() >= 10 {
            name[..10].to_string()
        } else {
            format!("{:<10}", name)
        };
        out.push_str(&padded_name);

        for j in 0..n {
            out.push_str(&format!("{:>10.6}", matrix.get(i, j)));
        }
        out.push('\n');
    }

    out
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::likelihood::search::MLSearchResult;
    use crate::parsimony::wagner::ParsimonyResult;
    use crate::tree::newick::parse_newick;
    use crate::tree::types::{DistanceMatrix, Tree};

    /// Helper to build a simple MLSearchResult for testing.
    fn make_ml_result() -> MLSearchResult {
        let tree = parse_newick("((A:0.1,B:0.2):0.05,C:0.3);").unwrap();
        MLSearchResult {
            tree,
            lnl: -15.123456,
            site_lnl: vec![-1.5, -2.0, -1.8, -1.9, -2.1, -1.7, -2.0, -2.12345],
        }
    }

    /// Helper to build a simple ParsimonyResult for testing.
    fn make_parsimony_result() -> ParsimonyResult {
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        ParsimonyResult {
            tree,
            score: 12,
            site_steps: vec![1, 0, 2, 1, 0, 1, 3, 1, 1, 2],
        }
    }

    /// Helper to build a simple tree for distance tests.
    fn make_distance_tree() -> Tree {
        parse_newick("((A:0.5,B:0.3):0.2,(C:0.4,D:0.6):0.1);").unwrap()
    }

    /// Helper to build a simple distance matrix.
    fn make_distance_matrix() -> DistanceMatrix {
        let mut matrix = DistanceMatrix::new(vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ]);
        matrix.set(0, 1, 0.5);
        matrix.set(0, 2, 1.2);
        matrix.set(0, 3, 1.4);
        matrix.set(1, 2, 1.0);
        matrix.set(1, 3, 1.3);
        matrix.set(2, 3, 0.6);
        matrix
    }

    // --- ML report tests ---

    #[test]
    fn test_ml_report_contains_model_name() {
        let result = make_ml_result();
        let report = format_ml_report(&result, "JC69");
        assert!(
            report.contains("JC69"),
            "ML report should contain model name"
        );
    }

    #[test]
    fn test_ml_report_contains_log_likelihood() {
        let result = make_ml_result();
        let report = format_ml_report(&result, "F84");
        assert!(
            report.contains("-15.123456"),
            "ML report should contain the log-likelihood value"
        );
    }

    #[test]
    fn test_ml_report_contains_header() {
        let result = make_ml_result();
        let report = format_ml_report(&result, "JC69");
        assert!(
            report.contains("Maximum Likelihood Analysis"),
            "ML report should contain analysis header"
        );
    }

    #[test]
    fn test_ml_report_contains_taxa_count() {
        let result = make_ml_result();
        let report = format_ml_report(&result, "JC69");
        assert!(
            report.contains("Number of taxa:     3"),
            "ML report should list taxa count"
        );
    }

    #[test]
    fn test_ml_report_contains_branch_table() {
        let result = make_ml_result();
        let report = format_ml_report(&result, "JC69");
        assert!(
            report.contains("Branch lengths:"),
            "ML report should include branch lengths section"
        );
        assert!(
            report.contains("0.100000"),
            "ML report should contain branch length 0.1"
        );
        assert!(
            report.contains("0.200000"),
            "ML report should contain branch length 0.2"
        );
    }

    #[test]
    fn test_ml_report_contains_site_summary() {
        let result = make_ml_result();
        let report = format_ml_report(&result, "JC69");
        assert!(
            report.contains("Per-site log-likelihood summary:"),
            "ML report should contain per-site summary"
        );
        assert!(
            report.contains("Number of sites:  8"),
            "ML report should report number of sites"
        );
    }

    #[test]
    fn test_ml_report_contains_leaf_names() {
        let result = make_ml_result();
        let report = format_ml_report(&result, "JC69");
        assert!(report.contains("A"), "ML report should contain leaf name A");
        assert!(report.contains("B"), "ML report should contain leaf name B");
        assert!(report.contains("C"), "ML report should contain leaf name C");
    }

    // --- Parsimony report tests ---

    #[test]
    fn test_parsimony_report_contains_score() {
        let result = make_parsimony_result();
        let report = format_parsimony_report(&result);
        assert!(
            report.contains("12"),
            "Parsimony report should contain the score"
        );
    }

    #[test]
    fn test_parsimony_report_contains_header() {
        let result = make_parsimony_result();
        let report = format_parsimony_report(&result);
        assert!(
            report.contains("Maximum Parsimony Analysis"),
            "Parsimony report should contain analysis header"
        );
    }

    #[test]
    fn test_parsimony_report_contains_step_summary() {
        let result = make_parsimony_result();
        let report = format_parsimony_report(&result);
        assert!(
            report.contains("Per-site step summary:"),
            "Parsimony report should contain step summary"
        );
        assert!(
            report.contains("Constant sites:"),
            "Parsimony report should list constant sites"
        );
        assert!(
            report.contains("Informative sites:"),
            "Parsimony report should list informative sites"
        );
    }

    #[test]
    fn test_parsimony_report_constant_sites_count() {
        let result = make_parsimony_result();
        let report = format_parsimony_report(&result);
        // site_steps has 2 zeros (constant sites)
        assert!(
            report.contains("Constant sites:       2"),
            "Parsimony report should correctly count constant sites"
        );
    }

    #[test]
    fn test_parsimony_report_contains_taxa() {
        let result = make_parsimony_result();
        let report = format_parsimony_report(&result);
        assert!(
            report.contains("Number of taxa:     4"),
            "Parsimony report should list taxa count"
        );
    }

    // --- Branch table tests ---

    #[test]
    fn test_branch_table_lists_all_branches() {
        let tree = parse_newick("((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.15);").unwrap();
        let table = format_branch_table(&tree);

        // Should have 6 branches (6 non-root nodes in a 7-node tree)
        let data_lines: Vec<&str> = table
            .lines()
            .skip(2) // skip header and separator lines
            .filter(|l| !l.is_empty())
            .collect();
        // The tree has 7 nodes: root + 2 internal + 4 leaves = 6 branches
        assert_eq!(
            data_lines.len(),
            6,
            "Branch table should list all 6 branches, got {}",
            data_lines.len()
        );
    }

    #[test]
    fn test_branch_table_contains_header() {
        let tree = parse_newick("(A:0.1,B:0.2);").unwrap();
        let table = format_branch_table(&tree);
        assert!(table.contains("From"));
        assert!(table.contains("To"));
        assert!(table.contains("Length"));
    }

    #[test]
    fn test_branch_table_no_branch_lengths() {
        let tree = parse_newick("(A,B,C);").unwrap();
        let table = format_branch_table(&tree);
        assert!(
            table.contains("n/a"),
            "Branch table should show 'n/a' for missing branch lengths"
        );
    }

    // --- Tree description tests ---

    #[test]
    fn test_tree_description_empty_tree() {
        let tree = Tree::new();
        let desc = format_tree_description(&tree);
        assert!(
            desc.contains("empty tree"),
            "Empty tree should be noted"
        );
    }

    #[test]
    fn test_tree_description_lists_taxa_count() {
        let tree = parse_newick("((A:0.1,B:0.2):0.05,C:0.3);").unwrap();
        let desc = format_tree_description(&tree);
        assert!(
            desc.contains("3 taxa"),
            "Description should mention 3 taxa"
        );
    }

    #[test]
    fn test_tree_description_lists_clades() {
        let tree = parse_newick("((A:0.1,B:0.2):0.05,C:0.3);").unwrap();
        let desc = format_tree_description(&tree);
        assert!(
            desc.contains("Clade"),
            "Description should list clades for internal nodes"
        );
        // The clade should contain A and B
        assert!(
            desc.contains("A") && desc.contains("B"),
            "Clade should list descendant taxa A and B"
        );
    }

    #[test]
    fn test_tree_description_has_root_label() {
        let tree = parse_newick("((A,B),C);").unwrap();
        let desc = format_tree_description(&tree);
        assert!(
            desc.contains("root"),
            "Description should mention the root"
        );
    }

    // --- Distance report tests ---

    #[test]
    fn test_distance_report_contains_method() {
        let tree = make_distance_tree();
        let matrix = make_distance_matrix();
        let report = format_distance_report(&tree, &matrix, "Neighbor-Joining");
        assert!(
            report.contains("Neighbor-Joining"),
            "Distance report should contain method name"
        );
    }

    #[test]
    fn test_distance_report_contains_matrix() {
        let tree = make_distance_tree();
        let matrix = make_distance_matrix();
        let report = format_distance_report(&tree, &matrix, "UPGMA");
        assert!(
            report.contains("Input distance matrix:"),
            "Distance report should include the distance matrix"
        );
        // Check a specific distance value
        assert!(
            report.contains("0.500000"),
            "Distance report matrix should include distance 0.5"
        );
    }

    #[test]
    fn test_distance_report_contains_taxa_count() {
        let tree = make_distance_tree();
        let matrix = make_distance_matrix();
        let report = format_distance_report(&tree, &matrix, "NJ");
        assert!(
            report.contains("Number of taxa:     4"),
            "Distance report should list taxa count"
        );
    }

    // --- Round-trip and integration tests ---

    #[test]
    fn test_ml_report_is_valid_utf8_string() {
        let result = make_ml_result();
        let report = format_ml_report(&result, "JC69");
        // Simply verifying it's a non-empty valid string
        assert!(!report.is_empty(), "ML report should be non-empty");
        // Check it's multi-line
        assert!(
            report.lines().count() > 10,
            "ML report should have multiple lines"
        );
    }

    #[test]
    fn test_parsimony_report_readable_format() {
        let result = make_parsimony_result();
        let report = format_parsimony_report(&result);
        // Ensure it has separator lines
        let dash_lines: Vec<&str> = report
            .lines()
            .filter(|l| l.chars().all(|c| c == '-') && l.len() > 10)
            .collect();
        assert!(
            dash_lines.len() >= 2,
            "Report should have separator lines for readability"
        );
    }

    #[test]
    fn test_format_branch_table_with_named_internal() {
        // Tree with labeled internal node
        let tree = parse_newick("((A:0.1,B:0.2)AB:0.05,C:0.3);").unwrap();
        let table = format_branch_table(&tree);
        assert!(
            table.contains("AB"),
            "Branch table should show internal node label 'AB'"
        );
    }

    #[test]
    fn test_site_lnl_summary_statistics() {
        let site_lnl = vec![-1.0, -2.0, -3.0, -4.0];
        let summary = format_site_lnl_summary(&site_lnl);
        // Mean should be -2.5
        assert!(
            summary.contains("-2.5"),
            "Summary should contain correct mean"
        );
        // Min should be -4.0
        assert!(
            summary.contains("-4.0"),
            "Summary should contain correct minimum"
        );
        // Max should be -1.0
        assert!(
            summary.contains("-1.0"),
            "Summary should contain correct maximum"
        );
    }

    #[test]
    fn test_distance_report_branch_table() {
        let tree = make_distance_tree();
        let matrix = make_distance_matrix();
        let report = format_distance_report(&tree, &matrix, "NJ");
        assert!(
            report.contains("Branch lengths:"),
            "Distance report should include branch lengths section"
        );
    }

    #[test]
    fn test_format_tree_description_four_taxa() {
        let tree = parse_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap();
        let desc = format_tree_description(&tree);
        assert!(desc.contains("4 taxa"));
        // Should have 2 clades (2 internal non-root nodes)
        let clade_count = desc.lines().filter(|l| l.starts_with("Clade")).count();
        assert_eq!(
            clade_count, 2,
            "Four-taxon tree should have 2 clades, got {}",
            clade_count
        );
    }

    #[test]
    fn test_ml_report_with_single_site() {
        let tree = parse_newick("(A:0.1,B:0.2);").unwrap();
        let result = MLSearchResult {
            tree,
            lnl: -1.386,
            site_lnl: vec![-1.386],
        };
        let report = format_ml_report(&result, "JC69");
        assert!(report.contains("Number of sites:  1"));
    }

    #[test]
    fn test_parsimony_report_zero_score() {
        let tree = parse_newick("((A,B),C);").unwrap();
        let result = ParsimonyResult {
            tree,
            score: 0,
            site_steps: vec![0, 0, 0, 0],
        };
        let report = format_parsimony_report(&result);
        assert!(report.contains("Parsimony score (total steps): 0"));
        assert!(report.contains("Constant sites:       4"));
    }
}

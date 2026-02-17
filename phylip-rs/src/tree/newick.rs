//! Newick format parser and writer.
//!
//! Implements reading and writing of the Newick (New Hampshire) tree format,
//! the standard text representation for phylogenetic trees. Supports branch
//! lengths, internal node labels, multifurcating trees, and various edge cases.
//!
//! Reference: <https://phylipweb.github.io/phylip/newicktree.html>

use std::fmt;

use super::types::Tree;

/// Errors that can occur during Newick parsing.
#[derive(Debug, Clone)]
pub enum NewickError {
    /// The input string is empty or contains only whitespace.
    EmptyInput,
    /// A closing parenthesis was found without a matching opening parenthesis.
    UnmatchedCloseParen(usize),
    /// An opening parenthesis was never closed.
    UnmatchedOpenParen(usize),
    /// The tree does not end with a semicolon.
    MissingSemicolon,
    /// An invalid branch length was encountered.
    InvalidBranchLength(String),
    /// Unexpected character encountered during parsing.
    UnexpectedCharacter(usize, char),
    /// General parse error with position and message.
    ParseError(usize, String),
}

impl fmt::Display for NewickError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NewickError::EmptyInput => write!(f, "empty input"),
            NewickError::UnmatchedCloseParen(pos) => {
                write!(f, "unmatched closing parenthesis at position {}", pos)
            }
            NewickError::UnmatchedOpenParen(pos) => {
                write!(f, "unmatched opening parenthesis at position {}", pos)
            }
            NewickError::MissingSemicolon => {
                write!(f, "tree does not end with semicolon")
            }
            NewickError::InvalidBranchLength(s) => {
                write!(f, "invalid branch length: '{}'", s)
            }
            NewickError::UnexpectedCharacter(pos, c) => {
                write!(f, "unexpected character '{}' at position {}", c, pos)
            }
            NewickError::ParseError(pos, msg) => {
                write!(f, "parse error at position {}: {}", pos, msg)
            }
        }
    }
}

impl std::error::Error for NewickError {}

/// Parse a Newick format string into a `Tree`.
///
/// The Newick format uses nested parentheses to represent phylogenetic tree
/// topology. Branch lengths are specified after a colon, and internal nodes
/// may optionally have labels.
///
/// # Examples
///
/// ```
/// use phylip_rs::tree::newick::parse_newick;
///
/// // Simple tree with branch lengths
/// let tree = parse_newick("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);").unwrap();
/// assert_eq!(tree.num_leaves(), 4);
///
/// // Tree without branch lengths
/// let tree = parse_newick("(A,B,(C,D));").unwrap();
/// assert_eq!(tree.num_leaves(), 4);
/// ```
pub fn parse_newick(input: &str) -> Result<Tree, NewickError> {
    let input = input.trim();
    if input.is_empty() {
        return Err(NewickError::EmptyInput);
    }

    // Must end with semicolon
    if !input.ends_with(';') {
        return Err(NewickError::MissingSemicolon);
    }

    // Strip trailing semicolon
    let input = input[..input.len() - 1].trim();
    if input.is_empty() {
        return Err(NewickError::EmptyInput);
    }

    let chars: Vec<char> = input.chars().collect();
    let mut tree = Tree::new();
    let mut pos = 0;

    let root = parse_subtree(&chars, &mut pos, &mut tree, None)?;
    tree.root = root;

    // Make sure we consumed all input
    // Skip any trailing whitespace
    while pos < chars.len() && chars[pos].is_whitespace() {
        pos += 1;
    }
    if pos < chars.len() {
        return Err(NewickError::UnexpectedCharacter(pos, chars[pos]));
    }

    Ok(tree)
}

/// Parse a subtree recursively.
///
/// A subtree is either:
/// - An internal node: `(child1,child2,...)label:length`
/// - A leaf node: `label:length`
fn parse_subtree(
    chars: &[char],
    pos: &mut usize,
    tree: &mut Tree,
    parent: Option<usize>,
) -> Result<usize, NewickError> {
    skip_whitespace(chars, pos);

    if *pos >= chars.len() {
        // Empty leaf node
        let node_id = tree.add_node(None, None, parent);
        return Ok(node_id);
    }

    if chars[*pos] == '(' {
        // Internal node
        let open_pos = *pos;
        *pos += 1; // consume '('

        // Create the internal node first (without label/branch_length yet)
        let node_id = tree.add_node(None, None, parent);

        // Parse children separated by commas
        loop {
            skip_whitespace(chars, pos);
            parse_subtree(chars, pos, tree, Some(node_id))?;
            skip_whitespace(chars, pos);

            if *pos >= chars.len() {
                return Err(NewickError::UnmatchedOpenParen(open_pos));
            }

            if chars[*pos] == ',' {
                *pos += 1; // consume ','
            } else if chars[*pos] == ')' {
                *pos += 1; // consume ')'
                break;
            } else {
                return Err(NewickError::UnexpectedCharacter(*pos, chars[*pos]));
            }
        }

        // Parse optional label and branch length for this internal node
        let (label, branch_length) = parse_label_and_length(chars, pos)?;
        tree.nodes[node_id].name = label;
        tree.nodes[node_id].branch_length = branch_length;

        Ok(node_id)
    } else {
        // Leaf node (or empty leaf if we hit , or ) immediately)
        let (label, branch_length) = parse_label_and_length(chars, pos)?;
        let node_id = tree.add_node(label, branch_length, parent);
        Ok(node_id)
    }
}

/// Parse an optional label followed by an optional `:branch_length`.
fn parse_label_and_length(
    chars: &[char],
    pos: &mut usize,
) -> Result<(Option<String>, Option<f64>), NewickError> {
    skip_whitespace(chars, pos);

    let label = parse_label(chars, pos);
    let branch_length = parse_branch_length(chars, pos)?;

    Ok((label, branch_length))
}

/// Parse a node label. Labels end at `:`, `,`, `)`, `(`, `;`, or end of input.
/// Underscores in labels are converted to spaces per the Newick convention.
fn parse_label(chars: &[char], pos: &mut usize) -> Option<String> {
    let mut label = String::new();

    while *pos < chars.len() {
        let c = chars[*pos];
        match c {
            ':' | ',' | ')' | '(' | ';' | '[' => break,
            '_' => {
                label.push(' ');
                *pos += 1;
            }
            c if c.is_whitespace() => {
                // Whitespace within a label context: stop if we already have
                // content (whitespace ends the label unless it's leading)
                if label.is_empty() {
                    *pos += 1;
                } else {
                    break;
                }
            }
            _ => {
                label.push(c);
                *pos += 1;
            }
        }
    }

    if label.is_empty() {
        None
    } else {
        Some(label)
    }
}

/// Parse an optional branch length (`:float_value`).
fn parse_branch_length(
    chars: &[char],
    pos: &mut usize,
) -> Result<Option<f64>, NewickError> {
    skip_whitespace(chars, pos);

    if *pos >= chars.len() || chars[*pos] != ':' {
        return Ok(None);
    }

    *pos += 1; // consume ':'
    skip_whitespace(chars, pos);

    let start = *pos;
    while *pos < chars.len() {
        match chars[*pos] {
            ',' | ')' | '(' | ';' | '[' => break,
            c if c.is_whitespace() => break,
            _ => *pos += 1,
        }
    }

    if *pos == start {
        // Colon with no number: treat as no branch length
        return Ok(None);
    }

    let num_str: String = chars[start..*pos].iter().collect();
    let num_str = num_str.trim();

    num_str
        .parse::<f64>()
        .map(Some)
        .map_err(|_| NewickError::InvalidBranchLength(num_str.to_string()))
}

/// Skip whitespace characters.
fn skip_whitespace(chars: &[char], pos: &mut usize) {
    while *pos < chars.len() && chars[*pos].is_whitespace() {
        *pos += 1;
    }
}

/// Write a `Tree` in Newick format.
///
/// Produces a Newick string ending with `;`. Branch lengths are included
/// when present, formatted with enough precision to avoid loss.
///
/// # Examples
///
/// ```
/// use phylip_rs::tree::newick::{parse_newick, write_newick};
///
/// let tree = parse_newick("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);").unwrap();
/// let newick = write_newick(&tree);
/// assert!(newick.ends_with(';'));
/// ```
pub fn write_newick(tree: &Tree) -> String {
    if tree.nodes.is_empty() {
        return ";".to_string();
    }

    let mut result = String::new();
    write_subtree(tree, tree.root, &mut result);
    result.push(';');
    result
}

/// Recursively write a subtree in Newick format.
fn write_subtree(tree: &Tree, node_id: usize, result: &mut String) {
    let node = &tree.nodes[node_id];

    if !node.children.is_empty() {
        result.push('(');
        for (i, &child_id) in node.children.iter().enumerate() {
            if i > 0 {
                result.push(',');
            }
            write_subtree(tree, child_id, result);
        }
        result.push(')');
    }

    // Write label (converting spaces back to underscores)
    if let Some(ref name) = node.name {
        let escaped = name.replace(' ', "_");
        result.push_str(&escaped);
    }

    // Write branch length
    if let Some(length) = node.branch_length {
        result.push(':');
        // Use enough precision to round-trip, but avoid unnecessary trailing zeros
        let formatted = format_branch_length(length);
        result.push_str(&formatted);
    }
}

/// Format a branch length value, avoiding unnecessary trailing zeros
/// but preserving enough precision for round-tripping.
fn format_branch_length(length: f64) -> String {
    // Use up to 10 decimal places, then strip trailing zeros
    let s = format!("{:.10}", length);
    let s = s.trim_end_matches('0');
    let s = s.trim_end_matches('.');
    // If the value was exactly zero, return "0"
    if s.is_empty() || s == "-" {
        "0".to_string()
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_tree() {
        let tree = parse_newick("(A,B);").unwrap();
        assert_eq!(tree.num_leaves(), 2);
        assert_eq!(tree.num_nodes(), 3);
        let names = tree.leaf_names();
        assert!(names.contains(&"A"));
        assert!(names.contains(&"B"));
    }

    #[test]
    fn test_tree_with_branch_lengths() {
        let tree = parse_newick("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);").unwrap();
        assert_eq!(tree.num_leaves(), 4);
        assert_eq!(tree.num_nodes(), 6); // 4 leaves + 1 internal + 1 root

        // Check that leaf names are present
        let names = tree.leaf_names();
        assert!(names.contains(&"A"));
        assert!(names.contains(&"B"));
        assert!(names.contains(&"C"));
        assert!(names.contains(&"D"));

        // Check branch lengths on leaves
        for node in tree.leaves() {
            assert!(node.branch_length.is_some());
        }
    }

    #[test]
    fn test_tree_without_branch_lengths() {
        let tree = parse_newick("(A,B,(C,D));").unwrap();
        assert_eq!(tree.num_leaves(), 4);
        for node in tree.leaves() {
            assert!(node.branch_length.is_none());
        }
    }

    #[test]
    fn test_fully_resolved_binary_tree() {
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        assert_eq!(tree.num_leaves(), 4);
        assert_eq!(tree.num_nodes(), 7); // 4 leaves + 2 internal + 1 root
    }

    #[test]
    fn test_internal_labels() {
        let tree = parse_newick("((A,B)AB,(C,D)CD)root;").unwrap();
        assert_eq!(tree.num_leaves(), 4);

        // Find the root node
        let root = &tree.nodes[tree.root];
        assert_eq!(root.name.as_deref(), Some("root"));

        // Find internal nodes with labels
        let internal_names: Vec<_> = tree
            .nodes
            .iter()
            .filter(|n| !n.is_leaf() && n.name.is_some())
            .map(|n| n.name.as_ref().unwrap().as_str())
            .collect();
        assert!(internal_names.contains(&"AB"));
        assert!(internal_names.contains(&"CD"));
        assert!(internal_names.contains(&"root"));
    }

    #[test]
    fn test_single_taxon() {
        let tree = parse_newick("A;").unwrap();
        assert_eq!(tree.num_leaves(), 1);
        assert_eq!(tree.num_nodes(), 1);
        assert_eq!(tree.nodes[0].name.as_deref(), Some("A"));
    }

    #[test]
    fn test_single_taxon_with_branch_length() {
        let tree = parse_newick("A:0.5;").unwrap();
        assert_eq!(tree.num_leaves(), 1);
        assert_eq!(tree.nodes[0].branch_length, Some(0.5));
    }

    #[test]
    fn test_whitespace_handling() {
        let tree = parse_newick("  ( A : 0.1 , B : 0.2 ) ; ").unwrap();
        assert_eq!(tree.num_leaves(), 2);
        let names = tree.leaf_names();
        assert!(names.contains(&"A"));
        assert!(names.contains(&"B"));
    }

    #[test]
    fn test_multiline() {
        let input = "(A:0.1,\n  B:0.2,\n  (C:0.3,\n   D:0.4):0.5);";
        let tree = parse_newick(input).unwrap();
        assert_eq!(tree.num_leaves(), 4);
    }

    #[test]
    fn test_underscore_to_space() {
        let tree = parse_newick("(Homo_sapiens,Pan_troglodytes);").unwrap();
        let names = tree.leaf_names();
        assert!(names.contains(&"Homo sapiens"));
        assert!(names.contains(&"Pan troglodytes"));
    }

    #[test]
    fn test_multifurcating() {
        let tree = parse_newick("(A,B,C,D,E);").unwrap();
        assert_eq!(tree.num_leaves(), 5);
        // Root should have 5 children
        assert_eq!(tree.nodes[tree.root].children.len(), 5);
    }

    #[test]
    fn test_empty_labels() {
        // Some leaves have empty labels (just commas)
        let tree = parse_newick("(,,(,));").unwrap();
        assert_eq!(tree.num_nodes(), 6); // 4 empty leaves + 1 internal + 1 root
    }

    #[test]
    fn test_deeply_nested() {
        let tree = parse_newick("((((A,B),C),D),E);").unwrap();
        assert_eq!(tree.num_leaves(), 5);
    }

    #[test]
    fn test_missing_semicolon() {
        let result = parse_newick("(A,B)");
        assert!(result.is_err());
        match result {
            Err(NewickError::MissingSemicolon) => {}
            other => panic!("expected MissingSemicolon, got {:?}", other),
        }
    }

    #[test]
    fn test_empty_input() {
        let result = parse_newick("");
        assert!(result.is_err());
        match result {
            Err(NewickError::EmptyInput) => {}
            other => panic!("expected EmptyInput, got {:?}", other),
        }
    }

    #[test]
    fn test_unmatched_open_paren() {
        let result = parse_newick("(A,B;");
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_branch_length() {
        let result = parse_newick("(A:abc,B:0.1);");
        assert!(result.is_err());
        match result {
            Err(NewickError::InvalidBranchLength(_)) => {}
            other => panic!("expected InvalidBranchLength, got {:?}", other),
        }
    }

    #[test]
    fn test_write_simple() {
        let tree = parse_newick("(A,B);").unwrap();
        let newick = write_newick(&tree);
        assert_eq!(newick, "(A,B);");
    }

    #[test]
    fn test_write_with_branch_lengths() {
        let tree = parse_newick("(A:0.1,B:0.2);").unwrap();
        let newick = write_newick(&tree);
        assert_eq!(newick, "(A:0.1,B:0.2);");
    }

    #[test]
    fn test_write_complex() {
        let tree = parse_newick("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);").unwrap();
        let newick = write_newick(&tree);
        assert_eq!(newick, "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);");
    }

    #[test]
    fn test_write_internal_labels() {
        let tree = parse_newick("((A,B)X,(C,D)Y)Z;").unwrap();
        let newick = write_newick(&tree);
        assert_eq!(newick, "((A,B)X,(C,D)Y)Z;");
    }

    #[test]
    fn test_roundtrip_simple() {
        let input = "(A,B,(C,D));";
        let tree = parse_newick(input).unwrap();
        let output = write_newick(&tree);
        let tree2 = parse_newick(&output).unwrap();
        assert_eq!(tree.num_leaves(), tree2.num_leaves());
        assert_eq!(tree.num_nodes(), tree2.num_nodes());
    }

    #[test]
    fn test_roundtrip_with_lengths() {
        let input = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);";
        let tree = parse_newick(input).unwrap();
        let output = write_newick(&tree);
        let tree2 = parse_newick(&output).unwrap();
        assert_eq!(tree.num_leaves(), tree2.num_leaves());
        assert_eq!(tree.num_nodes(), tree2.num_nodes());

        // Check that branch lengths survived the roundtrip
        for (n1, n2) in tree.nodes.iter().zip(tree2.nodes.iter()) {
            match (n1.branch_length, n2.branch_length) {
                (Some(l1), Some(l2)) => assert!((l1 - l2).abs() < 1e-10),
                (None, None) => {}
                _ => panic!("branch length mismatch"),
            }
        }
    }

    #[test]
    fn test_roundtrip_underscores() {
        let input = "(Homo_sapiens:0.1,Pan_troglodytes:0.2);";
        let tree = parse_newick(input).unwrap();
        // Names should have spaces internally
        let names = tree.leaf_names();
        assert!(names.contains(&"Homo sapiens"));
        // Writing should convert spaces back to underscores
        let output = write_newick(&tree);
        assert!(output.contains("Homo_sapiens"));
        assert!(output.contains("Pan_troglodytes"));
    }

    #[test]
    fn test_empty_tree() {
        let tree = Tree::new();
        let newick = write_newick(&tree);
        assert_eq!(newick, ";");
    }

    #[test]
    fn test_write_multifurcating() {
        let tree = parse_newick("(A,B,C,D);").unwrap();
        let newick = write_newick(&tree);
        assert_eq!(newick, "(A,B,C,D);");
    }

    #[test]
    fn test_branch_length_zero() {
        let tree = parse_newick("(A:0,B:0);").unwrap();
        assert_eq!(tree.nodes[1].branch_length, Some(0.0));
        assert_eq!(tree.nodes[2].branch_length, Some(0.0));
        let newick = write_newick(&tree);
        assert_eq!(newick, "(A:0,B:0);");
    }

    #[test]
    fn test_scientific_notation_branch_length() {
        let tree = parse_newick("(A:1e-5,B:2.5e-3);").unwrap();
        assert!((tree.nodes[1].branch_length.unwrap() - 1e-5).abs() < 1e-15);
        assert!((tree.nodes[2].branch_length.unwrap() - 2.5e-3).abs() < 1e-15);
    }

    #[test]
    fn test_negative_branch_length() {
        // Some tree formats include negative branch lengths (e.g., NJ trees)
        let tree = parse_newick("(A:-0.1,B:0.2);").unwrap();
        assert_eq!(tree.nodes[1].branch_length, Some(-0.1));
    }

    #[test]
    fn test_postorder_after_parse() {
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let order = tree.postorder();
        // All leaves should appear before their parent
        assert_eq!(order.len(), tree.num_nodes());
        // The root should be last
        assert_eq!(*order.last().unwrap(), tree.root);
    }

    #[test]
    fn test_preorder_after_parse() {
        let tree = parse_newick("((A,B),(C,D));").unwrap();
        let order = tree.preorder();
        assert_eq!(order.len(), tree.num_nodes());
        // The root should be first
        assert_eq!(order[0], tree.root);
    }

    #[test]
    fn test_large_tree() {
        // A tree with many taxa
        let taxa: Vec<String> = (0..100).map(|i| format!("T{}:0.{}", i, i + 1)).collect();
        let newick = format!("({});", taxa.join(","));
        let tree = parse_newick(&newick).unwrap();
        assert_eq!(tree.num_leaves(), 100);
    }

    #[test]
    fn test_display_newick_error() {
        let err = NewickError::EmptyInput;
        assert_eq!(format!("{}", err), "empty input");

        let err = NewickError::MissingSemicolon;
        assert_eq!(format!("{}", err), "tree does not end with semicolon");

        let err = NewickError::InvalidBranchLength("abc".to_string());
        assert_eq!(format!("{}", err), "invalid branch length: 'abc'");
    }
}

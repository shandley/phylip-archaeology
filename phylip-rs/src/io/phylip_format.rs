//! PHYLIP sequence format reader and writer.
//!
//! Supports reading PHYLIP-format sequence alignments in both sequential and
//! interleaved formats, with automatic format detection. Also supports reading
//! and writing PHYLIP distance matrix format (both full square and lower-triangular).

use std::fmt;

use crate::tree::types::{Alignment, Base, DistanceMatrix, Sequence};

/// Errors that can occur when parsing PHYLIP format files.
#[derive(Debug, Clone)]
pub enum PhylipFormatError {
    /// The input is empty.
    EmptyInput,
    /// The header line is malformed (expected "ntaxa nsites").
    InvalidHeader(String),
    /// A taxon name or sequence line could not be parsed.
    InvalidLine(usize, String),
    /// A character in a sequence is not a valid nucleotide.
    InvalidBase(usize, char),
    /// The number of taxa found doesn't match the header.
    TaxonCountMismatch { expected: usize, found: usize },
    /// Sequence lengths don't match the declared number of sites.
    SitesCountMismatch { expected: usize, found: usize, taxon: String },
    /// Sequences have different lengths.
    UnequalLengths,
    /// The distance matrix header is malformed.
    InvalidMatrixHeader(String),
    /// A distance value could not be parsed.
    InvalidDistance(usize, String),
    /// The distance matrix has the wrong number of entries.
    MatrixDimensionMismatch { row: usize, expected: usize, found: usize },
}

impl fmt::Display for PhylipFormatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PhylipFormatError::EmptyInput => write!(f, "empty input"),
            PhylipFormatError::InvalidHeader(s) => {
                write!(f, "invalid header line: '{}'", s)
            }
            PhylipFormatError::InvalidLine(line, msg) => {
                write!(f, "invalid line {}: {}", line, msg)
            }
            PhylipFormatError::InvalidBase(line, c) => {
                write!(f, "invalid base character '{}' at line {}", c, line)
            }
            PhylipFormatError::TaxonCountMismatch { expected, found } => {
                write!(
                    f,
                    "expected {} taxa but found {}",
                    expected, found
                )
            }
            PhylipFormatError::SitesCountMismatch {
                expected,
                found,
                taxon,
            } => {
                write!(
                    f,
                    "expected {} sites for taxon '{}' but found {}",
                    expected, taxon, found
                )
            }
            PhylipFormatError::UnequalLengths => {
                write!(f, "sequences have different lengths")
            }
            PhylipFormatError::InvalidMatrixHeader(s) => {
                write!(f, "invalid distance matrix header: '{}'", s)
            }
            PhylipFormatError::InvalidDistance(line, s) => {
                write!(f, "invalid distance value at line {}: '{}'", line, s)
            }
            PhylipFormatError::MatrixDimensionMismatch {
                row,
                expected,
                found,
            } => {
                write!(
                    f,
                    "row {} has {} distance values, expected {}",
                    row, found, expected
                )
            }
        }
    }
}

impl std::error::Error for PhylipFormatError {}

/// Parse a PHYLIP format sequence alignment.
///
/// Auto-detects whether the input is in sequential or interleaved format.
///
/// # PHYLIP sequential format
/// ```text
///  5  42
/// Turkey    AAGCTNGGGC ATTTCAGGGT GAGCCCGGGC AATACAGGGT AT
/// Salmo gairAAGCCTTGGC AGTGCAGGGT GAGCCGTGGC CGGGCACGGT AT
/// ...
/// ```
///
/// # PHYLIP interleaved format
/// ```text
///  5  42
/// Turkey    AAGCTNGGGC ATTTCAGGGT
/// Salmo gairAAGCCTTGGC AGTGCAGGGT
/// ...
///
/// GAGCCCGGGC AATACAGGGT AT
/// GAGCCGTGGC CGGGCACGGT AT
/// ...
/// ```
///
/// The first line contains the number of taxa and the number of sites.
/// Taxon names occupy the first 10 characters of their line.
pub fn read_phylip(input: &str) -> Result<Alignment, PhylipFormatError> {
    let input = input.trim();
    if input.is_empty() {
        return Err(PhylipFormatError::EmptyInput);
    }

    let lines: Vec<&str> = input.lines().collect();
    if lines.is_empty() {
        return Err(PhylipFormatError::EmptyInput);
    }

    // Parse header
    let header = lines[0].trim();
    let (ntaxa, nsites) = parse_header(header)?;

    if lines.len() < ntaxa + 1 {
        return Err(PhylipFormatError::TaxonCountMismatch {
            expected: ntaxa,
            found: lines.len() - 1,
        });
    }

    // Try to detect format: sequential vs interleaved.
    // In sequential format, each taxon gets enough lines to provide all sites.
    // In interleaved format, the first block has names + partial sequences,
    // then subsequent blocks have just sequence data.
    //
    // Heuristic: parse the first block (lines 1..ntaxa+1) as name+sequence.
    // If after parsing the first block, each sequence already has nsites bases,
    // it's sequential. Otherwise, it's interleaved.
    let first_block = parse_first_block(&lines[1..=ntaxa], ntaxa)?;

    // Count bases in the first block
    let first_block_len = count_bases_in_str(&first_block[0].1);

    if first_block_len >= nsites {
        // Sequential format (all bases on one line or already complete)
        parse_sequential(&lines, ntaxa, nsites)
    } else {
        // Interleaved format
        parse_interleaved(&lines, ntaxa, nsites)
    }
}

/// Parse the header line "ntaxa nsites".
fn parse_header(header: &str) -> Result<(usize, usize), PhylipFormatError> {
    let parts: Vec<&str> = header.split_whitespace().collect();
    if parts.len() < 2 {
        return Err(PhylipFormatError::InvalidHeader(header.to_string()));
    }

    let ntaxa = parts[0]
        .parse::<usize>()
        .map_err(|_| PhylipFormatError::InvalidHeader(header.to_string()))?;
    let nsites = parts[1]
        .parse::<usize>()
        .map_err(|_| PhylipFormatError::InvalidHeader(header.to_string()))?;

    Ok((ntaxa, nsites))
}

/// Parse the first block of lines into (name, sequence_text) pairs.
/// In PHYLIP format, the name occupies the first 10 characters.
fn parse_first_block(
    lines: &[&str],
    ntaxa: usize,
) -> Result<Vec<(String, String)>, PhylipFormatError> {
    let mut result = Vec::with_capacity(ntaxa);
    for (i, &line) in lines.iter().enumerate().take(ntaxa) {
        if line.is_empty() {
            return Err(PhylipFormatError::InvalidLine(
                i + 2,
                "empty line where taxon expected".to_string(),
            ));
        }

        let (name, seq_part) = split_name_sequence(line);
        if name.is_empty() {
            return Err(PhylipFormatError::InvalidLine(
                i + 2,
                "missing taxon name".to_string(),
            ));
        }
        result.push((name, seq_part));
    }
    Ok(result)
}

/// Split a PHYLIP line into name (first 10 chars, trimmed) and sequence portion.
fn split_name_sequence(line: &str) -> (String, String) {
    if line.len() >= 10 {
        let name = line[..10].trim().to_string();
        let seq = line[10..].to_string();
        (name, seq)
    } else {
        // Line is shorter than 10 chars: treat entire line as name
        (line.trim().to_string(), String::new())
    }
}

/// Count valid base characters in a string (ignoring whitespace and digits).
fn count_bases_in_str(s: &str) -> usize {
    s.chars()
        .filter(|c| !c.is_whitespace() && !c.is_ascii_digit())
        .filter(|c| Base::from_char(*c).is_some())
        .count()
}

/// Extract base characters from a string.
fn extract_bases(s: &str, line_num: usize) -> Result<Vec<Base>, PhylipFormatError> {
    let mut bases = Vec::new();
    for c in s.chars() {
        if c.is_whitespace() || c.is_ascii_digit() {
            continue;
        }
        match Base::from_char(c) {
            Some(b) => bases.push(b),
            None => return Err(PhylipFormatError::InvalidBase(line_num, c)),
        }
    }
    Ok(bases)
}

/// Parse PHYLIP sequential format.
fn parse_sequential(
    lines: &[&str],
    ntaxa: usize,
    nsites: usize,
) -> Result<Alignment, PhylipFormatError> {
    let mut sequences = Vec::with_capacity(ntaxa);

    let data_lines = &lines[1..];
    let mut line_idx = 0;

    for taxon_idx in 0..ntaxa {
        if line_idx >= data_lines.len() {
            return Err(PhylipFormatError::TaxonCountMismatch {
                expected: ntaxa,
                found: taxon_idx,
            });
        }

        let (name, seq_part) = split_name_sequence(data_lines[line_idx]);
        let mut bases = extract_bases(&seq_part, line_idx + 2)?;
        line_idx += 1;

        // Continue reading lines until we have nsites bases
        while bases.len() < nsites && line_idx < data_lines.len() {
            let line = data_lines[line_idx];
            // Continuation lines have no name prefix
            let more = extract_bases(line, line_idx + 2)?;
            if more.is_empty() {
                // Blank line separating taxa in sequential format
                line_idx += 1;
                continue;
            }
            bases.extend(more);
            line_idx += 1;
        }

        if bases.len() != nsites {
            return Err(PhylipFormatError::SitesCountMismatch {
                expected: nsites,
                found: bases.len(),
                taxon: name.clone(),
            });
        }

        sequences.push(Sequence::new(name, bases));
    }

    Alignment::new(sequences).map_err(|_| PhylipFormatError::UnequalLengths)
}

/// Parse PHYLIP interleaved format.
fn parse_interleaved(
    lines: &[&str],
    ntaxa: usize,
    nsites: usize,
) -> Result<Alignment, PhylipFormatError> {
    let data_lines = &lines[1..];

    // Parse first block: names + first chunk of sequence data
    let mut names: Vec<String> = Vec::with_capacity(ntaxa);
    let mut bases_per_taxon: Vec<Vec<Base>> = Vec::with_capacity(ntaxa);

    for i in 0..ntaxa {
        if i >= data_lines.len() {
            return Err(PhylipFormatError::TaxonCountMismatch {
                expected: ntaxa,
                found: i,
            });
        }
        let (name, seq_part) = split_name_sequence(data_lines[i]);
        let bases = extract_bases(&seq_part, i + 2)?;
        names.push(name);
        bases_per_taxon.push(bases);
    }

    // Parse subsequent blocks (sequence data only, no names)
    let mut line_idx = ntaxa;
    while bases_per_taxon[0].len() < nsites {
        // Skip blank lines between blocks
        while line_idx < data_lines.len() && data_lines[line_idx].trim().is_empty() {
            line_idx += 1;
        }

        if line_idx >= data_lines.len() {
            // Ran out of data lines before reaching nsites
            break;
        }

        // Read next block of ntaxa lines
        for taxon_idx in 0..ntaxa {
            if line_idx >= data_lines.len() {
                break;
            }
            let bases = extract_bases(data_lines[line_idx], line_idx + 2)?;
            bases_per_taxon[taxon_idx].extend(bases);
            line_idx += 1;
        }
    }

    // Validate lengths
    let mut sequences = Vec::with_capacity(ntaxa);
    for i in 0..ntaxa {
        if bases_per_taxon[i].len() != nsites {
            return Err(PhylipFormatError::SitesCountMismatch {
                expected: nsites,
                found: bases_per_taxon[i].len(),
                taxon: names[i].clone(),
            });
        }
        sequences.push(Sequence::new(names[i].clone(), bases_per_taxon[i].clone()));
    }

    Alignment::new(sequences).map_err(|_| PhylipFormatError::UnequalLengths)
}

/// Read a PHYLIP format distance matrix.
///
/// Supports both full square format and lower-triangular format.
///
/// # Full square format
/// ```text
/// 4
/// Alpha      0.0000  1.0000  2.0000  3.0000
/// Beta       1.0000  0.0000  2.0000  3.0000
/// Gamma      2.0000  2.0000  0.0000  3.0000
/// Delta      3.0000  3.0000  3.0000  0.0000
/// ```
///
/// # Lower-triangular format
/// ```text
/// 4
/// Alpha
/// Beta       1.0000
/// Gamma      2.0000  2.0000
/// Delta      3.0000  3.0000  3.0000
/// ```
pub fn read_distance_matrix(input: &str) -> Result<DistanceMatrix, PhylipFormatError> {
    let input = input.trim();
    if input.is_empty() {
        return Err(PhylipFormatError::EmptyInput);
    }

    let lines: Vec<&str> = input.lines().collect();
    if lines.is_empty() {
        return Err(PhylipFormatError::EmptyInput);
    }

    // First line: number of taxa
    let n = lines[0]
        .trim()
        .parse::<usize>()
        .map_err(|_| PhylipFormatError::InvalidMatrixHeader(lines[0].to_string()))?;

    if n == 0 {
        return Err(PhylipFormatError::InvalidMatrixHeader(
            "zero taxa".to_string(),
        ));
    }

    // Gather all data lines (skipping blank lines)
    let data_lines: Vec<&str> = lines[1..]
        .iter()
        .filter(|l| !l.trim().is_empty())
        .copied()
        .collect();

    // We may have multi-line entries for each taxon. Collect all content per taxon.
    // Strategy: parse greedily. The first token of each taxon's first line is the name.
    // Then distance values follow. A taxon might span multiple lines.

    // First, detect format: try parsing the first data line to see how many numbers it has.
    // In full square format, each taxon line has n distances (including self).
    // In lower-triangular, row i has i distances (row 0 has 0 distances, row 1 has 1, etc.).

    // Parse all taxon data - each taxon starts with a name (non-numeric first token)
    let mut names: Vec<String> = Vec::with_capacity(n);
    let mut all_distances: Vec<Vec<f64>> = Vec::with_capacity(n);

    let mut line_idx = 0;
    for taxon_idx in 0..n {
        if line_idx >= data_lines.len() {
            return Err(PhylipFormatError::TaxonCountMismatch {
                expected: n,
                found: taxon_idx,
            });
        }

        // Extract name from the first 10 characters (PHYLIP convention)
        let first_line = data_lines[line_idx];
        let (name, rest) = split_name_sequence(first_line);
        names.push(name);

        // Parse distances from the rest of this line
        let mut distances: Vec<f64> = Vec::new();
        parse_distances_from_str(&rest, line_idx + 2, &mut distances)?;
        line_idx += 1;

        // For continuation lines: keep reading until we have enough distances
        // We don't know yet if it's full or lower-tri, so read what's available
        // on this line. We'll handle continuation lines carefully.
        //
        // We need to be careful: a continuation line starts with distances
        // (no name prefix). We detect this by checking if the line's first
        // non-whitespace token parses as a number.
        while line_idx < data_lines.len() {
            let next_line = data_lines[line_idx];
            let first_token = next_line.trim().split_whitespace().next();
            if let Some(token) = first_token {
                if token.parse::<f64>().is_ok() {
                    // This is a continuation line (starts with a number)
                    parse_distances_from_str(next_line, line_idx + 2, &mut distances)?;
                    line_idx += 1;
                } else {
                    // Next taxon name
                    break;
                }
            } else {
                break;
            }
        }

        all_distances.push(distances);
    }

    // Detect format based on the number of distances in each row
    // Full square: each row has n distances
    // Lower-triangular: row 0 has 0, row 1 has 1, ..., row n-1 has n-1
    let is_full_square = if n > 0 && !all_distances.is_empty() {
        all_distances[0].len() == n
            || (n > 1 && all_distances.last().map_or(false, |d| d.len() == n))
    } else {
        false
    };

    let mut matrix = DistanceMatrix::new(names.clone());

    if is_full_square {
        // Full square format
        for i in 0..n {
            if all_distances[i].len() != n {
                return Err(PhylipFormatError::MatrixDimensionMismatch {
                    row: i,
                    expected: n,
                    found: all_distances[i].len(),
                });
            }
            for j in 0..n {
                if i != j {
                    matrix.set(i, j, all_distances[i][j]);
                }
            }
        }
    } else {
        // Lower-triangular format
        for i in 0..n {
            if all_distances[i].len() != i {
                return Err(PhylipFormatError::MatrixDimensionMismatch {
                    row: i,
                    expected: i,
                    found: all_distances[i].len(),
                });
            }
            for j in 0..i {
                matrix.set(i, j, all_distances[i][j]);
            }
        }
    }

    Ok(matrix)
}

/// Parse distance values from a string, appending to the provided vector.
fn parse_distances_from_str(
    s: &str,
    line_num: usize,
    distances: &mut Vec<f64>,
) -> Result<(), PhylipFormatError> {
    for token in s.split_whitespace() {
        match token.parse::<f64>() {
            Ok(d) => distances.push(d),
            Err(_) => {
                return Err(PhylipFormatError::InvalidDistance(
                    line_num,
                    token.to_string(),
                ))
            }
        }
    }
    Ok(())
}

/// Write a distance matrix in PHYLIP full square format.
///
/// Output format:
/// ```text
/// 4
/// Alpha      0.000000  1.000000  2.000000  3.000000
/// Beta       1.000000  0.000000  2.000000  3.000000
/// ...
/// ```
pub fn write_distance_matrix(matrix: &DistanceMatrix) -> String {
    let n = matrix.size();
    let mut result = String::new();

    // Header: number of taxa
    result.push_str(&format!("    {}\n", n));

    for i in 0..n {
        // Name padded to 10 characters
        let name = &matrix.names[i];
        let padded_name = if name.len() >= 10 {
            name[..10].to_string()
        } else {
            format!("{:<10}", name)
        };
        result.push_str(&padded_name);

        for j in 0..n {
            result.push_str(&format!("{:>10.6}", matrix.get(i, j)));
        }
        result.push('\n');
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- Sequence format tests ----

    #[test]
    fn test_read_sequential_simple() {
        let input = "\
 4  15
Alpha     AAGCTNGGGCATTTC
Beta      AAGCCTTGGCAGTGC
Gamma     TATTTAAAACCCGAG
Delta     AAGAATGGCCTTCAG";
        let alignment = read_phylip(input).unwrap();
        assert_eq!(alignment.ntaxa(), 4);
        assert_eq!(alignment.nsites(), 15);
        assert_eq!(alignment.taxon_name(0), "Alpha");
        assert_eq!(alignment.taxon_name(1), "Beta");
        assert_eq!(alignment.taxon_name(2), "Gamma");
        assert_eq!(alignment.taxon_name(3), "Delta");
    }

    #[test]
    fn test_read_sequential_with_spaces_in_sequence() {
        // PHYLIP allows spaces within sequence data (for readability)
        let input = "\
 3  20
Alpha     AAGCT NGGGG CATTT CAGGC
Beta      AAGCC TTGGC AGTGC AGGCC
Gamma     TATTT AAAAC CCGAG CAGGA";
        let alignment = read_phylip(input).unwrap();
        assert_eq!(alignment.ntaxa(), 3);
        assert_eq!(alignment.nsites(), 20);
    }

    #[test]
    fn test_read_interleaved() {
        let input = "\
 4  30
Alpha     AAGCTNGGGC ATTTCAGGGT
Beta      AAGCCTTGGC AGTGCAGGGT
Gamma     TATTTAAAAC CCGAGGGGGG
Delta     AAGAATGGCC TTCAGAGGGT

GAGCCCGGGC
GAGCCGTGGC
AAATTTCCCC
GAGCCCGGGC";
        let alignment = read_phylip(input).unwrap();
        assert_eq!(alignment.ntaxa(), 4);
        assert_eq!(alignment.nsites(), 30);
        assert_eq!(alignment.taxon_name(0), "Alpha");
    }

    #[test]
    fn test_read_interleaved_multiple_blocks() {
        let input = "\
 2  30
Taxon1    AAGCTNGGGC
Taxon2    AAGCCTTGGC

ATTTCAGGGT
AGTGCAGGGT

GAGCCCGGGC
GAGCCGTGGC";
        let alignment = read_phylip(input).unwrap();
        assert_eq!(alignment.ntaxa(), 2);
        assert_eq!(alignment.nsites(), 30);
    }

    #[test]
    fn test_read_phylip_empty() {
        let result = read_phylip("");
        assert!(result.is_err());
        match result {
            Err(PhylipFormatError::EmptyInput) => {}
            other => panic!("expected EmptyInput, got {:?}", other),
        }
    }

    #[test]
    fn test_read_phylip_invalid_header() {
        let result = read_phylip("not a valid header\n");
        assert!(result.is_err());
    }

    #[test]
    fn test_read_phylip_base_validation() {
        let alignment = read_phylip("\
 1  5
Taxon1    ACGTN").unwrap();
        assert_eq!(alignment.nsites(), 5);
        assert_eq!(alignment.base(0, 0), Base::A);
        assert_eq!(alignment.base(0, 4), Base::N);
    }

    #[test]
    fn test_read_phylip_gap_character() {
        let alignment = read_phylip("\
 1  5
Taxon1    ACG-N").unwrap();
        assert_eq!(alignment.base(0, 3), Base::Gap);
    }

    #[test]
    fn test_read_phylip_ten_char_name() {
        // PHYLIP names are exactly 10 characters
        let input = "\
 2  5
1234567890ACGTN
ABCDEFGHIJACGTN";
        let alignment = read_phylip(input).unwrap();
        assert_eq!(alignment.taxon_name(0), "1234567890");
        assert_eq!(alignment.taxon_name(1), "ABCDEFGHIJ");
    }

    #[test]
    fn test_classic_phylip_example() {
        // Based on the classic PHYLIP format example
        let input = "\
  5    42
Turkey    AAGCTNGGGC ATTTCAGGGT GAGCCCGGGC AATACAGGGT AT
Salmo gairAAGCCTTGGC AGTGCAGGGT GAGCCGTGGC CGGGCACGGT AT
H. SapiensACCGGTTGGC CGTTCAGGGT ACAGGTTGGC CGTTCAGGGT AA
Chimp     AAACCCTTGC CGTTACGCTT AAACCGAGGC CGGGACACTC AT
Gorilla   AAACCCTTGC CGGTACGCTT AAACCATTGC CGGTACGCTT AA";
        let alignment = read_phylip(input).unwrap();
        assert_eq!(alignment.ntaxa(), 5);
        assert_eq!(alignment.nsites(), 42);
        assert_eq!(alignment.taxon_name(0), "Turkey");
        assert_eq!(alignment.taxon_name(1), "Salmo gair");
        assert_eq!(alignment.taxon_name(2), "H. Sapiens");
    }

    // ---- Distance matrix tests ----

    #[test]
    fn test_read_distance_matrix_full_square() {
        let input = "\
    4
Alpha     0.0000  1.0000  2.0000  3.0000
Beta      1.0000  0.0000  2.5000  3.5000
Gamma     2.0000  2.5000  0.0000  4.0000
Delta     3.0000  3.5000  4.0000  0.0000";
        let matrix = read_distance_matrix(input).unwrap();
        assert_eq!(matrix.size(), 4);
        assert_eq!(matrix.names[0], "Alpha");
        assert_eq!(matrix.names[3], "Delta");
        assert!((matrix.get(0, 1) - 1.0).abs() < 1e-6);
        assert!((matrix.get(1, 2) - 2.5).abs() < 1e-6);
        assert!((matrix.get(0, 0)).abs() < 1e-6);
        // Symmetry
        assert!((matrix.get(0, 1) - matrix.get(1, 0)).abs() < 1e-6);
    }

    #[test]
    fn test_read_distance_matrix_lower_triangular() {
        let input = "\
    4
Alpha
Beta      1.0000
Gamma     2.0000  2.5000
Delta     3.0000  3.5000  4.0000";
        let matrix = read_distance_matrix(input).unwrap();
        assert_eq!(matrix.size(), 4);
        assert!((matrix.get(1, 0) - 1.0).abs() < 1e-6);
        assert!((matrix.get(2, 0) - 2.0).abs() < 1e-6);
        assert!((matrix.get(2, 1) - 2.5).abs() < 1e-6);
        assert!((matrix.get(3, 2) - 4.0).abs() < 1e-6);
        // Symmetry
        assert!((matrix.get(0, 1) - matrix.get(1, 0)).abs() < 1e-6);
    }

    #[test]
    fn test_read_distance_matrix_2x2() {
        let input = "\
2
A         0.0  0.5
B         0.5  0.0";
        let matrix = read_distance_matrix(input).unwrap();
        assert_eq!(matrix.size(), 2);
        assert!((matrix.get(0, 1) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_read_distance_matrix_empty() {
        let result = read_distance_matrix("");
        assert!(result.is_err());
    }

    #[test]
    fn test_write_distance_matrix() {
        let mut matrix = DistanceMatrix::new(vec![
            "Alpha".to_string(),
            "Beta".to_string(),
            "Gamma".to_string(),
        ]);
        matrix.set(0, 1, 1.0);
        matrix.set(0, 2, 2.0);
        matrix.set(1, 2, 3.0);

        let output = write_distance_matrix(&matrix);
        assert!(output.starts_with("    3\n"));
        assert!(output.contains("Alpha"));
        assert!(output.contains("Beta"));
        assert!(output.contains("Gamma"));

        // Verify it can be read back
        let matrix2 = read_distance_matrix(&output).unwrap();
        assert_eq!(matrix2.size(), 3);
        assert!((matrix2.get(0, 1) - 1.0).abs() < 1e-6);
        assert!((matrix2.get(0, 2) - 2.0).abs() < 1e-6);
        assert!((matrix2.get(1, 2) - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_write_read_roundtrip() {
        let mut matrix = DistanceMatrix::new(vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ]);
        matrix.set(0, 1, 0.5);
        matrix.set(0, 2, 1.0);
        matrix.set(0, 3, 1.5);
        matrix.set(1, 2, 0.8);
        matrix.set(1, 3, 1.2);
        matrix.set(2, 3, 0.3);

        let output = write_distance_matrix(&matrix);
        let matrix2 = read_distance_matrix(&output).unwrap();

        assert_eq!(matrix2.size(), 4);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (matrix.get(i, j) - matrix2.get(i, j)).abs() < 1e-6,
                    "mismatch at ({}, {}): {} vs {}",
                    i,
                    j,
                    matrix.get(i, j),
                    matrix2.get(i, j)
                );
            }
        }
    }

    #[test]
    fn test_display_phylip_format_error() {
        let err = PhylipFormatError::EmptyInput;
        assert_eq!(format!("{}", err), "empty input");

        let err = PhylipFormatError::InvalidHeader("bad".to_string());
        assert_eq!(format!("{}", err), "invalid header line: 'bad'");

        let err = PhylipFormatError::InvalidBase(5, 'Z');
        assert_eq!(format!("{}", err), "invalid base character 'Z' at line 5");
    }

    #[test]
    fn test_distance_matrix_names_preserved() {
        let input = "\
3
Human     0.0  0.5  1.2
Chimp     0.5  0.0  0.8
Gorilla   1.2  0.8  0.0";
        let matrix = read_distance_matrix(input).unwrap();
        assert_eq!(matrix.names[0], "Human");
        assert_eq!(matrix.names[1], "Chimp");
        assert_eq!(matrix.names[2], "Gorilla");
    }
}

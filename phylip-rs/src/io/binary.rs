//! Reader for binary (0/1) character data in PHYLIP format.
//!
//! Parses PHYLIP-formatted binary character matrices where each taxon is
//! scored for a set of characters as either ancestral (0) or derived (1).
//!
//! # Format
//!
//! The first line contains the number of taxa and number of characters.
//! Each subsequent line has a taxon name occupying the first 10 characters,
//! followed by a string of 0s and 1s (spaces allowed for readability).
//!
//! ```text
//! 5  10
//! Alpha     0110010001
//! Beta      0110010010
//! Gamma     1001101001
//! Delta     1001101010
//! Epsilon   0101010101
//! ```

use std::fmt;

use crate::compatibility::BinaryMatrix;

/// Errors that can occur when parsing binary PHYLIP format files.
#[derive(Debug, Clone)]
pub enum BinaryFormatError {
    /// The input is empty.
    EmptyInput,
    /// The header line is malformed (expected "ntaxa nchars").
    InvalidHeader(String),
    /// A taxon line could not be parsed.
    InvalidLine(usize, String),
    /// An invalid character was found (expected '0' or '1').
    InvalidCharacter(usize, char),
    /// The number of taxa found does not match the header.
    TaxonCountMismatch { expected: usize, found: usize },
    /// A taxon has the wrong number of characters.
    CharCountMismatch {
        expected: usize,
        found: usize,
        taxon: String,
    },
}

impl fmt::Display for BinaryFormatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BinaryFormatError::EmptyInput => write!(f, "empty input"),
            BinaryFormatError::InvalidHeader(s) => {
                write!(f, "invalid header line: '{}'", s)
            }
            BinaryFormatError::InvalidLine(line, msg) => {
                write!(f, "invalid line {}: {}", line, msg)
            }
            BinaryFormatError::InvalidCharacter(line, c) => {
                write!(f, "invalid character '{}' at line {} (expected 0 or 1)", c, line)
            }
            BinaryFormatError::TaxonCountMismatch { expected, found } => {
                write!(f, "expected {} taxa but found {}", expected, found)
            }
            BinaryFormatError::CharCountMismatch {
                expected,
                found,
                taxon,
            } => {
                write!(
                    f,
                    "expected {} characters for taxon '{}' but found {}",
                    expected, taxon, found
                )
            }
        }
    }
}

impl std::error::Error for BinaryFormatError {}

/// Parse a PHYLIP-format binary character matrix.
///
/// The first line contains the number of taxa and the number of characters.
/// Each subsequent line has a taxon name (first 10 characters) followed by
/// a sequence of `0` and `1` characters (spaces are ignored within the
/// character data).
///
/// # Examples
///
/// ```
/// use phylip_rs::io::binary::read_binary_phylip;
///
/// let input = "\
/// 3  4
/// Alpha     0110
/// Beta      0101
/// Gamma     1010";
///
/// let matrix = read_binary_phylip(input).unwrap();
/// assert_eq!(matrix.n_taxa(), 3);
/// assert_eq!(matrix.n_chars, 4);
/// ```
pub fn read_binary_phylip(input: &str) -> Result<BinaryMatrix, BinaryFormatError> {
    let input = input.trim();
    if input.is_empty() {
        return Err(BinaryFormatError::EmptyInput);
    }

    let lines: Vec<&str> = input.lines().collect();
    if lines.is_empty() {
        return Err(BinaryFormatError::EmptyInput);
    }

    // Parse header: "ntaxa nchars"
    let header = lines[0].trim();
    let parts: Vec<&str> = header.split_whitespace().collect();
    if parts.len() < 2 {
        return Err(BinaryFormatError::InvalidHeader(header.to_string()));
    }
    let ntaxa = parts[0]
        .parse::<usize>()
        .map_err(|_| BinaryFormatError::InvalidHeader(header.to_string()))?;
    let nchars = parts[1]
        .parse::<usize>()
        .map_err(|_| BinaryFormatError::InvalidHeader(header.to_string()))?;

    if ntaxa == 0 || nchars == 0 {
        return Err(BinaryFormatError::InvalidHeader(
            "ntaxa and nchars must be positive".to_string(),
        ));
    }

    // Parse data lines
    let data_lines: Vec<&str> = lines[1..]
        .iter()
        .filter(|l| !l.trim().is_empty())
        .copied()
        .collect();

    if data_lines.len() < ntaxa {
        return Err(BinaryFormatError::TaxonCountMismatch {
            expected: ntaxa,
            found: data_lines.len(),
        });
    }

    let mut taxa: Vec<String> = Vec::with_capacity(ntaxa);
    let mut characters: Vec<Vec<bool>> = Vec::with_capacity(ntaxa);

    for (i, &line) in data_lines.iter().enumerate().take(ntaxa) {
        let line_num = i + 2; // 1-indexed, header is line 1

        if line.is_empty() {
            return Err(BinaryFormatError::InvalidLine(
                line_num,
                "empty line where taxon expected".to_string(),
            ));
        }

        // Split name (first 10 chars) and character data
        let (name, char_part) = if line.len() >= 10 {
            let name = line[..10].trim().to_string();
            let data = &line[10..];
            (name, data)
        } else {
            // Short line: use entire line as name, no character data
            (line.trim().to_string(), "")
        };

        if name.is_empty() {
            return Err(BinaryFormatError::InvalidLine(
                line_num,
                "missing taxon name".to_string(),
            ));
        }

        // Parse character data
        let mut chars: Vec<bool> = Vec::with_capacity(nchars);
        for c in char_part.chars() {
            if c.is_whitespace() {
                continue;
            }
            match c {
                '0' => chars.push(false),
                '1' => chars.push(true),
                _ => return Err(BinaryFormatError::InvalidCharacter(line_num, c)),
            }
        }

        if chars.len() != nchars {
            return Err(BinaryFormatError::CharCountMismatch {
                expected: nchars,
                found: chars.len(),
                taxon: name,
            });
        }

        taxa.push(name);
        characters.push(chars);
    }

    BinaryMatrix::new(taxa, characters).map_err(|e| {
        BinaryFormatError::InvalidLine(0, e)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_binary_simple() {
        let input = "\
 3  4
Alpha     0110
Beta      0101
Gamma     1010";
        let matrix = read_binary_phylip(input).unwrap();
        assert_eq!(matrix.n_taxa(), 3);
        assert_eq!(matrix.n_chars, 4);
        assert_eq!(matrix.taxa[0], "Alpha");
        assert_eq!(matrix.taxa[1], "Beta");
        assert_eq!(matrix.taxa[2], "Gamma");
        // Alpha: 0110
        assert!(!matrix.characters[0][0]);
        assert!(matrix.characters[0][1]);
        assert!(matrix.characters[0][2]);
        assert!(!matrix.characters[0][3]);
    }

    #[test]
    fn test_read_binary_with_spaces() {
        let input = "\
 2  8
Alpha     0110 0001
Beta      0110 0010";
        let matrix = read_binary_phylip(input).unwrap();
        assert_eq!(matrix.n_taxa(), 2);
        assert_eq!(matrix.n_chars, 8);
    }

    #[test]
    fn test_read_binary_five_taxa() {
        let input = "\
5  10
Alpha     0110010001
Beta      0110010010
Gamma     1001101001
Delta     1001101010
Epsilon   0101010101";
        let matrix = read_binary_phylip(input).unwrap();
        assert_eq!(matrix.n_taxa(), 5);
        assert_eq!(matrix.n_chars, 10);
        assert_eq!(matrix.taxa[4], "Epsilon");
    }

    #[test]
    fn test_read_binary_empty() {
        let result = read_binary_phylip("");
        assert!(result.is_err());
        match result {
            Err(BinaryFormatError::EmptyInput) => {}
            other => panic!("expected EmptyInput, got {:?}", other),
        }
    }

    #[test]
    fn test_read_binary_invalid_header() {
        let result = read_binary_phylip("not valid\n");
        assert!(result.is_err());
    }

    #[test]
    fn test_read_binary_invalid_character() {
        let input = "\
 2  4
Alpha     012X
Beta      0101";
        let result = read_binary_phylip(input);
        assert!(result.is_err());
        match result {
            Err(BinaryFormatError::InvalidCharacter(_, c)) => {
                assert!(c == '2' || c == 'X');
            }
            other => panic!("expected InvalidCharacter, got {:?}", other),
        }
    }

    #[test]
    fn test_read_binary_wrong_char_count() {
        let input = "\
 2  4
Alpha     011
Beta      0101";
        let result = read_binary_phylip(input);
        assert!(result.is_err());
        match result {
            Err(BinaryFormatError::CharCountMismatch { expected: 4, found: 3, .. }) => {}
            other => panic!("expected CharCountMismatch, got {:?}", other),
        }
    }

    #[test]
    fn test_read_binary_too_few_taxa() {
        let input = "\
 3  4
Alpha     0110
Beta      0101";
        let result = read_binary_phylip(input);
        assert!(result.is_err());
        match result {
            Err(BinaryFormatError::TaxonCountMismatch { expected: 3, found: 2 }) => {}
            other => panic!("expected TaxonCountMismatch, got {:?}", other),
        }
    }

    #[test]
    fn test_read_binary_ten_char_names() {
        let input = "\
 2  4
1234567890 0110
ABCDEFGHIJ 0101";
        // Name is exactly 10 chars; the 11th char (space) starts character data
        let matrix = read_binary_phylip(input).unwrap();
        assert_eq!(matrix.taxa[0], "1234567890");
        assert_eq!(matrix.taxa[1], "ABCDEFGHIJ");
        assert_eq!(matrix.n_chars, 4);
    }

    #[test]
    fn test_read_binary_all_zeros() {
        let input = "\
 2  4
Alpha     0000
Beta      0000";
        let matrix = read_binary_phylip(input).unwrap();
        for t in 0..2 {
            for c in 0..4 {
                assert!(!matrix.characters[t][c]);
            }
        }
    }

    #[test]
    fn test_read_binary_all_ones() {
        let input = "\
 2  4
Alpha     1111
Beta      1111";
        let matrix = read_binary_phylip(input).unwrap();
        for t in 0..2 {
            for c in 0..4 {
                assert!(matrix.characters[t][c]);
            }
        }
    }

    #[test]
    fn test_display_binary_format_error() {
        let err = BinaryFormatError::EmptyInput;
        assert_eq!(format!("{}", err), "empty input");

        let err = BinaryFormatError::InvalidCharacter(3, 'X');
        assert_eq!(
            format!("{}", err),
            "invalid character 'X' at line 3 (expected 0 or 1)"
        );

        let err = BinaryFormatError::CharCountMismatch {
            expected: 10,
            found: 8,
            taxon: "Alpha".to_string(),
        };
        assert_eq!(
            format!("{}", err),
            "expected 10 characters for taxon 'Alpha' but found 8"
        );
    }
}

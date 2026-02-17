//! FASTA format reader.
//!
//! Parses FASTA-formatted sequence files into `Alignment` structures.
//! FASTA format uses `>` headers followed by sequence data on subsequent lines.

use std::fmt;

use crate::tree::types::{Alignment, Base, Sequence};

/// Errors that can occur when parsing FASTA format files.
#[derive(Debug, Clone)]
pub enum FastaError {
    /// The input is empty or contains no sequences.
    EmptyInput,
    /// The input does not start with a `>` header line.
    MissingHeader,
    /// A sequence header line is malformed.
    InvalidHeader(usize, String),
    /// An invalid character was found in a sequence.
    InvalidBase(usize, char),
    /// A header was found with no sequence data following it.
    EmptySequence(String),
    /// Sequences have different lengths (not a valid alignment).
    UnequalLengths,
}

impl fmt::Display for FastaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FastaError::EmptyInput => write!(f, "empty input"),
            FastaError::MissingHeader => {
                write!(f, "FASTA input does not start with '>' header")
            }
            FastaError::InvalidHeader(line, msg) => {
                write!(f, "invalid header at line {}: {}", line, msg)
            }
            FastaError::InvalidBase(line, c) => {
                write!(f, "invalid base character '{}' at line {}", c, line)
            }
            FastaError::EmptySequence(name) => {
                write!(f, "sequence '{}' has no data", name)
            }
            FastaError::UnequalLengths => {
                write!(f, "sequences have different lengths")
            }
        }
    }
}

impl std::error::Error for FastaError {}

/// Parse a FASTA format string into an `Alignment`.
///
/// FASTA format consists of one or more sequences, each starting with a
/// header line beginning with `>`, followed by one or more lines of
/// sequence data.
///
/// # Format
/// ```text
/// >Sequence1 optional description
/// AAGCTNGGGCATTTCAGGGT
/// GAGCCCGGGCAATACAGGGT
/// >Sequence2
/// AAGCCTTGGCAGTGCAGGGT
/// GAGCCGTGGCCGGGCACGGT
/// ```
///
/// The sequence name is taken from the header line (text after `>` up to
/// the first whitespace). The rest of the header line is treated as a
/// description and is discarded.
///
/// Sequences may span multiple lines. All sequences in the alignment must
/// have the same length.
///
/// # Examples
///
/// ```
/// use phylip_rs::io::fasta::read_fasta;
///
/// let input = ">Seq1\nACGT\n>Seq2\nACGT\n";
/// let alignment = read_fasta(input).unwrap();
/// assert_eq!(alignment.ntaxa(), 2);
/// assert_eq!(alignment.nsites(), 4);
/// ```
pub fn read_fasta(input: &str) -> Result<Alignment, FastaError> {
    let input = input.trim();
    if input.is_empty() {
        return Err(FastaError::EmptyInput);
    }

    let lines: Vec<&str> = input.lines().collect();
    if lines.is_empty() {
        return Err(FastaError::EmptyInput);
    }

    // The first non-empty, non-comment line should start with '>'
    let first_content = lines.iter().find(|l| {
        let trimmed = l.trim();
        !trimmed.is_empty() && !trimmed.starts_with('#')
    });
    match first_content {
        Some(line) if !line.trim().starts_with('>') => {
            return Err(FastaError::MissingHeader);
        }
        None => return Err(FastaError::EmptyInput),
        _ => {}
    }

    let mut sequences: Vec<Sequence> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_bases: Vec<Base> = Vec::new();

    for (line_idx, &line) in lines.iter().enumerate() {
        let trimmed = line.trim();
        let line_num = line_idx + 1;

        // Skip empty lines and comments
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed.starts_with('>') {
            // Save previous sequence if any
            if let Some(name) = current_name.take() {
                if current_bases.is_empty() {
                    return Err(FastaError::EmptySequence(name));
                }
                sequences.push(Sequence::new(name, current_bases));
                current_bases = Vec::new();
            }

            // Parse header line
            let header = &trimmed[1..]; // strip '>'
            let name = parse_fasta_name(header);
            if name.is_empty() {
                return Err(FastaError::InvalidHeader(
                    line_num,
                    "empty sequence name".to_string(),
                ));
            }
            current_name = Some(name);
        } else {
            // Sequence data line
            if current_name.is_none() {
                return Err(FastaError::MissingHeader);
            }

            for c in trimmed.chars() {
                if c.is_whitespace() || c.is_ascii_digit() {
                    continue;
                }
                match Base::from_char(c) {
                    Some(b) => current_bases.push(b),
                    None => return Err(FastaError::InvalidBase(line_num, c)),
                }
            }
        }
    }

    // Save the last sequence
    if let Some(name) = current_name.take() {
        if current_bases.is_empty() {
            return Err(FastaError::EmptySequence(name));
        }
        sequences.push(Sequence::new(name, current_bases));
    }

    if sequences.is_empty() {
        return Err(FastaError::EmptyInput);
    }

    Alignment::new(sequences).map_err(|_| FastaError::UnequalLengths)
}

/// Extract the sequence name from a FASTA header line (after the `>`).
/// Takes the first whitespace-delimited token as the name.
fn parse_fasta_name(header: &str) -> String {
    let trimmed = header.trim();
    match trimmed.split_whitespace().next() {
        Some(name) => name.to_string(),
        None => String::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::types::Base;

    #[test]
    fn test_simple_fasta() {
        let input = ">Seq1\nACGT\n>Seq2\nACGT\n";
        let alignment = read_fasta(input).unwrap();
        assert_eq!(alignment.ntaxa(), 2);
        assert_eq!(alignment.nsites(), 4);
        assert_eq!(alignment.taxon_name(0), "Seq1");
        assert_eq!(alignment.taxon_name(1), "Seq2");
    }

    #[test]
    fn test_multiline_sequences() {
        let input = ">Seq1\nACGT\nACGT\n>Seq2\nAAAA\nCCCC\n";
        let alignment = read_fasta(input).unwrap();
        assert_eq!(alignment.ntaxa(), 2);
        assert_eq!(alignment.nsites(), 8);
    }

    #[test]
    fn test_header_with_description() {
        let input = ">Seq1 this is a description\nACGT\n>Seq2 another description\nACGT\n";
        let alignment = read_fasta(input).unwrap();
        assert_eq!(alignment.taxon_name(0), "Seq1");
        assert_eq!(alignment.taxon_name(1), "Seq2");
    }

    #[test]
    fn test_various_bases() {
        let input = ">Seq1\nACGTNMRWSYKVHDB-\n>Seq2\nACGTNMRWSYKVHDB-\n";
        let alignment = read_fasta(input).unwrap();
        assert_eq!(alignment.nsites(), 16);
        assert_eq!(alignment.base(0, 0), Base::A);
        assert_eq!(alignment.base(0, 1), Base::C);
        assert_eq!(alignment.base(0, 2), Base::G);
        assert_eq!(alignment.base(0, 3), Base::T);
        assert_eq!(alignment.base(0, 4), Base::N);
        assert_eq!(alignment.base(0, 15), Base::Gap);
    }

    #[test]
    fn test_lowercase_sequences() {
        let input = ">Seq1\nacgt\n>Seq2\nacgt\n";
        let alignment = read_fasta(input).unwrap();
        assert_eq!(alignment.base(0, 0), Base::A);
        assert_eq!(alignment.base(0, 1), Base::C);
        assert_eq!(alignment.base(0, 2), Base::G);
        assert_eq!(alignment.base(0, 3), Base::T);
    }

    #[test]
    fn test_blank_lines_between_sequences() {
        let input = ">Seq1\nACGT\n\n>Seq2\nACGT\n";
        let alignment = read_fasta(input).unwrap();
        assert_eq!(alignment.ntaxa(), 2);
    }

    #[test]
    fn test_single_sequence() {
        let input = ">OnlyOne\nACGTACGT\n";
        let alignment = read_fasta(input).unwrap();
        assert_eq!(alignment.ntaxa(), 1);
        assert_eq!(alignment.nsites(), 8);
        assert_eq!(alignment.taxon_name(0), "OnlyOne");
    }

    #[test]
    fn test_many_sequences() {
        let mut input = String::new();
        for i in 0..50 {
            input.push_str(&format!(">Taxon{}\nACGTACGTAC\n", i));
        }
        let alignment = read_fasta(&input).unwrap();
        assert_eq!(alignment.ntaxa(), 50);
        assert_eq!(alignment.nsites(), 10);
    }

    #[test]
    fn test_empty_input() {
        let result = read_fasta("");
        assert!(result.is_err());
        match result {
            Err(FastaError::EmptyInput) => {}
            other => panic!("expected EmptyInput, got {:?}", other),
        }
    }

    #[test]
    fn test_missing_header() {
        let result = read_fasta("ACGT\nACGT\n");
        assert!(result.is_err());
        match result {
            Err(FastaError::MissingHeader) => {}
            other => panic!("expected MissingHeader, got {:?}", other),
        }
    }

    #[test]
    fn test_empty_sequence_name() {
        let result = read_fasta(">\nACGT\n");
        assert!(result.is_err());
        match result {
            Err(FastaError::InvalidHeader(_, _)) => {}
            other => panic!("expected InvalidHeader, got {:?}", other),
        }
    }

    #[test]
    fn test_empty_sequence_data() {
        let result = read_fasta(">Seq1\n>Seq2\nACGT\n");
        assert!(result.is_err());
        match result {
            Err(FastaError::EmptySequence(name)) => assert_eq!(name, "Seq1"),
            other => panic!("expected EmptySequence, got {:?}", other),
        }
    }

    #[test]
    fn test_unequal_lengths() {
        let result = read_fasta(">Seq1\nACGT\n>Seq2\nACG\n");
        assert!(result.is_err());
        match result {
            Err(FastaError::UnequalLengths) => {}
            other => panic!("expected UnequalLengths, got {:?}", other),
        }
    }

    #[test]
    fn test_invalid_base() {
        let result = read_fasta(">Seq1\nACG123T\n>Seq2\nACG1Z3T\n");
        // Digits are ignored, but Z is invalid
        assert!(result.is_err());
        match result {
            Err(FastaError::InvalidBase(_, 'Z')) => {}
            other => panic!("expected InvalidBase with 'Z', got {:?}", other),
        }
    }

    #[test]
    fn test_long_multiline_sequences() {
        let seq = "ACGT".repeat(250); // 1000 bases
        let chunk_size = 80;
        let mut fasta_seq = String::new();
        for chunk in seq.as_bytes().chunks(chunk_size) {
            fasta_seq.push_str(std::str::from_utf8(chunk).unwrap());
            fasta_seq.push('\n');
        }
        let input = format!(">Seq1\n{}>Seq2\n{}", fasta_seq, fasta_seq);
        let alignment = read_fasta(&input).unwrap();
        assert_eq!(alignment.ntaxa(), 2);
        assert_eq!(alignment.nsites(), 1000);
    }

    #[test]
    fn test_uracil() {
        // U should be treated as T (RNA sequences)
        let input = ">RNA1\nACGU\n>RNA2\nACGU\n";
        let alignment = read_fasta(input).unwrap();
        assert_eq!(alignment.base(0, 3), Base::T);
    }

    #[test]
    fn test_realistic_fasta() {
        let input = "\
>Homo_sapiens cytochrome b
ATGACCCCAATACGCAAAACTAACCCCCTAATAAAATTCATTAACCACTCATTCATCGAC
CTCCCAACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGC
>Pan_troglodytes cytochrome b
ATGACCCTAATACGCAAAACTAATCCCCTAATAAAATTCATTAATCACTCATTCATTGAC
CTCCCAACTCCATCCAACATCTCTGCATGATGAAACTTCGGCTCACTCCTAGGCACCTGC
>Gorilla_gorilla cytochrome b
ATGACTCCAATACGCAAAACTAACCCCCTAATAAAATTCATCAACCACTCATTCATTGAC
CTCCCAGCCCCATCTAACATCTCTGCATGATGAAACTTTGGCTCACTCCTCGGCGCCTGC";
        let alignment = read_fasta(input).unwrap();
        assert_eq!(alignment.ntaxa(), 3);
        assert_eq!(alignment.nsites(), 120);
        assert_eq!(alignment.taxon_name(0), "Homo_sapiens");
        assert_eq!(alignment.taxon_name(1), "Pan_troglodytes");
        assert_eq!(alignment.taxon_name(2), "Gorilla_gorilla");
    }

    #[test]
    fn test_display_fasta_error() {
        let err = FastaError::EmptyInput;
        assert_eq!(format!("{}", err), "empty input");

        let err = FastaError::MissingHeader;
        assert_eq!(
            format!("{}", err),
            "FASTA input does not start with '>' header"
        );

        let err = FastaError::InvalidBase(5, 'Z');
        assert_eq!(format!("{}", err), "invalid base character 'Z' at line 5");

        let err = FastaError::EmptySequence("test".to_string());
        assert_eq!(format!("{}", err), "sequence 'test' has no data");

        let err = FastaError::UnequalLengths;
        assert_eq!(format!("{}", err), "sequences have different lengths");
    }

    #[test]
    fn test_whitespace_only_input() {
        let result = read_fasta("   \n\n   \n");
        assert!(result.is_err());
    }

    #[test]
    fn test_comment_lines() {
        let input = "# This is a comment\n>Seq1\nACGT\n>Seq2\nACGT\n";
        let alignment = read_fasta(input).unwrap();
        assert_eq!(alignment.ntaxa(), 2);
    }
}

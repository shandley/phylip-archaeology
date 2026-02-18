//! Comparative methods for continuous character evolution.
//!
//! This module implements Felsenstein's phylogenetic comparative methods
//! for continuous character data:
//!
//! - `contrasts` - Phylogenetically independent contrasts (Felsenstein, 1985)
//! - `contml` - Maximum likelihood under Brownian motion (Felsenstein, 1973, 1981)
//!
//! These methods account for the non-independence of species data arising
//! from shared evolutionary history. Independent contrasts remove phylogenetic
//! correlation, while the Brownian motion ML approach directly models character
//! evolution as a random walk on the tree.

pub mod contrasts;
pub mod contml;

/// Continuous character data for multiple taxa.
///
/// Each taxon has measurements for one or more continuous characters
/// (e.g., body mass, limb length, metabolic rate). All taxa must have
/// the same number of characters.
#[derive(Debug, Clone)]
pub struct ContinuousData {
    /// Taxon names, one per row.
    pub taxa: Vec<String>,
    /// Character values: `values[taxon][character]` = measurement value.
    pub values: Vec<Vec<f64>>,
    /// Number of characters measured per taxon.
    pub n_chars: usize,
}

impl ContinuousData {
    /// Create a new ContinuousData from taxa names and values.
    ///
    /// Returns an error if the dimensions are inconsistent.
    pub fn new(taxa: Vec<String>, values: Vec<Vec<f64>>) -> Result<Self, String> {
        if taxa.is_empty() {
            return Err("no taxa provided".to_string());
        }
        if taxa.len() != values.len() {
            return Err(format!(
                "number of taxa ({}) does not match number of value rows ({})",
                taxa.len(),
                values.len()
            ));
        }
        let n_chars = values[0].len();
        if n_chars == 0 {
            return Err("no characters provided".to_string());
        }
        for (i, row) in values.iter().enumerate() {
            if row.len() != n_chars {
                return Err(format!(
                    "taxon {} has {} characters, expected {}",
                    taxa[i],
                    row.len(),
                    n_chars
                ));
            }
        }
        Ok(Self {
            taxa,
            values,
            n_chars,
        })
    }

    /// Number of taxa.
    pub fn n_taxa(&self) -> usize {
        self.taxa.len()
    }

    /// Get the value for a given taxon and character index.
    pub fn value(&self, taxon: usize, character: usize) -> f64 {
        self.values[taxon][character]
    }
}

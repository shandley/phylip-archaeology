//! Restriction site distance estimation (Nei & Li, 1979).
//!
//! Computes evolutionary distances from restriction site presence/absence
//! data, corresponding to PHYLIP's `restdist` program.
//!
//! Restriction enzymes cut DNA at specific recognition sequences (typically
//! 4, 6, or 8 base pairs). The presence or absence of restriction sites
//! across taxa provides information about nucleotide sequence divergence.
//!
//! # Distance formula
//!
//! Given two taxa i and j with restriction site presence/absence data:
//!
//! 1. Compute the proportion of shared sites:
//!    ```text
//!    S_hat = 2 * n_shared / (n_i + n_j)
//!    ```
//!    where `n_i` = number of sites present in taxon i, `n_j` = in taxon j,
//!    `n_shared` = number of sites present in both.
//!
//! 2. Estimate the per-nucleotide substitution rate:
//!    ```text
//!    p_hat = 1 - S_hat^(1/r)
//!    ```
//!    where `r` = length of the recognition sequence (number of bases).
//!
//! 3. Apply a Jukes-Cantor-like correction:
//!    ```text
//!    d = -(3/4) * ln(1 - (4/3) * p_hat)
//!    ```
//!
//! # References
//!
//! - Nei, M. & Li, W.-H. (1979). Mathematical model for studying genetic
//!   variation in terms of restriction endonucleases. *Proceedings of the
//!   National Academy of Sciences*, 76(10), 5269-5273.

use crate::tree::DistanceMatrix;

/// Restriction site data: presence/absence of restriction sites across taxa.
///
/// Each taxon has a boolean vector indicating which restriction sites are
/// present (true) or absent (false).
#[derive(Debug, Clone)]
pub struct RestrictionData {
    /// Taxon names.
    pub taxa: Vec<String>,
    /// sites\[taxon\]\[site\] = true if restriction site is present.
    pub sites: Vec<Vec<bool>>,
    /// Total number of restriction sites surveyed.
    pub n_sites: usize,
    /// Length of the recognition sequence in base pairs (typically 4, 6, or 8).
    pub site_length: usize,
}

impl RestrictionData {
    /// Create new restriction site data.
    ///
    /// # Arguments
    ///
    /// * `taxa` - Names of the taxa
    /// * `sites` - Presence/absence matrix: `sites[taxon][site]`
    /// * `site_length` - Length of the restriction enzyme recognition sequence
    ///
    /// # Errors
    ///
    /// Returns an error if the data is inconsistent (empty, unequal lengths,
    /// or invalid site length).
    pub fn new(
        taxa: Vec<String>,
        sites: Vec<Vec<bool>>,
        site_length: usize,
    ) -> Result<Self, String> {
        if taxa.is_empty() {
            return Err("no taxa provided".to_string());
        }
        if taxa.len() != sites.len() {
            return Err("number of taxa does not match number of site profiles".to_string());
        }
        if site_length == 0 {
            return Err("site length must be positive".to_string());
        }

        let n_sites = if sites.is_empty() {
            0
        } else {
            sites[0].len()
        };

        if sites.iter().any(|s| s.len() != n_sites) {
            return Err("all taxa must have the same number of sites".to_string());
        }

        Ok(Self {
            taxa,
            sites,
            n_sites,
            site_length,
        })
    }

    /// Number of taxa.
    pub fn ntaxa(&self) -> usize {
        self.taxa.len()
    }
}

/// Compute the Nei-Li distance between two taxa from restriction site data.
///
/// # Arguments
///
/// * `data` - The restriction site presence/absence data
/// * `taxon_i` - Index of the first taxon
/// * `taxon_j` - Index of the second taxon
///
/// # Returns
///
/// The estimated evolutionary distance, or an error if the distance cannot
/// be computed (e.g., no shared sites, or too divergent).
pub fn restriction_distance(
    data: &RestrictionData,
    taxon_i: usize,
    taxon_j: usize,
) -> Result<f64, String> {
    if taxon_i >= data.ntaxa() || taxon_j >= data.ntaxa() {
        return Err(format!(
            "taxon index out of range: i={}, j={}, ntaxa={}",
            taxon_i,
            taxon_j,
            data.ntaxa()
        ));
    }

    if taxon_i == taxon_j {
        return Ok(0.0);
    }

    let sites_i = &data.sites[taxon_i];
    let sites_j = &data.sites[taxon_j];

    // Count sites present in each taxon and shared between both
    let n_i: usize = sites_i.iter().filter(|&&s| s).count();
    let n_j: usize = sites_j.iter().filter(|&&s| s).count();
    let n_shared: usize = sites_i
        .iter()
        .zip(sites_j.iter())
        .filter(|(&a, &b)| a && b)
        .count();

    // Handle edge cases
    if n_i == 0 && n_j == 0 {
        // Both taxa have no sites at all -- undefined
        return Err(format!(
            "both taxa {} and {} have no restriction sites",
            taxon_i, taxon_j
        ));
    }

    if n_i + n_j == 0 {
        return Err("no restriction sites in either taxon".to_string());
    }

    // Identical site profiles
    if n_i == n_j && n_shared == n_i {
        return Ok(0.0);
    }

    // Proportion of shared sites
    let s_hat = 2.0 * n_shared as f64 / (n_i + n_j) as f64;

    if s_hat <= 0.0 {
        return Err(format!(
            "no shared restriction sites between taxa {} and {} (S_hat = 0); distance is infinite",
            taxon_i, taxon_j
        ));
    }

    let r = data.site_length as f64;

    // Estimate per-nucleotide substitution rate
    // p_hat = 1 - S_hat^(1/r)
    let p_hat = 1.0 - s_hat.powf(1.0 / r);

    // JC69-like correction
    // d = -(3/4) * ln(1 - (4/3) * p_hat)
    let inner = 1.0 - (4.0 / 3.0) * p_hat;

    if inner <= 0.0 {
        return Err(format!(
            "taxa {} and {} are too divergent for JC69 correction (1 - 4/3 * p_hat = {:.6})",
            taxon_i, taxon_j, inner
        ));
    }

    let d = -0.75 * inner.ln();

    if !d.is_finite() || d < 0.0 {
        return Err(format!(
            "computed distance between taxa {} and {} is not finite or negative (d = {:.6})",
            taxon_i, taxon_j, d
        ));
    }

    Ok(d)
}

/// Compute a full pairwise distance matrix from restriction site data.
///
/// # Returns
///
/// A symmetric `DistanceMatrix` with Nei-Li distances, or an error if
/// any pairwise distance computation fails.
pub fn compute_restriction_distance_matrix(
    data: &RestrictionData,
) -> Result<DistanceMatrix, String> {
    let n = data.ntaxa();
    let names: Vec<String> = data.taxa.clone();
    let mut matrix = DistanceMatrix::new(names);

    for i in 0..n {
        for j in (i + 1)..n {
            let d = restriction_distance(data, i, j)?;
            matrix.set(i, j, d);
        }
    }

    Ok(matrix)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: create RestrictionData from a list of (name, site_vector) pairs.
    fn make_restriction_data(
        taxa: &[(&str, &[bool])],
        site_length: usize,
    ) -> RestrictionData {
        let names: Vec<String> = taxa.iter().map(|(n, _)| n.to_string()).collect();
        let sites: Vec<Vec<bool>> = taxa.iter().map(|(_, s)| s.to_vec()).collect();
        RestrictionData::new(names, sites, site_length).unwrap()
    }

    // -----------------------------------------------------------------------
    // Identical site profiles: distance = 0
    // -----------------------------------------------------------------------

    #[test]
    fn test_identical_profiles() {
        let data = make_restriction_data(
            &[
                ("A", &[true, true, false, true, false]),
                ("B", &[true, true, false, true, false]),
            ],
            6,
        );
        let d = restriction_distance(&data, 0, 1).unwrap();
        assert!(
            d.abs() < 1e-10,
            "identical site profiles should have distance 0, got {}",
            d
        );
    }

    #[test]
    fn test_self_distance() {
        let data = make_restriction_data(
            &[("A", &[true, true, false, true, false])],
            6,
        );
        let d = restriction_distance(&data, 0, 0).unwrap();
        assert_eq!(d, 0.0);
    }

    // -----------------------------------------------------------------------
    // No shared sites: maximum/infinite distance
    // -----------------------------------------------------------------------

    #[test]
    fn test_no_shared_sites() {
        let data = make_restriction_data(
            &[
                ("A", &[true, true, false, false, false]),
                ("B", &[false, false, true, true, true]),
            ],
            6,
        );
        let result = restriction_distance(&data, 0, 1);
        assert!(
            result.is_err(),
            "no shared sites should return error (infinite distance)"
        );
    }

    // -----------------------------------------------------------------------
    // Known analytical examples
    // -----------------------------------------------------------------------

    #[test]
    fn test_known_distance_50_percent_shared() {
        // Taxon A: 4 sites present, Taxon B: 4 sites present, 2 shared
        // S_hat = 2 * 2 / (4 + 4) = 0.5
        let data = make_restriction_data(
            &[
                ("A", &[true, true, true, true, false, false]),
                ("B", &[true, true, false, false, true, true]),
            ],
            6,
        );
        let d = restriction_distance(&data, 0, 1).unwrap();

        // S_hat = 0.5
        let s_hat: f64 = 0.5;
        let r: f64 = 6.0;
        let p_hat = 1.0 - s_hat.powf(1.0 / r);
        let expected = -0.75 * (1.0 - (4.0 / 3.0) * p_hat).ln();

        assert!(
            (d - expected).abs() < 1e-10,
            "S_hat=0.5, r=6: expected {:.6}, got {:.6}",
            expected,
            d
        );
    }

    #[test]
    fn test_known_distance_75_percent_shared() {
        // Taxon A: 4 sites, Taxon B: 4 sites, 3 shared
        // S_hat = 2 * 3 / (4 + 4) = 0.75
        let data = make_restriction_data(
            &[
                ("A", &[true, true, true, true, false]),
                ("B", &[true, true, true, false, true]),
            ],
            6,
        );
        let d = restriction_distance(&data, 0, 1).unwrap();

        let s_hat: f64 = 0.75;
        let r: f64 = 6.0;
        let p_hat = 1.0 - s_hat.powf(1.0 / r);
        let expected = -0.75 * (1.0 - (4.0 / 3.0) * p_hat).ln();

        assert!(
            (d - expected).abs() < 1e-10,
            "S_hat=0.75, r=6: expected {:.6}, got {:.6}",
            expected,
            d
        );
    }

    #[test]
    fn test_known_distance_high_similarity() {
        // Nearly identical: 10 sites, 9 shared
        // A: 10 present, B: 10 present, 9 shared
        // S_hat = 2*9 / (10+10) = 0.9
        let mut a_sites = vec![true; 10];
        let mut b_sites = vec![true; 10];
        a_sites.push(false); // site 11: A absent, B present
        b_sites.push(true);
        let data = make_restriction_data(
            &[("A", &a_sites), ("B", &b_sites)],
            6,
        );
        let d = restriction_distance(&data, 0, 1).unwrap();

        let n_i: f64 = 10.0;
        let n_j: f64 = 11.0;
        let n_shared: f64 = 10.0;
        let s_hat: f64 = 2.0 * n_shared / (n_i + n_j);
        let r: f64 = 6.0;
        let p_hat: f64 = 1.0 - s_hat.powf(1.0 / r);
        let expected: f64 = -0.75 * (1.0 - (4.0 / 3.0) * p_hat).ln();

        assert!(
            (d - expected).abs() < 1e-10,
            "expected {:.6}, got {:.6}",
            expected,
            d
        );
    }

    // -----------------------------------------------------------------------
    // Different site lengths (4 vs 6 vs 8): affect distance scaling
    // -----------------------------------------------------------------------

    #[test]
    fn test_site_length_4_vs_6_vs_8() {
        // Same data, different site lengths
        // Longer recognition sequences -> larger p_hat -> larger distance
        let a_sites = vec![true, true, true, true, false, false, false, false];
        let b_sites = vec![true, true, false, false, true, true, false, false];
        // n_i = 4, n_j = 4, n_shared = 2
        // S_hat = 2*2/(4+4) = 0.5

        let data4 = make_restriction_data(&[("A", &a_sites), ("B", &b_sites)], 4);
        let data6 = make_restriction_data(&[("A", &a_sites), ("B", &b_sites)], 6);
        let data8 = make_restriction_data(&[("A", &a_sites), ("B", &b_sites)], 8);

        let d4 = restriction_distance(&data4, 0, 1).unwrap();
        let d6 = restriction_distance(&data6, 0, 1).unwrap();
        let d8 = restriction_distance(&data8, 0, 1).unwrap();

        // Shorter recognition sequence -> larger S_hat^(1/r) -> smaller p_hat -> smaller d
        // Wait, S_hat is the same. With larger r:
        // S_hat^(1/r) gets closer to 1 (for 0 < S_hat < 1)
        // So p_hat = 1 - S_hat^(1/r) gets smaller with larger r
        // So d gets smaller with larger r
        assert!(
            d4 > d6,
            "d(r=4)={:.6} should be > d(r=6)={:.6}",
            d4,
            d6
        );
        assert!(
            d6 > d8,
            "d(r=6)={:.6} should be > d(r=8)={:.6}",
            d6,
            d8
        );
    }

    #[test]
    fn test_site_length_analytical_r4() {
        // S_hat = 0.5, r = 4
        let data = make_restriction_data(
            &[
                ("A", &[true, true, false, false]),
                ("B", &[true, false, true, false]),
            ],
            4,
        );
        let d = restriction_distance(&data, 0, 1).unwrap();

        let s_hat: f64 = 0.5; // 2*1/(2+2) = 0.5
        let r: f64 = 4.0;
        let p_hat = 1.0 - s_hat.powf(1.0 / r);
        let expected = -0.75 * (1.0 - (4.0 / 3.0) * p_hat).ln();

        assert!(
            (d - expected).abs() < 1e-10,
            "r=4: expected {:.6}, got {:.6}",
            expected,
            d
        );
    }

    // -----------------------------------------------------------------------
    // Distance matrix computation
    // -----------------------------------------------------------------------

    #[test]
    fn test_distance_matrix() {
        let data = make_restriction_data(
            &[
                ("A", &[true, true, true, false, false]),
                ("B", &[true, true, false, true, false]),
                ("C", &[true, false, false, false, true]),
            ],
            6,
        );
        let matrix = compute_restriction_distance_matrix(&data).unwrap();

        assert_eq!(matrix.size(), 3);
        assert_eq!(matrix.names, vec!["A", "B", "C"]);

        // Symmetry
        assert_eq!(matrix.get(0, 1), matrix.get(1, 0));
        assert_eq!(matrix.get(0, 2), matrix.get(2, 0));
        assert_eq!(matrix.get(1, 2), matrix.get(2, 1));

        // Diagonal = 0
        assert_eq!(matrix.get(0, 0), 0.0);
        assert_eq!(matrix.get(1, 1), 0.0);
        assert_eq!(matrix.get(2, 2), 0.0);

        // All off-diagonal > 0
        assert!(matrix.get(0, 1) > 0.0);
        assert!(matrix.get(0, 2) > 0.0);
        assert!(matrix.get(1, 2) > 0.0);
    }

    #[test]
    fn test_distance_matrix_identical() {
        let data = make_restriction_data(
            &[
                ("A", &[true, true, false]),
                ("B", &[true, true, false]),
            ],
            6,
        );
        let matrix = compute_restriction_distance_matrix(&data).unwrap();
        assert!(
            matrix.get(0, 1).abs() < 1e-10,
            "identical profiles should have distance 0 in matrix"
        );
    }

    // -----------------------------------------------------------------------
    // Distance increases with fewer shared sites
    // -----------------------------------------------------------------------

    #[test]
    fn test_distance_increases_with_fewer_shared() {
        // More shared sites -> smaller distance
        let data_high = make_restriction_data(
            &[
                ("A", &[true, true, true, true, true]),
                ("B", &[true, true, true, true, false]),
            ],
            6,
        );
        let data_low = make_restriction_data(
            &[
                ("A", &[true, true, true, true, true]),
                ("B", &[true, true, false, false, false]),
            ],
            6,
        );

        let d_high = restriction_distance(&data_high, 0, 1).unwrap();
        let d_low = restriction_distance(&data_low, 0, 1).unwrap();

        assert!(
            d_high < d_low,
            "more shared sites ({:.6}) should give smaller distance than fewer shared ({:.6})",
            d_high,
            d_low
        );
    }

    // -----------------------------------------------------------------------
    // Edge cases
    // -----------------------------------------------------------------------

    #[test]
    fn test_both_no_sites_error() {
        let data = make_restriction_data(
            &[
                ("A", &[false, false, false]),
                ("B", &[false, false, false]),
            ],
            6,
        );
        let result = restriction_distance(&data, 0, 1);
        assert!(
            result.is_err(),
            "both taxa with no sites should return error"
        );
    }

    #[test]
    fn test_one_taxon_no_sites() {
        let data = make_restriction_data(
            &[
                ("A", &[true, true, true]),
                ("B", &[false, false, false]),
            ],
            6,
        );
        let result = restriction_distance(&data, 0, 1);
        // S_hat = 0 -> error (no shared sites)
        assert!(
            result.is_err(),
            "one taxon with no sites should return error (S_hat = 0)"
        );
    }

    #[test]
    fn test_out_of_range_index() {
        let data = make_restriction_data(
            &[("A", &[true, true]), ("B", &[true, false])],
            6,
        );
        let result = restriction_distance(&data, 0, 5);
        assert!(result.is_err());
    }

    #[test]
    fn test_restriction_data_validation() {
        // Empty taxa
        let result = RestrictionData::new(vec![], vec![], 6);
        assert!(result.is_err());

        // Mismatched taxa/sites count
        let result = RestrictionData::new(
            vec!["A".to_string()],
            vec![vec![true], vec![false]],
            6,
        );
        assert!(result.is_err());

        // Zero site length
        let result = RestrictionData::new(
            vec!["A".to_string()],
            vec![vec![true]],
            0,
        );
        assert!(result.is_err());

        // Unequal site vectors
        let result = RestrictionData::new(
            vec!["A".to_string(), "B".to_string()],
            vec![vec![true, false], vec![true]],
            6,
        );
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // All sites shared
    // -----------------------------------------------------------------------

    #[test]
    fn test_all_sites_shared() {
        // S_hat = 1.0 -> p_hat = 0 -> d = 0
        let data = make_restriction_data(
            &[
                ("A", &[true, true, true, true]),
                ("B", &[true, true, true, true]),
            ],
            6,
        );
        let d = restriction_distance(&data, 0, 1).unwrap();
        assert!(
            d.abs() < 1e-10,
            "all sites shared should give distance 0, got {}",
            d
        );
    }

    // -----------------------------------------------------------------------
    // Asymmetric site counts
    // -----------------------------------------------------------------------

    #[test]
    fn test_asymmetric_site_counts() {
        // A has 5 sites, B has 3 sites, 2 shared
        // S_hat = 2*2/(5+3) = 0.5
        let data = make_restriction_data(
            &[
                ("A", &[true, true, true, true, true, false, false]),
                ("B", &[true, true, false, false, false, true, false]),
            ],
            6,
        );
        let d = restriction_distance(&data, 0, 1).unwrap();

        let s_hat: f64 = 2.0 * 2.0 / (5.0 + 3.0);
        let r: f64 = 6.0;
        let p_hat: f64 = 1.0 - s_hat.powf(1.0 / r);
        let expected: f64 = -0.75 * (1.0 - (4.0 / 3.0) * p_hat).ln();

        assert!(
            (d - expected).abs() < 1e-10,
            "asymmetric: expected {:.6}, got {:.6}",
            expected,
            d
        );
    }
}

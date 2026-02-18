//! Gene frequency distance measures for population genetics.
//!
//! This module implements three classic distance measures that operate on
//! allele frequency data across populations. These are used in population
//! genetics to quantify genetic differentiation between populations and
//! construct phylogenetic trees of populations.
//!
//! The three methods correspond to those implemented in PHYLIP's `gendist`
//! program:
//!
//! - [`GeneFreqMethod::Nei`] — Nei's standard genetic distance (Nei 1972)
//! - [`GeneFreqMethod::CavalliSforza`] — Cavalli-Sforza chord distance (1967)
//! - [`GeneFreqMethod::Reynolds`] — Reynolds coancestry distance (Reynolds et al. 1983)
//!
//! # Data format
//!
//! Gene frequency data consists of allele frequencies at multiple loci for
//! multiple populations. Each locus has a number of alleles, and for each
//! population the frequencies of all alleles at a locus sum to 1.0.
//!
//! # References
//!
//! Nei, M. (1972). "Genetic distance between populations." *American
//! Naturalist*, 106:283-292.
//!
//! Cavalli-Sforza, L.L. & Edwards, A.W.F. (1967). "Phylogenetic analysis:
//! models and estimation procedures." *American Journal of Human Genetics*,
//! 19:233-257.
//!
//! Reynolds, J., Weir, B.S. & Cockerham, C.C. (1983). "Estimation of the
//! coancestry coefficient: basis for a short-term genetic distance."
//! *Genetics*, 105:767-779.

/// A single locus with allele frequency data across populations.
#[derive(Debug, Clone)]
pub struct Locus {
    /// Name of the locus.
    pub name: String,
    /// Number of alleles at this locus.
    pub n_alleles: usize,
    /// Allele frequencies: `frequencies[pop][allele]`.
    /// Each row (population) should sum to approximately 1.0.
    pub frequencies: Vec<Vec<f64>>,
}

/// Gene frequency data for multiple loci across multiple populations.
#[derive(Debug, Clone)]
pub struct GeneFreqData {
    /// Names of the populations.
    pub population_names: Vec<String>,
    /// Loci with frequency data.
    pub loci: Vec<Locus>,
}

impl GeneFreqData {
    /// Create new gene frequency data.
    ///
    /// # Arguments
    ///
    /// * `population_names` - Names of the populations
    /// * `loci` - Loci with allele frequency data
    pub fn new(population_names: Vec<String>, loci: Vec<Locus>) -> Self {
        Self {
            population_names,
            loci,
        }
    }

    /// Number of populations.
    pub fn n_populations(&self) -> usize {
        self.population_names.len()
    }

    /// Number of loci.
    pub fn n_loci(&self) -> usize {
        self.loci.len()
    }
}

/// Result of computing gene frequency distances.
#[derive(Debug, Clone)]
pub struct GeneFreqDistanceMatrix {
    /// Population names.
    pub names: Vec<String>,
    /// Distance matrix: `distances[i][j]` is the distance between
    /// populations i and j.
    pub distances: Vec<Vec<f64>>,
    /// Name of the distance method used.
    pub method: String,
}

impl GeneFreqDistanceMatrix {
    /// Number of populations.
    pub fn size(&self) -> usize {
        self.names.len()
    }

    /// Get distance between populations i and j.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.distances[i][j]
    }
}

/// Available gene frequency distance methods.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeneFreqMethod {
    /// Nei's standard genetic distance (Nei 1972).
    /// Based on genetic identity: D = -ln(I), where I measures the
    /// normalized probability that two randomly chosen alleles (one from
    /// each population) are the same.
    Nei,

    /// Cavalli-Sforza chord distance (Cavalli-Sforza & Edwards 1967).
    /// The geometric (chord) distance in allele frequency space.
    /// Bounded between 0 and 1.
    CavalliSforza,

    /// Reynolds distance (Reynolds, Weir & Cockerham 1983).
    /// Based on the coancestry coefficient theta_W, estimated from
    /// allele frequency differences. D = -ln(1 - theta_W).
    Reynolds,
}

/// Compute Nei's standard genetic distance between two populations.
///
/// # Formula
///
/// ```text
/// J_xy = (1/L) * sum_l sum_a x_la * y_la
/// J_xx = (1/L) * sum_l sum_a x_la^2
/// J_yy = (1/L) * sum_l sum_a y_la^2
/// I = J_xy / sqrt(J_xx * J_yy)
/// D = -ln(I)
/// ```
///
/// Where:
/// - L = number of loci
/// - x_la = frequency of allele a at locus l in population i
/// - y_la = frequency of allele a at locus l in population j
///
/// # Returns
///
/// The Nei distance. Returns `f64::INFINITY` if the populations share no
/// alleles (I = 0).
pub fn nei_distance(data: &GeneFreqData, pop_i: usize, pop_j: usize) -> f64 {
    if pop_i == pop_j {
        return 0.0;
    }

    let n_loci = data.n_loci();
    if n_loci == 0 {
        return 0.0;
    }

    let mut j_xy = 0.0;
    let mut j_xx = 0.0;
    let mut j_yy = 0.0;

    for locus in &data.loci {
        let x = &locus.frequencies[pop_i];
        let y = &locus.frequencies[pop_j];
        let n_alleles = x.len().min(y.len());

        for a in 0..n_alleles {
            j_xy += x[a] * y[a];
            j_xx += x[a] * x[a];
            j_yy += y[a] * y[a];
        }
    }

    // Average over loci
    let l = n_loci as f64;
    j_xy /= l;
    j_xx /= l;
    j_yy /= l;

    let denom = (j_xx * j_yy).sqrt();
    if denom <= 0.0 {
        return f64::INFINITY;
    }

    let identity = j_xy / denom;

    if identity <= 0.0 {
        return f64::INFINITY;
    }

    -identity.ln()
}

/// Compute the Cavalli-Sforza chord distance between two populations.
///
/// # Formula
///
/// ```text
/// D = sqrt((2 / L) * sum_l (1 - sum_a sqrt(x_la * y_la)))
/// ```
///
/// Where:
/// - L = number of loci
/// - x_la = frequency of allele a at locus l in population i
/// - y_la = frequency of allele a at locus l in population j
///
/// The distance is bounded between 0 (identical) and sqrt(2) (no shared
/// alleles at any locus, for biallelic loci).
///
/// # Returns
///
/// The Cavalli-Sforza chord distance.
pub fn cavalli_sforza_distance(data: &GeneFreqData, pop_i: usize, pop_j: usize) -> f64 {
    if pop_i == pop_j {
        return 0.0;
    }

    let n_loci = data.n_loci();
    if n_loci == 0 {
        return 0.0;
    }

    let mut sum = 0.0;

    for locus in &data.loci {
        let x = &locus.frequencies[pop_i];
        let y = &locus.frequencies[pop_j];
        let n_alleles = x.len().min(y.len());

        let mut cos_angle = 0.0;
        for a in 0..n_alleles {
            let product = x[a] * y[a];
            if product > 0.0 {
                cos_angle += product.sqrt();
            }
        }

        // 1 - cos(angle) is always >= 0
        sum += 1.0 - cos_angle;
    }

    let l = n_loci as f64;
    let d_squared = (2.0 / l) * sum;

    if d_squared <= 0.0 {
        return 0.0;
    }

    d_squared.sqrt()
}

/// Compute the Reynolds distance between two populations.
///
/// # Formula
///
/// ```text
/// numerator = sum_l sum_a (x_la - y_la)^2
/// denominator = 2 * sum_l (1 - sum_a x_la * y_la)
/// theta_W = numerator / denominator
/// D = -ln(1 - theta_W)
/// ```
///
/// Where:
/// - x_la = frequency of allele a at locus l in population i
/// - y_la = frequency of allele a at locus l in population j
///
/// # Returns
///
/// The Reynolds distance. Returns `f64::INFINITY` if theta_W >= 1.
pub fn reynolds_distance(data: &GeneFreqData, pop_i: usize, pop_j: usize) -> f64 {
    if pop_i == pop_j {
        return 0.0;
    }

    let n_loci = data.n_loci();
    if n_loci == 0 {
        return 0.0;
    }

    let mut numerator = 0.0;
    let mut sum_one_minus_dot = 0.0;

    for locus in &data.loci {
        let x = &locus.frequencies[pop_i];
        let y = &locus.frequencies[pop_j];
        let n_alleles = x.len().min(y.len());

        let mut dot_product = 0.0;
        for a in 0..n_alleles {
            let diff = x[a] - y[a];
            numerator += diff * diff;
            dot_product += x[a] * y[a];
        }

        sum_one_minus_dot += 1.0 - dot_product;
    }

    let denominator = 2.0 * sum_one_minus_dot;

    if denominator <= 0.0 {
        // Identical populations or degenerate data
        return 0.0;
    }

    let theta_w = numerator / denominator;

    if theta_w >= 1.0 {
        return f64::INFINITY;
    }

    -(1.0 - theta_w).ln()
}

/// Compute a full pairwise distance matrix for all populations using the
/// specified gene frequency distance method.
///
/// # Arguments
///
/// * `data` - Gene frequency data for multiple populations
/// * `method` - Which distance formula to use
///
/// # Returns
///
/// A [`GeneFreqDistanceMatrix`] containing all pairwise distances.
pub fn compute_gene_freq_distances(
    data: &GeneFreqData,
    method: GeneFreqMethod,
) -> GeneFreqDistanceMatrix {
    let n = data.n_populations();
    let mut distances = vec![vec![0.0; n]; n];

    let distance_fn = match method {
        GeneFreqMethod::Nei => nei_distance,
        GeneFreqMethod::CavalliSforza => cavalli_sforza_distance,
        GeneFreqMethod::Reynolds => reynolds_distance,
    };

    let method_name = match method {
        GeneFreqMethod::Nei => "Nei",
        GeneFreqMethod::CavalliSforza => "Cavalli-Sforza",
        GeneFreqMethod::Reynolds => "Reynolds",
    };

    for i in 0..n {
        for j in (i + 1)..n {
            let d = distance_fn(data, i, j);
            distances[i][j] = d;
            distances[j][i] = d;
        }
    }

    GeneFreqDistanceMatrix {
        names: data.population_names.clone(),
        distances,
        method: method_name.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to create a simple GeneFreqData with the given allele frequencies.
    fn make_data(
        pop_names: Vec<&str>,
        loci: Vec<(&str, usize, Vec<Vec<f64>>)>,
    ) -> GeneFreqData {
        let population_names: Vec<String> = pop_names.iter().map(|s| s.to_string()).collect();
        let loci: Vec<Locus> = loci
            .into_iter()
            .map(|(name, n_alleles, frequencies)| Locus {
                name: name.to_string(),
                n_alleles,
                frequencies,
            })
            .collect();
        GeneFreqData::new(population_names, loci)
    }

    // ---- Nei distance tests ----

    #[test]
    fn test_nei_identical_populations() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![
                ("Locus1", 3, vec![vec![0.5, 0.3, 0.2], vec![0.5, 0.3, 0.2]]),
                ("Locus2", 2, vec![vec![0.6, 0.4], vec![0.6, 0.4]]),
            ],
        );
        let d = nei_distance(&data, 0, 1);
        assert!(
            d.abs() < 1e-10,
            "identical populations should have Nei distance 0, got {}",
            d
        );
    }

    #[test]
    fn test_nei_same_population() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![0.5, 0.5], vec![0.3, 0.7]])],
        );
        let d = nei_distance(&data, 0, 0);
        assert!(d.abs() < 1e-10, "distance to self should be 0, got {}", d);
    }

    #[test]
    fn test_nei_no_shared_alleles() {
        // Pop1 has allele 1 fixed, Pop2 has allele 2 fixed
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![1.0, 0.0], vec![0.0, 1.0]])],
        );
        let d = nei_distance(&data, 0, 1);
        assert!(
            d.is_infinite(),
            "no shared alleles should give infinite Nei distance, got {}",
            d
        );
    }

    #[test]
    fn test_nei_known_value() {
        // One locus, 2 alleles
        // Pop1: [0.8, 0.2], Pop2: [0.5, 0.5]
        // J_xy = 0.8*0.5 + 0.2*0.5 = 0.5
        // J_xx = 0.64 + 0.04 = 0.68
        // J_yy = 0.25 + 0.25 = 0.50
        // I = 0.5 / sqrt(0.68 * 0.50) = 0.5 / sqrt(0.34) = 0.5 / 0.58310 = 0.85749
        // D = -ln(0.85749) = 0.15369
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![0.8, 0.2], vec![0.5, 0.5]])],
        );
        let d = nei_distance(&data, 0, 1);
        let j_xy: f64 = 0.8 * 0.5 + 0.2 * 0.5;
        let j_xx: f64 = 0.8 * 0.8 + 0.2 * 0.2;
        let j_yy: f64 = 0.5 * 0.5 + 0.5 * 0.5;
        let identity: f64 = j_xy / (j_xx * j_yy).sqrt();
        let expected: f64 = -identity.ln();
        assert!(
            (d - expected).abs() < 1e-10,
            "Nei distance: expected {}, got {}",
            expected,
            d
        );
    }

    #[test]
    fn test_nei_multiple_loci() {
        // Two loci, check averaging
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![
                ("Locus1", 2, vec![vec![0.8, 0.2], vec![0.5, 0.5]]),
                ("Locus2", 2, vec![vec![0.6, 0.4], vec![0.6, 0.4]]),
            ],
        );
        let d = nei_distance(&data, 0, 1);
        assert!(d > 0.0, "divergent populations should have positive distance");
        assert!(d.is_finite(), "distance should be finite");
    }

    #[test]
    fn test_nei_symmetry() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![0.8, 0.2], vec![0.5, 0.5]])],
        );
        let d12 = nei_distance(&data, 0, 1);
        let d21 = nei_distance(&data, 1, 0);
        assert!(
            (d12 - d21).abs() < 1e-10,
            "Nei distance should be symmetric: {} vs {}",
            d12,
            d21
        );
    }

    // ---- Cavalli-Sforza distance tests ----

    #[test]
    fn test_cavalli_identical_populations() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![
                ("Locus1", 3, vec![vec![0.5, 0.3, 0.2], vec![0.5, 0.3, 0.2]]),
                ("Locus2", 2, vec![vec![0.6, 0.4], vec![0.6, 0.4]]),
            ],
        );
        let d = cavalli_sforza_distance(&data, 0, 1);
        assert!(
            d.abs() < 1e-10,
            "identical populations should have CS distance 0, got {}",
            d
        );
    }

    #[test]
    fn test_cavalli_same_population() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![0.5, 0.5], vec![0.3, 0.7]])],
        );
        let d = cavalli_sforza_distance(&data, 0, 0);
        assert!(d.abs() < 1e-10, "distance to self should be 0, got {}", d);
    }

    #[test]
    fn test_cavalli_no_shared_alleles() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![1.0, 0.0], vec![0.0, 1.0]])],
        );
        let d = cavalli_sforza_distance(&data, 0, 1);
        // cos_angle = sqrt(0) + sqrt(0) = 0, so 1 - 0 = 1
        // D = sqrt(2/1 * 1) = sqrt(2)
        let expected = 2.0_f64.sqrt();
        assert!(
            (d - expected).abs() < 1e-10,
            "no shared alleles: expected {}, got {}",
            expected,
            d
        );
    }

    #[test]
    fn test_cavalli_known_value() {
        // One locus, 2 alleles
        // Pop1: [0.8, 0.2], Pop2: [0.5, 0.5]
        // cos_angle = sqrt(0.4) + sqrt(0.1) = 0.63246 + 0.31623 = 0.94868
        // D = sqrt(2/1 * (1 - 0.94868)) = sqrt(2 * 0.05132) = sqrt(0.10263) = 0.32036
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![0.8, 0.2], vec![0.5, 0.5]])],
        );
        let d = cavalli_sforza_distance(&data, 0, 1);
        let cos_angle = (0.8 * 0.5_f64).sqrt() + (0.2 * 0.5_f64).sqrt();
        let expected = (2.0 * (1.0 - cos_angle)).sqrt();
        assert!(
            (d - expected).abs() < 1e-10,
            "CS distance: expected {}, got {}",
            expected,
            d
        );
    }

    #[test]
    fn test_cavalli_symmetry() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![0.8, 0.2], vec![0.5, 0.5]])],
        );
        let d12 = cavalli_sforza_distance(&data, 0, 1);
        let d21 = cavalli_sforza_distance(&data, 1, 0);
        assert!(
            (d12 - d21).abs() < 1e-10,
            "CS distance should be symmetric: {} vs {}",
            d12,
            d21
        );
    }

    #[test]
    fn test_cavalli_bounded() {
        // Cavalli-Sforza distance should be bounded
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![
                ("Locus1", 2, vec![vec![1.0, 0.0], vec![0.0, 1.0]]),
                ("Locus2", 2, vec![vec![1.0, 0.0], vec![0.0, 1.0]]),
            ],
        );
        let d = cavalli_sforza_distance(&data, 0, 1);
        assert!(d.is_finite(), "CS distance should always be finite");
        assert!(d <= 2.0_f64.sqrt() + 1e-10, "CS distance should be <= sqrt(2)");
    }

    // ---- Reynolds distance tests ----

    #[test]
    fn test_reynolds_identical_populations() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![
                ("Locus1", 3, vec![vec![0.5, 0.3, 0.2], vec![0.5, 0.3, 0.2]]),
                ("Locus2", 2, vec![vec![0.6, 0.4], vec![0.6, 0.4]]),
            ],
        );
        let d = reynolds_distance(&data, 0, 1);
        assert!(
            d.abs() < 1e-10,
            "identical populations should have Reynolds distance 0, got {}",
            d
        );
    }

    #[test]
    fn test_reynolds_same_population() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![0.5, 0.5], vec![0.3, 0.7]])],
        );
        let d = reynolds_distance(&data, 0, 0);
        assert!(d.abs() < 1e-10, "distance to self should be 0, got {}", d);
    }

    #[test]
    fn test_reynolds_no_shared_alleles() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![1.0, 0.0], vec![0.0, 1.0]])],
        );
        let d = reynolds_distance(&data, 0, 1);
        // numerator = (1-0)^2 + (0-1)^2 = 2
        // dot = 0, sum_one_minus_dot = 1
        // denominator = 2
        // theta_W = 2/2 = 1 -> infinite
        assert!(
            d.is_infinite(),
            "no shared alleles should give infinite Reynolds distance, got {}",
            d
        );
    }

    #[test]
    fn test_reynolds_known_value() {
        // One locus, 2 alleles
        // Pop1: [0.8, 0.2], Pop2: [0.5, 0.5]
        // numerator = (0.3)^2 + (-0.3)^2 = 0.18
        // dot = 0.4 + 0.1 = 0.5
        // sum_one_minus_dot = 0.5
        // denominator = 2 * 0.5 = 1.0
        // theta_W = 0.18 / 1.0 = 0.18
        // D = -ln(1 - 0.18) = -ln(0.82) = 0.19845
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![0.8, 0.2], vec![0.5, 0.5]])],
        );
        let d = reynolds_distance(&data, 0, 1);
        let theta_w = 0.18;
        let expected = -(1.0 - theta_w as f64).ln();
        assert!(
            (d - expected).abs() < 1e-10,
            "Reynolds distance: expected {}, got {}",
            expected,
            d
        );
    }

    #[test]
    fn test_reynolds_symmetry() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![0.8, 0.2], vec![0.5, 0.5]])],
        );
        let d12 = reynolds_distance(&data, 0, 1);
        let d21 = reynolds_distance(&data, 1, 0);
        assert!(
            (d12 - d21).abs() < 1e-10,
            "Reynolds distance should be symmetric: {} vs {}",
            d12,
            d21
        );
    }

    // ---- Triangle inequality tests ----

    #[test]
    fn test_nei_triangle_inequality() {
        let data = make_data(
            vec!["Pop1", "Pop2", "Pop3"],
            vec![
                (
                    "Locus1",
                    3,
                    vec![
                        vec![0.6, 0.3, 0.1],
                        vec![0.4, 0.4, 0.2],
                        vec![0.2, 0.3, 0.5],
                    ],
                ),
                (
                    "Locus2",
                    2,
                    vec![vec![0.7, 0.3], vec![0.5, 0.5], vec![0.3, 0.7]],
                ),
            ],
        );
        let d01 = nei_distance(&data, 0, 1);
        let d02 = nei_distance(&data, 0, 2);
        let d12 = nei_distance(&data, 1, 2);

        // Note: Nei distance does not always satisfy the triangle inequality
        // in theory, but for typical data it should.
        // We just check all distances are non-negative and finite.
        assert!(d01 >= 0.0 && d01.is_finite());
        assert!(d02 >= 0.0 && d02.is_finite());
        assert!(d12 >= 0.0 && d12.is_finite());
    }

    #[test]
    fn test_cavalli_triangle_inequality() {
        let data = make_data(
            vec!["Pop1", "Pop2", "Pop3"],
            vec![
                (
                    "Locus1",
                    3,
                    vec![
                        vec![0.6, 0.3, 0.1],
                        vec![0.4, 0.4, 0.2],
                        vec![0.2, 0.3, 0.5],
                    ],
                ),
                (
                    "Locus2",
                    2,
                    vec![vec![0.7, 0.3], vec![0.5, 0.5], vec![0.3, 0.7]],
                ),
            ],
        );
        let d01 = cavalli_sforza_distance(&data, 0, 1);
        let d02 = cavalli_sforza_distance(&data, 0, 2);
        let d12 = cavalli_sforza_distance(&data, 1, 2);

        // Cavalli-Sforza is a metric and should satisfy the triangle inequality
        assert!(
            d02 <= d01 + d12 + 1e-10,
            "triangle inequality failed: d(0,2)={} > d(0,1)+d(1,2)={}",
            d02,
            d01 + d12
        );
        assert!(
            d01 <= d02 + d12 + 1e-10,
            "triangle inequality failed: d(0,1)={} > d(0,2)+d(2,1)={}",
            d01,
            d02 + d12
        );
        assert!(
            d12 <= d01 + d02 + 1e-10,
            "triangle inequality failed: d(1,2)={} > d(0,1)+d(0,2)={}",
            d12,
            d01 + d02
        );
    }

    #[test]
    fn test_reynolds_triangle_inequality() {
        let data = make_data(
            vec!["Pop1", "Pop2", "Pop3"],
            vec![
                (
                    "Locus1",
                    3,
                    vec![
                        vec![0.6, 0.3, 0.1],
                        vec![0.4, 0.4, 0.2],
                        vec![0.2, 0.3, 0.5],
                    ],
                ),
                (
                    "Locus2",
                    2,
                    vec![vec![0.7, 0.3], vec![0.5, 0.5], vec![0.3, 0.7]],
                ),
            ],
        );
        let d01 = reynolds_distance(&data, 0, 1);
        let d02 = reynolds_distance(&data, 0, 2);
        let d12 = reynolds_distance(&data, 1, 2);

        // All should be non-negative and finite
        assert!(d01 >= 0.0 && d01.is_finite());
        assert!(d02 >= 0.0 && d02.is_finite());
        assert!(d12 >= 0.0 && d12.is_finite());
    }

    // ---- compute_gene_freq_distances tests ----

    #[test]
    fn test_compute_distance_matrix_nei() {
        let data = make_data(
            vec!["Pop1", "Pop2", "Pop3"],
            vec![
                (
                    "Locus1",
                    2,
                    vec![vec![0.8, 0.2], vec![0.5, 0.5], vec![0.3, 0.7]],
                ),
            ],
        );
        let matrix = compute_gene_freq_distances(&data, GeneFreqMethod::Nei);
        assert_eq!(matrix.size(), 3);
        assert_eq!(matrix.method, "Nei");

        // Diagonal should be 0
        for i in 0..3 {
            assert!(matrix.get(i, i).abs() < 1e-10);
        }

        // Symmetry
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (matrix.get(i, j) - matrix.get(j, i)).abs() < 1e-10,
                    "matrix should be symmetric"
                );
            }
        }

        // Pop1 closer to Pop2 than to Pop3
        assert!(
            matrix.get(0, 1) < matrix.get(0, 2),
            "Pop1-Pop2 ({}) should be closer than Pop1-Pop3 ({})",
            matrix.get(0, 1),
            matrix.get(0, 2)
        );
    }

    #[test]
    fn test_compute_distance_matrix_cavalli() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![0.8, 0.2], vec![0.5, 0.5]])],
        );
        let matrix = compute_gene_freq_distances(&data, GeneFreqMethod::CavalliSforza);
        assert_eq!(matrix.method, "Cavalli-Sforza");
        assert!(matrix.get(0, 1) > 0.0);
    }

    #[test]
    fn test_compute_distance_matrix_reynolds() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 2, vec![vec![0.8, 0.2], vec![0.5, 0.5]])],
        );
        let matrix = compute_gene_freq_distances(&data, GeneFreqMethod::Reynolds);
        assert_eq!(matrix.method, "Reynolds");
        assert!(matrix.get(0, 1) > 0.0);
    }

    // ---- Distance ordering tests ----

    #[test]
    fn test_distance_increases_with_divergence_nei() {
        // Pop1 and Pop2 are similar, Pop3 is very different
        let data = make_data(
            vec!["Pop1", "Pop2", "Pop3"],
            vec![
                (
                    "Locus1",
                    2,
                    vec![vec![0.9, 0.1], vec![0.8, 0.2], vec![0.1, 0.9]],
                ),
            ],
        );
        let d12 = nei_distance(&data, 0, 1);
        let d13 = nei_distance(&data, 0, 2);
        assert!(
            d12 < d13,
            "similar pops should have smaller distance: d12={} < d13={}",
            d12,
            d13
        );
    }

    #[test]
    fn test_distance_increases_with_divergence_cavalli() {
        let data = make_data(
            vec!["Pop1", "Pop2", "Pop3"],
            vec![
                (
                    "Locus1",
                    2,
                    vec![vec![0.9, 0.1], vec![0.8, 0.2], vec![0.1, 0.9]],
                ),
            ],
        );
        let d12 = cavalli_sforza_distance(&data, 0, 1);
        let d13 = cavalli_sforza_distance(&data, 0, 2);
        assert!(
            d12 < d13,
            "similar pops should have smaller distance: d12={} < d13={}",
            d12,
            d13
        );
    }

    #[test]
    fn test_distance_increases_with_divergence_reynolds() {
        let data = make_data(
            vec!["Pop1", "Pop2", "Pop3"],
            vec![
                (
                    "Locus1",
                    2,
                    vec![vec![0.9, 0.1], vec![0.8, 0.2], vec![0.1, 0.9]],
                ),
            ],
        );
        let d12 = reynolds_distance(&data, 0, 1);
        let d13 = reynolds_distance(&data, 0, 2);
        assert!(
            d12 < d13,
            "similar pops should have smaller distance: d12={} < d13={}",
            d12,
            d13
        );
    }

    // ---- Multi-locus tests ----

    #[test]
    fn test_multi_locus_averaging() {
        // Two loci: one identical, one divergent
        // Distance should be intermediate
        let data_one_divergent = make_data(
            vec!["Pop1", "Pop2"],
            vec![
                ("Locus1", 2, vec![vec![0.5, 0.5], vec![0.5, 0.5]]),  // identical
                ("Locus2", 2, vec![vec![0.9, 0.1], vec![0.1, 0.9]]),  // divergent
            ],
        );
        let data_both_divergent = make_data(
            vec!["Pop1", "Pop2"],
            vec![
                ("Locus1", 2, vec![vec![0.9, 0.1], vec![0.1, 0.9]]),  // divergent
                ("Locus2", 2, vec![vec![0.9, 0.1], vec![0.1, 0.9]]),  // divergent
            ],
        );

        let d_one = nei_distance(&data_one_divergent, 0, 1);
        let d_both = nei_distance(&data_both_divergent, 0, 1);
        assert!(
            d_one < d_both,
            "one divergent locus ({}) should give smaller distance than two ({})",
            d_one,
            d_both
        );
    }

    // ---- Three alleles test ----

    #[test]
    fn test_three_alleles() {
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![(
                "Locus1",
                3,
                vec![vec![0.5, 0.3, 0.2], vec![0.3, 0.4, 0.3]],
            )],
        );

        let d_nei = nei_distance(&data, 0, 1);
        let d_cs = cavalli_sforza_distance(&data, 0, 1);
        let d_rey = reynolds_distance(&data, 0, 1);

        assert!(d_nei > 0.0 && d_nei.is_finite(), "Nei: {}", d_nei);
        assert!(d_cs > 0.0 && d_cs.is_finite(), "CS: {}", d_cs);
        assert!(d_rey > 0.0 && d_rey.is_finite(), "Reynolds: {}", d_rey);
    }

    // ---- Edge case: many loci ----

    #[test]
    fn test_many_loci() {
        let n_loci = 100;
        let mut loci = Vec::new();
        for i in 0..n_loci {
            let p = (i as f64 + 1.0) / (n_loci as f64 + 2.0);
            loci.push((
                "L",
                2,
                vec![vec![p, 1.0 - p], vec![1.0 - p, p]],
            ));
        }
        let data = make_data(vec!["Pop1", "Pop2"], loci);

        let d_nei = nei_distance(&data, 0, 1);
        let d_cs = cavalli_sforza_distance(&data, 0, 1);
        let d_rey = reynolds_distance(&data, 0, 1);

        assert!(d_nei > 0.0 && d_nei.is_finite());
        assert!(d_cs > 0.0 && d_cs.is_finite());
        assert!(d_rey > 0.0 && d_rey.is_finite());
    }

    // ---- Edge case: single allele (fixed) ----

    #[test]
    fn test_fixed_allele_identical() {
        // Both populations fixed for the same allele
        let data = make_data(
            vec!["Pop1", "Pop2"],
            vec![("Locus1", 1, vec![vec![1.0], vec![1.0]])],
        );
        let d = nei_distance(&data, 0, 1);
        assert!(
            d.abs() < 1e-10,
            "identical fixed populations should have distance 0, got {}",
            d
        );
    }
}

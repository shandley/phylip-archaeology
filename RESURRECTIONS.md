# RESURRECTIONS.md — Five Lost Phylogenetics Tools, Resurrected

## Overview

Felsenstein's [software catalog](https://evolution.genetics.washington.edu/phylip/software.html) once listed over 400 phylogenetics tools. Many have vanished: websites gone dark, downloads broken, source code lost. These tools embodied important algorithmic ideas that risk being forgotten.

We resurrected five of them as modern Rust reimplementations within `phylip-rs`, preserving their core algorithms with zero external dependencies. Each implementation was validated with comprehensive test suites and documented with full algorithmic detail.

| Tool | Author(s) | Year | Algorithm | Lines | Tests |
|------|-----------|------|-----------|-------|-------|
| **TipDate** | Rambaut | 2000 | Serial-sample molecular clock ML | 1,307 | 17 |
| **scaleboot/AU** | Shimodaira | 2002 | Approximately Unbiased bootstrap test | 1,237 | 20 |
| **DIVA** | Ronquist | 1997 | Dispersal-Vicariance biogeographic analysis | 1,695 | 49 |
| **TREEMAP** | Page | 1994 | Host-parasite tree reconciliation | 955 | 18 |
| **PLATO** | Grassly & Holmes | 1997 | Sliding-window recombination detection | 1,106 | 26 |
| **Total** | | | | **6,300** | **130** |

---

## 1. TipDate — Serial-Sample Molecular Clock

**Original**: TipDate v1.2 by Andrew Rambaut (University of Oxford, 2000). The original Java application is no longer available from its original distribution site.

**What it does**: Estimates the rate of molecular evolution from sequences sampled at different calendar dates. When viral isolates or ancient DNA samples span years or decades, the sampling times act as calibration points that anchor the molecular clock without requiring fossil evidence.

**Algorithm**:
1. Root the tree and assign known sampling dates to tips.
2. Initialize the substitution rate from root-to-tip regression (a linear fit of genetic divergence against sampling date).
3. Set branch lengths as `r * |t_parent - t_child|`, where `r` is the rate and times are estimated for internal nodes.
4. Optimize `r` and internal node dates using golden section search to maximize the log-likelihood under Felsenstein's pruning algorithm.
5. Compare three nested models via likelihood ratio test:
   - **Model 0**: Free branch lengths (no clock)
   - **Model 1**: Strict molecular clock (all tips contemporaneous)
   - **Model 2**: Dated-tip clock (SRDT) — the TipDate model

**Key insight**: The LRT between Model 2 and Model 1 tests whether the temporal signal in the data is statistically significant. A significant result means the sequences contain enough clock-like signal to estimate divergence times from tip dates alone.

**Implementation**: `phylip-rs/src/likelihood/tipdate.rs` (1,307 lines)

**Public API**:
- `tipdate_optimize()` — Estimate rate and node dates
- `tipdate_likelihood()` — Compute likelihood for a given rate
- `tipdate_lrt()` — Likelihood ratio test (dated-tip vs. strict clock)

**Tests validated**:
- Two-taxon analytical case (known branch length = rate * time difference)
- Five-taxon heterochronous tree with temporal signal
- Contemporaneous tips reduce to standard molecular clock
- Root date precedes all tip dates
- Temporal ordering of node dates respected
- LRT properties (non-negative statistic, proper df)

**Why it matters**: TipDate's approach became the foundation for Bayesian methods like BEAST. The core idea — that sampling times provide free calibration — revolutionized viral phylodynamics and ancient DNA studies.

**Reference**: Rambaut, A. (2000). Estimating the rate of molecular evolution: incorporating non-contemporaneous sequences into maximum likelihood phylogenies. *Bioinformatics*, 16(4), 395-399.

---

## 2. scaleboot / AU Test — Approximately Unbiased Hypothesis Testing

**Original**: scaleboot R package by Hidetoshi Shimodaira (Tokyo Institute of Technology / Osaka University, 2002). The R package exists but the underlying algorithm is complex enough that standalone implementations are rare; the original CONSEL C implementation has become difficult to build.

**What it does**: Tests whether a set of candidate phylogenetic trees can be statistically distinguished. Standard bootstrap proportions (BP) are biased when the set of trees being compared was selected *because* they were good — a form of selection bias. The AU test corrects this bias through multiscale bootstrap resampling.

**Algorithm**:
1. Compute per-site log-likelihoods for each candidate tree.
2. For each scale factor `r` in {0.5, 0.6, ..., 1.4}, resample `floor(r * n)` sites with replacement and record which tree "wins" (highest total log-likelihood).
3. For each tree `k`, compute the bootstrap proportion `BP_k(r)` at each scale.
4. Transform to rejection z-values: `z_k(r) = Phi^{-1}(1 - BP_k(r))`.
5. Fit the linear model `z_k(r) = d1_k * sqrt(r) + d2_k / sqrt(r)` by least squares.
6. The AU p-value is `p_AU_k = 1 - Phi(d1_k + d2_k)`.

The two-parameter model captures both the "bias" component (d2, which corrects for selection) and the "variance" component (d1, which reflects genuine signal). When d2 = 0, the AU test reduces to the standard BP.

**Also implements**:
- **KH test** (Kishino-Hasegawa 1989): Pairwise comparison of two trees using the variance of site log-likelihood differences.
- **SH test** (Shimodaira-Hasegawa 1999): Conservative multiple-comparison correction using the maximum of centered test statistics.
- **Confidence set construction**: The set of trees not rejected at a given significance level.

**Implementation**: `phylip-rs/src/bootstrap/au_test.rs` (1,237 lines)

**Public API**:
- `au_test()` — Full AU test from site log-likelihood matrix
- `au_test_from_trees()` — Convenience wrapper: compute site log-likelihoods from trees, alignment, and model, then run AU test

**Mathematical functions implemented from scratch**:
- `normal_cdf()` — Abramowitz & Stegun rational approximation
- `normal_quantile()` — Acklam (2004) rational approximation with Newton-Raphson refinement

**Tests validated**:
- Best tree receives highest AU p-value, worst receives lowest
- AU p-values sum approximately to the number of trees
- BP at scale r=1 matches standard bootstrap proportion
- KH test: best tree p-value = 1.0 (by definition)
- SH test more conservative than KH test
- Confidence set includes all trees not rejected at alpha
- Reproducibility with fixed random seed
- Normal CDF/quantile round-trip accuracy

**Why it matters**: The AU test solved a fundamental statistical problem in phylogenetics: how to perform valid hypothesis testing when the hypotheses (tree topologies) were themselves discovered from the data. It remains the gold standard for tree topology testing.

**References**:
- Shimodaira, H. (2002). An approximately unbiased test of phylogenetic tree selection. *Systematic Biology*, 51(3), 492-508.
- Shimodaira, H. & Hasegawa, M. (1999). Multiple comparisons of log-likelihoods with applications to phylogenetic inference. *Molecular Biology and Evolution*, 16, 1114-1116.
- Kishino, H. & Hasegawa, M. (1989). Evaluation of the maximum likelihood estimate of the evolutionary tree topologies. *Journal of Molecular Evolution*, 29, 170-179.

---

## 3. DIVA — Dispersal-Vicariance Analysis

**Original**: DIVA v1.1 by Fredrik Ronquist (University of Uppsala, 1996-1997). The original DOS executable and its website are no longer available.

**What it does**: Reconstructs ancestral geographic distributions on a phylogenetic tree. Given the current distributions of species across geographic areas, DIVA infers where their ancestors lived and what biogeographic events (dispersal, vicariance, extinction) shaped the present pattern.

**Algorithm**:
DIVA uses bottom-up dynamic programming on area bit-vectors:

1. **Leaf nodes**: The observed distribution has cost 0; all others cost infinity.
2. **Internal nodes**: For each candidate ancestral area set *A*, enumerate all ways to partition *A* between two daughter lineages:
   - **Vicariance** (cost 0): Split *A* into disjoint subsets — one daughter gets areas {1,2}, the other gets {3,4}.
   - **Duplication** (cost 0): Both daughters inherit the full set *A* (within-area speciation).
   - **Dispersal** (cost 1 per area gained): Bridge the gap between inherited and optimal daughter distribution.
   - **Extinction** (cost 1 per area lost): Account for range contraction.
3. **Root**: Select the distribution with minimum total cost. Multiple equally optimal solutions may exist.

**Key data structure**: `AreaSet` — a `u32` bit-vector where bit *i* represents presence in area *i*. Supports up to 20 areas. Set operations (union, intersection, difference, subset test) are single bitwise instructions.

**Complexity**: O(3^k) per node for *k* areas (from enumerating all bipartitions of each area set via the binomial theorem). Practical for k <= 10-12.

**Implementation**: `phylip-rs/src/biogeography/diva.rs` (1,663 lines) + `phylip-rs/src/biogeography/mod.rs` (32 lines)

**Public API**:
- `diva_optimize()` — Run DIVA with default event costs
- `diva_optimize_with_costs()` — Run DIVA with custom event costs
- `AreaSet` — Bit-vector type with full set algebra
- `DivaResult` — Optimal ancestral distributions and total cost

**Tests validated**:
- Perfect vicariance scenarios (2, 3, 4 areas) yield cost 0
- Dispersal-required scenarios yield correct non-zero costs
- Max range constraint limits ancestral area set sizes
- Gondwanan vicariance scenario (continental fragmentation)
- Island colonization scenario (stepping-stone dispersal)
- Symmetric trees with symmetric area assignments
- Error handling (missing leaf areas, non-binary trees, invalid constraints)
- AreaSet operations (union, intersection, difference, subset, iteration)

**Why it matters**: DIVA formalized the intuition that vicariance (geographic splitting) is the "null" mode of speciation, and that dispersal and extinction require explanation. It bridged the gap between narrative biogeography and quantitative cladistic methods, and its event-based framework influenced all subsequent parametric biogeographic methods (DEC, BioGeoBEARS).

**Reference**: Ronquist, F. (1997). Dispersal-vicariance analysis: a new approach to the quantification of historical biogeography. *Systematic Biology*, 46:195-203.

---

## 4. TREEMAP — Host-Parasite Tree Reconciliation

**Original**: TREEMAP v1.0 by Rod Page (University of Glasgow, 1994). Originally a Macintosh Classic application; the original binary and its distribution are no longer available, though later versions (TREEMAP 2, 3) existed in various states.

**What it does**: Maps a parasite (or symbiont) phylogeny onto a host phylogeny to reconstruct the coevolutionary history. Given leaf-level associations (which parasite infects which host), the algorithm determines the minimum-cost mapping of internal parasite nodes to host nodes, classifying each divergence event.

**Algorithm**:
LCA-based reconciliation by postorder traversal of the parasite tree:

1. **Leaf mapping**: Each parasite leaf maps to its known host via the association `phi`.
2. **Internal nodes**: For parasite node `p` with children `p1`, `p2`:
   - Compute `sigma(p) = LCA_H(sigma(p1), sigma(p2))` — the least common ancestor in the host tree.
   - **Cospeciation**: If `sigma(p1)` and `sigma(p2)` descend from different children of `sigma(p)` — host and parasite speciated together.
   - **Duplication**: If both map to the same subtree of `sigma(p)` — parasite speciated within a single host lineage.
3. **Sorting events**: Between `sigma(p)` and `sigma(child)`, the number of "missing" host speciation events is `depth(sigma(child)) - depth(sigma(p)) - 1`. Each represents a sorting (lineage loss) event.
4. **Host-switching**: Detected when the parasite mapping is inconsistent with strict vertical inheritance.

**Event costs** (Page 1994 defaults):

| Event | Cost | Interpretation |
|-------|------|----------------|
| Cospeciation | 0 | Congruent speciation — the expected pattern |
| Duplication | 1 | Intra-host speciation — requires explanation |
| Sorting | 1 | Lineage loss — parasite failed to persist |
| Host-switch | 2 | Lateral transfer — most costly event |

**Implementation**: `phylip-rs/src/reconciliation/treemap.rs` (903 lines) + `phylip-rs/src/reconciliation/mod.rs` (52 lines)

**Public API**:
- `reconcile()` — Perform full reconciliation analysis
- `ReconciliationResult` — Mapping, events, costs, event counts
- `EventCosts` — Configurable event cost weights
- `lca()` — Least common ancestor computation

**Tests validated**:
- Perfect cospeciation (3 and 4 taxa) — all events are cospeciations, cost = 0
- Single duplication event — detected and costed correctly
- All-duplication scenario — parasite tree fully incongruent with host
- Sorting events counted correctly based on depth differences
- Custom event costs applied correctly
- Mismatched tree sizes handled
- Error cases (empty tree, missing leaf mapping, host not found)
- LCA correctness (siblings, across subtrees, leaf-ancestor pairs)

**Why it matters**: TREEMAP pioneered the quantitative analysis of coevolution by treating it as a tree-mapping problem. The event-based framework (cospeciation, duplication, sorting, host-switching) became the standard vocabulary for studying host-parasite, gene-species, and area-phylogeny associations. Its ideas directly influenced Jane, Notung, and other reconciliation tools.

**References**:
- Page, R.D.M. (1994). Maps between trees and cladistic analysis of historical associations among genes, organisms, and areas. *Systematic Biology*, 43:58-77.
- Page, R.D.M. (1994). Parallel phylogenies: reconstructing the history of host-parasite assemblages. *Cladistics*, 10:155-173.

---

## 5. PLATO — Partial Likelihoods Assessed Through Optimisation

**Original**: PLATO by Nick Grassly & Eddie Holmes (University of Oxford, 1997). The original C implementation is no longer available from its distribution site.

**What it does**: Scans a multiple sequence alignment with a sliding window to detect regions where the phylogenetic signal deviates from the overall tree. Anomalous regions may indicate recombination breakpoints, gene conversion, variation in selective pressure, or alignment errors.

**Algorithm**:
1. Fit a substitution model and tree to the full alignment (using Felsenstein's pruning algorithm).
2. Compute per-site log-likelihoods under the fitted model.
3. Slide a window of size W with step size S across the alignment. For each window position j, compute:
   - Partial log-likelihood: `PLL(j) = sum of site_lnL[j..j+W]`
   - Expected PLL under homogeneity: `E[PLL] = (W/N) * total_lnL`
   - Deviation: `delta(j) = PLL(j) - E[PLL]`
4. Compute z-scores from the distribution of delta values across all windows.
5. Flag windows where `|z| > threshold` as anomalous (default threshold = 3.0).

**Parametric bootstrap for formal testing**:
1. Simulate `B` alignments under the fitted null model on the fitted tree.
2. For each simulated alignment, run the same sliding-window scan and record `min(delta)`.
3. The p-value is the proportion of simulated `min(delta)` values at least as extreme as the observed `min(delta)`.

**Sequence simulation**: Implements full Markov chain simulation along the tree:
- Draws root state from equilibrium frequencies
- Simulates substitutions along each branch using transition probability matrices
- Supports JC69, K2P, F81, F84 models via the `SubstitutionModel` trait

**Implementation**: `phylip-rs/src/likelihood/plato.rs` (1,106 lines)

**Public API**:
- `plato_scan()` — Run the sliding-window scan
- `plato_parametric_bootstrap()` — Formal significance testing via parametric bootstrap
- `simulate_alignment()` — Simulate sequences under a model on a tree
- `PlatoResult` — Window positions, PLLs, z-scores, anomalous windows
- `BootstrapResult` — Bootstrap p-value and null distribution

**Tests validated**:
- Homogeneous alignment produces no anomalous windows
- Concatenated alignment (two different trees) detects breakpoint region
- Window positions computed correctly for various step sizes
- Z-scores have mean ~0 under the null model
- Mean window PLL close to expected value
- Parametric bootstrap returns correct number of replicates
- Bootstrap p-value ~1.0 under the null (no recombination)
- Sequence simulation produces valid bases at correct frequencies
- Short branches produce low divergence (validates Markov simulation)
- Error handling (invalid window size, step size, taxa mismatch)

**Why it matters**: PLATO was one of the first tools to use likelihood-based methods for recombination detection. Its sliding-window approach influenced later methods like GARD and PhyML-based scanning approaches. The parametric bootstrap framework provided proper statistical testing for what had previously been ad hoc visual inspection of likelihood surfaces.

**Reference**: Grassly, N. C. & Holmes, E. C. (1997). A likelihood method for the detection of selection and recombination using nucleotide sequences. *Molecular Biology and Evolution*, 14(3), 239-247.

---

## Cross-Cutting Observations

### Shared Mathematical Infrastructure

All five tools reuse the same zero-dependency mathematical primitives:

- **Felsenstein pruning algorithm** (TipDate, PLATO, AU test): The workhorse of phylogenetic likelihood, computing P(data|tree, model) in O(n * s * k^2) time.
- **Golden section search** (TipDate): One-dimensional optimization without derivatives.
- **Normal CDF/quantile** (AU test): Implemented from scratch using Abramowitz-Stegun and Acklam approximations.
- **Chi-squared survival function** (TipDate): For likelihood ratio test p-values.
- **Dynamic programming on trees** (DIVA, TREEMAP): Bottom-up postorder traversal with memoization.
- **Bit-vector set operations** (DIVA): Efficient representation of combinatorial state spaces.
- **Markov chain simulation** (PLATO): Generating synthetic data under known models.
- **Pseudorandom number generation** (AU test, PLATO): LCG-based RNG for reproducible resampling.

### Algorithmic Patterns

Three recurring patterns emerged across these tools:

1. **Likelihood as a lens**: TipDate, PLATO, and the AU test all use log-likelihoods not as final answers but as diagnostic tools — to test temporal signal, spatial heterogeneity, and tree selection bias.

2. **Event-based reconstruction**: DIVA and TREEMAP both classify evolutionary history as a sequence of discrete events (dispersal/vicariance, cospeciation/duplication) scored by parsimony-like costs. This "event vocabulary" approach independently emerged in biogeography and coevolution.

3. **Bootstrap resampling as null distribution**: The AU test and PLATO both generate null distributions through resampling — multiscale bootstrap for tree comparison, parametric bootstrap for recombination detection. The shared insight is that analytical null distributions are often inadequate for complex phylogenetic hypotheses.

### Why These Tools Were Lost

Each tool disappeared for different reasons:
- **TipDate**: Superseded by BEAST/BEAST2, which generalized its approach to Bayesian inference.
- **scaleboot/CONSEL**: The R package survived, but the standalone C implementation (CONSEL) became difficult to build; the algorithm's complexity discouraged reimplementation.
- **DIVA**: Platform-dependent (DOS), superseded by model-based methods (DEC, BioGeoBEARS).
- **TREEMAP**: Macintosh Classic application, lost to platform obsolescence.
- **PLATO**: Superseded by more sophisticated methods (GARD, 3SEQ), original website went dark.

In every case, the *algorithm* remained valuable even as the *software* disappeared. These reimplementations preserve the algorithmic ideas in a form that will outlast any particular platform or distribution mechanism.

---

## Updated Project Statistics

With the five resurrected tools, `phylip-rs` now contains:

| Metric | Before | After | Change |
|--------|--------|-------|--------|
| Lines of Rust | 35,805 | 42,105 | +6,300 |
| Source files | 58 | 65 | +7 |
| Unit tests | 934 | 1,062 | +128 |
| Doc tests | 25 | 30 | +5 |
| Validation tests | 91 | 91 | — |
| **Total tests** | **1,050** | **1,150** | **+100** |
| PHYLIP programs covered | 29/36 | 29/36 | — |
| Resurrected tools | 0 | 5 | +5 |
| New modules | — | 2 | +2 (`biogeography`, `reconciliation`) |
| External dependencies | 0 | 0 | — |

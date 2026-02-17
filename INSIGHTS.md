# Hidden Insights in PHYLIP: What Felsenstein Understood That We're Losing

*A technical analysis of the algorithms in Joe Felsenstein's PHYLIP software, reimplemented in Rust as part of the PHYLIP Archaeology project.*

## 1. The Pruning Algorithm Is a General-Purpose Computational Pattern

The most important algorithm in PHYLIP is Felsenstein's pruning algorithm (1981). On the surface, it computes the likelihood of DNA sequences on a phylogenetic tree. But its true nature is far more general.

The core recurrence is:

```
L_node(b) = PRODUCT over children c of
              SUM over states j of P_c(b, j) * L_c(j)
```

This is a **message-passing algorithm on a tree-structured graphical model** — a special case of the sum-product algorithm (belief propagation). Felsenstein published this in 1981, **seven years** before Judea Pearl formalized belief propagation (1988) and over a decade before the machine learning community embraced variational inference on graphical models.

The algorithm is not limited to 4 nucleotide states. It works on any discrete state space: 20 amino acids, 61 codons, or arbitrary categorical variables. It is not limited to molecular evolution. It works on **any tree-structured process** where:

- Observed data sits at the leaves
- Hidden states exist at internal nodes
- Transitions between states follow a known probabilistic model
- The tree topology and branch lengths are given

### Modern applications beyond phylogenetics

| Domain | Leaves (observed) | States | Tree structure |
|--------|-------------------|--------|----------------|
| Tumor phylogenetics | Single-cell genotypes | Mutation states | Clonal lineage tree |
| Cell lineage tracing | CRISPR barcode states | Barcode characters | Cell division tree |
| Language evolution | Modern word forms | Proto-language states | Language family tree |
| Gene expression evolution | Expression levels (binned) | Expression categories | Species tree |
| Viral quasispecies | Sampled sequences | Nucleotide states | Transmission tree |
| Cultural evolution | Present/absent traits | Cultural trait states | Population tree |

Most modern tools in these domains reinvent the pruning algorithm (often poorly). Felsenstein had the clean, numerically stable version — with underflow prevention via log-scaling — in 1981.

### Implementation

See [`phylip-rs/src/likelihood/pruning.rs`](phylip-rs/src/likelihood/pruning.rs). Key design features:

- **Underflow prevention**: When conditional likelihoods become extremely small (deep trees, many taxa), values are rescaled and a log-scaling factor is accumulated per site. This avoids floating-point underflow without expensive log-space arithmetic.
- **Model-agnostic**: The algorithm takes a `SubstitutionModel` trait object, decoupling the tree traversal from the specific probability model.
- **Postorder traversal**: The algorithm processes nodes from tips to root, ensuring children are computed before parents — a natural fit for any bottom-up tree computation.

---

## 2. The F84 Model: Exploiting Symmetry to Avoid Matrix Exponentiation

The F84 substitution model ([`phylip-rs/src/likelihood/models.rs`](phylip-rs/src/likelihood/models.rs)) reveals a deep insight about parameterization.

Most modern phylogenetics software computes transition probabilities via general matrix exponentiation: P(t) = exp(Qt), which requires eigendecomposition of the rate matrix Q. This is numerically fragile for ill-conditioned matrices and computationally expensive.

Felsenstein took a different approach. He observed that the purine/pyrimidine biochemical distinction creates a **block symmetry** in the rate matrix. By decomposing transitions into three components — identity, within-class (transition), and between-class (transversion) — each with its own exponential decay:

```
P(i->j; t) = delta(i,j) * exp1
            + same_class(i,j) * (pi_j / pi_class) * (exp2 - exp1)
            + pi_j * (1 - exp2)
```

he obtained **closed-form transition probabilities** that require only two `exp()` calls instead of a full eigendecomposition. This is mathematically equivalent to what physicists do with Lie algebra decompositions, but Felsenstein derived it from biological reasoning.

### Why this matters

For models with biological symmetry — and many biologically meaningful models *do* have symmetry — the Felsenstein approach of exploiting structure analytically is:

1. **Faster**: Two `exp()` calls vs. eigendecomposition + matrix multiply
2. **More numerically stable**: No ill-conditioned eigenvector matrices
3. **More interpretable**: Parameters map directly to biological quantities (transition/transversion ratio, base frequencies)

This principle applies to codon models with synonymous/nonsynonymous structure, amino acid models with chemical groupings, and epigenetic state models with methylation/demethylation asymmetry.

---

## 3. Site-Pattern Compression: A Data Structure Insight Worth Stealing

[`phylip-rs/src/likelihood/optimized.rs`](phylip-rs/src/likelihood/optimized.rs) implements `CompressedAlignment`: group identical alignment columns and weight them.

For a 1000-taxon, 10,000-site alignment, there might be only 500-2000 unique column patterns. Since the pruning algorithm is O(sites × nodes × states²), replacing "sites" with "unique patterns" gives a **5-10x speedup with zero approximation error**.

### The broader principle

This is the same insight as:
- **Suffix arrays** in genome assembly
- **Feature hashing** in machine learning
- **Run-length encoding** in data compression
- **Column-oriented storage** in databases (Parquet, Arrow)

The general rule: **before doing expensive per-column computation on categorical data, compress identical columns**.

### Where modern bioinformatics should use this but doesn't

| Data type | Column = | Typical redundancy |
|-----------|----------|--------------------|
| Variant call matrices (VCF) | Genotype pattern across samples | High (many monomorphic sites) |
| Microbiome OTU tables | Abundance pattern across samples | Moderate |
| Single-cell count matrices | Expression pattern across cells | Very high (many zeros) |
| Methylation arrays | Methylation state across CpG sites | High |
| GWAS genotype matrices | SNP genotype pattern | Moderate to high |

Tools like IQ-TREE and RAxML adopted pattern compression for phylogenetics. But most per-site statistical tests (GWAS, selection scans, differential methylation analysis) still process every column independently.

---

## 4. The Bootstrap: Weight Vectors, Not Physical Resampling

Felsenstein's bootstrap (1985) is widely used but often misunderstood. The key insight in the PHYLIP implementation is that **you never need to physically resample the data**.

A bootstrap replicate is fully described by a **weight vector** — an integer for each site saying how many times it was sampled. A site with weight 3 contributes 3× to the likelihood. Combined with site-pattern compression, this means:

1. Generate a weight vector (microseconds)
2. Re-compress with weights (fast)
3. Recompute likelihood (same cost as one evaluation)

No data copying, no new alignment objects.

### Block bootstrap and spatial correlation

The implementation also includes a **circular block bootstrap**: sample contiguous blocks of sites, wrapping around at sequence boundaries. Felsenstein recognized that adjacent sites in a gene may not be independent due to linkage and secondary structure.

This is the same principle as the **moving block bootstrap** in time series econometrics (Künsch 1989), but Felsenstein applied it to biological sequences in 1985.

**Where this matters today**: Most single-cell and spatial omics methods don't account for spatial or genomic correlation in their confidence estimates. The block bootstrap is the simplest correct approach, and it's been available since 1985.

---

## 5. Discrete Gamma Rates: A General Mixture Model Framework

Yang's (1994) discrete gamma model, implemented in [`phylip-rs/src/likelihood/gamma.rs`](phylip-rs/src/likelihood/gamma.rs), accounts for among-site rate variation. The underlying insight is Felsenstein's: **a single substitution rate across all sites is biologically absurd**.

But look at what the computation actually is:

```
L_site = (1/k) * sum_{c=1}^{k} L(site | rate_c)
```

This is a **mixture model** where the site likelihood is a weighted average over hidden rate categories. The log-sum-exp computation used for numerical stability:

```rust
let max_lnl = site_lnls.iter().copied()
    .fold(f64::NEG_INFINITY, f64::max);
let sum_exp: f64 = site_lnls.iter()
    .map(|&lnl| (lnl - max_lnl).exp()).sum();
```

is **identical to the computation in**:
- Gaussian mixture models (EM algorithm, E-step)
- Hidden Markov models (forward algorithm)
- Variational autoencoders (ELBO computation)

Felsenstein and Yang were doing mixture modeling with log-sum-exp stabilization decades before deep learning made it standard practice.

### The framework generalizes

The gamma rate heterogeneity framework accommodates **any source of hidden heterogeneity** in site-wise analysis:

| Source of heterogeneity | What varies | Alpha interpretation |
|-------------------------|-------------|---------------------|
| Functional constraint | Substitution rate | Low alpha = many constrained sites |
| Selection pressure | dN/dS ratio | Gamma-distributed omega |
| Methylation dynamics | CpG mutation rate | Rate variation across CpG sites |
| Chromatin accessibility | Mutation rate | Open vs. closed chromatin |
| RNA secondary structure | Evolutionary rate | Paired vs. unpaired sites |

### Mathematics from first principles

The implementation builds the entire gamma distribution machinery from scratch — Lanczos approximation for ln(Gamma), series and continued fraction representations for the incomplete gamma function, Halley's method for inversion, and Abramowitz & Stegun's rational approximation for the normal quantile function. This is a complete mathematical toolkit implemented in ~300 lines of Rust with no dependencies.

---

## 6. Fitch Parsimony: Bitwise Computation on Trees

The Fitch algorithm in [`phylip-rs/src/parsimony/wagner.rs`](phylip-rs/src/parsimony/wagner.rs) represents nucleotide state sets as bit masks:

```rust
pub struct StateSet(u8);  // Bit 0=A, 1=C, 2=G, 3=T, 4=gap

pub fn intersection(self, other: StateSet) -> StateSet {
    StateSet(self.0 & other.0)
}
pub fn union(self, other: StateSet) -> StateSet {
    StateSet(self.0 | other.0)
}
```

The entire parsimony score computation reduces to: if AND is non-zero, take AND (free); else take OR (cost +1). This computes minimum character state changes on the whole tree in O(sites × nodes) time — **linear in both dimensions**.

### Where parsimony beats likelihood

Modern genomics has largely abandoned parsimony for probabilistic methods. But for several applications, parsimony gives exact answers to the right question:

- **Somatic mutation counting in tumor phylogenies**: You want minimum changes, not maximum likelihood
- **CRISPR barcode lineage tracing**: Character-based data, not sequence substitution
- **Pangenomics**: Gene presence/absence on a tree — minimum gain/loss events
- **Antimicrobial resistance tracking**: Minimum acquisition/loss events for resistance genes
- **Ancestral state reconstruction**: When the tree is known and you want the most parsimonious explanation

The bitwise representation also scales naturally:
- `u64` → 64-state alphabets
- SIMD (`__m256i`) → 256 sites processed per instruction
- GPU → millions of sites in parallel

---

## 7. Model Selection: The Discipline of Counting Parameters

[`phylip-rs/src/likelihood/model_selection.rs`](phylip-rs/src/likelihood/model_selection.rs) implements AIC, BIC, and AICc — information criteria that penalize model complexity.

The critical function is parameter counting:

| Model | Substitution params | Branch length params | Total |
|-------|-------------------|---------------------|-------|
| JC69 | 0 | 2m - 3 | 2m - 3 |
| F84 | 4 (3 freq + 1 ts/tv) | 2m - 3 | 2m + 1 |
| + Gamma | +1 (alpha) | — | +1 |
| + Invariant sites | +1 (p_inv) | — | +1 |

where m = number of taxa.

The insight: **every free parameter costs you roughly 1-2 log-likelihood units of penalty**. This is Occam's razor made quantitative. A more complex model must improve the likelihood by at least this much per parameter to be justified.

Modern genomics routinely fits models with thousands of parameters (regularized regression, neural networks) but rarely asks: "How many of these parameters are actually justified by the data?" The PHYLIP tradition of explicit parameter counting and model comparison is a discipline the field could use more of.

---

## 8. What Felsenstein Understood That We're Losing

The deepest insight in this codebase isn't any single algorithm. It's the **relationship between the algorithms**.

Felsenstein understood that parsimony, distance methods, and maximum likelihood are all views of the same underlying problem — inferring evolutionary history from observed data. He proved when each one works and when it fails:

- **The Felsenstein zone** (1978): When two long branches are separated by a short internal branch, parsimony converges on the **wrong** tree with probability 1 as you add data. Likelihood gets it right.
- **Consistency**: ML is statistically consistent (converges to the truth with enough data). Parsimony is not, in general.
- **Efficiency**: NJ on good distances gives a reasonable topology in seconds. Use it as a starting point for ML, which refines it in minutes.

### Three things being lost to time

**1. The discipline of derivation from first principles.**
Our `gamma.rs` implements the Lanczos approximation, continued fractions, and Halley's method from scratch. Modern practitioners install scipy and call `gammainc()`. They can't debug it, modify it, or adapt it to new problems.

**2. The understanding that statistical phylogenetics is really about continuous-time Markov chains on trees.**
This abstraction applies to any tree-structured stochastic process — gene regulation cascades, language evolution, cultural transmission, epidemic spread. But it's taught as "a phylogenetics thing" rather than "a fundamental computational framework."

**3. The habit of testing against known analytical results.**
Our tests compare against hand-calculated values: JC69 transition probabilities from the formula, Yang's published gamma rates for alpha=0.5, the exponential CDF as a special case of the incomplete gamma function. This grounds code in theory. Most modern bioinformatics tests are regression tests — they verify the output hasn't changed, not that it's correct.

---

## Summary

| Algorithm | PHYLIP year | "Rediscovered" as | Modern application |
|-----------|-------------|-------------------|-------------------|
| Pruning algorithm | 1981 | Belief propagation (1988) | Any tree-structured inference |
| F84 closed-form P(t) | 1984 | Lie algebra decomposition | Structured substitution models |
| Site-pattern compression | ~1980s | Column-oriented databases | Any per-column computation |
| Bootstrap weights | 1985 | Weighted resampling | Spatial/genomic correlation |
| Block bootstrap | 1985 | Moving block bootstrap (1989) | Correlated site data |
| Discrete gamma rates | 1994 | Mixture models / log-sum-exp | Hidden heterogeneity in omics |
| Fitch parsimony | 1971 | Bitwise set operations | Character-based lineage tracing |
| Model selection (AIC/BIC) | ~1980s | Regularization theory | Parameter counting discipline |

Felsenstein built an intellectual framework, not just software. The algorithms are timeless. The software is archaeology. The *thinking* is what's worth preserving.

---

*This analysis is part of the [PHYLIP Archaeology](https://github.com/shandley/phylip-archaeology) project, which preserves and modernizes Joe Felsenstein's PHYLIP phylogenetics software as a pure-Rust implementation with zero external dependencies.*

## References

- Felsenstein, J. (1978). Cases in which parsimony or compatibility methods will be positively misleading. *Systematic Zoology*, 27, 401-410.
- Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. *Journal of Molecular Evolution*, 17, 368-376.
- Felsenstein, J. (1984). A likelihood approach to character weighting and what it tells us about parsimony and compatibility. *Biological Journal of the Linnean Society*, 16, 183-196.
- Felsenstein, J. (1985). Confidence limits on phylogenies: an approach using the bootstrap. *Evolution*, 39, 783-791.
- Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
- Fitch, W.M. (1971). Toward defining the course of evolution. *Systematic Zoology*, 20, 406-416.
- Pearl, J. (1988). *Probabilistic Reasoning in Intelligent Systems*. Morgan Kaufmann.
- Saitou, N. & Nei, M. (1987). The neighbor-joining method. *Molecular Biology and Evolution*, 4, 406-425.
- Yang, Z. (1994). Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites. *Journal of Molecular Evolution*, 39, 306-314.

# Reflection: What We Built and What We Learned

*Notes on the PHYLIP Archaeology project — what emerged from reimplementing 45 years of algorithmic thinking in a modern language.*

---

## The Project

We set out to do archaeology — to excavate, preserve, and understand one of the most important codebases in the history of computational biology. Joe Felsenstein's PHYLIP, first released in 1980, contains algorithms that remain the mathematical backbone of modern phylogenetics. But like all software from that era, it is at risk: link rot, bitwise decay, and the steady loss of institutional knowledge as the researchers who understood these methods retire or move on.

What emerged was something more than preservation. By reimplementing every algorithm from first principles — every mathematical function, every matrix decomposition, every statistical test — we were forced to understand them deeply. And in that understanding, patterns emerged that connect PHYLIP's 1980s algorithms to problems being solved (or poorly solved) across modern computational biology, machine learning, and data science.

## By the Numbers

| Metric | Value |
|--------|-------|
| Lines of Rust | **20,749** |
| Source files | **35** |
| Unit tests | **561** |
| Doc tests | **17** |
| Total tests | **578** |
| Compiler warnings | **0** |
| External dependencies | **0** |
| CLI commands | **5** |
| Interactive demonstrations | **2** |
| Substitution models | **4** (JC69, F84, Poisson, WAG) |

The zero-dependency constraint was deliberate and consequential. Every mathematical function — the Lanczos approximation for the gamma function, continued fractions for the incomplete gamma integral, Halley's method for numerical inversion, Abramowitz & Stegun's rational approximation for the normal quantile — is implemented from scratch. The code is its own textbook. You can read `gamma.rs` and learn how the incomplete gamma function works, because there is nowhere else to hide.

## What We Implemented

### Core Algorithms

**Maximum Likelihood** — Felsenstein's pruning algorithm (1981), the foundational method. A postorder tree traversal that computes conditional likelihoods at each node, with underflow prevention via per-site log-scaling factors. Decoupled from specific substitution models via a trait interface. ML tree search via stepwise addition with SPR rearrangement. Branch length optimization by Newton-Raphson with numerical derivatives. NNI (Nearest-Neighbor Interchange) for local tree refinement.

**Substitution Models** — JC69 (equal rates, analytical P(t)), F84 (unequal frequencies with purine/pyrimidine symmetry, closed-form P(t) requiring only two exponentials instead of eigendecomposition), Poisson model for protein sequences, and the WAG empirical amino acid rate matrix. A 20-state pruning algorithm for protein data.

**Parsimony** — Fitch algorithm (1971) with bitwise state set operations: intersection for free steps, union for costly ones. Wagner parsimony tree search with SPR. Ancestral state reconstruction via the Fitch preorder pass.

**Distance Methods** — Neighbor-Joining (Saitou & Nei 1987) with the Q-criterion for partner selection. Fitch-Margoliash weighted least squares. ML pairwise distances via Newton-Raphson optimization of single-pair likelihoods.

**Statistical Support** — Bootstrap resampling using weight vectors (no data copying). Block bootstrap for spatially correlated sites. Delete-fraction jackknife. Full bootstrap + ML integration: replicate ML searches with consensus support values. Consensus trees: strict, majority-rule, extended majority-rule, and threshold methods, all built on bipartition (split) representations.

**Rate Heterogeneity** — Yang's (1994) discrete gamma model with k rate categories. Complete gamma distribution machinery from first principles. Alpha parameter optimization via golden section search. Integration with the pruning algorithm as a mixture model with log-sum-exp stabilization.

**Model Selection** — AIC, BIC, and AICc with proper parameter counting. Akaike weights for model averaging. The discipline of asking: "How many parameters are justified by the data?"

**I/O** — PHYLIP interleaved and sequential format parser. FASTA parser (DNA and protein). Newick tree format reader and writer. PHYLIP-style formatted output reports. Command-line interface with five analysis modes.

### Interactive Demonstrations

**The Felsenstein Zone** — A simulation that makes Felsenstein's 1978 result visceral. Four taxa on a tree with two long branches (branch length 0.80) separated by a short internal branch (0.01). DNA sequences simulated along this tree, then all three possible 4-taxon topologies evaluated under parsimony and ML. The result is dramatic:

| Sites | Parsimony correct | ML correct |
|-------|-------------------|------------|
| 100 | 8% | 88% |
| 500 | 2% | 100% |
| 1,000 | 0% | 100% |
| 5,000 | 0% | 100% |
| 10,000 | 0% | 100% |

More data makes parsimony *more wrong*. This is statistical inconsistency made concrete.

**Language Evolution** — The exact same pruning algorithm that was designed for DNA applied to human language data: cognate class assignments for 37 vocabulary items across English, German, French, Italian, Spanish, and Portuguese. Not a single line of algorithm code changes. The system correctly identifies the known language family tree (Germanic vs. Romance, with correct subgroupings) and reconstructs "proto-language" states at the root. This demonstrates that Felsenstein's algorithm is a general-purpose inference engine for discrete states on trees — not just a "DNA algorithm."

## What We Discovered

The deepest insights came not from any single algorithm but from seeing them together — as facets of a coherent intellectual framework that Felsenstein built over decades.

### 1. The Pruning Algorithm Anticipated Belief Propagation

Felsenstein published the pruning algorithm in 1981. Judea Pearl formalized belief propagation — the general algorithm for inference on graphical models — in 1988. The pruning algorithm is a special case: message passing on a tree-structured graphical model with discrete states. Felsenstein had the clean, numerically stable, practically useful version seven years before the general theory existed.

### 2. Biological Symmetry Eliminates Computational Overhead

The F84 model exploits purine/pyrimidine biochemical symmetry to decompose transition probabilities into three analytical components, each with its own exponential decay. Two `exp()` calls replace a full matrix eigendecomposition. This is faster, more stable, and more interpretable. The principle — exploit domain-specific structure rather than reaching for general numerical machinery — is underappreciated in modern bioinformatics.

### 3. Site-Pattern Compression Is a Universal Data Structure Insight

Group identical alignment columns and weight them. For a 1000-taxon, 10,000-site alignment with 1,500 unique patterns, this gives a 6.7x speedup with zero approximation error. The principle applies everywhere categorical columns are processed independently: VCF genotype matrices, OTU tables, single-cell count matrices, methylation arrays. Most tools outside phylogenetics still don't do this.

### 4. Bootstrap Weight Vectors Eliminate Data Copying

A bootstrap replicate is fully described by a weight vector — integers saying how many times each site was sampled. Combined with pattern compression, this means bootstrap replicates cost almost nothing: generate weights (microseconds), recompress (fast), evaluate (same cost as one analysis). No data duplication.

### 5. Discrete Gamma Rates Are a Mixture Model Framework

The log-sum-exp computation for integrating over rate categories is identical to what modern deep learning uses for mixture models, variational inference, and the forward algorithm in HMMs. Felsenstein and Yang were doing this in 1994, decades before it became standard practice in machine learning.

### 6. Parsimony Has Irreplaceable Applications

Despite being statistically inconsistent in general, parsimony answers exactly the right question for several modern applications: minimum somatic mutation counts in tumor phylogenies, CRISPR barcode character changes in lineage tracing, gene gain/loss events in pangenomics. The bitwise Fitch implementation scales naturally to SIMD and GPU parallelism.

### 7. Parameter Counting Is a Lost Discipline

PHYLIP's model selection computes AIC and BIC with explicit parameter counts. Every free parameter costs roughly 1-2 log-likelihood units of penalty. Modern genomics routinely fits models with thousands of parameters but rarely asks whether they are justified. The PHYLIP tradition of quantitative Occam's razor deserves revival.

## What This Means

Three things are being lost to time.

**The discipline of derivation from first principles.** Our `gamma.rs` implements the Lanczos approximation, continued fractions, and Halley's method from scratch. Modern practitioners call `scipy.stats.gamma.ppf()`. They cannot debug it, modify it, or adapt it to new problems. When the library changes behavior, they have no recourse.

**The understanding that statistical phylogenetics is really about continuous-time Markov chains on trees.** This abstraction applies to any tree-structured stochastic process — gene regulation cascades, language evolution, cultural transmission, epidemic spread. But it is taught as "a phylogenetics thing" rather than "a fundamental computational framework."

**The habit of testing against known analytical results.** Our tests compare against hand-calculated values: JC69 transition probabilities from the formula, Yang's published gamma rates for alpha = 0.5, the exponential CDF as a special case of the incomplete gamma function. This grounds code in theory. Most modern bioinformatics tests are regression tests — they verify the output hasn't changed, not that it's correct.

## The Repository

The complete project is at **https://github.com/shandley/phylip-archaeology**.

What started as archaeology became a tribute — not to old software, but to the kind of deep, principled thinking that built a field. The algorithms are timeless. The software is archaeology. The *thinking* is what is worth preserving.

---

*This reflection is part of the [PHYLIP Archaeology](https://github.com/shandley/phylip-archaeology) project.*

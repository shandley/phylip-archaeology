# LLM-Assisted Code Archaeology Recovers the Algorithmic Legacy of PHYLIP

**Scott A. Handley**

Department of Pathology and Immunology, Washington University School of Medicine, St. Louis, MO, USA

---

## Abstract

Scientific software encodes algorithmic knowledge often absent from published papers, existing only in aging source code. We demonstrate that large language models can serve as archaeological instruments, reading and resurrecting algorithms from legacy codebases. Applying this to PHYLIP — the most influential phylogenetics package, first released in 1980 — we identified 20 cross-disciplinary connections linking phylogenetics to electrical engineering, information geometry, algebraic statistics, and coding theory, including the recognition that Felsenstein's 1981 pruning algorithm is mathematically identical to belief propagation, formalized seven years later. The recovered algorithms were validated through a 36,745-line Rust reimplementation with 959 tests and zero dependencies. Benchmarking against modern tools shows comparable results on small datasets and quantifies the search heuristic gap accumulated over four decades. Analysis of Felsenstein's software catalog revealed that 23 of 407 phylogenetics tools have been permanently lost. LLM-assisted code archaeology offers a scalable approach to preserving computational knowledge across sciences.

---

## Introduction

Scientific software has a half-life. As programming languages evolve, operating systems change, and maintainers retire, the algorithms embedded in research code become increasingly inaccessible [1]. Unlike mathematical proofs, which persist in print, computational methods are often documented only in their implementations — implementations that compile on systems no longer in use, written in idioms no longer taught, by researchers no longer available to explain their design decisions. The result is a growing corpus of computational knowledge that, while technically still extant in version-controlled repositories and archived websites, is effectively lost — comprehensible only to specialists who may themselves be retiring.

This problem is particularly acute in fields where foundational software was written decades ago by individual researchers who embedded mathematical innovations in their code without always publishing them as separate contributions. Phylogenetics — the inference of evolutionary relationships from molecular data — is an exemplary case. The field's most influential software package, PHYLIP (PHYLogeny Inference Package), was first released in October 1980 by Joe Felsenstein at the University of Washington [2]. Distributed initially on magnetic tapes in an era before the IBM PC, PHYLIP grew to encompass 36 programs implementing every major class of phylogenetic algorithm: maximum likelihood, parsimony, distance methods, bootstrapping, and consensus trees. Felsenstein's papers have accumulated over 162,000 citations, with his bootstrap paper alone exceeding 32,000 citations — placing it among the 100 most-cited papers in all of science [3,4].

Despite this influence, PHYLIP's algorithmic content has never been systematically extracted and documented. Felsenstein's textbook [16] describes the statistical theory, but the source code contains implementation-level details — connections between phylogenetic algorithms and circuit theory, information geometry, algebraic statistics, and error-correcting codes — that have rarely been systematically documented. The code also contains numerical techniques (underflow prevention via log-scaling, site-pattern compression, closed-form transition probabilities from symmetry exploitation) that represent practical wisdom about scientific computation. As PHYLIP's C codebase ages and the era of its development recedes, these insights risk being lost.

We propose **LLM-assisted code archaeology** as a methodology for recovering algorithmic knowledge from legacy scientific software. Using a large language model (Claude, Anthropic) in an iterative human-AI collaboration, we read, analyzed, and reimplemented PHYLIP's algorithms in a modern language (Rust), producing a validated codebase of 36,745 lines with 959 tests and zero external dependencies. The process revealed 20 cross-disciplinary connections between PHYLIP's algorithms and techniques in other fields — some previously documented but underappreciated, others apparently novel (Table 1, Supplementary Note 1). We benchmarked the reimplementation against three contemporary maximum likelihood tools to quantify what four decades of heuristic engineering have added beyond Felsenstein's original algorithm designs. In parallel, we analyzed Felsenstein's hand-curated software catalog of 407 phylogenetics tools, documenting the preservation status of an entire discipline's computational infrastructure and finding that 23 tools have been permanently lost.

---

## Results

### LLM-assisted code archaeology

We developed an iterative workflow for extracting algorithmic knowledge from legacy source code using a large language model (**Fig. 1**). The process consisted of five phases: (1) the LLM reads and annotates original source code, identifying algorithms and data structures; (2) algorithms are extracted and expressed in mathematical notation, independent of their implementation language; (3) extracted algorithms are cross-referenced against the published literature to identify whether the insight has been previously documented; (4) each algorithm is reimplemented from the mathematical description in a modern language (Rust), with the zero-dependency constraint forcing implementation of all mathematical primitives from first principles; and (5) the reimplementation is validated against known analytical results and hand-calculated values.

Human domain expertise was essential at every stage: directing the LLM to algorithmically interesting code sections, recognizing cross-disciplinary connections, designing validation tests against known analytical results rather than regression tests, and correcting the LLM when it misunderstood biological context. The LLM's primary contribution was its ability to read and comprehend large volumes of C code rapidly, synthesize algorithmic descriptions across multiple source files, and produce working implementations that could be iteratively refined. When the LLM encountered PHYLIP's independent contrasts implementation, for example, it correctly extracted the mathematical recurrence and produced a working Rust implementation, but the insight that the variance formula is the parallel resistor equation required human recognition of the circuit theory connection.

The complete reimplementation — phylip-rs — covers 29 of PHYLIP's 36 programs (the seven excluded programs — drawgram, drawtree, dnamove, dolmove, move, retree, and factor — are interactive graphical or editing tools with no novel algorithmic content). The codebase comprises 60 source files organized into 11 modules: likelihood inference (7,831 lines), substitution models (6,905 lines), parsimony methods (5,884 lines), distance methods (3,010 lines), input/output (2,582 lines), tree operations (1,969 lines), comparative methods (1,927 lines), compatibility analysis (1,758 lines), bootstrapping (1,355 lines), phylogenetic invariants (830 lines), and consensus trees (814 lines). All 959 tests pass with zero compiler warnings. The zero-dependency constraint — no external crates, no linear algebra libraries — was a deliberate design choice: it forced the reimplementation of all mathematical primitives (gamma function, matrix eigendecomposition, Newton-Raphson optimization, pseudorandom number generation) from first principles, ensuring that every algorithm is fully transparent and auditable.

### Cross-disciplinary algorithmic connections

The code archaeology process identified 20 connections between PHYLIP's algorithms and techniques in fields outside phylogenetics (Table 1; Supplementary Note 1). These range from apparently novel equivalences to well-documented results that nonetheless remain underappreciated in computational biology. We highlight four case studies that illustrate the methodology's capacity to surface and validate cross-disciplinary connections.

**Felsenstein's pruning algorithm is belief propagation.** Felsenstein's pruning algorithm [5] computes the likelihood of sequence data on a phylogenetic tree via a postorder traversal, computing conditional likelihoods at each internal node as the product over children of summed transition probabilities weighted by child likelihoods. The core recurrence — L_node(b) = PRODUCT_c SUM_j P_c(b,j) * L_c(j) — is a message-passing algorithm on a tree-structured graphical model, mathematically identical to the sum-product algorithm (belief propagation) that Judea Pearl formalized in 1988 [6]. Felsenstein published this algorithm in 1981, seven years before Pearl's book (**Fig. 2a**). The equivalence has been noted in the graphical models literature [30], but remains underappreciated in computational biology, where the pruning algorithm is typically presented as a phylogenetics-specific technique rather than a general inference method. The algorithm works on any discrete state space and any tree-structured process where observed data sits at leaves, hidden states exist at internal nodes, and transitions follow a known probabilistic model. Modern applications include tumor phylogenetics, cell lineage tracing via CRISPR barcodes, and language evolution — domains where practitioners frequently reinvent the pruning algorithm, often without the numerical stability features (log-scaling for underflow prevention) that Felsenstein included in 1981.

**Independent contrasts are Kirchhoff's circuit laws.** Felsenstein's independent contrasts algorithm [7] computes phylogenetic regression in O(n) time instead of the naive O(n^3) matrix inversion, via a postorder traversal that propagates weighted averages upward through the tree. The core variance formula — v = (v_L * v_R) / (v_L + v_R) — is the parallel resistor formula from electrical engineering. This is not a loose analogy: it is a formal mathematical equivalence (**Fig. 2b**). Branch lengths in the phylogenetic tree correspond exactly to resistances in a circuit, tip trait values to boundary voltages, contrasts at internal nodes to currents, and the weighted average propagation to Kirchhoff's voltage law. Numerical verification confirms the equivalence to eight decimal places. Felsenstein derived this from statistical considerations in 1985; the connection to circuit theory was implicit in the mathematics but never stated.

**The genetic code step matrix connects parsimony to coding theory.** PHYLIP's protein parsimony implementation constructs a 20x20 cost matrix from the genetic code, computing the minimum number of nucleotide substitutions required to convert each amino acid to every other (**Fig. 2c**). This Sankoff step matrix [21] implicitly encodes the genetic code's error-tolerance properties — a connection to coding theory that was present in PHYLIP's source code but not discussed in its documentation. The biological conclusion that the genetic code is optimized for error tolerance was independently established by Freeland and Hurst [8], who used a similar randomization approach (z-score = -2.76 against random codes). Our reimplementation of the step matrix construction and randomization test reproduces their finding, demonstrating that PHYLIP's weighted parsimony machinery contained the computational infrastructure for this analysis decades before it was published.

**LogDet distance: compositional robustness from determinant factorization.** The LogDet distance [9,10] isolates pure evolutionary signal by exploiting a determinant factorization: the divergence matrix F between two sequences factors as F = diag(pi) x P(t), and taking the log-determinant cancels the base frequency terms entirely (**Fig. 2d**). This property was described in the original publications, but the reimplementation process clarifies *why* it works in a way that the original mathematical treatments leave implicit: the log-determinant decomposition separates the composition-dependent terms from the rate-dependent terms algebraically, not just statistically. When base composition varies across lineages — common in bacterial genomes (GC content ranges from 25% to 75%), mitochondrial DNA, and ancient DNA — standard distances (JC69, K2P) produce biased estimates. LogDet gives correct distances regardless of compositional heterogeneity. Our reimplementation computes the 4x4 determinant via Laplace cofactor expansion from first principles, requiring no linear algebra library — demonstrating that even numerically sophisticated algorithms can be made fully self-contained and transparent.

**Table 1. Twenty cross-disciplinary connections identified in PHYLIP source code.**

| # | Algorithm | Year | Cross-disciplinary connection | Published? |
|---|-----------|------|-------------------------------|------------|
| 1 | Pruning algorithm | 1981 | Belief propagation (Pearl, 1988) | Algorithm yes; connection noted [30] |
| 2 | F84 closed-form P(t) | 1984 | Lie algebra decomposition | Implementation only |
| 3 | Site-pattern compression | ~1980s | Column-oriented databases | Implementation only |
| 4 | Bootstrap weight vectors | 1985 | Weighted resampling | Algorithm yes; generalization no |
| 5 | Discrete gamma rates | 1994 | Mixture models / log-sum-exp | Algorithm yes; connection no |
| 6 | Fitch parsimony | 1971 | Bitwise set operations | Algorithm yes; implementation no |
| 7 | Model selection (AIC/BIC) | ~1980s | Regularization theory | Well-documented in statistics |
| 8 | First-principles derivation | — | Computational self-sufficiency | Methodology observation |
| 9 | Independent contrasts | 1985 | Kirchhoff's circuit laws | Algorithm yes; connection no |
| 10 | Contml stereographic projection | 1973 | Hellinger embedding / information geometry | Algorithm yes; connection no |
| 11 | Hendy-Penny supplement bound | 1982 | Dual decomposition / Lagrangian relaxation | Bound yes; connection no |
| 12 | Dollo parsimony | 1977 | Min-cut on trees | Algorithm yes; connection no |
| 13 | LogDet distance | 1994 | Determinantal factorization | Algorithm yes; robustness implicit |
| 14 | Kitsch scrunch | ~1980s | Pool-adjacent-violators / isotonic regression | Implementation only |
| 15 | Clique analysis | 1986 | Bron-Kerbosch / splits equivalence | Algorithm yes; equivalence no |
| 16 | Lake's invariants | 1987 | Algebraic statistics / variety ideals | Algorithm yes; connection no |
| 17 | O(n) Brownian ML | 1973 | Tree-structured Gaussian processes | Algorithm yes; connection no |
| 18 | Felsenstein-Churchill HMM | 1996 | Baum-Welch for rate correlation | Published; connection well-known |
| 19 | Protein Sankoff | ~1990s | Genetic code optimization / coding theory | Implementation only |
| 20 | Score-ordered B&B | 1982 | Greedy-guided exact search (A*, alpha-beta) | Algorithm yes; connection no |

The "Published?" column indicates whether the algorithm was published in a paper (most were) versus whether the specific cross-disciplinary connection we identified was previously documented. Entries range from apparently novel cross-disciplinary connections (e.g., independent contrasts / circuit theory) through underappreciated but documented equivalences (e.g., pruning / belief propagation [30]) to well-established techniques included because their presence in PHYLIP illustrates the breadth of algorithmic content in a single codebase. Original algorithm descriptions: Fitch parsimony [20], Sankoff [21], Hendy-Penny branch-and-bound [22], Brownian ML on continuous characters [23], Felsenstein-Churchill HMM [24], Bron-Kerbosch clique-finding [25], discrete gamma rates [17]. Full mathematical descriptions and numerical demonstrations for all 20 entries are provided in Supplementary Note 1.

### Software catalog preservation

Felsenstein maintained a hand-curated catalog of phylogenetics software — effectively a package manager for an entire scientific discipline, predating SourceForge (1999), Bioconductor (2001), and GitHub (2008). We analyzed all 407 tools documented in this catalog, extracting metadata from the original HTML pages and cross-referencing each tool's URL against the Wayback Machine (**Fig. 3**).

Of the 407 tools, 196 (48.2%) are archived but no longer actively maintained, 137 (33.7%) appear dormant with websites still accessible but no evidence of recent development activity, 34 (8.4%) are completely unreachable with their original URLs returning errors, and 40 (9.8%) have unknown status due to missing or ambiguous URL information. The Wayback Machine has preserved snapshots of 203 tools, but 23 tools have no archived copy anywhere and may be permanently lost.

Tool publication peaked during 2000-2004 with 89 new tools, coinciding with the genomics revolution and the rise of Bayesian phylogenetics (**Fig. 3a**). The catalog reveals a clear methodological paradigm shift over time (**Fig. 3c**): parsimony and distance methods dominated the 1980s-1990s, maximum likelihood methods rose sharply in the late 1990s following improvements in computational power and heuristic search strategies, and Bayesian methods emerged after 2000 with MrBayes and BEAST. By the 2005-2009 era, maximum likelihood had become the dominant paradigm, reflecting a field-wide shift from algorithmic simplicity to statistical rigor. Programming language preferences evolved from C dominance (90 tools) through Java (60 tools) to Python and R in recent years, tracking broader trends in scientific computing. The catalog functions as a longitudinal record of an entire discipline's computational infrastructure — the only such record that exists for phylogenetics.

### Benchmarking against modern tools

To quantify the gap between PHYLIP's original algorithm designs and modern implementations, we benchmarked phylip-rs against three contemporary tools: IQ-TREE 3 [11], RAxML-NG [12], and VeryFastTree [13]. We simulated 36 datasets under the JC69 model (10-500 taxa, 500-5,000 sites, 3 replicates per condition) and ran all tools single-threaded for fair comparison (**Fig. 4**). To ensure scoring consistency, all inferred trees were evaluated under JC69 using phylip-rs; we validated this scorer against IQ-TREE's internal likelihood calculation and confirmed agreement to four decimal places on shared topologies (Supplementary Note 4).

On small datasets (10-20 taxa), phylip-rs often found trees with log-likelihoods comparable to those of modern tools. In some replicates all five tools converged to the same optimum (e.g., lnL = -35,848.85 for all tools on one 20-taxon, 5,000-site dataset; normalized Robinson-Foulds distance [29] = 0.059 for all tools). However, performance was variable: across 15 datasets at 10-20 taxa, phylip-rs ML matched IQ-TREE's optimum (within 0.01 lnL units) in 4 cases, was within 1 lnL unit in 4 cases, and found substantially worse trees (gaps of 1-81 lnL units) in 7 cases. The largest gaps occurred on 20-taxon datasets with 5,000 sites, where the combinatorial search space exceeds phylip-rs's simple NNI strategy. This demonstrates that Felsenstein's likelihood computation is mathematically sound — the same scoring function reaches the same answer — but the tree search strategy is where modern heuristics diverge.

At 50 taxa, the search gap widens: phylip-rs ML required 234 seconds versus 14 seconds for IQ-TREE (17x slower), finding trees 1-16 log-likelihood units worse. Beyond 100 taxa, phylip-rs ML could not complete within the 600-second timeout, while IQ-TREE and VeryFastTree continued to scale. At the largest dataset size (500 taxa), only VeryFastTree (14 seconds) and phylip-rs NJ (180 seconds) completed; IQ-TREE and RAxML-NG also exceeded the timeout at this scale. This pattern quantifies four decades of heuristic engineering: the likelihood functions are mathematically equivalent, but modern tools employ SIMD vectorization, sophisticated tree rearrangement strategies (lazy SPR, stochastic NNI), and site-pattern-aware likelihood computation that PHYLIP's original designs lack.

Phylip-rs's neighbor-joining implementation [18] completed all 36 datasets with competitive accuracy, demonstrating that deterministic algorithms transfer directly across implementations. At matched dataset sizes (10-20 taxa), phylip-rs used 2.7-4.4 MB versus 54-55 MB for IQ-TREE and 19-20 MB for RAxML-NG, though IQ-TREE's memory footprint includes substantial runtime overhead from its execution environment. VeryFastTree, optimized for speed with approximate methods, completed all datasets in a median of 1.3 seconds and never timed out.

### Validated reimplementation

The phylip-rs reimplementation serves as both a validation artifact and an educational resource. The zero-dependency constraint — no external crates, no linear algebra libraries, no random number generator crates — forced the implementation of all mathematical primitives from first principles: the gamma function via Lanczos approximation, matrix exponentiation via eigendecomposition, Newton-Raphson optimization, continued fraction evaluation, and pseudorandom number generation. This constraint ensures that every algorithm is fully transparent and self-contained.

The 959 tests validate against known analytical results rather than regression tests. For example: the JC69 distance for two sequences differing at exactly 25% of sites should be 0.3041 (the Jukes-Cantor correction); the likelihood of a star tree under JC69 has a closed-form solution; neighbor-joining on a distance matrix computed from a known tree should recover that tree exactly; and the sum of squared standardized independent contrasts should equal the generalized least squares solution computed via direct matrix inversion. This approach catches algorithmic errors that regression tests miss, because it tests the mathematics rather than the implementation's consistency with itself.

Ten interactive demonstrations, distributed as compilable examples, illustrate cross-disciplinary applications: the Felsenstein zone [19] (statistical consistency of ML versus parsimony), Kirchhoff contrasts (circuit theory equivalence), genetic code optimization (coding theory), compositional bias correction (LogDet), clock-constrained inference (isotonic regression via the pool-adjacent-violators algorithm), stereographic projection [27,28] (information geometry on allele frequency simplices), Dollo parsimony (combinatorial optimization), branch-and-bound pruning [22] (dual decomposition bounds), Lake's invariants [26] (algebraic statistics on site-pattern frequencies), and language evolution (the pruning algorithm applied to Indo-European cognate data). Each demonstration is self-contained, serving both as validation and as an educational resource for researchers learning these algorithms.

---

## Discussion

LLM-assisted code archaeology offers a scalable methodology for recovering algorithmic knowledge from legacy scientific software. Our application to PHYLIP demonstrates that a single codebase can contain connections spanning multiple disciplines — connections that, while sometimes documented individually, have not been systematically cataloged or validated. The methodology is not specific to phylogenetics: any field with foundational software written decades ago by individual researchers — BLAST and HMMER in sequence analysis [14], PAML in molecular evolution [15], CESM in climate modeling, ROOT and Geant4 in particle physics — could benefit from systematic code archaeology. The common thread is software written by mathematically sophisticated researchers who embedded algorithmic innovations in their implementations without always publishing them separately.

The human-AI collaboration was essential and its dynamics are instructive for future applications. The LLM excelled at reading and comprehending large volumes of C code, identifying algorithmic patterns across multiple source files, and producing working reimplementations from mathematical descriptions. It was particularly effective at the mechanical aspects of code archaeology: tracing data flow across functions, identifying variable correspondences between mathematical notation and code, and generating well-structured implementations. But the cross-disciplinary connections — recognizing the circuit theory equivalence in independent contrasts, the coding theory implication of the protein step matrix, or the information geometry underlying the stereographic projection — required human domain expertise. The LLM served as a powerful archaeological instrument; the human researcher provided the intellectual framework for interpreting what was found.

The LLM also made characteristic errors that are instructive for future applications of this methodology. It occasionally conflated mathematical similarity with mathematical identity, proposing connections that were analogies rather than equivalences — for example, suggesting that two algorithms were "identical" when they shared structural features but operated on different mathematical objects. It sometimes struggled with the biological motivation for algorithmic choices — for instance, why Felsenstein chose the F84 model's purine/pyrimidine parameterization over the more general GTR, a decision rooted in biochemical constraints that the LLM could not independently assess. It produced implementations that were syntactically correct but numerically unstable, requiring human intervention to add log-scaling for underflow prevention or to handle edge cases in matrix decomposition. And it required careful prompting to distinguish insights genuinely absent from the literature from results that were simply not immediately apparent in the code. These failure modes suggest that LLM-assisted code archaeology is best viewed as an augmentation of expert analysis, not a replacement — the LLM amplifies the reach of a domain expert but cannot substitute for domain knowledge.

Our analysis of Felsenstein's software catalog reveals a broader preservation crisis in scientific computing. The permanent loss of 23 tools with no archived copy represents algorithmic knowledge that may be irrecoverable. The catalog itself — a hand-curated resource maintained by a single researcher over two decades — is a form of scientific infrastructure that has no natural successor. Unlike code repositories, which are increasingly archived by platforms like Software Heritage and GitHub Arctic Code Vault, hand-curated catalogs that contextualize software within a field's intellectual history have no systematic preservation mechanism. Notably, the dead tools are not uniformly distributed: tools from the 1990s-2000s peak era are most at risk, as their university-hosted web pages are decommissioned faster than personal research pages from the pre-web era (which were often preserved in print documentation). The pattern is concerning because this peak era coincides with the rise of maximum likelihood and Bayesian methods that now dominate the field — meaning the most methodologically relevant historical software is the most endangered. As the generation of researchers who wrote foundational bioinformatics software approaches retirement, the urgency of preservation increases.

Several limitations should be noted. We examined a single codebase (PHYLIP), and the methodology's effectiveness on other legacy software remains to be demonstrated. The seven excluded interactive programs (drawing and editing tools) represent a limitation of our Rust reimplementation, though these programs contain no novel algorithms. The benchmarking comparison uses simulated data under the same model (JC69) employed by phylip-rs, which represents an idealized scenario; real-world datasets with model misspecification would likely show larger performance gaps, and phylip-rs's limited model repertoire (JC69, F84) means it cannot be directly compared against modern tools on empirical datasets where GTR+G is standard. We did not compare against the original PHYLIP C programs, which would provide a direct validation of reimplementation fidelity; our validation instead relies on agreement with known analytical results. Some of the cross-disciplinary connections we highlight (notably the pruning/belief propagation equivalence) have been noted in the graphical models literature [30], and our contribution is to systematize and validate these connections rather than to claim priority for all of them. Finally, our catalog analysis is limited to Felsenstein's curated entries and does not cover phylogenetics tools developed after the catalog's last major update.

The algorithms Felsenstein embedded in PHYLIP are not historical curiosities. Belief propagation on trees powers modern probabilistic programming frameworks [6]. The parallel resistor formula for variance propagation appears in spatial statistics, sensor networks, and Gaussian process inference on tree-structured domains. The Sankoff algorithm with weighted costs is used in ancestral genome reconstruction, morphological character evolution, and natural language processing. LogDet distances are increasingly relevant as metagenomic studies encounter organisms with extreme compositional biases. Algebraic invariants are a growth area in computational algebra and model identifiability [26]. By excavating these connections from aging C code and expressing them in a modern, validated, zero-dependency implementation, we aim to ensure that the next generation of computational scientists can build on foundations that might otherwise be lost to bit rot and institutional memory.

---

## Methods

### LLM interaction

All code archaeology was performed using Claude (Anthropic; claude-sonnet-4-20250514 and claude-opus-4-20250514 models) via the Claude Code command-line interface. The workflow was iterative: the LLM read original PHYLIP C source files, produced algorithm descriptions, generated Rust implementations, and refined them based on test failures and human feedback. The human researcher directed which code sections to examine, designed validation tests based on domain knowledge, identified cross-disciplinary connections, and corrected biological misunderstandings. In accordance with Nature Methods policy, the LLM is not listed as an author; its substantial contribution to code analysis and implementation is disclosed here.

### Reimplementation architecture

phylip-rs is implemented in Rust (edition 2021) with zero external dependencies (no crates). The architecture uses trait-based polymorphism: `SubstitutionModel` for pluggable substitution models (JC69, K2P, F81, F84, Poisson, WAG) and `ParsimonyScorer` for pluggable parsimony methods (Fitch, Wagner, Dollo, Camin-Sokal, Sankoff). All mathematical functions are implemented from first principles, including the gamma function (Lanczos approximation), matrix eigendecomposition, Newton-Raphson optimization, and the Mersenne Twister pseudorandom number generator.

### Validation methodology

Tests validate against known analytical results rather than regression tests. Examples include: verifying JC69 distances against the closed-form correction formula; checking that neighbor-joining on distances computed from a known tree recovers that tree; confirming that the pruning algorithm likelihood matches hand-calculated values for small trees; and verifying that independent contrasts produce identical results to direct matrix inversion for the variance-covariance computation.

### Benchmark setup

Simulated datasets were generated using a self-contained JC69 sequence simulator (no external tools). Random binary trees were constructed via iterative taxon addition with branch lengths drawn from an exponential distribution, rescaled so mean root-to-tip distance equals 0.1 substitutions per site. The dataset matrix comprised 12 conditions (10-500 taxa, 500-5,000 sites) with 3 replicates each (36 total datasets). All tools were run single-threaded on Apple Silicon (M4): phylip-rs ML and NJ, IQ-TREE 3.0.1 (`-m JC -T 1`), RAxML-NG 1.2.2 (`--model JC --threads 1`), and VeryFastTree 4.0.5 (`-nt -nocat -threads 1`; VeryFastTree does not support JC69 directly and uses its default nucleotide model without rate categories). Wall time was measured via `time.perf_counter()`, peak memory via GNU time (`gtime`). Log-likelihood was evaluated under JC69 using phylip-rs for all inferred trees; we validated this scorer against IQ-TREE's internal likelihood calculator and confirmed agreement to four decimal places on shared topologies, confirming that scoring differences reflect tree search quality rather than evaluator bias (Supplementary Note 4). Topology accuracy was measured as normalized Robinson-Foulds distance to the true tree. The default timeout was 600 seconds (10 minutes) for all tools, with a reduced timeout of 120 seconds (2 minutes) for phylip-rs ML on datasets with 200 or more taxa, where exhaustive NNI search over the full topology space is impractical; this asymmetry means phylip-rs ML is given less time on larger datasets, and the timeout results should be interpreted accordingly.

### Catalog analysis

The 407 tools in Felsenstein's software catalog were extracted from HTML snapshots of the catalog pages (migrated to `phylipweb.github.io` in 2023). Metadata (release year, author, programming language, categories) was extracted via multi-pass regular expression matching against full tool descriptions. Release years were obtained for 261 tools (64%) through a combination of hand-curated dates for well-known tools, pattern matching for version strings and publication years, and extraction of the earliest year mentioned in each description. Tool preservation status was assessed by checking original URLs and cross-referencing against the Wayback Machine.

### Code and data availability

phylip-rs source code, benchmark scripts, simulated datasets, and catalog analysis code are available at https://github.com/shandley/phylip-archaeology. The benchmark results (180 runs) and all generated figures are included in the repository. phylip-rs requires only the Rust compiler (no external dependencies) and builds on all major platforms.

---

## References

1. Hinsen, K. The approximation tower in computational science: why testing scientific software is difficult. *Computing in Science & Engineering* **17**, 72-77 (2015).
2. Felsenstein, J. PHYLIP — Phylogeny Inference Package (Version 3.2). *Cladistics* **5**, 164-166 (1989).
3. Felsenstein, J. Confidence limits on phylogenies: an approach using the bootstrap. *Evolution* **39**, 783-791 (1985).
4. Van Noorden, R., Maher, B. & Nuzzo, R. The top 100 papers. *Nature* **514**, 550-553 (2014).
5. Felsenstein, J. Evolutionary trees from DNA sequences: a maximum likelihood approach. *J. Mol. Evol.* **17**, 368-376 (1981).
6. Pearl, J. *Probabilistic Reasoning in Intelligent Systems* (Morgan Kaufmann, 1988).
7. Felsenstein, J. Phylogenies and the comparative method. *Am. Nat.* **125**, 1-15 (1985).
8. Freeland, S. J. & Hurst, L. D. The genetic code is one in a million. *J. Mol. Evol.* **47**, 238-248 (1998).
9. Lockhart, P. J., Steel, M. A., Hendy, M. D. & Penny, D. Recovering evolutionary trees under a more realistic model of sequence evolution. *Mol. Biol. Evol.* **11**, 605-612 (1994).
10. Lake, J. A. Reconstructing evolutionary trees from DNA and protein sequences: paralinear distances. *Proc. Natl Acad. Sci. USA* **91**, 1455-1459 (1994).
11. Wong, T. K. F. et al. IQ-TREE 3: phylogenomic inference software using complex evolutionary models. *Preprint at* https://doi.org/10.32942/X2P62N (2025).
12. Kozlov, A. M., Darriba, D., Flouri, T., Morel, B. & Stamatakis, A. RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics* **35**, 4453-4455 (2019).
13. Piñeiro, C., Abuín, J. M. & Pichel, J. C. Very Fast Tree: speeding up the estimation of phylogenies for large alignments through parallelization and vectorization strategies. *Bioinformatics* **36**, 4658-4659 (2020).
14. Altschul, S. F. et al. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Res.* **25**, 3389-3402 (1997).
15. Yang, Z. PAML 4: phylogenetic analysis by maximum likelihood. *Mol. Biol. Evol.* **24**, 1586-1591 (2007).
16. Felsenstein, J. *Inferring Phylogenies* (Sinauer Associates, 2004).
17. Yang, Z. Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods. *J. Mol. Evol.* **39**, 306-314 (1994).
18. Saitou, N. & Nei, M. The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Mol. Biol. Evol.* **4**, 406-425 (1987).
19. Felsenstein, J. Cases in which parsimony or compatibility methods will be positively misleading. *Syst. Zool.* **27**, 401-410 (1978).
20. Fitch, W. M. Toward defining the course of evolution: minimum change for a specific tree topology. *Syst. Zool.* **20**, 406-416 (1971).
21. Sankoff, D. Minimal mutation trees of sequences. *SIAM J. Appl. Math.* **28**, 35-42 (1975).
22. Hendy, M. D. & Penny, D. Branch and bound algorithms to determine minimal evolutionary trees. *Math. Biosci.* **59**, 277-290 (1982).
23. Felsenstein, J. Maximum-likelihood estimation of evolutionary trees from continuous characters. *Am. J. Hum. Genet.* **25**, 471-492 (1973).
24. Felsenstein, J. & Churchill, G. A. A hidden Markov model approach to variation among sites in rate of evolution. *Mol. Biol. Evol.* **13**, 93-104 (1996).
25. Bron, C. & Kerbosch, J. Algorithm 457: finding all cliques of an undirected graph. *Commun. ACM* **16**, 575-577 (1973).
26. Lake, J. A. A rate-independent technique for analysis of nucleic acid sequences: evolutionary parsimony. *Mol. Biol. Evol.* **4**, 167-191 (1987).
27. Amari, S. *Differential-Geometrical Methods in Statistics* (Springer, 1985).
28. Cavalli-Sforza, L. L. & Edwards, A. W. F. Phylogenetic analysis: models and estimation procedures. *Am. J. Hum. Genet.* **19**, 233-257 (1967).
29. Robinson, D. F. & Foulds, L. R. Comparison of phylogenetic trees. *Math. Biosci.* **53**, 131-147 (1981).
30. Koller, D. & Friedman, N. *Probabilistic Graphical Models: Principles and Techniques* (MIT Press, 2009).

---

## Figure Legends

**Figure 1. LLM-assisted code archaeology workflow.** Schematic of the iterative process for extracting algorithmic knowledge from legacy scientific software. The workflow proceeds through five phases: source code reading and annotation by the LLM; algorithm extraction into mathematical notation; cross-referencing against published literature; reimplementation in a modern language with validation; and identification of cross-disciplinary connections. Human domain expertise directs the process at every stage.

**Figure 2. Cross-disciplinary connections recovered from PHYLIP source code.** (a) Felsenstein's pruning algorithm (1981) is mathematically identical to belief propagation, formalized by Pearl seven years later (1988). The tree diagram shows the postorder message-passing computation. (b) Independent contrasts variance propagation follows the parallel resistor formula from electrical circuit theory. Branch lengths correspond to resistances; the weighted average at each node follows Kirchhoff's voltage law. Numerical equivalence is exact to eight decimal places. (c) The genetic code step matrix: a 20x20 heatmap showing minimum nucleotide substitution costs between amino acids, derived from codon assignments. The z-score of -2.76 against 10,000 random codes confirms the genetic code is optimized for error tolerance. (d) LogDet distance isolates evolutionary signal by factoring out base composition bias via the determinant of the divergence matrix.

**Figure 3. Preservation status of Felsenstein's phylogenetics software catalog.** (a) Timeline of tool publications by five-year bins, stacked by methodological category, showing the shift from parsimony to maximum likelihood to Bayesian methods. Peak development occurred in 2000-2004 (89 tools). (b) Preservation status of all 407 tools: 196 archived, 137 dormant, 34 dead, 40 unknown. Twenty-three tools have no Wayback Machine archive and may be permanently lost. (c) Methodological paradigm shift: proportion of tools by category across eras.

**Figure 4. Benchmarking phylip-rs against modern phylogenetics tools.** (a) Log-likelihood (scored under JC69) versus wall time for each tool on each dataset. phylip-rs ML overlaps with modern tools on small datasets but diverges at larger scales. (b) Normalized Robinson-Foulds distance to the true tree by dataset size, showing comparable topology accuracy across tools where phylip-rs completes. (c) Scaling behavior (wall time versus number of taxa, log-log scale, 1,000-site datasets) with O(n^2) and O(n^3) reference lines. phylip-rs ML scales steeply beyond 50 taxa; NJ and VeryFastTree scale favorably to 500 taxa.

---

## Acknowledgments

We thank Joe Felsenstein for creating PHYLIP and maintaining the software catalog that made this work possible. His contributions to computational biology — the pruning algorithm, the bootstrap, independent contrasts, and the Newick format — are the foundations on which modern phylogenetics rests.

## Competing Interests

The author declares no competing interests.

## Funding

No specific funding was received for this work.

---

*Correspondence: shandley@wustl.edu*

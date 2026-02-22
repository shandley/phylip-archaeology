# Algorithms Outlive Implementations: Resurrecting Five Lost Phylogenetics Tools

**Scott A. Handley**

Department of Pathology and Immunology, Washington University School of Medicine, St. Louis, MO, USA

---

## Abstract

Scientific software encodes algorithmic knowledge that exists nowhere else — not in papers, not in textbooks, only in implementations that compile on systems no longer in use. When software disappears, the algorithms it embodied risk disappearing with it. We resurrected five phylogenetics tools — TipDate, the AU test, DIVA, TREEMAP, and PLATO — that were listed in Felsenstein's software catalog but whose original implementations have become unavailable. Each tool introduced an algorithmic idea that shaped its subfield and influenced successor methods still in wide use today, yet each died for a different reason: platform obsolescence, supersession, build system decay, or website failure. We reimplemented all five as modern Rust code within a shared zero-dependency library (6,300 lines, 128 tests), preserving their core algorithms in a form designed to outlast any particular platform. The resurrections illustrate a general pattern: algorithms are intellectual heritage, but without active preservation, the implementations that carry them have a half-life measured in years.

---

## Introduction

Software has a half-life. Programming languages fall out of use, operating systems drop backward compatibility, websites go dark, and the researchers who maintained foundational tools retire. The algorithms embedded in this software — often documented only in their implementations — become inaccessible. Unlike mathematical proofs, which persist in print, computational methods are tied to the platforms that run them. When the platform dies, the method dies with it.

This problem is quantifiable. Joe Felsenstein's hand-curated catalog of phylogenetics software — effectively a package manager for an entire scientific discipline, predating SourceForge, Bioconductor, and GitHub — documents 407 tools developed between the 1970s and the 2010s [1]. Our analysis of this catalog found that 196 tools (48%) are archived but unmaintained, 34 (8%) are completely unreachable, and 23 have no archived copy anywhere — not even in the Wayback Machine [2]. These 23 tools may be permanently lost.

The catalog also reveals when tools are most vulnerable. Publication peaked during 2000-2004 with 89 new tools, coinciding with the genomics revolution. Tools from this era — written for platforms that predate modern version control and package management, but too recent to have been preserved in print documentation — are disproportionately at risk. They occupy a preservation gap: too old for GitHub, too new for the filing cabinet.

Yet the algorithms in these tools are not historical curiosities. They introduced ideas that became foundational: tip-dating for viral phylodynamics, multiscale bootstrap for hypothesis testing, event-based biogeographic reconstruction, cophylogenetic reconciliation analysis, and likelihood-based recombination detection. In every case, the algorithmic idea survived — absorbed into successor methods — while the original implementation vanished. The question is whether reimplementation can serve as a form of algorithmic conservation, preserving not just the idea but the implementation details that papers leave out.

As part of a broader project applying LLM-assisted code archaeology to Felsenstein's PHYLIP package [2], we identified five lost tools from the software catalog whose algorithms remained influential but whose implementations had become unavailable. We reimplemented all five within phylip-rs, a zero-dependency Rust library covering 29 of PHYLIP's 36 programs [2]. Here we describe what each tool did, why it mattered, how it was lost, and what we recovered.

---

## Five Algorithms, Five Extinctions, Five Resurrections

### TipDate: sampling times as free calibration (Rambaut, 2000)

Before TipDate [3], molecular clock calibration required fossils or biogeographic events — external information that was unavailable for many organisms, particularly viruses. Andrew Rambaut recognized that when sequences are sampled at different calendar dates, the sampling times themselves serve as calibration points. The algorithm roots a tree, assigns known dates to tips, and optimizes a substitution rate and internal node dates by maximum likelihood, comparing three nested models via likelihood ratio test: free branch lengths, a strict clock, and a dated-tip (SRDT) clock. A significant LRT between the SRDT and strict clock models indicates that the sequences contain measurable temporal signal.

This idea became the foundation for Bayesian molecular dating. BEAST [4] and BEAST2 [5] generalized TipDate's maximum likelihood framework to full Bayesian inference with relaxed clocks, coalescent priors, and model averaging — but the core insight (sampling times are free calibrations) originated here. TipDate's original Java application is no longer available from its distribution site.

**Resurrected**: 1,307 lines of Rust, 17 tests. Implements golden section search optimization, root-to-tip regression initialization, and likelihood ratio testing against the strict molecular clock. Validated on two-taxon analytical cases and five-taxon heterochronous trees with known temporal signal.

### AU test: correcting selection bias in tree comparison (Shimodaira, 2002)

Standard bootstrap proportions are biased when the candidate trees being compared were selected because they fit the data well — a form of selection bias that inflates confidence in the best tree. Hidetoshi Shimodaira's Approximately Unbiased (AU) test [6] corrects this bias through multiscale bootstrap resampling: resample at multiple scale factors (50%-140% of original sites), fit a two-parameter model to rejection z-values across scales, and extract a bias-corrected p-value. The two parameters capture variance (genuine signal) and bias (selection artifact) separately.

The AU test remains the gold standard for phylogenetic hypothesis testing, cited in virtually every study that compares alternative tree topologies. The R package scaleboot survives, but the standalone C implementation (CONSEL [7]) has become difficult to compile on modern systems, and the algorithm's complexity has discouraged independent reimplementation.

**Resurrected**: 1,237 lines of Rust, 20 tests. Implements the full AU test plus the KH test [8] and SH test [9] for comparison, confidence set construction, and from-scratch normal CDF (Abramowitz-Stegun approximation) and quantile (Acklam approximation with Newton-Raphson refinement) functions. Validated on best-tree identification, p-value summation properties, and round-trip CDF/quantile accuracy.

### DIVA: vicariance as the null hypothesis (Ronquist, 1997)

Before DIVA [10], historical biogeography was largely narrative — researchers told stories about how organisms came to occupy their present ranges. Fredrik Ronquist formalized the intuition that vicariance (geographic splitting during speciation) is the expected mode of speciation, while dispersal and extinction require explanation. The algorithm uses bottom-up dynamic programming on area bit-vectors: for each internal node, it enumerates all ways to partition an ancestral area set between daughter lineages, scoring vicariance at zero cost and charging one unit per dispersal or extinction event.

DIVA's event-cost framework directly influenced parametric biogeographic methods (DEC [11], BioGeoBEARS [12]) that now dominate the field. The original DOS executable and its website are no longer available.

**Resurrected**: 1,695 lines of Rust, 47 tests. Implements bit-vector area sets (u32, supporting up to 32 areas), subset enumeration, disjoint bipartition generation, configurable event costs, and maximum range constraints. Validated on perfect vicariance scenarios, Gondwanan continental fragmentation, and island stepping-stone colonization.

### TREEMAP: coevolution as a tree-mapping problem (Page, 1994)

Rod Page pioneered the quantitative analysis of coevolution by treating it as a reconciliation problem [13,14]: given a host phylogeny, a parasite phylogeny, and leaf-level associations (which parasite infects which host), find the minimum-cost mapping of parasite divergence events to host tree nodes. Each event is classified as cospeciation (host and parasite speciated together, cost 0), duplication (parasite speciated within a host lineage, cost 1), sorting (parasite lineage lost, cost 1), or host-switching (lateral transfer, cost 2).

This event vocabulary became the standard framework for studying host-parasite associations, gene-species reconciliation, and area-phylogeny relationships, directly influencing tools like Jane [15] and Notung [16]. TREEMAP was originally a Macintosh Classic application; the original binary is no longer available.

**Resurrected**: 955 lines of Rust, 18 tests. Implements LCA-based reconciliation by postorder traversal, configurable event costs, sorting event detection via depth differences, and per-node event classification. Validated on perfect cospeciation, single and complete duplication scenarios, and custom cost weightings.

### PLATO: likelihood scanning for mosaic genomes (Grassly & Holmes, 1997)

Recombination creates mosaic sequences where different genomic regions have different evolutionary histories — a violation of the assumption underlying every phylogenetic method that all sites share a single tree. Nick Grassly and Eddie Holmes developed PLATO [17] to detect such regions by sliding a window across an alignment and computing partial log-likelihoods under a fitted tree. Regions where the local likelihood deviates significantly from the genome-wide expectation are flagged as potentially recombinant, with formal significance assessed via parametric bootstrap.

PLATO was among the first likelihood-based recombination detection tools, influencing GARD [18] and other scanning approaches. Its original C implementation and website are no longer available.

**Resurrected**: 1,106 lines of Rust, 26 tests. Implements sliding-window likelihood scanning, z-score anomaly detection, full Markov chain sequence simulation along trees for parametric bootstrap, and configurable window size, step size, and threshold. Validated on homogeneous alignments (no false positives), concatenated alignments with known breakpoints (correct detection), and parametric bootstrap p-value calibration.

---

## Cross-Cutting Patterns

Three observations emerged from the five resurrections.

**Shared mathematical infrastructure.** All five tools build on a remarkably small set of primitives: Felsenstein's pruning algorithm for likelihood computation (used by TipDate, the AU test, and PLATO), dynamic programming on tree structures (DIVA, TREEMAP), and bootstrap resampling for null distribution construction (AU test, PLATO). Reimplementing them within a shared library revealed these connections explicitly. TipDate and PLATO, for instance, appear to address unrelated questions (molecular dating versus recombination detection), but both reduce to optimizing parameters of a likelihood function computed by the same pruning algorithm. The shared infrastructure — 42,105 lines covering 29 PHYLIP programs and these 5 resurrected tools, with zero external dependencies — means that each resurrected tool benefits from validated implementations of substitution models, tree I/O, and numerical optimization that it need not reimplement.

**A taxonomy of software death.** The five tools died for five distinct reasons: platform obsolescence (TREEMAP on Mac Classic, DIVA on DOS), supersession by a more general method (TipDate by BEAST), build system decay (CONSEL's C code no longer compiles cleanly), and hosting failure (PLATO's website went dark). These causes are not independent — tools on obsolete platforms are also less likely to be ported, and superseded tools are less likely to have their hosting maintained. But the taxonomy suggests that different preservation strategies are needed: platform-independent reimplementation addresses obsolescence, while archiving addresses hosting failure. The most vulnerable tools are those from 1994-2002 — old enough that their platforms are gone, but young enough that they were never distributed in print.

**What papers leave out.** In every case, the published paper described the algorithm at a level sufficient for reimplementation, but omitted implementation details that matter in practice: how to initialize the optimization (TipDate uses root-to-tip regression), how to handle numerical edge cases (the AU test requires careful treatment of z-values near zero), how to enumerate combinatorial structures efficiently (DIVA's subset enumeration via bitmask tricks). These details — the craft knowledge of scientific programming — are precisely what is lost when software disappears. Reimplementation from the paper alone recovers the algorithm; reimplementation with access to the original code (when available) recovers the craft.

---

## Conclusion

We resurrected five phylogenetics tools — 6,300 lines of Rust, 128 tests, zero external dependencies — preserving algorithms that shaped viral phylodynamics, statistical hypothesis testing, historical biogeography, coevolutionary analysis, and recombination detection. Each tool's algorithmic idea survived its implementation, living on in successor methods. But the implementation details — the numerical tricks, the edge case handling, the design decisions — were at risk of permanent loss.

Algorithms are intellectual heritage. Like biological species, they can go extinct when their habitat (the platform that runs them) disappears. Reimplementation in a modern, dependency-free, well-tested language is a form of conservation biology for computational ideas. The five tools resurrected here are available as part of phylip-rs at https://github.com/shandley/phylip-archaeology.

---

## Acknowledgments

We thank Joe Felsenstein for creating and maintaining the phylogenetics software catalog that made this work possible. The original algorithm designers — Andrew Rambaut, Hidetoshi Shimodaira, Fredrik Ronquist, Rod Page, Nick Grassly, and Eddie Holmes — created the ideas we sought to preserve. All reimplementation was performed using Claude (Anthropic) via the Claude Code command-line interface in an iterative human-AI collaboration; in accordance with journal policy, the LLM is not listed as an author but its contribution is disclosed here.

## Competing Interests

The author declares no competing interests.

## Funding

No specific funding was received for this work.

---

*Correspondence: shandley@wustl.edu*

---

## References

1. Felsenstein J. PHYLIP — Phylogeny Inference Package (Version 3.2). Cladistics. 1989;5:164-166.
2. Handley SA. LLM-assisted code archaeology recovers the algorithmic legacy of PHYLIP. [Companion manuscript].
3. Rambaut A. Estimating the rate of molecular evolution: incorporating non-contemporaneous sequences into maximum likelihood phylogenies. Bioinformatics. 2000;16(4):395-399.
4. Drummond AJ, Rambaut A. BEAST: Bayesian evolutionary analysis by sampling trees. BMC Evol Biol. 2007;7:214.
5. Bouckaert R, Heled J, Kuhnert D, Vaughan T, Wu CH, Xie D, et al. BEAST 2: a software platform for Bayesian evolutionary analysis. PLoS Comput Biol. 2014;10(4):e1003537.
6. Shimodaira H. An approximately unbiased test of phylogenetic tree selection. Syst Biol. 2002;51(3):492-508.
7. Shimodaira H, Hasegawa M. CONSEL: for assessing the confidence of phylogenetic tree selection. Bioinformatics. 2001;17(12):1246-1247.
8. Kishino H, Hasegawa M. Evaluation of the maximum likelihood estimate of the evolutionary tree topologies from DNA sequence data, and the branching order in Hominoidea. J Mol Evol. 1989;29:170-179.
9. Shimodaira H, Hasegawa M. Multiple comparisons of log-likelihoods with applications to phylogenetic inference. Mol Biol Evol. 1999;16:1114-1116.
10. Ronquist F. Dispersal-vicariance analysis: a new approach to the quantification of historical biogeography. Syst Biol. 1997;46:195-203.
11. Ree RH, Smith SA. Maximum likelihood inference of geographic range evolution by dispersal, local extinction, and cladogenesis. Syst Biol. 2008;57(1):4-14.
12. Matzke NJ. BioGeoBEARS: BioGeography with Bayesian (and likelihood) evolutionary analysis in R scripts. R package version 0.2.1. 2013.
13. Page RDM. Maps between trees and cladistic analysis of historical associations among genes, organisms, and areas. Syst Biol. 1994;43:58-77.
14. Page RDM. Parallel phylogenies: reconstructing the history of host-parasite assemblages. Cladistics. 1994;10:155-173.
15. Conow C, Fielder D, Ovadia Y, Libeskind-Hadas R. Jane: a new tool for the cophylogeny reconstruction problem. Algorithms Mol Biol. 2010;5:16.
16. Durand D, Halldorsson BV, Vernot B. A hybrid micro-macroevolutionary approach to gene tree reconstruction. J Comput Biol. 2006;13(2):320-335.
17. Grassly NC, Holmes EC. A likelihood method for the detection of selection and recombination using nucleotide sequences. Mol Biol Evol. 1997;14(3):239-247.
18. Kosakovsky Pond SL, Posada D, Gravenor MB, Woelk CH, Frost SDW. GARD: a genetic algorithm for recombination detection. Bioinformatics. 2006;22(24):3096-3098.

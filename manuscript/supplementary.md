# Supplementary Information

## LLM-Assisted Code Archaeology Recovers the Algorithmic Legacy of PHYLIP

Scott A. Handley

---

## Supplementary Note 1: Complete Algorithmic Case Studies

The 20 algorithmic insights recovered from PHYLIP source code, with full mathematical descriptions, implementation details, cross-disciplinary connections, and numerical demonstrations.

*Content: Full text of INSIGHTS.md (964 lines), covering:*

| # | Algorithm | Cross-disciplinary Connection |
|---|-----------|------------------------------|
| 1 | Pruning algorithm | Belief propagation (Pearl, 1988) |
| 2 | F84 closed-form P(t) | Lie algebra decomposition |
| 3 | Site-pattern compression | Column-oriented databases |
| 4 | Bootstrap weight vectors | Weighted resampling |
| 5 | Discrete gamma rates | Mixture models / log-sum-exp |
| 6 | Fitch parsimony | Bitwise set operations |
| 7 | Model selection (AIC/BIC) | Regularization theory |
| 8 | First-principles derivation | Lost computational discipline |
| 9 | Independent contrasts | Kirchhoff's circuit laws |
| 10 | Contml stereographic projection | Hellinger embedding / information geometry |
| 11 | Hendy-Penny supplement bound | Dual decomposition / Lagrangian relaxation |
| 12 | Dollo parsimony | Min-cut on trees |
| 13 | LogDet distance | Determinantal factorization |
| 14 | Kitsch scrunch | Pool-adjacent-violators / isotonic regression |
| 15 | Clique analysis | Bron-Kerbosch / Splits Equivalence |
| 16 | Lake's invariants | Algebraic statistics / variety ideals |
| 17 | O(n) Brownian ML | Tree-structured Gaussian processes |
| 18 | Felsenstein-Churchill HMM | Baum-Welch for rate correlation |
| 19 | Protein Sankoff | Weighted dynamic programming / genetic code optimization |
| 20 | Score-ordered B&B | Greedy-guided exact search (A*, alpha-beta) |

---

## Supplementary Table 1: Software Catalog

Complete listing of all 407 phylogenetics tools from Felsenstein's catalog with preservation status, release year (where known), programming language, methodological categories, and URL status.

*Content: Tabulated from catalog/analysis/tools_enriched.json*

---

## Supplementary Table 2: Benchmark Results

Complete results for all 180 benchmark runs (36 datasets x 5 tools), including wall time, peak memory, log-likelihood (scored under JC69), Robinson-Foulds distance to true tree, and timeout status.

*Content: benchmarks/results/benchmark_results.csv*

---

## Supplementary Table 3: phylip-rs Module Summary

| Module | Files | Lines | Tests | Description |
|--------|-------|-------|-------|-------------|
| likelihood/ | 10 | 7,831 | 205 | Pruning algorithm, ML search, NNI, gamma rates, model selection, clock ML |
| models/ | 10 | 6,905 | 219 | JC69, K2P, F81, F84, protein (WAG/Poisson), LogDet, gene frequency distances |
| parsimony/ | 8 | 5,884 | 134 | Wagner, Dollo, Camin-Sokal, branch-and-bound, Sankoff, protein parsimony |
| distance/ | 6 | 3,010 | 71 | Neighbor-joining, UPGMA, Fitch-Margoliash, Kitsch, ML pairwise distances |
| io/ | 5 | 2,582 | 80 | PHYLIP format, FASTA, binary data, output reports |
| tree/ | 5 | 1,969 | 70 | Tree types, Newick I/O, splits (bitvectors), Robinson-Foulds distance |
| comparative/ | 3 | 1,927 | 34 | Independent contrasts, Brownian motion ML (contml) |
| compatibility/ | 3 | 1,758 | 37 | Clique analysis (Bron-Kerbosch), DNA compatibility |
| bootstrap/ | 3 | 1,355 | 41 | Resampling, consensus trees, ML bootstrap |
| invariants/ | 2 | 830 | 18 | Lake's and Cavender's phylogenetic invariants |
| consensus/ | 1 | 814 | 25 | Strict/majority-rule/extended consensus |
| **Total** | **60** | **36,745** | **961** | |

---

## Supplementary Note 2: Interactive Demonstration Outputs

Output from 10 compilable example programs demonstrating cross-disciplinary applications of PHYLIP algorithms. Each demonstration is self-contained and requires only the Rust compiler to build and run.

1. **felsenstein_zone** — Statistical consistency: parsimony accuracy drops to 0% as data increases (long branch attraction), while ML accuracy rises to 100%
2. **kirchhoff_contrasts** — Circuit theory: variance propagation matches parallel resistor formula to 8 decimal places
3. **genetic_code_distances** — Coding theory: z-score = -2.76 vs 10,000 random codes
4. **compositional_bias** — LogDet correctly recovers tree under compositional heterogeneity where JC69 fails
5. **clock_constraints** — Isotonic regression: pool-adjacent-violators enforces ultrametric constraint on ML tree
6. **chord_geometry** — Information geometry: Cavalli-Sforza chord = Euclidean distance after Hellinger embedding
7. **dollo_gain_loss** — Combinatorial optimization: Dollo parsimony infers gene gain/loss history
8. **supplement_bound** — Dual decomposition: Hendy-Penny bound prunes 99.7% of search space
9. **lake_invariants** — Algebraic statistics: polynomial invariants distinguish topologies from site-pattern frequencies
10. **language_evolution** — Pruning algorithm applied to Indo-European cognate data recovers known language family tree

---

## Supplementary Note 3: Dataset Generation

Simulated datasets were generated using a self-contained JC69 sequence simulator with no external dependencies. Random binary trees were constructed via iterative taxon addition with branch lengths drawn from an exponential distribution, rescaled so that the mean root-to-tip distance equals 0.1 substitutions per site. Sequences were evolved down the tree using JC69 transition probabilities: P(same) = 1/4 + 3/4 * exp(-4t/3), P(diff) = 1/4 - 1/4 * exp(-4t/3). Each replicate used a deterministic seed derived from the condition parameters for full reproducibility.

| Taxa | Sites | Replicates | Purpose |
|------|-------|------------|---------|
| 10 | 500, 1000 | 3 | Small baseline |
| 20 | 500, 1000, 5000 | 3 | Medium |
| 50 | 1000, 5000 | 3 | Scaling transition |
| 100 | 1000, 5000 | 3 | Large |
| 200 | 1000, 5000 | 3 | Scalability probe |
| 500 | 1000 | 3 | Boundary test |

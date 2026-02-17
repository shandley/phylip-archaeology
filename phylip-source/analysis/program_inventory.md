# PHYLIP Program Inventory

**Version:** 3.698 (September 2019)
**Author:** Joseph Felsenstein, University of Washington
**Source location:** `phylip-source/src/download/phylip-3.698.zip` (extracted to `/tmp/phylip-3.698-extracted/phylip-3.698/src/`)
**License:** BSD 2-Clause

## Overview

PHYLIP (Phylogeny Inference Package) contains **35 executable programs** built from **49 C source files** (35 program files + 14 library/support files). The package has been distributed since 1980 and can infer phylogenies by parsimony, compatibility, distance matrix methods, and maximum likelihood. It handles nucleotide sequences, protein sequences, gene frequencies, restriction sites, distance matrices, discrete characters, and continuous characters.

**Total source code:** 75,820 lines of C across 49 `.c` files, plus 1,691 lines across 13 `.h` header files.

---

## Source File Architecture

### Library/Support Files (shared code, no main function)

These files provide shared functionality used by multiple programs:

| File | Lines | Bytes | Description |
|------|------:|------:|-------------|
| `phylip.c` / `phylip.h` | 3,207 + 720 | 87,953 + 25,924 | Core library: I/O, memory, random numbers, tree data structures, common utilities |
| `seq.c` / `seq.h` | 4,155 + 234 | 118,377 + 9,108 | Molecular sequence support: DNA/protein input, site handling, tree operations for sequence data |
| `disc.c` / `disc.h` | 926 + 109 | 27,161 + 4,017 | Discrete character support: input, tree operations, output for 0/1 character data |
| `discrete.c` / `discrete.h` | 3,109 + 170 | 89,827 + 7,114 | Extended discrete character support: multistate characters, tree operations |
| `dollo.c` / `dollo.h` | 394 + 53 | 12,353 + 2,252 | Dollo parsimony support: character state reconstruction, tree scoring |
| `wagner.c` / `wagner.h` | 516 + 56 | 15,772 + 2,504 | Wagner parsimony support: mixed method character handling, tree scoring |
| `dist.c` / `dist.h` | 561 + 59 | 14,647 + 2,296 | Distance matrix support: I/O, tree operations for distance data |
| `cont.c` / `cont.h` | 269 + 45 | 8,430 + 1,890 | Continuous character support: tree allocation, view operations |
| `cons.c` / `cons.h` | 1,557 + 59 | 42,094 + 1,606 | Consensus tree support: partition computation, tree reconstruction |
| `draw.c` / `draw.h` | 3,318 + 107 | 98,030 + 3,709 | Tree drawing core: graphics primitives, output formats (PostScript, SVG, etc.) |
| `draw2.c` | 1,527 | 48,930 | Tree drawing extensions: splines, curves, BMP output |
| `mlclock.c` / `mlclock.h` | 584 + 39 | 15,961 + 1,382 | Molecular clock support: node time estimation for clock-constrained ML |
| `moves.c` / `moves.h` | 297 + 34 | 7,722 + 1,049 | Interactive tree manipulation: cursor movement, display utilities |
| `printree.c` / `printree.h` | 188 + 6 | 4,557 + 156 | Tree printing for ML clock programs: ASCII tree display |

---

## Program Inventory by Category

### 1. DNA/RNA Sequence Programs

Programs that analyze nucleotide sequence data.

| Program | Source File | Lines | Bytes | Algorithm | Description |
|---------|-----------|------:|------:|-----------|-------------|
| **Dnapars** | `dnapars.c` | 1,663 | 51,531 | Wagner parsimony (heuristic) | Unrooted parsimony for DNA sequences. Uses Fitch (1971) method to count base changes. Treats gaps as a fifth state. Supports transversion parsimony, ancestral state reconstruction, character weights, and branch length inference. |
| **Dnacomp** | `dnacomp.c` | 1,179 | 34,315 | Compatibility (heuristic) | DNA compatibility method. Finds the largest set of sites that are compatible with the same tree. Particularly appropriate when rates vary greatly among sites. |
| **Dnapenny** | `dnapenny.c` | 819 | 24,556 | Parsimony (branch-and-bound) | Exact branch-and-bound search for all most parsimonious trees for nucleic acid sequences. Practical for up to ~10-11 species. |
| **Dnamove** | `dnamove.c` | 2,363 | 62,473 | Parsimony (interactive) | Interactive DNA parsimony. Allows manual tree construction and evaluation by parsimony and compatibility criteria, with display of reconstructed ancestral bases. |
| **Dnaml** | `dnaml.c` | 2,619 | 76,609 | Maximum likelihood (heuristic) | DNA maximum likelihood. Implements the HKY85/F84 substitution model with unequal base frequencies and different transition/transversion rates. Supports gamma-distributed rates, hidden Markov model of rate variation, and invariant sites. Uses Felsenstein-Churchill (1996) algorithm. |
| **Dnamlk** | `dnamlk.c` | 2,252 | 66,399 | Maximum likelihood with clock (heuristic) | DNA ML with molecular clock constraint. Same model as Dnaml but enforces ultrametric tree (all tips equidistant from root). |
| **Dnadist** | `dnadist.c` | 1,335 | 39,064 | Distance computation | Computes pairwise distances from DNA sequences. Models: Jukes-Cantor, Kimura 2-parameter, F84, LogDet. Supports gamma and gamma+invariant sites rate correction. |
| **Dnainvar** | `dnainvar.c` | 828 | 26,653 | Phylogenetic invariants | For exactly 4 species: computes Lake's and Cavender's phylogenetic invariants to test alternative tree topologies. Also tabulates nucleotide pattern frequencies. |

**Dependencies:** All DNA programs link against `seq.c` and `phylip.c`. Dnamlk also uses `mlclock.c` and `printree.c`. Dnamove also uses `moves.c`.

---

### 2. Protein Sequence Programs

Programs that analyze amino acid sequence data.

| Program | Source File | Lines | Bytes | Algorithm | Description |
|---------|-----------|------:|------:|-----------|-------------|
| **Protpars** | `protpars.c` | 1,962 | 55,652 | Parsimony (heuristic) | Protein parsimony. Counts only nucleotide changes that alter the amino acid (assumes silent changes occur more easily). Uses genetic code to determine which amino acid changes require minimum nucleotide substitutions. |
| **Proml** | `proml.c` | 3,246 | 108,014 | Maximum likelihood (heuristic) | Protein maximum likelihood. Uses Jones-Taylor-Thornton (JTT) or Dayhoff PAM probability models. Supports gamma-distributed rates, HMM rate variation, and invariant sites. |
| **Promlk** | `promlk.c` | 2,998 | 99,714 | Maximum likelihood with clock (heuristic) | Protein ML with molecular clock. Same models as Proml but enforces ultrametric tree. |
| **Protdist** | `protdist.c` | 2,006 | 71,189 | Distance computation | Computes protein distances. Models: JTT, PMB (Henikoff/Tillier), Dayhoff PAM, Kimura (1983) approximation, categories model. Supports gamma and gamma+invariant sites correction. |

**Dependencies:** All protein programs link against `seq.c` and `phylip.c`. Promlk also uses `mlclock.c` and `printree.c`.

---

### 3. Restriction Sites Programs

Programs for restriction enzyme site data.

| Program | Source File | Lines | Bytes | Algorithm | Description |
|---------|-----------|------:|------:|-----------|-------------|
| **Restml** | `restml.c` | 2,528 | 67,819 | Maximum likelihood (heuristic) | Restriction sites ML using the Jukes-Cantor symmetrical model of nucleotide change. Very slow. Does not distinguish transitions from transversions. |
| **Restdist** | `restdist.c` | 690 | 18,196 | Distance computation | Computes distances from restriction sites or restriction fragments data. Based on Nei and Li (1979) with modifications. Can also handle RAPD and AFLP data. |

**Dependencies:** Link against `seq.c` and `phylip.c`.

---

### 4. Distance Matrix Programs

Programs that work with pairwise distance matrices.

| Program | Source File | Lines | Bytes | Algorithm | Description |
|---------|-----------|------:|------:|-----------|-------------|
| **Fitch** | `fitch.c` | 1,203 | 32,449 | Fitch-Margoliash / Least Squares / Minimum Evolution (heuristic) | Distance-based tree estimation under the additive tree model. Implements Fitch-Margoliash criterion, least squares, and minimum evolution. No evolutionary clock assumed. Supports negative branch lengths optionally. |
| **Kitsch** | `kitsch.c` | 1,028 | 29,724 | Fitch-Margoliash with clock (heuristic) | Same as Fitch but with evolutionary clock constraint (ultrametric tree). All tip species assumed contemporaneous. |
| **Neighbor** | `neighbor.c` | 629 | 17,365 | Neighbor-Joining / UPGMA | Implements Saitou and Nei (1987) Neighbor-Joining method and UPGMA clustering. NJ does not assume a clock; UPGMA does. Written by Mary Kuhner and Jon Yamato. |

**Dependencies:** Link against `dist.c` and `phylip.c`.

---

### 5. Gene Frequency / Continuous Character Programs

Programs for allele frequency and quantitative trait data.

| Program | Source File | Lines | Bytes | Algorithm | Description |
|---------|-----------|------:|------:|-----------|-------------|
| **Contml** | `contml.c` | 1,560 | 42,769 | Maximum likelihood (heuristic) | ML for gene frequencies and continuous characters under Brownian motion / genetic drift model. Based on Edwards and Cavalli-Sforza (1964) model. Does not assume a molecular clock. |
| **Gendist** | `gendist.c` | 443 | 11,418 | Distance computation | Computes genetic distances from gene frequency data. Three measures: Nei's genetic distance (1972), Cavalli-Sforza chord measure (1967), Reynolds-Weir-Cockerham distance (1983). |
| **Contrast** | `contrast.c` | 964 | 30,226 | Comparative method | Computes Felsenstein's (1985) phylogenetically independent contrasts. Reads a tree and continuous character data, produces independent contrasts for multivariate statistical analysis. Can correct for within-species variation. |

**Dependencies:** Link against `cont.c` and `phylip.c`.

---

### 6. Discrete Character (0,1) Programs

Programs for binary character data (morphological, restriction site presence/absence).

| Program | Source File | Lines | Bytes | Algorithm | Description |
|---------|-----------|------:|------:|-----------|-------------|
| **Pars** | `pars.c` | 1,649 | 52,193 | Wagner parsimony, multistate (heuristic) | General parsimony for discrete characters with up to 8 states. Wagner criterion (minimum changes). Supports multifurcations, ancestral states, character weights, branch lengths. |
| **Mix** | `mix.c` | 1,180 | 33,350 | Wagner / Camin-Sokal parsimony (heuristic) | Mixed method parsimony for 0/1 characters. Each character can independently use Wagner (reversible) or Camin-Sokal (irreversible: only 0->1) parsimony. Defaults to Wagner. |
| **Move** | `move.c` | 1,655 | 44,024 | Wagner / Camin-Sokal parsimony (interactive) | Interactive version of Mix. Manual tree construction and evaluation with display of reconstructed states. |
| **Penny** | `penny.c` | 843 | 25,353 | Wagner / Camin-Sokal parsimony (branch-and-bound) | Exact search for all most parsimonious trees for 0/1 characters using branch-and-bound. Supports Wagner, Camin-Sokal, and mixed criteria. Practical for ~10-11 species. |
| **Dollop** | `dollop.c` | 1,023 | 30,379 | Dollo / Polymorphism parsimony (heuristic) | Dollo parsimony: assumes complex features are gained once and lost multiple times (irreversible gain). Also implements polymorphism parsimony. |
| **Dolmove** | `dolmove.c` | 1,602 | 42,005 | Dollo / Polymorphism parsimony (interactive) | Interactive version of Dollop. Manual tree construction with Dollo or polymorphism criteria. |
| **Dolpenny** | `dolpenny.c` | 736 | 21,861 | Dollo / Polymorphism parsimony (branch-and-bound) | Exact branch-and-bound search for Dollo or polymorphism parsimony criteria. Practical for ~10-11 species. |
| **Clique** | `clique.c` | 1,532 | 41,006 | Compatibility / Clique (heuristic) | Finds the largest clique of mutually compatible 0/1 characters and the trees they imply. Based on Le Quesne (1969), Estabrook et al. (1976). |

**Dependencies:** Mix, Move, Penny link against `disc.c`, `wagner.c`, `phylip.c`. Move also uses `moves.c`. Dollop, Dolmove, Dolpenny link against `disc.c`, `dollo.c`, `phylip.c`. Dolmove also uses `moves.c`. Pars links against `discrete.c`, `phylip.c`. Clique links against `disc.c`, `phylip.c`.

---

### 7. Tree Utility Programs

Programs for tree visualization, manipulation, consensus, and comparison.

| Program | Source File | Lines | Bytes | Algorithm | Description |
|---------|-----------|------:|------:|-----------|-------------|
| **Drawgram** | `drawgram.c` | 2,179 | 66,877 | Tree plotting | Plots rooted phylogenies: cladograms, phenograms, circular trees. Supports PostScript, SVG, PCL, BMP, and other output formats. Interactive preview. |
| **Drawtree** | `drawtree.c` | 3,151 | 97,134 | Tree plotting | Plots unrooted tree diagrams. Same output format support as Drawgram. Uses equal-daylight algorithm for branch angle optimization. |
| **Consense** | `consense.c` | 443 | 11,868 | Consensus tree | Computes majority-rule consensus trees (Ml methods of Margush and McMorris, 1981). Includes strict consensus. Reads multiple trees from a tree file. |
| **Treedist** | `treedist.c` | 1,298 | 41,609 | Tree distance | Computes distances between trees: Branch Score Distance (Kuhner and Felsenstein, 1994) using branch lengths, and Symmetric Difference (Robinson and Foulds, 1981) based on topology. |
| **Retree** | `retree.c` | 3,329 | 86,050 | Tree editing (interactive) | Interactive tree editor. Reads, displays, and allows rearrangement, rerooting, and output of trees. Can midpoint-root trees. Does not evaluate any optimality criterion. |

**Dependencies:** Drawgram and Drawtree link against `draw.c`, `draw2.c`, `phylip.c`. Consense and Treedist link against `cons.c`, `phylip.c`. Retree links against `moves.c`, `phylip.c`.

---

### 8. Data Transformation Programs

Programs that transform data for use by other programs.

| Program | Source File | Lines | Bytes | Algorithm | Description |
|---------|-----------|------:|------:|-----------|-------------|
| **Seqboot** | `seqboot.c` | 1,683 | 48,077 | Bootstrap / Jackknife / Permutation | General resampling tool. Creates multiple datasets by bootstrap, delete-half jackknife, or permutation. Works with molecular sequences, restriction sites, gene frequencies, and discrete characters. Enables statistical confidence assessment via bootstrap proportions. |
| **Factor** | `factor.c` | 594 | 17,577 | Character factoring | Converts multistate characters into binary (0,1) characters for use with discrete character programs. Also provides a way to delete characters and recode data. |

**Dependencies:** Seqboot links against `seq.c`, `phylip.c`. Factor links against `phylip.c`.

---

### 9. Supplementary Programs (not in main distribution)

| Program | Source File | Lines | Bytes | Algorithm | Description |
|---------|-----------|------:|------:|-----------|-------------|
| **Threshml** | `threshml.c` (in `download/threshml/`) | 1,714 | 52,130 | Threshold ML (MCMC) | Maximum likelihood for threshold characters (discrete characters modeled as underlying continuous liabilities). Uses MCMC to estimate covariance matrices. Not part of the standard PHYLIP build. |

---

## Complete Source File Summary Table

### Program Source Files (with main function)

| # | Program Name | File | Lines | Bytes | Data Type | Algorithm Category | Library Dependencies |
|---|-------------|------|------:|------:|-----------|-------------------|---------------------|
| 1 | Clique | `clique.c` | 1,532 | 41,006 | Discrete (0,1) | Compatibility | disc, phylip |
| 2 | Consense | `consense.c` | 443 | 11,868 | Trees | Consensus | cons, phylip |
| 3 | Contml | `contml.c` | 1,560 | 42,769 | Gene freq / Continuous | Maximum Likelihood | cont, phylip |
| 4 | Contrast | `contrast.c` | 964 | 30,226 | Continuous | Comparative Method | cont, phylip |
| 5 | Dnacomp | `dnacomp.c` | 1,179 | 34,315 | DNA | Compatibility | seq, phylip |
| 6 | Dnadist | `dnadist.c` | 1,335 | 39,064 | DNA | Distance Computation | seq, phylip |
| 7 | Dnainvar | `dnainvar.c` | 828 | 26,653 | DNA | Invariants | seq, phylip |
| 8 | Dnaml | `dnaml.c` | 2,619 | 76,609 | DNA | Maximum Likelihood | seq, phylip |
| 9 | Dnamlk | `dnamlk.c` | 2,252 | 66,399 | DNA | ML + Clock | seq, mlclock, printree, phylip |
| 10 | Dnamove | `dnamove.c` | 2,363 | 62,473 | DNA | Parsimony (interactive) | seq, moves, phylip |
| 11 | Dnapars | `dnapars.c` | 1,663 | 51,531 | DNA | Parsimony | seq, phylip |
| 12 | Dnapenny | `dnapenny.c` | 819 | 24,556 | DNA | Parsimony (B&B) | seq, phylip |
| 13 | Dollop | `dollop.c` | 1,023 | 30,379 | Discrete (0,1) | Dollo Parsimony | disc, dollo, phylip |
| 14 | Dolmove | `dolmove.c` | 1,602 | 42,005 | Discrete (0,1) | Dollo Parsimony (interactive) | disc, dollo, moves, phylip |
| 15 | Dolpenny | `dolpenny.c` | 736 | 21,861 | Discrete (0,1) | Dollo Parsimony (B&B) | disc, dollo, phylip |
| 16 | Drawgram | `drawgram.c` | 2,179 | 66,877 | Trees | Visualization | draw, draw2, phylip |
| 17 | Drawtree | `drawtree.c` | 3,151 | 97,134 | Trees | Visualization | draw, draw2, phylip |
| 18 | Factor | `factor.c` | 594 | 17,577 | Discrete (multistate) | Data Transformation | phylip |
| 19 | Fitch | `fitch.c` | 1,203 | 32,449 | Distance matrix | Least Squares | dist, phylip |
| 20 | Gendist | `gendist.c` | 443 | 11,418 | Gene frequencies | Distance Computation | phylip |
| 21 | Kitsch | `kitsch.c` | 1,028 | 29,724 | Distance matrix | Least Squares + Clock | dist, phylip |
| 22 | Mix | `mix.c` | 1,180 | 33,350 | Discrete (0,1) | Wagner/Camin-Sokal Parsimony | disc, wagner, phylip |
| 23 | Move | `move.c` | 1,655 | 44,024 | Discrete (0,1) | Wagner/Camin-Sokal (interactive) | disc, wagner, moves, phylip |
| 24 | Neighbor | `neighbor.c` | 629 | 17,365 | Distance matrix | Neighbor-Joining / UPGMA | dist, phylip |
| 25 | Pars | `pars.c` | 1,649 | 52,193 | Discrete (multistate) | Wagner Parsimony | discrete, phylip |
| 26 | Penny | `penny.c` | 843 | 25,353 | Discrete (0,1) | Wagner/Camin-Sokal (B&B) | disc, wagner, phylip |
| 27 | Proml | `proml.c` | 3,246 | 108,014 | Protein | Maximum Likelihood | seq, phylip |
| 28 | Promlk | `promlk.c` | 2,998 | 99,714 | Protein | ML + Clock | seq, mlclock, printree, phylip |
| 29 | Protdist | `protdist.c` | 2,006 | 71,189 | Protein | Distance Computation | seq, phylip |
| 30 | Protpars | `protpars.c` | 1,962 | 55,652 | Protein | Parsimony | seq, phylip |
| 31 | Restdist | `restdist.c` | 690 | 18,196 | Restriction sites | Distance Computation | seq, phylip |
| 32 | Restml | `restml.c` | 2,528 | 67,819 | Restriction sites | Maximum Likelihood | seq, phylip |
| 33 | Retree | `retree.c` | 3,329 | 86,050 | Trees | Tree Editing | moves, phylip |
| 34 | Seqboot | `seqboot.c` | 1,683 | 48,077 | All sequence types | Resampling | seq, phylip |
| 35 | Treedist | `treedist.c` | 1,298 | 41,609 | Trees | Tree Comparison | cons, phylip |

### Library/Support Files (no main function)

| # | File | Header | Lines (.c) | Lines (.h) | Bytes (.c) | Bytes (.h) | Purpose |
|---|------|--------|----------:|----------:|----------:|----------:|---------|
| 1 | `phylip.c` | `phylip.h` | 3,207 | 720 | 87,953 | 25,924 | Core library |
| 2 | `seq.c` | `seq.h` | 4,155 | 234 | 118,377 | 9,108 | Sequence data operations |
| 3 | `disc.c` | `disc.h` | 926 | 109 | 27,161 | 4,017 | Discrete character operations |
| 4 | `discrete.c` | `discrete.h` | 3,109 | 170 | 89,827 | 7,114 | Multistate discrete operations |
| 5 | `dollo.c` | `dollo.h` | 394 | 53 | 12,353 | 2,252 | Dollo parsimony operations |
| 6 | `wagner.c` | `wagner.h` | 516 | 56 | 15,772 | 2,504 | Wagner parsimony operations |
| 7 | `dist.c` | `dist.h` | 561 | 59 | 14,647 | 2,296 | Distance matrix operations |
| 8 | `cont.c` | `cont.h` | 269 | 45 | 8,430 | 1,890 | Continuous character operations |
| 9 | `cons.c` | `cons.h` | 1,557 | 59 | 42,094 | 1,606 | Consensus tree operations |
| 10 | `draw.c` | `draw.h` | 3,318 | 107 | 98,030 | 3,709 | Drawing core |
| 11 | `draw2.c` | -- | 1,527 | -- | 48,930 | -- | Drawing extensions |
| 12 | `mlclock.c` | `mlclock.h` | 584 | 39 | 15,961 | 1,382 | Molecular clock support |
| 13 | `moves.c` | `moves.h` | 297 | 34 | 7,722 | 1,049 | Interactive tree moves |
| 14 | `printree.c` | `printree.h` | 188 | 6 | 4,557 | 156 | ML tree printing |

---

## Programs Sorted by Code Complexity (lines of code, descending)

| Rank | Program | Lines | Category |
|-----:|---------|------:|----------|
| 1 | Retree | 3,329 | Tree Editing |
| 2 | Proml | 3,246 | Protein ML |
| 3 | Drawtree | 3,151 | Visualization |
| 4 | Promlk | 2,998 | Protein ML+Clock |
| 5 | Dnaml | 2,619 | DNA ML |
| 6 | Restml | 2,528 | Restriction Sites ML |
| 7 | Dnamove | 2,363 | DNA Parsimony (interactive) |
| 8 | Dnamlk | 2,252 | DNA ML+Clock |
| 9 | Drawgram | 2,179 | Visualization |
| 10 | Protdist | 2,006 | Protein Distance |
| 11 | Protpars | 1,962 | Protein Parsimony |
| 12 | Seqboot | 1,683 | Resampling |
| 13 | Dnapars | 1,663 | DNA Parsimony |
| 14 | Move | 1,655 | Discrete Parsimony (interactive) |
| 15 | Pars | 1,649 | Discrete Parsimony |
| 16 | Dolmove | 1,602 | Dollo Parsimony (interactive) |
| 17 | Contml | 1,560 | Gene Freq ML |
| 18 | Clique | 1,532 | Compatibility |
| 19 | Dnadist | 1,335 | DNA Distance |
| 20 | Treedist | 1,298 | Tree Comparison |
| 21 | Fitch | 1,203 | Distance Least Squares |
| 22 | Mix | 1,180 | Wagner/Camin-Sokal Parsimony |
| 23 | Dnacomp | 1,179 | DNA Compatibility |
| 24 | Kitsch | 1,028 | Distance LS+Clock |
| 25 | Dollop | 1,023 | Dollo Parsimony |
| 26 | Contrast | 964 | Comparative Method |
| 27 | Penny | 843 | Parsimony (B&B) |
| 28 | Dnainvar | 828 | Invariants |
| 29 | Dnapenny | 819 | DNA Parsimony (B&B) |
| 30 | Dolpenny | 736 | Dollo Parsimony (B&B) |
| 31 | Restdist | 690 | Restriction Distance |
| 32 | Neighbor | 629 | Neighbor-Joining/UPGMA |
| 33 | Factor | 594 | Data Transformation |
| 34 | Consense | 443 | Consensus Trees |
| 35 | Gendist | 443 | Genetic Distance |

---

## Data Type Cross-Reference

| Data Type | Programs |
|-----------|----------|
| **DNA sequences** | Dnapars, Dnacomp, Dnapenny, Dnamove, Dnaml, Dnamlk, Dnadist, Dnainvar |
| **Protein sequences** | Protpars, Proml, Promlk, Protdist |
| **Restriction sites** | Restml, Restdist |
| **Distance matrices** | Fitch, Kitsch, Neighbor |
| **Gene frequencies** | Contml, Gendist |
| **Continuous characters** | Contml, Contrast |
| **Discrete characters (0,1)** | Mix, Move, Penny, Dollop, Dolmove, Dolpenny, Clique |
| **Discrete characters (multistate)** | Pars, Factor |
| **Trees** | Drawgram, Drawtree, Consense, Treedist, Retree |
| **All data types** | Seqboot (resampling) |

---

## Build System Notes

- Build is via `Makefile.unx` (Unix), `Makefile.osx` (macOS), or `Makefile.cyg` (Cygwin)
- Compiler: GCC by default
- All programs link against `phylip.o` and the math library (`-lm`)
- No external library dependencies beyond libc and libm
- Drawing programs (Drawgram, Drawtree) also built as shared libraries (`libdrawgram.so`, `libdrawtree.so`) for Java GUI interfaces
- Java GUI interfaces for Drawgram and Drawtree are provided as JAR files
- Font files (`font1` through `font6`) are required at runtime for drawing programs

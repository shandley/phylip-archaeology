# PHYLIP Archaeology

**Preserving, mining, and modernizing the algorithmic legacy of phylogenetics**

---

> *"In 1980, before the World Wide Web, before GenBank went public, before most biologists
> had ever touched a command line, Joe Felsenstein released PHYLIP -- and quietly changed
> the course of evolutionary biology."*

## Mission

This project is an archaeological expedition into one of the most important archives in
the history of bioinformatics: **PHYLIP** (PHYLogeny Inference Package) and Joe
Felsenstein's comprehensive [catalog of 407 phylogenetics software packages](https://phylipweb.github.io/phylip/software.html).

Our goals:

1. **Preserve** the source code, algorithms, and software catalog before links rot and
   history is lost
2. **Mine** the codebase for algorithms and computational ideas that shaped modern
   phylogenetics
3. **Modernize** the most valuable algorithms with clean, safe Rust implementations
4. **Validate** the reimplementation against PHYLIP's original C executables
5. **Honor** the extraordinary contributions of Joe Felsenstein, who built the
   computational foundation of an entire scientific discipline

## Why This Matters

PHYLIP was first released in **1980**. It is one of the oldest and most influential
bioinformatics software packages ever created. The algorithms it implements --
Felsenstein's pruning algorithm for maximum likelihood, neighbor-joining, bootstrap
resampling, parsimony methods, and many more -- remain the mathematical backbone of
modern tools like RAxML, IQ-TREE, BEAST, and MrBayes.

Many of the ideas in PHYLIP anticipated developments in machine learning by decades.
The pruning algorithm (1981) is a special case of belief propagation, published seven
years before Pearl's formalization. The discrete gamma rate model is a mixture model
with log-sum-exp stabilization, standard in deep learning but pioneered here in 1994.
Site-pattern compression, closed-form matrix decomposition, weight-based bootstrapping
-- these techniques remain relevant far beyond phylogenetics.

## Key Findings

**20 cross-disciplinary connections** recovered from PHYLIP's source code, linking
phylogenetics to electrical engineering (independent contrasts = Kirchhoff's circuit
laws), information geometry (chord distance on hyperspheres), algebraic statistics
(Lake's invariants), coding theory (genetic code optimization), and probabilistic
graphical models (pruning = belief propagation). See [INSIGHTS.md](INSIGHTS.md).

**No bugs found in PHYLIP.** Direct comparison of 27 PHYLIP programs against our
reimplementation found zero mathematical errors in the original C code. Every
difference was attributable to deliberate design choices, different search heuristics,
or different model implementations. The one algorithmic discrepancy was in *our* code
(Dollo parsimony scoring), not PHYLIP's. A codebase written primarily by one researcher
over four decades proved mathematically correct across all algorithmic programs tested.

**23 tools permanently lost.** Analysis of Felsenstein's hand-curated catalog of 407
phylogenetics tools found that 23 have no archived copy anywhere and may be
irrecoverable.

## phylip-rs: The Rust Reimplementation

The heart of this project is `phylip-rs` -- a complete, modern Rust reimplementation
of PHYLIP's core algorithms. **Zero external dependencies.** Every mathematical
function -- the gamma function, matrix exponentiation, Newton-Raphson optimization,
continued fractions -- is implemented from first principles using only `std`.

### By the Numbers

| Metric | Value |
|--------|-------|
| Lines of Rust | **35,805** |
| Source files | **58** |
| Unit tests | **934** |
| Doc tests | **25** |
| Validation tests | **91** |
| Total tests | **1,050** |
| Compiler warnings | **0** |
| External dependencies | **0** |
| PHYLIP programs covered | **29/36** |
| PHYLIP programs compared | **27** |
| CLI commands | **9** |
| Interactive demonstrations | **10** |

### Algorithms Implemented

**Maximum Likelihood**
- Felsenstein's pruning algorithm (1981) -- the foundational ML phylogenetics algorithm
- ML tree search via stepwise addition + SPR rearrangement
- NNI (Nearest-Neighbor Interchange) tree refinement
- Branch length optimization via Newton-Raphson with numerical derivatives
- Transition/transversion ratio estimation (counting and ML methods)
- Discrete gamma site rate heterogeneity (Yang 1994)
- Model selection via AIC, BIC, and AICc with Akaike weights
- Optimized engine with site-pattern compression and transition matrix caching
- Clock-constrained ML (dnamlk) -- height-parameterized branch lengths

**Substitution Models & Distance Formulas**
- JC69 (Jukes-Cantor 1969) -- equal rates
- F84 (Felsenstein 1984) -- unequal frequencies, ts/tv distinction
- Poisson and WAG models for protein sequences
- 20-state pruning algorithm for amino acid data
- LogDet/Paralinear distance -- compositionally robust, 4x4 determinant from first principles
- Protein distances: Kimura (1983), Poisson, PAM/Dayhoff
- Restriction site distances: Nei-Li (1979)
- Gene frequency distances: Nei's genetic distance, Cavalli-Sforza chord, Reynolds

**Parsimony**
- Fitch algorithm (1971) with bitwise state set operations
- Wagner parsimony tree search via stepwise addition + SPR
- Ancestral state reconstruction (Fitch preorder pass)
- Dollo parsimony -- derived state arises once, losses free (dollop)
- Camin-Sokal irreversible parsimony -- 0->1 only (mix)
- Branch-and-bound exact search -- Hendy-Penny algorithm, guaranteed optimal (dnapenny)
- Multistate parsimony -- up to 32 states with Sankoff weighted step matrix (pars)
- Protein parsimony -- genetic code step matrix, 20-state Sankoff (protpars)
- Pluggable ParsimonyScorer trait for custom scoring criteria

**Distance Methods**
- Neighbor-Joining (Saitou & Nei 1987)
- UPGMA (Sokal & Michener 1958) -- ultrametric clustering
- Fitch-Margoliash weighted least squares
- Kitsch -- clock-constrained Fitch-Margoliash (ultrametric least squares)
- ML pairwise distances via Newton-Raphson optimization

**Tree Comparison**
- Robinson-Foulds distance (symmetric difference of bipartitions)
- Normalized Robinson-Foulds distance
- Branch Score Distance (Kuhner-Felsenstein) -- Euclidean distance in split space

**Compatibility & Invariants**
- Character compatibility analysis with pairwise compatibility testing
- Maximum clique finding via Bron-Kerbosch with pivoting (clique)
- DNA compatibility search -- maximize sites without homoplasy (dnacomp)
- Lake's phylogenetic invariants for 4-taxon problems (dnainvar)
- Cavender's invariants

**Comparative Methods**
- Felsenstein's independent contrasts (1985) -- standardized contrasts on trees (contrast)
- PIC correlation testing between continuous traits
- Brownian motion ML -- contrasts-based O(n) likelihood (contml)
- ML tree search for continuous character data

**Statistical Support**
- Bootstrap resampling (Felsenstein 1985)
- Block bootstrap for correlated sites
- Delete-fraction jackknife
- Bootstrap + ML integration (replicate ML searches with support values)
- Consensus trees: strict, majority-rule, extended majority-rule, threshold

**I/O and Interface**
- PHYLIP interleaved/sequential format parser
- FASTA parser (DNA and protein)
- Binary (0/1) character data parser
- Newick tree format reader/writer
- PHYLIP-style output reports
- Command-line interface with 9 analysis commands

### Quick Start

```bash
# Build the library and CLI
cd phylip-rs
cargo build --release

# Run maximum likelihood analysis
cargo run --release -- ml --input alignment.fasta --model f84

# Run parsimony analysis
cargo run --release -- parsimony --input alignment.fasta

# Run distance-based analysis (Neighbor-Joining)
cargo run --release -- distance --input alignment.fasta --method nj

# Run bootstrap analysis (100 replicates)
cargo run --release -- bootstrap --input alignment.fasta --replicates 100

# Run all 1,050 tests
cargo test
```

### Architecture

```
phylip-rs/src/
  tree/          Core types (Tree, Alignment, Base), Newick I/O,
                 splits, Robinson-Foulds & Branch Score distances
  models/        Substitution models (JC69, F84, protein/WAG),
                 LogDet, protein distances, restriction, gene freq
  likelihood/    Pruning algorithm, ML search, NNI, gamma rates,
                 ts/tv estimation, model selection, clock ML
  parsimony/     Wagner, Dollo, Camin-Sokal, branch-and-bound,
                 multistate (Sankoff), protein parsimony
  distance/      NJ, UPGMA, Fitch-Margoliash, Kitsch, ML distances
  bootstrap/     Resampling, consensus trees, ML bootstrap
  consensus/     Strict/majority-rule/extended consensus
  compatibility/ Clique analysis (Bron-Kerbosch), DNA compatibility
  comparative/   Independent contrasts, Brownian motion ML
  invariants/    Lake's & Cavender's phylogenetic invariants
  io/            PHYLIP format, FASTA, binary data, output reports
  main.rs        CLI binary with 9 analysis commands
```

## Validation

phylip-rs is validated through four complementary strategies and **1,050 tests**:

1. **Analytical tests** (32) -- mathematical formulas verified against hand-calculated values
2. **Classic dataset tests** (15) -- published results from Saitou & Nei, Felsenstein, Kimura, Lake, Cavender reproduced
3. **PHYLIP C comparison** (33) -- direct comparison against PHYLIP 3.697 C executables on identical inputs across 27 programs
4. **Medium-scale integration** (11) -- statistical properties verified on simulated datasets of 8-100 taxa

The PHYLIP comparison tests use hardcoded reference values (run without PHYLIP installed) with optional live execution when PHYLIP binaries are available:

```bash
# Run standard tests (no external dependencies)
cargo test -p phylip-rs

# Download PHYLIP 3.697 and run comparison tests
cd validation && bash setup.sh
PHYLIP_EXE_DIR=validation/phylip-3.697/exe cargo test -p phylip-rs --test validation_phylip -- --ignored
```

See [validation/VALIDATION_REPORT.md](validation/VALIDATION_REPORT.md) for the full report covering all 33 comparison tests with inputs, commands, reference values, tolerances, and known differences.

## Benchmarking

We benchmarked phylip-rs against IQ-TREE 3, RAxML-NG, and VeryFastTree on 36 simulated
datasets (10-500 taxa, 500-5,000 sites) under JC69. Key findings:

- **Small datasets (10-20 taxa)**: phylip-rs ML is competitive, sometimes finding the same optimum as modern tools
- **50+ taxa**: modern tools are 10-20x faster with better optima, reflecting four decades of search heuristic engineering
- **NJ transfers perfectly**: phylip-rs NJ completed all 36 datasets with competitive accuracy
- **The scoring function is identical**: all trees evaluated under JC69 matched IQ-TREE to 4 decimal places

See [benchmarks/](benchmarks/) for scripts, data, and results (180 benchmark runs).

## Interactive Demonstrations

Ten executable examples demonstrate algorithmic insights, each proving a specific claim
with concrete numerical results. Run any with `cargo run --release --example <name>`.

| Example | What It Demonstrates |
|---------|---------------------|
| `felsenstein_zone` | Parsimony converges on the wrong tree; ML gets it right |
| `language_evolution` | DNA pruning algorithm applied to human languages -- same code |
| `compositional_bias` | LogDet recovers true tree under GC-bias where JC69 fails |
| `kirchhoff_contrasts` | Independent contrasts = circuit theory (match to 8+ decimals) |
| `supplement_bound` | Branch-and-bound examines ~0.05% of tree space for 9 taxa |
| `dollo_gain_loss` | Dollo vs Fitch: same topology, different biological stories |
| `chord_geometry` | Chord distance is a literal chord on a hypersphere |
| `clock_constraints` | Kitsch outperforms UPGMA when the molecular clock is violated |
| `lake_invariants` | Polynomial invariants identify true topology at 500+ sites |
| `genetic_code_distances` | The genetic code is optimized for error tolerance (z = -2.76) |

## Software Catalog

Felsenstein maintained a hand-curated catalog of phylogenetics software -- effectively
a package manager for an entire scientific discipline. We analyzed all 407 tools,
documenting their preservation status:

- **196** (48%) archived but unmaintained
- **137** (34%) dormant (websites accessible, no recent development)
- **34** (8%) completely unreachable
- **23** tools have no Wayback Machine archive and may be permanently lost

Browse the catalog: **[Interactive Explorer](https://shandley.github.io/phylip-archaeology/)**

## Project Structure

```
phylip-archaeology/
├── README.md              # This file
├── LICENSE                # MIT License
├── INSIGHTS.md            # 20 cross-disciplinary algorithmic connections
├── TRIBUTE.md             # Historical narrative of Felsenstein's contributions
├── REFLECTION.md          # What we built and what we learned
├── Cargo.toml             # Workspace root (phylip-rs + phylip-wasm)
├── phylip-rs/             # Rust reimplementation (35,805 lines, 1,050 tests)
│   ├── src/               # Library and CLI source (58 files)
│   └── examples/          # 10 interactive demonstrations
├── manuscript/            # Nature Methods article draft
│   ├── manuscript.md      # Full manuscript (~6,000 words)
│   ├── supplementary.md   # Supplementary materials
│   └── figures/           # Figures 1-2 with generation scripts
├── validation/            # PHYLIP 3.697 comparison infrastructure
│   ├── VALIDATION_REPORT.md  # Detailed report (33 tests, 27 programs)
│   ├── README.md          # Quick start guide
│   └── setup.sh           # Downloads and compiles PHYLIP 3.697
├── benchmarks/            # Performance benchmarking pipeline
│   ├── generate_data.py   # JC69 dataset simulator
│   ├── run_benchmarks.py  # Multi-tool benchmark runner
│   ├── results/           # benchmark_results.csv (180 runs)
│   └── figures/           # Figure 4 (PDF/PNG)
├── catalog/               # Software catalog preservation (407 tools)
│   └── analysis/          # Scraping, enrichment, figures
├── docs/                  # Interactive WASM demo (GitHub Pages)
├── phylip-wasm/           # WebAssembly bindings for browser demo
├── phylip-source/         # PHYLIP C source code archive and analysis
├── algorithms/            # Extracted algorithm documentation
└── timeline/              # Historical data and visualizations
```

## Key Principles

- **Fidelity first**: Preserve original algorithms exactly before modernizing
- **Zero dependencies**: The code is its own textbook -- every function from first principles
- **Validation**: 1,050 tests including direct comparison against 27 PHYLIP 3.697 C programs
- **Attribution**: Every algorithm traces back to its originator and key papers
- **Accessibility**: Clear documentation for both historians and practitioners
- **Respect**: This is archaeology, not criticism -- honor the constraints of the era

## References

- Felsenstein, J. (1978). Cases in which parsimony or compatibility methods will be
  positively misleading. *Systematic Zoology*, 27, 401-410.
- Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood
  approach. *Journal of Molecular Evolution*, 17, 368-376.
- Felsenstein, J. (1985). Confidence limits on phylogenies: an approach using the
  bootstrap. *Evolution*, 39, 783-791.
- Felsenstein, J. (1989). PHYLIP - Phylogeny Inference Package (Version 3.2).
  *Cladistics*, 5, 164-166.
- Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
- Yang, Z. (1994). Maximum likelihood phylogenetic estimation from DNA sequences with
  variable rates over sites. *Journal of Molecular Evolution*, 39, 306-314.
- Saitou, N. & Nei, M. (1987). The neighbor-joining method. *Molecular Biology and
  Evolution*, 4, 406-425.
- PHYLIP home page: https://phylipweb.github.io/phylip/
- PHYLIP source: https://github.com/phylipweb/phylip

## License

This project is released under the [MIT License](LICENSE).

The original PHYLIP source code has its own open-source license (since v3.696).
See the PHYLIP repository for details.

# PHYLIP Archaeology

**Preserving, mining, and modernizing the algorithmic legacy of phylogenetics**

---

> *"In 1980, before the World Wide Web, before GenBank went public, before most biologists
> had ever touched a command line, Joe Felsenstein released PHYLIP -- and quietly changed
> the course of evolutionary biology."*

## Mission

This project is an archaeological expedition into one of the most important archives in
the history of bioinformatics: **PHYLIP** (PHYLogeny Inference Package) and Joe
Felsenstein's comprehensive [catalog of 392+ phylogenetics software packages](https://phylipweb.github.io/phylip/software.html).

Our goals:

1. **Preserve** the source code, algorithms, and software catalog before links rot and
   history is lost
2. **Mine** the codebase for algorithms and computational ideas that shaped modern
   phylogenetics
3. **Modernize** the most valuable algorithms with clean, safe Rust implementations
4. **Honor** the extraordinary contributions of Joe Felsenstein, who built the
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

See [INSIGHTS.md](INSIGHTS.md) for a deep analysis of what Felsenstein understood that
the field is at risk of losing.

See [TRIBUTE.md](TRIBUTE.md) for a full historical narrative of Felsenstein's
contributions.

## phylip-rs: The Rust Reimplementation

The heart of this project is `phylip-rs` -- a complete, modern Rust reimplementation
of PHYLIP's core algorithms. **Zero external dependencies.** Every mathematical
function -- the gamma function, matrix exponentiation, Newton-Raphson optimization,
continued fractions -- is implemented from first principles using only `std`.

### By the Numbers

| Metric | Value |
|--------|-------|
| Lines of Rust | **20,749** |
| Source files | **35** |
| Unit tests | **561** |
| Doc tests | **17** |
| Compiler warnings | **0** |
| External dependencies | **0** |

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

**Substitution Models**
- JC69 (Jukes-Cantor 1969) -- equal rates
- F84 (Felsenstein 1984) -- unequal frequencies, ts/tv distinction
- Poisson and WAG models for protein sequences
- 20-state pruning algorithm for amino acid data

**Parsimony**
- Fitch algorithm (1971) with bitwise state set operations
- Wagner parsimony tree search via stepwise addition + SPR
- Ancestral state reconstruction (Fitch preorder pass)

**Distance Methods**
- Neighbor-Joining (Saitou & Nei 1987)
- Fitch-Margoliash weighted least squares
- ML pairwise distances via Newton-Raphson optimization

**Statistical Support**
- Bootstrap resampling (Felsenstein 1985)
- Block bootstrap for correlated sites
- Delete-fraction jackknife
- Bootstrap + ML integration (replicate ML searches with support values)
- Consensus trees: strict, majority-rule, extended majority-rule, threshold

**I/O and Interface**
- PHYLIP interleaved/sequential format parser
- FASTA parser (DNA and protein)
- Newick tree format reader/writer
- PHYLIP-style output reports
- Command-line interface with 5 analysis commands

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
```

### Architecture

```
phylip-rs/src/
  tree/          Core types (Tree, Alignment, Base), Newick I/O
  models/        Substitution models (JC69, F84, protein/WAG)
  likelihood/    Pruning algorithm, ML search, NNI, gamma rates,
                 ts/tv estimation, model selection, optimized engine
  parsimony/     Wagner parsimony tree search
  distance/      NJ, Fitch-Margoliash, ML pairwise distances
  bootstrap/     Resampling, consensus trees, ML bootstrap
  consensus/     Strict/majority-rule/extended consensus
  io/            PHYLIP format parser, FASTA, output reports
  main.rs        CLI binary with 5 analysis commands
```

## Interactive Demonstrations

Two executable examples demonstrate insights from the algorithms in contexts beyond
traditional phylogenetics.

### The Felsenstein Zone: When More Data Makes You More Wrong

```bash
cargo run --release --example felsenstein_zone
```

Simulates Felsenstein's famous 1978 result: maximum parsimony converges on the
**wrong** tree with increasing confidence as you add data, while maximum likelihood
correctly recovers the truth. The simulation generates DNA sequences along a tree
with two long branches (long branch attraction), then evaluates all three possible
4-taxon topologies under both methods.

```
  Sites | Parsimony correct |    ML correct | Parsimony picks T3 (wrong)
--------|-------------------|---------------|---------------------------
    100 |                8% |           88% |                        92%
    500 |                2% |          100% |                        98%
   1000 |                0% |          100% |                       100%
   5000 |                0% |          100% |                       100%
  10000 |                0% |          100% |                       100%
```

### Language Evolution: DNA Code Analyzes Human Languages

```bash
cargo run --example language_evolution
```

Applies the **exact same pruning algorithm** -- designed for DNA -- to linguistic
data: cognate class assignments for 37 vocabulary items across English, German,
French, Italian, Spanish, and Portuguese. Not a single line of code changes. The
algorithm correctly identifies the known language family tree and reconstructs
"proto-language" states at the root. Demonstrates that Felsenstein's algorithm is a
general-purpose inference engine for discrete states on trees, not just a "DNA
algorithm."

## Project Structure

```
phylip-archaeology/
├── INSIGHTS.md            # Deep analysis of PHYLIP's algorithmic insights
├── TRIBUTE.md             # Historical narrative of Felsenstein's contributions
├── catalog/               # Software catalog preservation (392+ tools)
├── phylip-source/         # PHYLIP C source code archive and analysis
├── algorithms/            # Extracted algorithm documentation
├── phylip-rs/             # Modern Rust reimplementation (20,749 lines)
│   ├── src/               # Library and CLI source
│   └── examples/          # Interactive demonstrations
└── timeline/              # Historical data and visualizations
```

## Key Principles

- **Fidelity first**: Preserve original algorithms exactly before modernizing
- **Zero dependencies**: The code is its own textbook -- every function from first principles
- **Validation**: 578 tests verify correctness against known analytical results
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

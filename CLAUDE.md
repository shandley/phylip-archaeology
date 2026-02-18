# CLAUDE.md — Project Context for Claude Code

## Project Overview

**PHYLIP Archaeology** — Preserving, mining, and modernizing the algorithmic legacy of Joe Felsenstein's PHYLIP (PHYLogeny Inference Package, first released 1980). The heart of the project is `phylip-rs`, a modern Rust reimplementation of PHYLIP's core algorithms with **zero external dependencies**.

## Repository Structure

```
phylip-archaeology/
├── CLAUDE.md              # This file — project context
├── README.md              # Project overview and documentation
├── INSIGHTS.md            # Deep analysis of PHYLIP's algorithmic insights
├── TRIBUTE.md             # Historical narrative of Felsenstein's contributions
├── REFLECTION.md          # What we built and what we learned
├── Cargo.toml             # Workspace root
├── catalog/               # Software catalog preservation (392+ tools)
├── phylip-source/         # PHYLIP C source code archive and analysis
├── algorithms/            # Extracted algorithm documentation
├── phylip-rs/             # Modern Rust reimplementation
│   ├── Cargo.toml
│   ├── src/               # Library and CLI source (58 files)
│   └── examples/          # Interactive demonstrations (2)
└── timeline/              # Historical data and visualizations
```

## Current Stats

| Metric | Value |
|--------|-------|
| Lines of Rust | 35,805 |
| Source files | 58 |
| Unit tests | 934 |
| Doc tests | 25 |
| Total tests | 959 |
| Compiler warnings | 0 |
| External dependencies | 0 |
| PHYLIP programs covered | ~30/36 |
| CLI commands | 9 |

## phylip-rs Module Map

```
phylip-rs/src/
  tree/           Core types (Tree, Alignment, Base), Newick I/O,
                  splits (bitvectors), Robinson-Foulds & Branch Score distances
  models/         Substitution models (JC69, K2P, F81, F84, protein/WAG),
                  LogDet, protein distances, restriction site, gene frequency
  likelihood/     Pruning algorithm, ML search, NNI, gamma rates,
                  ts/tv estimation, model selection, clock-constrained ML
  parsimony/      Wagner, Dollo, Camin-Sokal, branch-and-bound,
                  multistate (Sankoff), protein parsimony, ParsimonyScorer trait
  distance/       NJ, UPGMA, Fitch-Margoliash, Kitsch, ML pairwise distances
  bootstrap/      Resampling, consensus trees, ML bootstrap
  consensus/      Strict/majority-rule/extended consensus
  compatibility/  Clique analysis (Bron-Kerbosch), DNA compatibility
  comparative/    Independent contrasts, Brownian motion ML (contml)
  invariants/     Lake's & Cavender's phylogenetic invariants
  io/             PHYLIP format, FASTA, binary data, output reports
  main.rs         CLI binary with 9 analysis commands
  lib.rs          Library root
```

## Build & Test

```bash
cd phylip-rs
cargo build                    # Build library and CLI
cargo test                     # Run all 959 tests
cargo build --examples         # Build interactive demonstrations
cargo run --release -- --help  # CLI usage
```

## Key Design Principles

- **Zero dependencies** — every mathematical function from first principles (gamma function, matrix exponentiation, Newton-Raphson, continued fractions, etc.)
- **Fidelity** — algorithms match original PHYLIP behavior, validated against known analytical results
- **Trait-based extensibility** — SubstitutionModel trait, ParsimonyScorer trait for pluggable algorithms
- **Test-grounded** — tests compare against hand-calculated values and published results, not just regression

## Git State

- **Branch:** main
- **Unpushed:** 1 commit (algorithm expansion from ~12 to ~30 PHYLIP programs, +15,299 lines)
- **Remote:** origin (GitHub: shandley/phylip-archaeology)

## 6 Programs NOT Covered (interactive/drawing tools, not algorithmic)

drawgram, drawtree, dnamove, dolmove, move, retree, factor

## Previously Discussed Next Steps (Not Yet Started)

These were proposed before the algorithm expansion and should be re-evaluated:

1. **Academic paper** — "PHYLIP Archaeology: Rediscovering the Algorithmic Foundations of Phylogenetics" targeting Molecular Ecology or Bioinformatics Application Note
2. **WASM demo** — Browser-based interactive Felsenstein Zone visualization
3. **Performance benchmarking** — Compare phylip-rs against IQ-TREE, RAxML-NG, PAUP* on standard datasets
4. **Software catalog preservation** — Scrape and archive the 392+ entries from Felsenstein's software catalog before links rot
5. **Cross-domain demonstrations** — More examples applying PHYLIP algorithms outside phylogenetics (tumor phylogenetics, cultural evolution)
6. **Community engagement** — Share with Joe Felsenstein, phylogenetics community

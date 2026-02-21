# CLAUDE.md — Project Context for Claude Code

## Project Overview

**PHYLIP Archaeology** — Preserving, mining, and modernizing the algorithmic legacy of Joe Felsenstein's PHYLIP (PHYLogeny Inference Package, first released 1980). The heart of the project is `phylip-rs`, a modern Rust reimplementation of PHYLIP's core algorithms with **zero external dependencies**.

## Repository Structure

```
phylip-archaeology/
├── CLAUDE.md              # This file — project context
├── README.md              # Project overview and documentation
├── INSIGHTS.md            # Deep analysis of 20 algorithmic case studies
├── TRIBUTE.md             # Historical narrative of Felsenstein's contributions
├── REFLECTION.md          # What we built and what we learned
├── PAPER_PLAN.md           # Nature Methods paper planning document
├── Cargo.toml             # Workspace root
├── phylip-rs/             # Modern Rust reimplementation
│   ├── Cargo.toml
│   ├── src/               # Library and CLI source (57 files)
│   └── examples/          # 10 interactive demonstrations
├── manuscript/            # Nature Methods Article draft
│   ├── manuscript.md      # Full manuscript (~3,000 words + Methods)
│   └── supplementary.md   # Supplementary material outline
├── benchmarks/            # Performance benchmarking pipeline
│   ├── generate_data.py   # JC69 dataset simulator
│   ├── run_benchmarks.py  # Multi-tool benchmark runner
│   ├── plot_figure5.py    # Figure generation script
│   ├── results/           # benchmark_results.csv (180 runs)
│   └── figures/           # Figure 5 (PDF/PNG)
├── catalog/               # Software catalog preservation (407 tools)
│   └── analysis/          # Scraping, enrichment, Figure 4
├── validation/            # Validation infrastructure (PHYLIP 3.697 comparison)
│   ├── VALIDATION_REPORT.md  # Living validation report (Supplementary Note 5)
│   ├── README.md          # Quick start and test overview
│   └── setup.sh           # Downloads and compiles PHYLIP 3.697
├── phylip-source/         # PHYLIP C source code archive and analysis
├── algorithms/            # Extracted algorithm documentation
└── timeline/              # Historical data and visualizations
```

## Current Stats

| Metric | Value |
|--------|-------|
| Lines of Rust | 35,805 |
| Source files | 58 |
| Unit tests | 934 |
| Doc tests | 25 |
| Validation tests | 91 |
| Total tests | 1,050 |
| Compiler warnings | 0 |
| External dependencies | 0 |
| PHYLIP programs covered | ~30/36 |
| CLI commands | 9 |
| Interactive demonstrations | 10 |
| Algorithmic case studies | 20 |
| Software catalog tools | 407 |
| Benchmark datasets | 36 |

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
cargo test                     # Run all 1,050 tests
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
- **Remote:** origin (GitHub: shandley/phylip-archaeology)

## 6 Programs NOT Covered (interactive/drawing tools, not algorithmic)

drawgram, drawtree, dnamove, dolmove, move, retree, factor

## Project Status

### Completed
- **Rust reimplementation** — 35,805 lines, 1,050 tests, ~30/36 PHYLIP programs, zero dependencies
- **Validation suite** — 91 validation tests across 4 strategies (analytical, classic datasets, PHYLIP 3.697 comparison, medium-scale integration)
- **20 algorithmic case studies** — Cross-disciplinary insights documented in INSIGHTS.md
- **10 interactive demonstrations** — Compilable examples in phylip-rs/examples/
- **Software catalog analysis** — 407 tools scraped, enriched, analyzed; Figure 4 generated
- **Performance benchmarking** — 180 runs (36 datasets x 5 tools); Figure 5 generated
- **Manuscript draft** — Nature Methods Article in manuscript/manuscript.md (~3,000 words + Methods)

### Remaining
- **Figures 1-2** — Workflow schematic and algorithmic discovery panels (require manual/programmatic creation)
- **Manuscript revision** — Polish, expand to full word budget, finalize references
- **Community engagement** — Share with Joe Felsenstein, phylogenetics community
- **Submission** — Final polish, cover letter, Nature Methods submission

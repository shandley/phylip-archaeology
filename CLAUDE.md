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
├── LICENSE                # MIT License
├── Cargo.toml             # Workspace root (phylip-rs + phylip-wasm)
├── .github/workflows/     # CI: test (ubuntu + macos), format check
├── phylip-rs/             # Modern Rust reimplementation
│   ├── Cargo.toml
│   ├── src/               # Library and CLI source (58 files)
│   ├── examples/          # 10 interactive demonstrations
│   └── tests/             # 4 validation test files (91 tests)
├── manuscript/            # Nature Methods Article draft
│   ├── manuscript.md      # Full manuscript (~6,000 words)
│   ├── supplementary.md   # Supplementary materials
│   └── figures/           # Figures 1-2 with generation scripts
├── validation/            # PHYLIP 3.697 comparison infrastructure
│   ├── VALIDATION_REPORT.md  # Detailed report (33 tests, 27 programs)
│   ├── README.md          # Quick start and test overview
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
└── timeline/              # Historical data (events.json)
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
| PHYLIP programs covered | 29/36 |
| PHYLIP programs compared | 27 |
| CLI commands | 9 |
| Interactive demonstrations | 10 |
| Algorithmic case studies | 20 |
| Software catalog tools | 407 |
| Benchmark datasets | 36 |

## Validation Test Breakdown

| File | Category | Tests |
|------|----------|-------|
| validation_analytical.rs | Analytical formulas | 32 |
| validation_classics.rs | Classic published datasets | 15 |
| validation_phylip.rs | PHYLIP 3.697 C comparison | 33 |
| validation_medium.rs | Medium-scale integration | 11 |

Run validation tests:
```bash
cargo test -p phylip-rs --test validation_analytical --test validation_classics --test validation_medium --test validation_phylip
```

Run PHYLIP live comparison (requires PHYLIP 3.697 binaries):
```bash
cd validation && bash setup.sh
PHYLIP_EXE_DIR=validation/phylip-3.697/exe cargo test -p phylip-rs --test validation_phylip -- --ignored
```

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
cargo clippy -- -D warnings    # Lint check
```

## Key Design Principles

- **Zero dependencies** — every mathematical function from first principles (gamma function, matrix exponentiation, Newton-Raphson, continued fractions, etc.)
- **Fidelity** — algorithms match original PHYLIP behavior, validated against PHYLIP 3.697 C executables
- **Trait-based extensibility** — SubstitutionModel trait, ParsimonyScorer trait for pluggable algorithms
- **Test-grounded** — tests compare against hand-calculated values, published results, and PHYLIP C output

## Key Finding

**No bugs found in PHYLIP.** Direct comparison of 27 programs against our reimplementation found zero mathematical errors in the original C code. Every difference was attributable to deliberate design choices, different search heuristics, or different model implementations. The one algorithmic discrepancy was in *our* code (Dollo parsimony scoring), not PHYLIP's.

## Known Limitations

- **Dollo parsimony (DolloScorer)** — Uses upward-only (postorder) pass; PHYLIP uses two-pass (upward + downward correction placing gain at MRCA). Our scorer can overcount losses. Documented in `parsimony/dollo.rs`.
- **7 programs not covered** — drawgram, drawtree, dnamove, dolmove, move, retree, factor (interactive/drawing tools, not algorithmic)

## Git State

- **Branch:** main
- **Remote:** origin (GitHub: shandley/phylip-archaeology)
- **Release:** v0.1.0
- **CI:** GitHub Actions (test on ubuntu + macos, format check)

## Project Status

### Completed
- **Rust reimplementation** — 35,805 lines, 1,050 tests, 29/36 PHYLIP programs, zero dependencies
- **Validation suite** — 91 validation tests across 4 strategies (analytical, classic datasets, PHYLIP 3.697 comparison, medium-scale integration)
- **PHYLIP 3.697 comparison** — 33 tests across 27 programs, zero bugs found in original C code
- **20 algorithmic case studies** — Cross-disciplinary insights documented in INSIGHTS.md
- **10 interactive demonstrations** — Compilable examples in phylip-rs/examples/
- **Software catalog analysis** — 407 tools scraped, enriched, analyzed; interactive explorer at GitHub Pages
- **Performance benchmarking** — 180 runs (36 datasets x 5 tools)
- **Figures 1-4** — Workflow schematic, algorithmic discovery, catalog analysis, benchmarks
- **Manuscript** — Nature Methods Article (~6,000 words) with supplementary materials
- **CI/CD** — GitHub Actions: test (ubuntu + macos), format check
- **Release** — v0.1.0 tagged and published on GitHub
- **Repository cleanup** — Dev artifacts, empty scaffolding, lock files removed

### Remaining
- **Manuscript revision** — Final copy-editing, cover letter
- **Community engagement** — Share with Joe Felsenstein, phylogenetics community
- **Submission** — Nature Methods submission
- **Optional: Fix Dollo scoring** — Implement two-pass algorithm to match PHYLIP's behavior

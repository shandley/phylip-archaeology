# Algorithm Extraction

This directory documents the core algorithms found in PHYLIP and the broader
phylogenetics software ecosystem. Each algorithm entry includes:

- Mathematical foundation and derivation
- Pseudocode (language-agnostic)
- Original C implementation (extracted from PHYLIP source, annotated)
- Key references and papers
- Historical context
- Relationship to other algorithms

## Priority Algorithms

These algorithms will be documented and reimplemented first, based on their
foundational importance to the field:

| # | Algorithm | Origin | PHYLIP Program |
|---|-----------|--------|----------------|
| 1 | Felsenstein pruning | Felsenstein 1981 | `dnaml`, `proml` |
| 2 | Neighbor-joining | Saitou & Nei 1987 | `neighbor` |
| 3 | UPGMA | Sokal & Michener 1958 | `neighbor` |
| 4 | Fitch-Margoliash | Fitch & Margoliash 1967 | `fitch` |
| 5 | Wagner parsimony | Kluge & Farris 1969 | `dnapars`, `pars` |
| 6 | Bootstrap resampling | Felsenstein 1985 | `seqboot` |
| 7 | Consensus trees | Adams 1972; Margush & McMorris 1981 | `consense` |
| 8 | JC69 model | Jukes & Cantor 1969 | `dnadist` |
| 9 | K2P model | Kimura 1980 | `dnadist` |
| 10 | F81 model | Felsenstein 1981 | `dnaml` |
| 11 | F84 model | Felsenstein 1984 | `dnaml`, `dnadist` |

## Directory Structure

Each algorithm has its own directory under `entries/`:

```
entries/
├── felsenstein-pruning/
│   ├── README.md           # Full documentation
│   ├── original_source.c   # Extracted from PHYLIP
│   └── references.bib      # Key papers
├── neighbor-joining/
├── upgma/
└── ...
```

## Index

The `index.json` file provides a machine-readable index of all documented
algorithms, their relationships, and links to implementations.

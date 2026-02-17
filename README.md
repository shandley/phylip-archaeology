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

Felsenstein's software catalog page was essentially a **hand-curated package manager for
an entire scientific discipline**, maintained years before GitHub, Bioconductor, or
conda-forge existed. It cataloged 392+ phylogenetics tools and likely guided thousands
of researchers to the software they needed to complete their work.

Many of those tools have disappeared from the internet. University FTP servers have been
decommissioned. Personal web pages have gone dark. This project aims to preserve what
remains and extract the algorithmic wisdom embedded in these programs.

See [TRIBUTE.md](TRIBUTE.md) for a full historical narrative of Felsenstein's
contributions.

## Project Structure

```
phylip-archaeology/
├── catalog/               # Software catalog preservation (392+ tools as structured data)
├── phylip-source/         # PHYLIP C source code archive and analysis
├── algorithms/            # Extracted algorithm documentation (math, pseudocode, history)
├── phylip-rs/             # Modern Rust reimplementations of core algorithms
├── timeline/              # Structured historical data and visualizations
├── notebooks/             # Exploratory analysis notebooks
└── .github/workflows/     # CI for Rust crate and periodic link checking
```

### catalog/

Structured preservation of Felsenstein's software catalog. Every tool indexed with
metadata: name, author, status (alive/dead/archived), categories, algorithms
implemented, languages, platforms, citations, and influence lineage.

### phylip-source/

Archaeological analysis of PHYLIP's C source code. Program inventory, algorithm
mapping, and annotated analysis of key files like `dnaml.c` (maximum likelihood),
`neighbor.c` (neighbor-joining), `dnapars.c` (parsimony), and `seqboot.c` (bootstrap).

### algorithms/

Formalized documentation of each core algorithm: mathematical foundations, pseudocode,
original C implementation (extracted and annotated), key references, and historical
context. This is the bridge between the original code and the modern reimplementations.

### phylip-rs/

A Rust crate (`phylip-rs`) providing modern, safe, well-tested implementations of the
core phylogenetic algorithms. Every implementation is validated against original PHYLIP
output to ensure fidelity.

### timeline/

Structured data tracking the history of PHYLIP, Felsenstein's contributions, and the
evolution of the phylogenetics software ecosystem. Includes visualization scripts.

## Roadmap

| Phase | Goal | Status |
|-------|------|--------|
| 1. Foundation | Repo structure, README, tribute, timeline | In Progress |
| 2. Catalog Preservation | Scrape, structure, and assess all 392+ tools | Planned |
| 3. Source Archaeology | Deep analysis of PHYLIP's C source code | Planned |
| 4. Algorithm Extraction | Formalize core algorithms with math and pseudocode | Planned |
| 5. Rust Reimplementation | Modern `phylip-rs` crate with validation tests | Planned |
| 6. Analysis & Visualization | Ecosystem graphs, timelines, benchmarks | Planned |
| 7. Publication & Outreach | Polish, publish crate, community engagement | Planned |

## Key Principles

- **Fidelity first**: Preserve original algorithms exactly before modernizing
- **Validation**: Every Rust implementation must reproduce original PHYLIP output
- **Attribution**: Every algorithm traces back to its originator and key papers
- **Accessibility**: Clear documentation for both historians and practitioners
- **Respect**: This is archaeology, not criticism -- honor the constraints of the era

## References

- Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood
  approach. *Journal of Molecular Evolution*, 17, 368-376.
- Felsenstein, J. (1985). Confidence limits on phylogenies: an approach using the
  bootstrap. *Evolution*, 39, 783-791.
- Felsenstein, J. (1989). PHYLIP - Phylogeny Inference Package (Version 3.2).
  *Cladistics*, 5, 164-166.
- Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
- PHYLIP home page: https://phylipweb.github.io/phylip/
- PHYLIP source: https://github.com/phylipweb/phylip

## License

This project is released under the [MIT License](LICENSE).

The original PHYLIP source code has its own open-source license (since v3.696).
See the PHYLIP repository for details.

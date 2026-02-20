# Paper Plan: Nature Methods

**Working Title:** "LLM-Assisted Archaeological Bioinformatics: Recovering the Algorithmic Legacy of PHYLIP"

**Target Journal:** Nature Methods (Article — not Analysis, which is for benchmarking existing tools)

**Status:** Manuscript drafted (manuscript/manuscript.md)

**Last Updated:** 2026-02-20

---

## 1. Central Thesis

Large language models can serve as archaeological instruments for legacy scientific software — reading, understanding, and resurrecting algorithms buried in aging codebases. We demonstrate this by applying LLM-assisted code archaeology to PHYLIP, one of the most influential bioinformatics packages ever created (1980–present), recovering algorithmic insights that were never formally published, validating them through modern reimplementation, and preserving a decaying software catalog of 392+ phylogenetics tools.

## 2. Why Nature Methods

Nature Methods publishes Articles introducing new computational methodologies for biological research. (Note: "Analysis" articles are specifically for benchmarking existing tools, which is not our primary contribution.) This paper introduces **LLM-assisted code archaeology** as a methodology, demonstrated on phylogenetics but applicable to any field with legacy scientific software. The contributions are:

- A new methodology (LLM as archaeological tool for scientific code)
- Concrete algorithmic discoveries recovered from legacy code
- A validated, zero-dependency reimplementation as proof of understanding
- A preservation framework for endangered scientific software catalogs

The audience extends beyond phylogenetics to anyone maintaining or inheriting legacy computational infrastructure in the sciences.

## 3. Core Contributions (Four Pillars)

### Pillar 1: LLM as Archaeological Instrument
- **Novelty:** No prior systematic use of LLMs to excavate, analyze, and reimplement legacy scientific software
- **Methodology:** Describe the workflow — code reading, algorithm extraction, cross-referencing with literature, reimplementation, validation
- **Generalizability:** Applicable to BLAST, HMMER, PAML, climate models (CESM, WRF), particle physics (ROOT, Geant4), structural biology (CNS, X-PLOR)
- **Honest assessment:** What did the LLM get right? What required human domain expertise to guide? Where did it struggle?
- **Process documentation:** How many iterations, what kinds of errors, how the human-AI collaboration worked in practice

### Pillar 2: Algorithmic Resurrection
- **20 case studies** documenting algorithms recovered from PHYLIP source code
- **Key discoveries** (insights never published in Felsenstein's papers or textbook):
  - Independent contrasts = Kirchhoff's circuit laws (contrasts variance = effective resistance)
  - Cavalli-Sforza chord distance as literal chord on a hypersphere (information geometry connection)
  - Hendy-Penny supplement bound as dual decomposition (placed vs. unplaced taxa)
  - Dollo parsimony as minimum edge-cut / max-flow problem
  - LogDet as determinant-based marginal frequency correction
  - The genetic code as an error-correcting code (z-score = -2.76 vs random codes)
  - Kitsch node-height optimization as isotonic regression on tree partial orders
  - Bron-Kerbosch clique finding for character compatibility (graph theory ↔ phylogenetics)
  - Lake's invariants as algebraic statistics (polynomial varieties)
  - Felsenstein-Churchill HMM for autocorrelated rates
- **Cross-disciplinary connections:** Each algorithm connects to a different field (electrical engineering, information geometry, combinatorial optimization, algebraic geometry, coding theory)

### Pillar 3: Code Preservation and Cataloging
- **Felsenstein's software catalog:** 407 phylogenetics tools curated since the 1980s
- **Preservation status:** 196 archived (48%), 137 dormant (34%), 34 dead (8%), 40 unknown (10%); 23 permanently lost
- **Archival approach:** Systematic scraping, Wayback Machine cross-referencing, metadata extraction
- **The broader problem:** Scientific software as cultural heritage — what happens when maintainers retire?
- **STATUS: COMPLETE** — Figure 4 generated (catalog/analysis/figures/)

### Pillar 4: Historical Perspective
- Felsenstein's pruning algorithm (1981) anticipated belief propagation (Pearl, 1988) by 7 years
- Discrete gamma rates (Yang, 1994, building on Felsenstein's framework) anticipated mixture models used in deep learning
- Site-pattern compression anticipated the "sufficient statistics" approach in modern ML
- Bootstrap resampling for phylogenetics (1985) was one of the earliest applications in biology
- This is not hagiography — it's correcting the historical record of computational ideas

## 4. Paper Structure (Nature Methods Analysis format)

### Abstract (150 words)
- Problem: Legacy scientific software contains algorithmic knowledge at risk of being lost
- Approach: LLM-assisted code archaeology applied to PHYLIP
- Results: 20 algorithmic insights recovered, validated Rust reimplementation (35,805 lines, 959 tests), connections to 10+ fields outside phylogenetics
- Impact: New methodology for computational science preservation

### Introduction (~500 words)
- The preservation crisis: scientific software ages, maintainers retire, algorithms are lost
- PHYLIP as case study: 45+ years old, most-cited phylogenetics software, algorithms underpin all modern tools
- The opportunity: LLMs can read and understand legacy code in ways that scale beyond human effort
- What we did and what we found

### Results (~2000 words)

**The LLM archaeology workflow**
- How we used Claude to read, analyze, and reimplement PHYLIP's C codebase
- The iterative process: code reading → algorithm extraction → cross-referencing → reimplementation → validation
- Quantitative assessment of the process (accuracy, iteration count, types of errors)

**Algorithmic discoveries**
- Select 4-6 of the strongest case studies for the main text
- Emphasis on insights NOT found in any published paper (only discoverable from code)
- Concrete numerical demonstrations proving each claim
- Candidates for main text:
  1. Contrasts = Kirchhoff (beautiful, easy to explain, provably identical)
  2. Genetic code step matrix optimization (z = -2.76, visually striking)
  3. LogDet compositional robustness (practical importance, clean demonstration)
  4. Branch-and-bound supplement bound (elegant pruning, combinatorial insight)
  5. Felsenstein Zone (the most famous result, good for accessibility)

**Validated reimplementation**
- 35,805 lines of Rust, 959 tests, zero dependencies
- ~30/36 PHYLIP programs covered
- Performance characteristics vs. original C code
- What validation against known analytical results looks like

**Software catalog preservation**
- State of the 392+ entries: how many alive, how many dead, how many partially archived
- What was recovered, what was lost
- Metadata analysis: publication dates, citation patterns, technology evolution

**Benchmarking against modern tools**
- phylip-rs vs IQ-TREE 3, RAxML-NG, VeryFastTree on 36 simulated JC69 datasets (180 runs)
- Small datasets (10-20 taxa): phylip-rs finds identical optima as modern tools
- Large (100+ taxa): phylip-rs ML times out — quantifies the 40-year engineering gap
- phylip-rs uses 18x less memory than IQ-TREE (3.9 MB vs 54.5 MB)
- Not claiming performance parity — measuring what 40 years of heuristic engineering buys
- **STATUS: COMPLETE** — Figure 5 generated (benchmarks/figures/)

### Discussion (~500 words)
- LLM-assisted code archaeology as a general methodology
- What this means for scientific software preservation
- Limitations: what the LLM missed, what required human expertise
- Call to action: other legacy codebases that need this treatment
- The responsibility to preserve computational heritage

### Methods (~1000 words, can be longer)
- Detailed LLM workflow (model, prompts, iteration strategy)
- Reimplementation architecture
- Validation methodology
- Benchmarking setup
- Catalog scraping and analysis

### Figures (6 max for Analysis)

1. **Overview figure** — The LLM archaeology workflow (code → analysis → insight → reimplementation → validation). Schematic.

2. **Algorithmic discovery panel** — 2-3 side-by-side demonstrations:
   - (a) Kirchhoff contrasts: tree with circuit diagram overlay, numerical equivalence
   - (b) Genetic code heatmap: 20x20 step matrix, z-score distribution vs random codes
   - (c) Felsenstein Zone: parsimony vs ML accuracy curves crossing

3. **Cross-disciplinary connections** — Network/Sankey diagram showing PHYLIP algorithms connecting to other fields (electrical engineering, information geometry, algebraic geometry, coding theory, combinatorial optimization, etc.)

4. **Preservation status** — Software catalog analysis:
   - (a) Timeline: when the 392+ tools were published
   - (b) Link status: alive, dead, partially archived (stacked bar or waffle chart)
   - (c) Citation decay or technology evolution

5. **Benchmarking** — Performance comparison panel:
   - (a) ML search: phylip-rs vs IQ-TREE vs RAxML on datasets of increasing size
   - (b) NJ/distance: scaling behavior
   - (c) Accuracy on simulated data where truth is known

6. **Reimplementation validation** — Test coverage and validation methodology:
   - (a) Module coverage map
   - (b) Example validation: phylip-rs output vs PHYLIP original output vs analytical result

### Supplementary Material
- All 20 case studies (full INSIGHTS.md content)
- Complete software catalog with link status
- Benchmark raw data and scripts
- All 10 interactive demonstrations with output

## 5. Additional Work Needed

### High Priority (Required for Paper)

#### A. Software Catalog Analysis — COMPLETE
- [x] Scrape all 407 entries from Felsenstein's catalog page
- [x] Check each link: alive, dead, redirected, archived (Wayback Machine)
- [x] Extract metadata: name, author, year, language, method type, citation count
- [x] Categorize by method type (ML, parsimony, distance, Bayesian, etc.)
- [x] Analyze temporal patterns: when were tools created, when did they die?
- [x] Identify the "lost" tools — 23 permanently lost
- [x] Create visualization of the catalog's health (Figure 4: 3-panel)

#### B. Benchmarking — COMPLETE
- [x] Select standard benchmark datasets: 36 datasets (10/20/50/100/200/500 taxa, 500-5000 sites)
- [x] Benchmark against: IQ-TREE 3, RAxML-NG, VeryFastTree (plus phylip-rs ML and NJ)
- [x] Metrics: wall time, memory (gtime), log-likelihood (scored under JC69), Robinson-Foulds to true tree
- [x] Run on Apple Silicon (M4), single-threaded, 600s timeout
- [x] All 180 runs complete (benchmarks/results/benchmark_results.csv)
- [x] Figure 5 generated (benchmarks/figures/)

#### C. LLM Process Documentation
- [ ] Review conversation logs from this project
- [ ] Document the iterative process: how many rounds, what kinds of corrections
- [ ] Categorize LLM contributions: code understanding, algorithm extraction, implementation, debugging, insight generation
- [ ] Identify failures: where did the LLM misunderstand the code or algorithms?
- [ ] Quantify: lines of code written per session, test pass rate over time
- [ ] Document the human role: domain expertise, validation, direction-setting

#### D. Validation Against Original PHYLIP
- [ ] Compile original PHYLIP C code
- [ ] Run both implementations on identical inputs
- [ ] Compare outputs: tree topologies, branch lengths, scores, bootstrap values
- [ ] Document any discrepancies and explain them

### Medium Priority (Strengthens Paper)

#### E. Additional Cross-Domain Demonstrations
- [ ] Tumor phylogenetics: apply to single-cell mutation data
- [ ] Cultural evolution: language family reconstruction (extends existing demo)
- [ ] Epidemiology: viral phylogenetics on public SARS-CoV-2 data
- [ ] Document analysis: stemmatic analysis of manuscript traditions

#### F. Community Feedback
- [ ] Share with Joe Felsenstein (if appropriate at this stage)
- [ ] Get feedback from phylogenetics community (Twitter/Mastodon, phylogenetics Slack)
- [ ] Identify potential reviewers and their likely concerns

### Lower Priority (Nice to Have)

#### G. WASM Interactive Demo
- [ ] Browser-based Felsenstein Zone visualization
- [ ] Could be linked from the paper as supplementary material
- [ ] Compile phylip-rs core to WASM

#### H. Teaching Materials
- [ ] Jupyter-style walkthrough of each algorithm
- [ ] Could accompany the paper as educational supplement

## 6. Figures — Detailed Planning

### Figure 1: The LLM Code Archaeology Workflow
- **Type:** Schematic / flow diagram
- **Content:**
  - Input: Legacy codebase (PHYLIP C, 1980-2024)
  - Step 1: LLM reads and annotates source code
  - Step 2: Algorithm extraction and cross-referencing with literature
  - Step 3: Modern reimplementation (Rust)
  - Step 4: Validation (tests, analytical results, comparison with original)
  - Step 5: Insight generation (cross-disciplinary connections)
  - Output: Preserved algorithms, new insights, validated code
- **Style:** Clean, Nature Methods style, horizontal flow

### Figure 2: Algorithmic Discoveries (Multi-panel)
- Panel (a): Kirchhoff contrasts — phylogenetic tree overlaid with circuit diagram, table showing variance = resistance to 8 decimal places
- Panel (b): Genetic code step matrix — 20x20 heatmap, histogram of real code z-score vs 1000 random codes
- Panel (c): Felsenstein Zone — parsimony vs ML accuracy curves, showing convergence to wrong answer

### Figure 3: Cross-Disciplinary Connection Map
- **Type:** Network diagram or chord diagram
- **Content:** PHYLIP algorithms in center, connected to fields they anticipate or relate to
- Connections: belief propagation, circuit theory, information geometry, algebraic geometry, error-correcting codes, isotonic regression, max-flow/min-cut, HMMs

### Figure 4: Software Catalog Preservation
- Panel (a): Timeline of tool publications (histogram by decade)
- Panel (b): Link status (alive/dead/archived) — stacked bar chart
- Panel (c): Category breakdown (ML, parsimony, Bayesian, distance, etc.)

### Figure 5: Benchmarking
- Panel (a): Log-likelihood vs wall time scatter (phylip-rs, IQ-TREE, RAxML-NG)
- Panel (b): Robinson-Foulds accuracy on simulated data
- Panel (c): Scaling behavior (time vs number of taxa)

### Figure 6: Validation
- Panel (a): Test coverage by module (bar chart)
- Panel (b): phylip-rs vs original PHYLIP output comparison (specific example)
- Panel (c): Accuracy on known analytical test cases

## 7. Author and Contributor Considerations

### Potential authors
- **Scott Handley** — PI, conceived project, directed LLM archaeology, domain expertise
- **Claude (Anthropic)** — LLM contributor. Nature Methods policy on AI authorship should be reviewed. Most likely acknowledged rather than listed as author, per current journal policies.

### Potential collaborators to invite
- **Joe Felsenstein** — Original PHYLIP author. Could be invited to write a companion commentary or be consulted for accuracy. His endorsement would be enormously valuable.
- **Anthropic collaborator?** — Someone from Anthropic interested in LLM-for-science applications

### Acknowledgments
- Joe Felsenstein for creating PHYLIP and maintaining the software catalog
- The phylogenetics community
- Anthropic for Claude

## 8. Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| "Just a reimplementation" criticism | Emphasize the *methodology* (LLM archaeology) and *discoveries* (insights not in any paper), not just the code |
| AI authorship controversy | Follow Nature Methods policy exactly; be transparent about LLM role |
| Benchmarking shows poor performance | Frame as pedagogical, not competitive; emphasize zero-dependency constraint |
| Felsenstein objects | Engage early, frame as tribute, give him veto on historical claims |
| Reviewers want more than phylogenetics | Include cross-domain demonstrations and discuss generalizability |
| "Why Rust?" question | Zero-dependency constraint forces understanding of every algorithm; Rust's type system catches errors; performance is adequate |

## 9. Timeline

| Phase | Tasks | Status |
|-------|-------|--------|
| **Phase 1: Additional analysis** | Software catalog scraping, benchmarking | COMPLETE |
| **Phase 2: Writing** | Draft manuscript | COMPLETE (first draft) |
| **Phase 3: Figures** | Create all figures | Figures 3-4 complete; Figures 1-2 needed |
| **Phase 4: Revision** | Polish manuscript, expand, finalize references | Next |
| **Phase 5: Review** | Internal review, feedback from Felsenstein (if invited) | Pending |
| **Phase 6: Submission** | Final polish, cover letter, submission | Pending |

## 10. Key References to Cite

- Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
- Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. *J Mol Evol*, 17, 368-376.
- Felsenstein, J. (1985). Confidence limits on phylogenies: an approach using the bootstrap. *Evolution*, 39, 783-791.
- Pearl, J. (1988). *Probabilistic Reasoning in Intelligent Systems*. Morgan Kaufmann.
- Yang, Z. (1994). Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites. *J Mol Evol*, 39, 306-314.
- Freeland, S.J. & Hurst, L.D. (1998). The genetic code is one in a million. *J Mol Evol*, 47, 238-248.
- Stamatakis, A. (2014). RAxML version 8. *Bioinformatics*, 30, 1312-1313.
- Nguyen, L.-T. et al. (2015). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. *Mol Biol Evol*, 32, 268-274.
- [Nature Methods AI policy — check current version]
- [Papers on scientific software preservation — TBD]
- [Papers on LLMs for code understanding — TBD]

## 11. Open Questions

1. **What is Nature Methods' current policy on AI/LLM contributions?** Need to check before submission.
2. **Should we contact Felsenstein before or after drafting?** Before seems more respectful.
3. **How much of the LLM process can we reconstruct from logs?** Need to check what's available.
4. **Can we compile and run the original PHYLIP C code on modern systems?** Need to test.
5. **What standard benchmark datasets should we use?** TreeBASE? Simulated? Both?
6. **Is there prior work on "code archaeology" we should cite?** Software heritage, mining software repositories.
7. **Should we release phylip-rs as a crate on crates.io?** Increases impact but also invites performance criticism.

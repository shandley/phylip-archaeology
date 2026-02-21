# Validation: phylip-rs vs PHYLIP C

This directory contains infrastructure for validating phylip-rs against the
original PHYLIP C implementation (v3.697) by Joe Felsenstein.

For the full validation report with results, tolerances, known differences,
and reproduction instructions, see **[VALIDATION_REPORT.md](VALIDATION_REPORT.md)**.

## Quick Start

```bash
# 1. Download and compile PHYLIP 3.697
cd validation
bash setup.sh

# 2. Run PHYLIP comparison tests (30 tests)
PHYLIP_EXE_DIR=validation/phylip-3.697/exe cargo test -p phylip-rs --test validation_phylip -- --ignored

# 3. Run all other validation tests (58 tests, no external dependencies)
cargo test -p phylip-rs --test validation_analytical --test validation_classics --test validation_medium
```

## Programs Compared

| PHYLIP Program | phylip-rs Module | What's Compared |
|---|---|---|
| dnadist (JC69) | models::jc69 | Distance matrix values (tol: 1e-3) |
| dnadist (K2P) | models::k2p | Distance ranking preservation |
| neighbor (NJ) | distance::neighbor_joining | Topology (RF=0) + branch lengths (5%) |
| neighbor (UPGMA) | distance::upgma | Topology + ultrametric property |
| fitch | distance::fitch_margoliash | WLS score + topology |
| kitsch | distance::kitsch | WLS score + ultrametric property |
| dnapars | parsimony::wagner | Score (exact) + topology |
| dnapenny | parsimony::branch_and_bound | Score (exact, guaranteed optimal) |
| dnaml | likelihood::pruning | Log-likelihood (same range) |
| dnamlk | likelihood::clock | Clock lnL + ultrametric property |
| dnacomp | compatibility::dna_compat | Compatible sites (±1) |
| dnainvar | invariants | Lake's + Cavender's invariant values |
| protdist | models::protein_distances | Protein distances (tol: 0.05) |
| protpars | parsimony::protein_parsimony | Score (exact) |
| clique | compatibility::clique | Clique size + tree |
| dollop | parsimony::dollo | Dollo score (heuristic) |
| mix | parsimony::wagner | Binary Wagner score (exact) |
| penny | parsimony::branch_and_bound | Binary B&B score (exact) |
| pars | parsimony::multistate | Multistate score (±1) |
| gendist | models::gene_freq | Nei distances (tol: 0.01) |
| restdist | models::restriction | Nei-Li distances (tol: 0.02) |
| contml | comparative::contml | Brownian ML lnL |
| contrast | comparative::contrasts | PIC correlations (tol: 0.15) |
| seqboot+consense | bootstrap+consensus | Pipeline validation |
| treedist | tree::distances | Robinson-Foulds distance (exact) |

# SEQBOOT Analysis - Bootstrap and Jackknife Resampling

**Source file**: `phylip-3.698/src/seqboot.c` (1683 lines)
**Dependencies**: `phylip.h`, `seq.h`

## Overview

SEQBOOT generates resampled data sets for assessing phylogenetic confidence.
It implements several resampling strategies:

1. **Bootstrap** (Felsenstein 1985): Sample sites with replacement
2. **Jackknife** (delete-half or delete-fraction): Sample sites without replacement
3. **Permutation (ILD test)**: Permute character order across all species
4. **Lockhart permutation**: Permute characters independently within each species
5. **Rewrite**: Format conversion (PHYLIP <-> NEXUS <-> XML) without resampling

The bootstrap is the most widely used method for assessing phylogenetic support.
The typical workflow is: SEQBOOT -> tree inference (e.g., DNAML) on each
replicate -> CONSENSE to summarize support. Bootstrap proportions approximate
the probability that a clade is real.

SEQBOOT supports multiple data types (DNA sequences, morphological characters,
restriction sites, gene frequencies) and multiple output formats (PHYLIP,
NEXUS, XML).

## Key Data Structures

### Data Type Enums (lines 35-41)
```c
typedef enum { seqs, morphology, restsites, genefreqs } datatype;
typedef enum { dna, rna, protein } seqtype;
```

### Original Data Storage (lines 121-127)
```c
Char **nodep;      // nodep[spp][sites] - molecular/morphological data
double **nodef;    // nodef[spp][sites] - gene frequency data
Char *factor;      // factor[sites] - factors file content
long *factorr;     // factorr[sites] - factor group assignments [1..groups]
long *alleles;     // alleles[loci] - number of alleles per locus
```

### Weight and Mapping Arrays (lines 117, 131-143)
```c
steptr oldweight;     // original site weights
steptr weight;        // current (resampled) weights
steptr where;         // where[loci] -> first site of each group
steptr how_many;      // how_many[loci] -> number of sites per group

// After removing zero-weight sites:
long *newwhere;       // newwhere[newgroups] -> site index
long *newhowmany;     // newhowmany[newgroups] -> group size

// After bootstrap resampling:
long *newerfactor;    // newerfactor[newersites] -> group assignment
long *newerwhere;     // newerwhere[newergroups] -> site index
long *newerhowmany;   // newerhowmany[newergroups] -> group size
long **charorder;     // charorder[spp][newergroups] - permutation per species
long **sppord;        // sppord[newergroups][spp] - species permutation per site
```

The three-level mapping (original -> new -> newer) handles:
1. Original sites with their weights
2. Sites after removing zero-weight entries (`new*`)
3. Sites after bootstrap resampling (`newer*`)

### Key Configuration Variables (lines 78-101)
- `bootstrap`, `jackknife`, `permute`, `ild`, `lockhart`, `rewrite`:
  mutually exclusive resampling method flags
- `fracsample`: sampling fraction (1.0 for bootstrap, 0.5 for jackknife)
- `regular`: use default sampling fraction
- `blocksize`: block size for block bootstrap (default 1)
- `reps`: number of replicates (default 100)
- `xml`, `nexus`: output format flags
- `justwts`: output weights only (no data rewrite)
- `factors`: use factors file for grouping characters

### Random Number State (line 104)
```c
longer seed;  // random number seed (array of longs)
```
Uses PHYLIP's built-in linear congruential generator via `randum(seed)`.

## Core Algorithms

### 1. bootweights() - Generate Resampled Weights (lines 1069-1151)

This is the central function. It generates a new set of weights representing
one bootstrap/jackknife replicate.

**Bootstrap** (lines 1105-1115):
```c
blocks = fracsample * newgroups / blocksize;
for (i = 1; i <= blocks; i++) {
    j = (long)(newgroups * randum(seed)) + 1;  // random start position
    for (k = 0; k < blocksize; k++) {
        weight[j - 1]++;
        j++;
        if (j > newgroups) j = 1;  // wrap around for block bootstrap
    }
}
```

Standard bootstrap (`blocksize=1`): Draw `newgroups` sites with replacement.
Each site's weight is the number of times it was sampled (Poisson-like
distribution with mean 1).

Block bootstrap (`blocksize>1`): Draw `newgroups/blocksize` blocks of
contiguous sites, each of length `blocksize`. Blocks wrap around at the
sequence end. This preserves local sequence correlation.

Partial bootstrap (`fracsample<1.0`): Sample fewer than n sites, producing
smaller resampled data sets.

**Jackknife** (lines 1078-1101):
```c
q = (long)(newgroups * fracsample + 0.5);  // number of sites to keep
r = newgroups;
p = q / r;
for (i = 0; i < newgroups; i++) {
    if (randum(seed) < p) {
        weight[i]++;
        q--;
    }
    r--;
    p = q / r;  // update probability to maintain exact count
}
```

This uses a streaming selection algorithm (similar to reservoir sampling)
that selects exactly `q` of `n` sites without replacement. Each site's
weight is 0 or 1. The probability is updated after each decision to ensure
the exact target count is achieved.

**Permutation** (lines 1102-1104): All weights are 1 (no resampling of
characters). Species order is permuted instead (see `sppermute()`).

**Rewrite** (lines 1116-1118): All weights are 1. Only format conversion.

After weight generation, the function computes `newergroups` and `newersites`
(lines 1121-1126) and allocates the `newer*` arrays accordingly.

### 2. sppermute() / charpermute() - Permutation Tests (lines 1167-1176)

**sppermute()** (line 1168): Permutes the species order for a given character
group. Used for the ILD (Incongruence Length Difference) test: if characters
from different partitions are phylogenetically congruent, permuting species
assignments should not improve the fit.

**charpermute()** (line 1174): Permutes the character order for a given
species. Used for Lockhart's test of compositional heterogeneity.

Both use Fisher-Yates shuffle via `permute_vec()` (lines 1154-1164):
```c
for (i = 1; i < n; i++) {
    k = (long)((i+1) * randum(seed));
    swap(a[i], a[k]);
}
```

### 3. writedata() - Output Resampled Data (lines 1179-1330)

Writes one resampled data set in the chosen output format:

1. **Header**: Writes species count and site count. For NEXUS, writes
   the full NEXUS header with `#NEXUS`, `BEGIN DATA`, `DIMENSIONS`,
   and `FORMAT` blocks.

2. **Data**: For each species, writes the resampled sequence by looking up
   characters via the `charorder` and `newerwhere` mappings:
   ```c
   charstate = nodep[sppord[charorder[j][k]][j] - 1]
                    [newerwhere[charorder[j][k]] + n2];
   ```
   This double indirection handles both character resampling (via
   `charorder`) and species permutation (via `sppord`).

3. **Format specifics**:
   - PHYLIP: Standard interleaved or sequential format
   - NEXUS: Full NEXUS block with MATRIX
   - XML: `<alignment>` with `<sequence>` elements

### 4. writeweights() - Output Bootstrap Weights (lines 1333-1373)

When `justwts` is true, outputs weights instead of full data sets. This is
more space-efficient for large datasets. The weights are encoded as single
characters: '0'-'9' for weights 0-9, 'A'-'Z' for weights 10-35.

### 5. bootwrite() - Main Resampling Loop (lines 1549-1601)

Orchestrates the complete resampling process:

```c
for (rr = 1; rr <= reps; rr++) {
    bootweights();                    // generate weights
    initialize charorder[i][j] = j;   // identity permutation

    if (ild)
        charpermute(0, newergroups);  // permute chars (same for all spp)
    if (lockhart)
        for each species: charpermute(i, newergroups);  // independent permutes

    if (!justwts || permute)
        writedata();                  // write resampled data
    if (justwts)
        writeweights();              // write weights only

    if (categories) writecategories();
    if (factors) writefactors();
    if (mixture) writeauxdata(mixdata, outmixfile);
    if (ancvar) writeauxdata(ancdata, outancfile);
}
```

### 6. inputoptions() - Factor Group Handling (lines 566-631)

Sets up the factor-to-site mapping. For gene frequency data, groups are
defined by loci. For other data with a factors file, groups are defined by
the factors file. Sites within the same factor group are always resampled
together (they cannot be independently bootstrapped).

This is critical for:
- Gene frequency data (alleles at one locus must stay together)
- Restriction site data with multiple enzymes
- Morphological characters with linked multi-state codings

## Data Type Support

### Sequences (seqs)
Standard DNA/RNA/protein sequences. Each character is independently
resampleable. Supports IUPAC ambiguity codes.

### Morphology
Discrete morphological characters. Can use factors file to group linked
characters. Supports mixture and ancestor files.

### Restriction Sites (restsites)
Restriction enzyme presence/absence data. Enzyme count may be in the input
file (`enzymes` flag).

### Gene Frequencies (genefreqs)
Allele frequency data. Sites are grouped by locus (via `alleles[]` array).
Resampling operates at the locus level, keeping all alleles at a locus
together.

## Block Bootstrap

The block bootstrap (`blocksize > 1`) resamples contiguous blocks of sites
rather than individual sites. This is appropriate when there is serial
correlation along the sequence (e.g., codon structure, secondary structure).
Blocks wrap around from the end to the beginning of the sequence (line 1112:
`if (j > newgroups) j = 1`), treating the sequence as circular.

## Memory Management

The three-tier allocation scheme handles dynamic sizing:

1. **allocrest()** (lines 903-919): Allocates arrays sized by original
   `sites` and `loci`. Called once.

2. **allocnew()** (lines 942-950): Allocates arrays sized by `newgroups`
   (after removing zero-weight sites). Called once.

3. **allocnewer()** (lines 964-1007): Allocates/reallocates arrays sized
   by `newergroups` and `newersites` (after bootstrapping). Called per
   replicate, but uses static variables to avoid unnecessary reallocation
   when sizes don't increase.

## I/O

**Input**: PHYLIP format data file (sequences, morphology, restriction
sites, or gene frequencies). Optional weights, categories, mixture,
ancestors, and factors files.

**Output**: Multiple resampled data sets concatenated in a single file,
or just weights if `justwts` is true. Supported formats:
- PHYLIP interleaved/sequential (default)
- NEXUS (with proper header)
- XML (for sequences only)

## Complexity

- **Time**: O(reps * spp * sites) total. Each replicate requires O(sites)
  for weight generation and O(spp * sites) for data output.
- **Space**: O(spp * sites) for the original data, plus O(sites) for
  weight and mapping arrays.

## Modernization Notes for Rust Reimplementation

1. **Use a proper RNG**: Replace PHYLIP's custom LCG with a well-tested
   RNG from the `rand` crate (e.g., `rand::rngs::StdRng` for
   reproducibility with a seed).

2. **Resampling as iterators**: Express bootstrap/jackknife as iterators
   over weight vectors:
   ```rust
   fn bootstrap_weights(n: usize, rng: &mut impl Rng) -> Vec<usize> {
       let mut weights = vec![0; n];
       for _ in 0..n {
           weights[rng.gen_range(0..n)] += 1;
       }
       weights
   }
   ```

3. **Output format trait**:
   ```rust
   trait SequenceWriter {
       fn write_header(&mut self, spp: usize, sites: usize);
       fn write_sequence(&mut self, name: &str, seq: &[u8]);
       fn write_footer(&mut self);
   }
   ```
   Implement for PHYLIP, NEXUS, and XML.

4. **Eliminate three-tier mapping**: The original->new->newer mapping is
   confusing. Use a single `ResampledDataset` struct that handles the
   complete mapping from original data to output.

5. **Parallel replicate generation**: Each replicate is independent.
   Generate replicates in parallel using `rayon`, then serialize output
   (or write to separate files).

6. **Factor groups as types**: Model factor groups as a proper type:
   ```rust
   struct FactorGroups {
       group_of_site: Vec<usize>,    // site -> group
       sites_in_group: Vec<Range<usize>>, // group -> site range
   }
   ```

7. **Streaming output**: For large datasets, don't store all replicates
   in memory. Write each replicate to the output file as it's generated.

8. **Config struct**: Replace the 20+ boolean flags with a proper
   configuration enum:
   ```rust
   enum ResamplingMethod {
       Bootstrap { block_size: usize, fraction: f64 },
       Jackknife { fraction: f64 },
       Permute,
       ILD,
       Lockhart,
       Rewrite,
   }
   ```

9. **Validation**: Add input validation (e.g., gene frequencies sum to
   <= 1.0, sequence characters are valid) using Rust's type system rather
   than runtime checks scattered throughout the code.

10. **Test with known seeds**: The original code's random number generator
    is deterministic given a seed. Preserve this property for regression
    testing, but also support modern RNGs for production use.

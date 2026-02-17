# PHYLIP Algorithm Map

**Version:** 3.698 (September 2019)
**Purpose:** Maps each algorithm to the programs that implement it, with references to specific source files and key functions. Designed to guide algorithm extraction and Rust reimplementation.

---

## Table of Contents

1. [Maximum Likelihood Algorithms](#1-maximum-likelihood-algorithms)
2. [Parsimony Algorithms](#2-parsimony-algorithms)
3. [Distance Matrix Algorithms](#3-distance-matrix-algorithms)
4. [Distance Computation Algorithms](#4-distance-computation-algorithms)
5. [Compatibility Algorithms](#5-compatibility-algorithms)
6. [Consensus and Tree Comparison Algorithms](#6-consensus-and-tree-comparison-algorithms)
7. [Resampling and Data Transformation Algorithms](#7-resampling-and-data-transformation-algorithms)
8. [Tree Visualization Algorithms](#8-tree-visualization-algorithms)
9. [Interactive Tree Manipulation](#9-interactive-tree-manipulation)
10. [Comparative Methods](#10-comparative-methods)
11. [Phylogenetic Invariants](#11-phylogenetic-invariants)
12. [Core Infrastructure](#12-core-infrastructure)
13. [Algorithm-to-Program Cross-Reference](#13-algorithm-to-program-cross-reference)
14. [Dependency Graph](#14-dependency-graph)
15. [Reimplementation Priority Recommendations](#15-reimplementation-priority-recommendations)

---

## 1. Maximum Likelihood Algorithms

### 1.1 DNA Maximum Likelihood (F84/HKY model)

**Programs:** Dnaml, Dnamlk

**Algorithm:** Computes the likelihood of a phylogenetic tree given DNA sequence data, using the F84 substitution model (equivalent to HKY85 in the limit). Uses the "pruning algorithm" (Felsenstein, 1981) for likelihood computation and the algorithm of Felsenstein and Churchill (1996) for speed optimization via conditional likelihood storage at internal nodes.

**Substitution Model Features:**
- Unequal base frequencies (freqa, freqc, freqg, freqt)
- Transition/transversion rate ratio
- Gamma-distributed rate variation across sites
- Hidden Markov Model of rate variation
- Invariant sites proportion
- Multiple rate categories

**Source:** `dnaml.c` (2,619 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~1800 | Main tree search driver |
| `evaluate()` | ~800 | Compute log-likelihood of current tree |
| `nuview()` | ~700 | Compute conditional likelihoods at internal node (pruning algorithm core) |
| `insert_()` | ~1100 | Insert taxon into tree |
| `addtraverse()` | ~1200 | Try all insertion points for a taxon |
| `rearrange()` | ~1300 | Local rearrangement (SPR) |
| `globrearrange()` | ~1400 | Global rearrangement search |
| `smooth()` | ~950 | Optimize branch lengths iteratively |
| `update()` | ~900 | Optimize single branch length |
| `getoptions()` | ~100 | Parse user menu selections |
| `getbasefreqs()` | in `seq.c` | Compute empirical base frequencies |
| `makevalues()` | in `seq.c` | Initialize tip conditional likelihoods |

**Source (clock version):** `dnamlk.c` (2,252 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~1600 | Main tree search with clock constraint |
| `evaluate()` | ~700 | Log-likelihood with clock-constrained branches |
| `nuview()` | ~600 | Conditional likelihood computation |
| `setuptree()` | ~300 | Initialize clock-constrained tree |
| `cur_node_eval()` | in `mlclock.c` | Evaluate node time under clock |
| `evaluate_tyme()` | in `mlclock.c` | Optimize node time |
| `mlk_printree()` | in `printree.c` | Display clock tree |

**Library dependencies:** `seq.c`, `phylip.c`, `mlclock.c` (clock version), `printree.c` (clock version)

---

### 1.2 Protein Maximum Likelihood (JTT/PAM/PMB models)

**Programs:** Proml, Promlk

**Algorithm:** Maximum likelihood for amino acid sequences. Uses the same pruning algorithm and tree search strategy as Dnaml, but with 20x20 amino acid substitution rate matrices instead of 4x4 nucleotide matrices.

**Substitution Models:**
- Jones-Taylor-Thornton (JTT) model
- Dayhoff PAM model
- PMB (Henikoff/Tillier) model
- Gamma-distributed rate variation
- HMM rate variation
- Invariant sites

**Source:** `proml.c` (3,246 lines -- the largest program in PHYLIP)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~2400 | Main tree search driver |
| `evaluate()` | ~1100 | Log-likelihood computation for protein data |
| `nuview()` | ~900 | Conditional likelihood at internal node (20-state) |
| `insert_()` | ~1500 | Insert taxon into tree |
| `addtraverse()` | ~1600 | Traverse insertion points |
| `rearrange()` | ~1700 | Local SPR rearrangement |
| `globrearrange()` | ~1800 | Global rearrangement |
| `smooth()` | ~1300 | Branch length optimization |
| `update()` | ~1200 | Single branch optimization |
| `getbasefreqs()` | ~700 | Set up amino acid frequencies from model |
| `transition()` | ~400 | Set up transition probability matrix from model |

**Source (clock version):** `promlk.c` (2,998 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~2200 | Tree search with clock |
| `evaluate()` | ~900 | Clock-constrained likelihood |
| `nuview()` | ~800 | Conditional likelihood computation |

**Library dependencies:** `seq.c`, `phylip.c`, `mlclock.c` (clock), `printree.c` (clock)

---

### 1.3 Restriction Sites Maximum Likelihood

**Program:** Restml

**Algorithm:** Maximum likelihood for restriction site presence/absence data, using the Jukes-Cantor model. Does not distinguish transitions from transversions. Very computationally intensive.

**Source:** `restml.c` (2,528 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~1800 | Main tree search |
| `evaluate()` | ~800 | Log-likelihood for restriction sites |
| `nuview()` | ~600 | Conditional likelihood for site data |
| `insert_()` | ~1000 | Taxon insertion |
| `addtraverse()` | ~1100 | Traverse insertion points |
| `rearrange()` | ~1200 | Local rearrangement |
| `globrearrange()` | ~1300 | Global rearrangement |
| `smooth()` | ~900 | Branch length optimization |

**Library dependencies:** `seq.c`, `phylip.c`

---

### 1.4 Continuous Character / Gene Frequency Maximum Likelihood

**Program:** Contml

**Algorithm:** Restricted maximum likelihood (REML) under the Brownian motion model for continuous characters, or the genetic drift model of Edwards and Cavalli-Sforza (1964) for gene frequencies. Assumes characters evolve independently at equal rates.

**Source:** `contml.c` (1,560 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~1362 | Main tree search |
| `evaluate()` | ~639 | Log-likelihood under Brownian motion / drift |
| `sumlikely()` | ~588 | Sum likelihoods across alleles/characters |
| `distance()` | ~665 | Compute expected distance between nodes |
| `makedists()` | ~691 | Compute all pairwise expected distances |
| `makebigv()` | ~704 | Compute "big V" variance parameter |
| `nuview()` | ~772 | Conditional computation at internal node |
| `insert_()` | ~831 | Insert taxon |
| `addtraverse()` | ~931 | Traverse insertion points |
| `rearrange()` | ~1022 | Local rearrangement |
| `globrearrange()` | ~958 | Global rearrangement |
| `smooth()` | ~821 | Branch length optimization |
| `transformgfs()` | ~485 | Transform gene frequencies for analysis |

**Library dependencies:** `cont.c`, `phylip.c`

---

## 2. Parsimony Algorithms

### 2.1 DNA Parsimony (Wagner criterion)

**Programs:** Dnapars (heuristic), Dnapenny (branch-and-bound), Dnamove (interactive), Dnacomp (compatibility variant)

**Algorithm:** Counts the minimum number of nucleotide substitutions needed to explain sequence data on a given tree using the Fitch (1971) algorithm. Each site is treated independently. Gaps can be treated as a fifth state. Transversion parsimony is optionally available.

#### Dnapars (heuristic search)
**Source:** `dnapars.c` (1,663 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~1200 | Main heuristic search using stepwise addition + rearrangement |
| `evaluate()` | ~650 | Count total steps on tree |
| `fillin()` | in `seq.c` | Fitch algorithm: compute intersection/union of state sets at internal nodes |
| `sumnsteps()` | in `seq.c` | Sum number of steps across sites |
| `insert_()` | ~800 | Insert taxon |
| `addtraverse()` | ~850 | Try all insertion points |
| `rearrange()` | ~900 | Local SPR rearrangement |
| `globrearrange()` | ~950 | Global SPR rearrangement |
| `add()` | in `seq.c` | Add branch to tree (topology change) |
| `re_move()` | in `seq.c` | Remove branch from tree |
| `tryadd()` | ~750 | Try adding at a specific location |
| `describe()` | ~1050 | Output tree description with reconstructed states |

#### Dnapenny (branch-and-bound)
**Source:** `dnapenny.c` (819 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~600 | Branch-and-bound search driver |
| `evaluate()` | ~350 | Score partial tree for bounding |
| `addpreorder()` | ~500 | Add species in preorder traversal with bounding |
| `add()` | in `seq.c` | Topology manipulation |

#### Dnamove (interactive)
**Source:** `dnamove.c` (2,363 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~1800 | Interactive loop: display tree, accept commands |
| `evaluate()` | ~800 | Score current tree |
| `add()` | in `seq.c` | Add branch |
| `re_move()` | in `seq.c` | Remove branch |

#### Dnacomp (compatibility)
**Source:** `dnacomp.c` (1,179 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~850 | Heuristic search for most compatible tree |
| `evaluate()` | ~500 | Count number of compatible sites |
| `insert_()` | ~600 | Insert taxon |
| `addtraverse()` | ~650 | Try insertion points |
| `rearrange()` | ~700 | Local rearrangement |

**Library dependencies:** `seq.c`, `phylip.c`. Dnamove also uses `moves.c`.

---

### 2.2 Protein Parsimony

**Program:** Protpars

**Algorithm:** Parsimony for protein sequences. Counts only nucleotide changes that alter the amino acid, using a cost matrix derived from the genetic code. Intermediate between the Eck-Dayhoff (1966) approach (any amino acid change costs 1) and the Fitch (1971) approach (count nucleotide changes).

**Source:** `protpars.c` (1,962 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~1400 | Main heuristic search |
| `evaluate()` | ~700 | Score tree by protein parsimony |
| `fillin()` | in `seq.c` | State set computation at internal nodes |
| `insert_()` | ~900 | Insert taxon |
| `addtraverse()` | ~950 | Traverse insertion points |
| `rearrange()` | ~1050 | Local rearrangement |
| `globrearrange()` | ~1100 | Global rearrangement |
| `protpreorder()` | ~800 | Preorder traversal for ancestral state reconstruction |

**Library dependencies:** `seq.c`, `phylip.c`

---

### 2.3 Discrete Character Parsimony (Wagner / Camin-Sokal)

**Programs:** Mix (heuristic), Penny (branch-and-bound), Move (interactive), Pars (multistate heuristic)

**Algorithm:**
- **Wagner parsimony:** Changes between states 0 and 1 are both allowed and equally weighted. Finds trees requiring the fewest total changes.
- **Camin-Sokal parsimony:** Only changes 0->1 are allowed (irreversible). Rooted trees only.
- **Mixed:** Each character can independently use either criterion.

#### Mix (heuristic)
**Source:** `mix.c` (1,180 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~850 | Heuristic search with stepwise addition + SPR |
| `evaluate()` | ~500 | Count steps under Wagner/Camin-Sokal |
| `fillin()` | in `wagner.c` | State set computation for 0/1 characters |
| `postorder()` | in `wagner.c` | Postorder traversal for step counting |
| `add()` | in `disc.c` | Add branch to tree |
| `re_move()` | in `disc.c` | Remove branch from tree |

#### Penny (branch-and-bound)
**Source:** `penny.c` (843 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~600 | Branch-and-bound exact search |
| `evaluate()` | ~350 | Score for bounding |
| `addpreorder()` | ~500 | Add with bound checking |

#### Move (interactive)
**Source:** `move.c` (1,655 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~1200 | Interactive loop |
| `evaluate()` | ~600 | Score current tree |

#### Pars (multistate)
**Source:** `pars.c` (1,649 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~1200 | Heuristic search for multistate characters (up to 8 states) |
| `evaluate()` | ~600 | Count steps for multistate Wagner parsimony |
| `fillin()` | in `discrete.c` | Multistate Fitch algorithm |
| `sumnsteps()` | in `discrete.c` | Sum steps across characters |
| `add()` | in `discrete.c` | Add branch |
| `re_move()` | in `discrete.c` | Remove branch |

**Library dependencies:** Mix, Penny use `disc.c`, `wagner.c`, `phylip.c`. Move also uses `moves.c`. Pars uses `discrete.c`, `phylip.c`.

---

### 2.4 Dollo / Polymorphism Parsimony

**Programs:** Dollop (heuristic), Dolpenny (branch-and-bound), Dolmove (interactive)

**Algorithm:**
- **Dollo parsimony:** A complex feature can arise only once but can be lost multiple times. Named after Louis Dollo. The tree is rooted. (Le Quesne, 1974; Farris, 1977).
- **Polymorphism parsimony:** A generalization where the ancestral species may have been polymorphic (had both states present), with the two types of events being gain of polymorphism and loss of one state.

#### Dollop (heuristic)
**Source:** `dollop.c` (1,023 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~750 | Heuristic search for Dollo/polymorphism parsimony |
| `evaluate()` | ~400 | Score under Dollo criterion |
| `fillin()` | in `dollo.c` | Dollo state reconstruction at internal nodes |
| `postorder()` | in `dollo.c` | Postorder traversal |
| `correct()` | in `dollo.c` | Correct counts for Dollo criterion |
| `add()` | in `disc.c` | Add branch |
| `re_move()` | in `disc.c` | Remove branch |

#### Dolpenny (branch-and-bound)
**Source:** `dolpenny.c` (736 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~550 | Branch-and-bound exact search |
| `evaluate()` | ~300 | Score for bounding |
| `addpreorder()` | ~450 | Add with bounding |

#### Dolmove (interactive)
**Source:** `dolmove.c` (1,602 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~1200 | Interactive loop |
| `evaluate()` | ~500 | Score current tree |

**Library dependencies:** `disc.c`, `dollo.c`, `phylip.c`. Dolmove also uses `moves.c`.

---

## 3. Distance Matrix Algorithms

### 3.1 Fitch-Margoliash / Least Squares

**Programs:** Fitch, Kitsch

**Algorithm:** Fits an additive tree to a distance matrix by minimizing a weighted least-squares criterion. The Fitch-Margoliash criterion weights residuals by 1/D^2 (where D is the observed distance). Other weightings are available. Branch lengths are iteratively optimized.

**Fitch** (no clock): `fitch.c` (1,203 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~850 | Main search: stepwise addition + rearrangement |
| `evaluate()` | ~450 | Compute weighted sum of squares |
| `insert_()` | ~550 | Insert taxon and optimize branch lengths |
| `addtraverse()` | ~600 | Try all insertion points |
| `rearrange()` | ~700 | Local SPR rearrangement |
| `globrearrange()` | ~750 | Global rearrangement |
| `smooth()` | ~400 | Iterative branch length optimization |
| `update()` | ~350 | Optimize single branch length |
| `fillnamedist()` | ~300 | Read distance matrix |

**Kitsch** (with clock): `kitsch.c` (1,028 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~750 | Search with clock constraint |
| `evaluate()` | ~400 | FM criterion with ultrametric constraint |
| `insert_()` | ~500 | Insert with clock |
| `addtraverse()` | ~550 | Traverse insertion points |
| `rearrange()` | ~650 | Local rearrangement |

**Minimum Evolution option:** Both Fitch and Kitsch can use the Minimum Evolution criterion, which minimizes the sum of branch lengths rather than the least squares criterion.

**Library dependencies:** `dist.c`, `phylip.c`

---

### 3.2 Neighbor-Joining / UPGMA

**Program:** Neighbor

**Algorithm:**
- **Neighbor-Joining (Saitou and Nei, 1987):** Agglomerative clustering that does not assume an evolutionary clock. At each step, joins the pair of nodes that minimizes the total branch length. Produces an unrooted tree.
- **UPGMA:** Unweighted Pair Group Method with Arithmetic Mean. Agglomerative clustering assuming an evolutionary clock. Joins the pair with smallest average distance. Produces a rooted tree.

**Source:** `neighbor.c` (629 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~400 | Main clustering loop |
| `jointree()` | ~350 | NJ: join pair of nodes that minimizes tree length |
| `choose_best_pair()` | ~300 | Select pair to join next |
| `setup_tree()` | ~200 | Initialize star topology |
| `describe()` | ~500 | Output tree with branch lengths |

**Library dependencies:** `dist.c`, `phylip.c`

---

## 4. Distance Computation Algorithms

### 4.1 DNA Distance Models

**Program:** Dnadist

**Algorithm:** Computes pairwise evolutionary distances from aligned DNA sequences using various substitution models.

**Models implemented:**
- **Jukes-Cantor (1969):** Equal base frequencies, equal substitution rates. Distance = -3/4 * ln(1 - 4/3 * p), where p = observed proportion of differences.
- **Kimura 2-parameter (1980):** Equal base frequencies, different transition/transversion rates.
- **F84 (Felsenstein, 1984):** Unequal base frequencies, transition/transversion rate ratio. Equivalent model to what Dnaml uses.
- **LogDet (Barry and Hartigan, 1987; Lake, 1994; Lockhart et al., 1994):** Based on the log of the determinant of the divergence matrix. Consistent even when base composition varies among lineages.

**Rate correction options:**
- Gamma-distributed rates across sites
- Gamma + invariant sites
- Multiple rate categories

**Source:** `dnadist.c` (1,335 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `makedists()` | ~800 | Main distance computation loop |
| `makev()` | ~700 | Compute distance for one pair under chosen model |
| `getbasefreqs()` | in `seq.c` | Compute/set base frequencies |
| `empiricalfreqs()` | in `seq.c` | Estimate frequencies from data |
| `getoptions()` | ~100 | Model selection menu |
| `inittable()` | ~600 | Initialize rate category table |

**Key variables:**
- `jukes`, `kimura`, `f84`, `logdet`: Model selection booleans
- `ttratio`: Transition/transversion ratio
- `freqa`, `freqc`, `freqg`, `freqt`: Base frequencies
- `gama`, `invar`: Rate variation flags
- `cvi`: Coefficient of variation for gamma distribution
- `invarfrac`: Proportion of invariant sites

**Library dependencies:** `seq.c`, `phylip.c`

---

### 4.2 Protein Distance Models

**Program:** Protdist

**Algorithm:** Computes pairwise evolutionary distances from aligned protein sequences.

**Models implemented:**
- **JTT (Jones, Taylor, Thornton, 1992):** Empirical amino acid substitution matrix derived from protein databases.
- **PMB (Henikoff and Tillier):** Another empirical matrix.
- **Dayhoff PAM (1978):** The original PAM (Point Accepted Mutation) matrix.
- **Kimura (1983):** Approximation to PAM distances based on proportion of different amino acids.
- **Categories model:** Groups amino acids into categories and models changes between/within categories.

**Source:** `protdist.c` (2,006 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `makedists()` | ~1947 | ML protein distance computation |
| `makev()` | ~1600 | Compute distance for one pair |
| `transition()` | ~400 | Set up 20x20 transition probability matrix |
| `maketrans()` | ~500 | Build eigenvalue decomposition of rate matrix |
| `qreigen()` | ~1628 | QR eigenvector/eigenvalue decomposition for symmetric matrix |
| `getoptions()` | ~100 | Model selection menu |
| `reallocategories()` | ~200 | Set up amino acid category assignments |

**Amino acid model data:**
- Lines ~100-400 contain the JTT, PAM, and PMB rate matrices as hard-coded double arrays
- Line ~320: dcmut version of PAM model

**Library dependencies:** `seq.c`, `phylip.c`

---

### 4.3 Restriction Sites/Fragments Distance

**Program:** Restdist

**Algorithm:** Computes distances from restriction site or restriction fragment data using the method of Nei and Li (1979) with modifications. Can also handle RAPD and AFLP data.

**Source:** `restdist.c` (690 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `makedists()` | ~400 | Main distance computation |
| `makev()` | ~350 | Compute one pairwise distance |
| `getoptions()` | ~100 | Options menu |

**Library dependencies:** `seq.c`, `phylip.c`

---

### 4.4 Genetic Distance

**Program:** Gendist

**Algorithm:** Computes genetic distances from allele frequency data.

**Distance measures:**
- **Nei's genetic distance (1972):** Based on an infinite-isoalleles neutral mutation model. D = -ln(I), where I is Nei's identity.
- **Cavalli-Sforza chord measure (1967):** Geometric distance on a hypersphere. Appropriate for pure genetic drift (no mutation).
- **Reynolds, Weir, and Cockerham (1983):** Also for pure genetic drift. Based on coancestry coefficients.

**Source:** `gendist.c` (443 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `makedists()` | ~250 | Main distance computation |
| `getoptions()` | ~60 | Distance measure selection |
| `inputdata()` | ~150 | Read gene frequency data |

**Library dependencies:** `phylip.c`

---

## 5. Compatibility Algorithms

### 5.1 Clique Compatibility

**Program:** Clique

**Algorithm:** Finds the largest clique of mutually compatible two-state characters and the trees they imply. Two characters are compatible if there exists a tree on which both can evolve without homoplasy. Based on Le Quesne (1969), Estabrook, Johnson, and McMorris (1976a, 1976b).

**Source:** `clique.c` (1,532 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `Compatible()` | ~571 | Test if two characters are compatible |
| `SetUp()` | ~626 | Build compatibility matrix |
| `GetMaxCliques()` | ~1425 | Find all maximum cliques |
| `Gen1()` / `Gen2()` | ~693, ~1368 | Generate cliques by recursive search |
| `Intersect()` | ~669 | Set intersection for clique building |
| `CountStates()` | ~679 | Count states in character |
| `reconstruct()` | ~1020 | Reconstruct tree from clique |
| `DoAll()` | ~1318 | Process all cliques found |

**Library dependencies:** `disc.c`, `phylip.c`

---

### 5.2 DNA Compatibility

**Program:** Dnacomp

**Algorithm:** Compatibility method for DNA sequences. Evaluates each tree topology by counting how many sites are compatible (can evolve without homoplasy on that tree). Searches for the tree with the most compatible sites.

**Source:** `dnacomp.c` (1,179 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `evaluate()` | ~500 | Count compatible sites |
| `maketree()` | ~850 | Heuristic search |
| `insert_()` | ~600 | Insert taxon |
| `addtraverse()` | ~650 | Try insertion points |
| `rearrange()` | ~700 | Local rearrangement |

**Library dependencies:** `seq.c`, `phylip.c`

---

## 6. Consensus and Tree Comparison Algorithms

### 6.1 Consensus Trees

**Program:** Consense

**Algorithm:** Majority-rule consensus tree (Margush and McMorris, 1981). Also supports strict consensus. Each bipartition (split) is included in the consensus tree if it appears in more than the threshold proportion of input trees (default 50% for majority rule, 100% for strict).

**Core logic is in the library file `cons.c` (1,557 lines).**

**Source:** `consense.c` (443 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `getoptions()` | ~64 | Set consensus type and threshold |
| `count_siblings()` | ~250 | Count children at each node |
| `treeout()` | ~276 | Write consensus tree in Newick format |

**Key functions in `cons.c`:**
| Function | Line | Purpose |
|----------|------|---------|
| `enterpartition()` | ~200 | Record a bipartition from an input tree |
| `censor()` | ~269 | Remove partitions below threshold |
| `compress()` | ~290 | Compress partition table |
| `sort()` | ~319 | Sort partitions by frequency |
| `compatible()` | ~349 | Check if two partitions are compatible |
| `eliminate()` | ~383 | Eliminate incompatible partitions |
| `reconstruct()` | ~598 | Build consensus tree from retained partitions |
| `bigsubset()` | ~459 | Find subset relationships between partitions |
| `recontraverse()` | ~502 | Recursive tree building |
| `printset()` | ~416 | Print partition sets |
| `printree()` | ~741 | Display consensus tree |

**Library dependencies:** `cons.c`, `phylip.c`

---

### 6.2 Tree Distance Computation

**Program:** Treedist

**Algorithm:** Computes distances between trees read from an input file.

**Distance measures:**
- **Symmetric Difference (Robinson and Foulds, 1981):** Counts bipartitions present in one tree but not the other. Purely topological, ignores branch lengths.
- **Branch Score Distance (Kuhner and Felsenstein, 1994):** Uses branch lengths. Sum of squared differences of branch lengths for matching bipartitions, plus squared lengths of unmatched bipartitions.

**Source:** `treedist.c` (1,298 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~900 | Main computation loop over tree pairs |
| `compute_distances()` | ~700 | Compute both distance measures |

**Shares extensive code with `cons.c` for bipartition handling.**

**Library dependencies:** `cons.c`, `phylip.c`

---

## 7. Resampling and Data Transformation Algorithms

### 7.1 Bootstrap / Jackknife / Permutation

**Program:** Seqboot

**Algorithm:** Creates multiple resampled datasets for statistical confidence assessment.

**Resampling methods:**
- **Bootstrap:** Resample sites/characters with replacement. The standard method for assessing phylogenetic confidence (Felsenstein, 1985a).
- **Delete-half jackknife:** Randomly delete half the sites without replacement.
- **Permutation:** Randomly permute species assignments at each site. For testing null hypotheses.

**Data types supported:** Molecular sequences, restriction sites, gene frequencies, discrete characters.

**Source:** `seqboot.c` (1,683 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~1200 | Main resampling loop (create nreps datasets) |
| `bootweights()` | ~900 | Generate bootstrap weights for sites |
| `permute()` | ~950 | Permutation resampling |
| `seqboot_inputdata()` | ~500 | Read original data set |
| `getoptions()` | ~100 | Select resampling method and parameters |

**Key type definitions:**
```c
typedef enum { seqs, morphology, restsites, genefreqs } datatype;
typedef enum { dna, rna, protein } seqtype;
```

**Library dependencies:** `seq.c`, `phylip.c`

---

### 7.2 Character Factoring

**Program:** Factor

**Algorithm:** Converts multistate characters into binary (0,1) characters using a specified character-state tree (a graph showing which states can transform into which). Each binary factor represents one edge in the character-state tree.

**Source:** `factor.c` (594 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~400 | Process character-state trees |
| `getoptions()` | ~100 | Options |
| `inputoptions()` | ~200 | Read factor specifications |

**Library dependencies:** `phylip.c`

---

## 8. Tree Visualization Algorithms

### 8.1 Rooted Tree Plotting

**Program:** Drawgram

**Algorithm:** Plots rooted phylogenies in various styles: cladogram, phenogram, curvogram, eurogram, swoopogram, circular (radial) tree. Uses recursive coordinate assignment and various output drivers.

**Source:** `drawgram.c` (2,179 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `plottree()` | ~1500 | Main tree plotting driver |
| `calculate()` | ~1200 | Compute node coordinates |
| `drawline()` | ~800 | Draw a branch line |
| `getoptions()` | ~100 | Interactive option menu |

**Library dependencies:** `draw.c`, `draw2.c`, `phylip.c`

---

### 8.2 Unrooted Tree Plotting

**Program:** Drawtree

**Algorithm:** Plots unrooted tree diagrams using the equal-daylight algorithm for optimal branch angle assignment. Supports various output formats.

**Source:** `drawtree.c` (3,151 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `plottree()` | ~2000 | Main plotting driver |
| `calculate()` | ~1500 | Compute coordinates with equal-daylight |
| `improvtree()` | ~1200 | Iteratively improve branch angles |
| `getoptions()` | ~100 | Interactive option menu |

**Output formats (in `draw.c` and `draw2.c`):** PostScript, SVG, PCL (HP LaserJet), PICT, BMP, Xbm, FIG, VRML, POV-Ray, Tektronix, X Windows, REGIS.

**Library dependencies:** `draw.c`, `draw2.c`, `phylip.c`

---

## 9. Interactive Tree Manipulation

### 9.1 Tree Editing

**Program:** Retree

**Algorithm:** Interactive tree editor. Allows the user to:
- Rearrange tree topology (prune and regraft)
- Reroot the tree
- Rotate branches
- Midpoint-root
- Rename taxa
- Write modified tree to file

Does **not** evaluate any optimality criterion -- purely for manipulation.

**Source:** `retree.c` (3,329 lines -- second largest program)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~2500 | Interactive command loop |
| `reroot()` | ~1500 | Reroot tree at specified node |
| `rearrange()` | ~1800 | Prune and regraft branch |
| `coordinates()` | ~1200 | Compute display coordinates |
| `drawline()` | ~1300 | Draw ASCII tree |
| `printree()` | ~1400 | Display current tree |
| `treeout()` | ~2000 | Write tree in Newick format |

**Library dependencies:** `moves.c`, `phylip.c`

---

## 10. Comparative Methods

### 10.1 Independent Contrasts

**Program:** Contrast

**Algorithm:** Felsenstein's (1985d) phylogenetically independent contrasts. Given a phylogenetic tree and continuous character data at the tips, computes standardized contrasts between sister taxa that are statistically independent under a Brownian motion model of evolution. These contrasts can then be used in standard statistical analyses (regression, correlation, ANOVA) without violating the assumption of independent observations.

**Additional features:**
- Covariance/correlation/regression between characters
- Correction for within-species sampling variation (when multiple individuals are measured)
- Log-likelihood computation under multivariate Brownian motion

**Source:** `contrast.c` (964 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `makecontrasts()` | ~400 | Compute independent contrasts at each internal node |
| `contbetween()` | ~300 | Between-species contrasts |
| `contwithin()` | ~250 | Within-species variation correction |
| `nuview()` | ~350 | Compute ancestral values at internal nodes |
| `writecontrasts()` | ~450 | Output contrasts |
| `regressions()` | ~500 | Compute regressions between characters |
| `logdet()` | ~550 | Log-determinant of covariance matrix |
| `invert()` | ~570 | Matrix inversion |
| `emiterate()` | ~650 | EM algorithm for maximum likelihood estimation of covariance components |
| `initcovars()` | ~600 | Initialize covariance matrices |
| `newcovars()` | ~630 | Update covariances in EM iteration |

**Library dependencies:** `cont.c`, `phylip.c`

---

## 11. Phylogenetic Invariants

### 11.1 Lake's and Cavender's Invariants

**Program:** Dnainvar

**Algorithm:** For exactly 4 species, computes phylogenetic invariants that can distinguish between the three possible unrooted tree topologies. Lake (1987) called his version "evolutionary parsimony." Cavender and Felsenstein (1987) developed a related approach. These invariants are linear functions of nucleotide pattern frequencies that are expected to be zero under one topology and nonzero under others.

**Source:** `dnainvar.c` (828 lines)
| Function | Line | Purpose |
|----------|------|---------|
| `maketree()` | ~550 | Main computation |
| `makeinv()` | ~350 | Compute invariant values |
| `getbasefreqs()` | in `seq.c` | Compute base frequencies |
| `getnums()` | ~200 | Count nucleotide pattern frequencies |
| `getpatterns()` | ~250 | Tabulate all 256 possible site patterns |

**Library dependencies:** `seq.c`, `phylip.c`

---

## 12. Core Infrastructure

### 12.1 Core Library (`phylip.c` / `phylip.h`)

The foundation of all PHYLIP programs. **3,207 lines.**

**Key functional areas:**

| Area | Key Functions | Purpose |
|------|--------------|---------|
| **I/O** | `openfile()`, `getstryng()`, `scan_eoln()`, `eoff()`, `eoln()` | File handling, input parsing |
| **Random numbers** | `randum()`, `randumize()`, `normrand()` | Pseudorandom number generation (linear congruential) |
| **Tree I/O** | `treeread()`, `treeout()` | Read/write Newick format trees |
| **Tree data structures** | `node`, `tree`, `pointarray` | Fundamental tree node types |
| **Memory** | `gnu()`, `chuck()`, `alloctree()`, `freetree()` | Node allocation/deallocation linked lists |
| **Rate categories** | `initlaguerrecat()`, `initgammacat()`, `lgr()` | Gamma distribution discretization |
| **User interaction** | `getoptions()`, `initseed()`, `initjumble()`, `initoutgroup()` | Menu-driven option handling |
| **Math** | `logfac()`, `glaguerre()`, `hermite()` | Special mathematical functions |
| **Bootstrap** | `initseed()`, `initjumble()` | Random seed management |

**Key data types (from `phylip.h`):**
```c
typedef struct node {        /* tree node */
  struct node *next, *back;  /* ring of nodes at fork, back pointer */
  long index;                /* node number */
  double v;                  /* branch length */
  double tyme;               /* node time (for clock models) */
  ...
} node;

typedef struct tree {
  node **nodep;              /* array of node pointers */
  node *root;
  double likelihood;
  ...
} tree;
```

---

### 12.2 Sequence Library (`seq.c` / `seq.h`)

Shared code for all molecular sequence programs. **4,155 lines.**

| Area | Key Functions | Purpose |
|------|--------------|---------|
| **Data input** | `inputdata()` | Read aligned sequences |
| **Frequencies** | `getbasefreqs()`, `empiricalfreqs()` | Compute/set nucleotide/amino acid frequencies |
| **Site handling** | `sitesort()`, `sitecombine()`, `sitescrunch()` | Compress sites by identical patterns (major speedup) |
| **Tip initialization** | `makevalues()`, `makevalues2()`, `alloctip()` | Set up conditional likelihoods at tips |
| **Fitch algorithm** | `fillin()`, `sumnsteps()` | Core parsimony operations |
| **Tree operations** | `add()`, `re_move()`, `reroot()` | Tree topology manipulation |
| **Memory** | `allocx()`, `prot_allocx()`, `alloctemp()` | Allocate likelihood arrays |
| **Traversals** | `preorder()`, `postorder()` | Tree traversal helpers |

---

### 12.3 Discrete Character Library (`disc.c` / `disc.h`)

For 0/1 character programs. **926 lines.**

| Function | Purpose |
|----------|---------|
| `inputdata()` | Read 0/1 character matrix |
| `inputancestors()` | Read ancestral state assignments |
| `add()`, `add2()`, `add3()` | Tree topology operations (different node types) |
| `re_move()`, `re_move2()`, `re_move3()` | Remove branch operations |
| `coordinates()` | Compute display coordinates |
| `treeout()` | Write Newick tree |
| `standev()` | Statistical test for tree comparison |
| `guesstates()` | Estimate ancestral states |

---

### 12.4 Multistate Discrete Library (`discrete.c` / `discrete.h`)

For multistate character programs (Pars). **3,109 lines.**

Same general structure as `disc.c` but handles up to 8 character states with packed bit representations.

---

### 12.5 Drawing Library (`draw.c` / `draw2.c` / `draw.h`)

For Drawgram and Drawtree. **3,318 + 1,527 = 4,845 lines total.**

| Area | Key Functions | Purpose |
|------|--------------|---------|
| **PostScript** | `postscript_header()`, `plot()` | PostScript output |
| **SVG** | SVG output functions | Scalable Vector Graphics |
| **Splines** | `splyne()`, `swoopspline()`, `curvespline()` | Curved branch drawing |
| **BMP** | `write_bmp_header()`, `write_full_pic()` | Bitmap output |
| **Coordinates** | `computeAngle()` | Geometric computations |
| **Setup** | `setupgraphics()` | Initialize graphics system |

---

## 13. Algorithm-to-Program Cross-Reference

| Algorithm | Programs | Search Strategy |
|-----------|----------|-----------------|
| **Maximum Likelihood (DNA)** | Dnaml, Dnamlk | Heuristic (stepwise addition + SPR) |
| **Maximum Likelihood (Protein)** | Proml, Promlk | Heuristic (stepwise addition + SPR) |
| **Maximum Likelihood (Restriction sites)** | Restml | Heuristic (stepwise addition + SPR) |
| **Maximum Likelihood (Gene freq/Continuous)** | Contml | Heuristic (stepwise addition + SPR) |
| **Wagner Parsimony (DNA)** | Dnapars, Dnapenny, Dnamove | Heuristic, Branch-and-bound, Interactive |
| **Wagner Parsimony (Protein)** | Protpars | Heuristic |
| **Wagner Parsimony (Discrete 0/1)** | Mix, Penny, Move | Heuristic, Branch-and-bound, Interactive |
| **Wagner Parsimony (Multistate)** | Pars | Heuristic |
| **Camin-Sokal Parsimony** | Mix, Penny, Move | Heuristic, Branch-and-bound, Interactive |
| **Dollo Parsimony** | Dollop, Dolpenny, Dolmove | Heuristic, Branch-and-bound, Interactive |
| **Polymorphism Parsimony** | Dollop, Dolpenny, Dolmove | Heuristic, Branch-and-bound, Interactive |
| **Compatibility (DNA)** | Dnacomp | Heuristic |
| **Compatibility (Discrete)** | Clique | Heuristic (clique finding) |
| **Fitch-Margoliash / Least Squares** | Fitch, Kitsch | Heuristic |
| **Minimum Evolution** | Fitch, Kitsch | Heuristic |
| **Neighbor-Joining** | Neighbor | Agglomerative |
| **UPGMA** | Neighbor | Agglomerative |
| **DNA Distance (JC/K2P/F84/LogDet)** | Dnadist | Direct computation |
| **Protein Distance (JTT/PAM/PMB/Kimura)** | Protdist | Direct computation |
| **Restriction Distance** | Restdist | Direct computation |
| **Genetic Distance (Nei/Cavalli-Sforza/Reynolds)** | Gendist | Direct computation |
| **Phylogenetic Invariants** | Dnainvar | Direct computation |
| **Independent Contrasts** | Contrast | Direct computation |
| **Bootstrap/Jackknife/Permutation** | Seqboot | Resampling |
| **Character Factoring** | Factor | Data transformation |
| **Consensus Trees (Majority-rule)** | Consense | Bipartition counting |
| **Tree Distances (RF/Branch Score)** | Treedist | Bipartition comparison |
| **Tree Plotting (Rooted)** | Drawgram | Coordinate assignment |
| **Tree Plotting (Unrooted)** | Drawtree | Equal-daylight algorithm |
| **Tree Editing** | Retree | Interactive |

---

## 14. Dependency Graph

```
                           phylip.c / phylip.h
                          /    |    |    |    \
                         /     |    |    |     \
                        /      |    |    |      \
                    seq.c   disc.c cont.c dist.c cons.c   draw.c/draw2.c   moves.c
                   / | \     / | \    |     |      |        |      |          |
                  /  |  \   /  |  \   |     |      |        |      |          |
  [DNA programs]    |  [Protein]  |   |     |      |        |      |    [Interactive]
  dnapars          |  protpars   |   |     |      |        |      |     dnamove
  dnacomp          |  proml     dollo.c  |   |    consense |      |     move
  dnapenny         |  promlk    |   | cont.c |   treedist  |      |     dolmove
  dnaml            |  protdist  |   |    |   |             |      |     retree
  dnamlk           |            |   | contml |           drawgram  |
  dnadist        discrete.c  dollop |   |    |           drawtree  |
  dnainvar         |        dolpenny| contrast               |
  restml          pars     dolmove  |   |                     |
  restdist         |          |     |   |                     |
  seqboot          |   wagner.c  gendist                     |
                   |     / | \                                |
                   |    /  |  \                            mlclock.c
                   |  mix penny                            printree.c
                   |  move                                    |
                   |                                       dnamlk
                   |                                       promlk
                factor
```

---

## 15. Reimplementation Priority Recommendations

Based on usage frequency, algorithmic importance, and code complexity:

### Tier 1: Core Infrastructure (implement first)
1. **`phylip.c`/`phylip.h`** -- All programs depend on this. Tree data structures, I/O, random numbers.
2. **`seq.c`/`seq.h`** -- Used by 16 programs. Sequence input, site compression, Fitch algorithm.
3. **Newick tree I/O** (in `phylip.c`) -- Reading and writing tree files.

### Tier 2: Most-Used Programs
4. **Dnaml** -- The most popular program. DNA maximum likelihood.
5. **Dnapars** -- DNA parsimony, widely used.
6. **Neighbor** -- Neighbor-Joining, very commonly used and algorithmically simple.
7. **Dnadist** -- DNA distance computation, feeds into Neighbor/Fitch/Kitsch.
8. **Seqboot** -- Bootstrap resampling, essential for confidence assessment.
9. **Consense** -- Consensus trees, pairs with Seqboot.
10. **Fitch** -- Distance-based tree estimation.

### Tier 3: Important Secondary Programs
11. **Proml** -- Protein ML (largest program, 3,246 lines).
12. **Protdist** -- Protein distances.
13. **Protpars** -- Protein parsimony.
14. **Dnamlk** -- DNA ML with clock.
15. **Promlk** -- Protein ML with clock.

### Tier 4: Specialized Programs
16. **Pars** -- Multistate parsimony.
17. **Kitsch** -- Distance with clock.
18. **Mix/Penny** -- Discrete character parsimony.
19. **Dollop/Dolpenny** -- Dollo parsimony.
20. **Contrast** -- Independent contrasts.
21. **Contml** -- Gene frequency ML.
22. **Gendist** -- Genetic distances.

### Tier 5: Visualization and Interactive (may not need Rust reimplementation)
23. **Drawgram/Drawtree** -- Tree visualization (consider using existing Rust plotting libraries).
24. **Retree** -- Interactive editor (consider TUI library).
25. **Dnamove/Move/Dolmove** -- Interactive parsimony (niche use).

### Tier 6: Niche/Legacy Programs
26. **Clique** -- Compatibility method (rarely used today).
27. **Dnacomp** -- DNA compatibility.
28. **Dnainvar** -- Invariants (4 species only).
29. **Factor** -- Character factoring.
30. **Restdist/Restml** -- Restriction site methods (declining use).
31. **Dnapenny/Dolpenny/Penny** -- Branch-and-bound (limited to small datasets).

---

## Shared Algorithm Patterns

Several algorithm patterns are replicated across many programs, making them good candidates for shared Rust traits/generics:

1. **Stepwise Addition + SPR Search:** Used by Dnaml, Proml, Dnapars, Protpars, Fitch, Mix, Dollop, Dnacomp, Contml, Restml. The pattern is:
   - Build initial 3-taxon tree
   - Add taxa one at a time, trying all insertion points
   - After all taxa added, do local SPR rearrangements
   - Optionally do global rearrangements
   - Optionally repeat with different random addition orders (jumble)

2. **Branch-and-Bound Search:** Used by Dnapenny, Penny, Dolpenny. Pattern:
   - Build partial tree with first few taxa
   - Recursively try adding remaining taxa at all positions
   - Prune search when partial tree score exceeds best known

3. **Fitch/Pruning Algorithm:** Used by all parsimony and ML programs. Pattern:
   - Postorder traversal computing state sets (parsimony) or conditional likelihoods (ML) at internal nodes
   - Preorder traversal for ancestral state reconstruction

4. **Distance Matrix I/O:** Used by Fitch, Kitsch, Neighbor. Shared in `dist.c`.

5. **Bipartition/Split Operations:** Used by Consense, Treedist. Shared in `cons.c`.

6. **Interactive Menu System:** All programs share a common pattern for user interaction, defined through phylip.c.

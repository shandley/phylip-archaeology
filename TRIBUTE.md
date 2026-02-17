# A Tribute to Joe Felsenstein and the PHYLIP Legacy

*An appreciation of one of the most consequential contributions in the
history of computational biology*

---

## The Man

Joseph "Joe" Felsenstein was born on May 9, 1942, in Philadelphia,
Pennsylvania. He earned his B.S. in Zoology (Honors) from the University of
Wisconsin-Madison in 1964, studying under the legendary population geneticist
James F. Crow. He completed his Ph.D. in Zoology at the University of Chicago
in 1968 under Richard Lewontin, one of the founders of molecular population
genetics. After a postdoctoral fellowship at the Institute of Animal Genetics
in Edinburgh, he joined the University of Washington in 1967, where he would
spend his entire career — rising from Assistant Professor of Genetics to
Professor of Genome Sciences, with adjunct appointments in Statistics, Zoology,
Biology, and Computer Science.

In a remarkable family connection, Joe's younger brother **Lee Felsenstein** was
a pioneering personal computer designer — creator of the Sol-20 and the Osborne
1 (the first mass-produced portable computer), founding member of the Homebrew
Computer Club, and co-founder of Community Memory, the first public computerized
bulletin board system. Early versions of PHYLIP were developed on machines
designed by Joe's own brother, a small piece of history that ties the birth of
personal computing to the birth of computational phylogenetics.

Joe Felsenstein's honors tell part of his story: elected to the American Academy
of Arts and Sciences (1992), President of the Society for the Study of Evolution
(1993), the Sewall Wright Award (1993), elected to the National Academy of
Sciences (1999), the Weldon Memorial Prize from Oxford (2000), the Darwin-Wallace
Medal from the Linnean Society of London (2009), the International Prize for
Biology from the Japan Society for the Promotion of Science (2013), and the
Mendel Medal from the U.K. Genetics Society (2026). A moth species, *Ufeus
felsensteini*, was named in his honor. His Google Scholar profile shows over
**162,000 citations**. But the deeper story is in the tools and ideas he gave
to the world.

## 1980: PHYLIP is Born

In October 1980, Joe Felsenstein released version 1.0 of PHYLIP — the PHYLogeny
Inference Package. To appreciate what this meant, consider the world it entered.

The IBM PC would not be released for another year. The World Wide Web would not
exist for another decade. GenBank had not yet opened to the public. Most
biologists had never used a computer for anything. DNA sequencing was a
painstaking manual process involving radioactive labeling and gel
electrophoresis. The entire GenBank database, when it launched in 1982, would
contain fewer sequences than a single modern sequencing run produces in minutes.

In this environment, Felsenstein wrote a package of programs — initially in
Pascal, later rewritten in C — that could infer evolutionary trees from
molecular data. The first version was distributed on magnetic tapes, physically
mailed to researchers who requested them. This was open-source software before
the term "open source" existed, free scientific software in an era when the
concept was still radical.

Version 1.0 focused on parsimony and distance methods. The earliest surviving
version in archives is 1.7, and a reading of its "bugs" file reveals that only
minor fixes separated it from the original release — evidence that the initial
design was remarkably solid.

## The Algorithms That Changed a Field

Felsenstein's contributions were not merely practical — they were foundational
to the mathematics of evolutionary biology.

### The Felsenstein Zone (1978)

Before any of his most famous work, Felsenstein made a contribution that would
reshape the entire debate over phylogenetic methods. In 1978, he published
"Cases in which parsimony or compatibility methods will be positively
misleading" in *Systematic Zoology*. The paper demonstrated that maximum
parsimony — then the dominant method — could be **statistically inconsistent**:
it could converge on the *wrong* tree topology even with infinite data.

This occurs in what is now called the **Felsenstein zone**: tree topologies with
short internal branches and long external branches, where convergent mutations
on long branches are mistaken for shared ancestry. The phenomenon, known as
**long branch attraction**, remains one of the most important concepts in
phylogenetics. This paper was a key motivator for the adoption of model-based
methods (maximum likelihood, Bayesian inference) that could handle such cases
correctly.

### The Pruning Algorithm (1981)

In 1981, Felsenstein published "Evolutionary trees from DNA sequences: a maximum
likelihood approach" in the *Journal of Molecular Evolution*. This paper
introduced what is now universally known as **Felsenstein's pruning algorithm**
— a dynamic programming method for computing the likelihood of observing
sequence data given a phylogenetic tree.

The algorithm traverses a tree from tips to root, computing conditional
likelihoods at each internal node. Its elegance lies in its efficiency: it
reduces what would be an exponential computation to one that is linear in the
number of taxa. The paper has accumulated over **9,000 citations** and is the
second most cited paper in the *Journal of Molecular Evolution*. This single
algorithm is the mathematical backbone of every modern maximum likelihood
phylogenetics program — RAxML, IQ-TREE, PhyML, GARLI, and many more. Without
it, genome-scale phylogenetics would be computationally intractable.

(Felsenstein had actually proposed the approach earlier, in two 1973 papers on
maximum likelihood for discrete and continuous characters, but the 1981 DNA
sequence paper is the one that made it practical and widely adopted.)

### The Phylogenetic Bootstrap (1985)

In 1985, Felsenstein published "Confidence limits on phylogenies: an approach
using the bootstrap" in *Evolution*. Before this paper, phylogenetic trees were
presented as point estimates with no measure of uncertainty. Felsenstein adapted
Efron's bootstrap — a statistical resampling technique — to the problem of
phylogenetic confidence.

The method is breathtakingly simple: resample columns from a sequence alignment
with replacement, reconstruct the tree, repeat hundreds or thousands of times,
and report how often each branch appears. This gave the field its first widely
adopted measure of phylogenetic support.

The paper has accumulated over **32,000 citations**, placing it in the **top 100
most cited scientific papers of all time** — across all fields of science. It
still receives approximately 2,000 citations per year. Felsenstein himself has
noted it is "the most cited paper ever produced at my university." Today,
virtually every published phylogenetic tree includes bootstrap values.

### Phylogenetically Independent Contrasts (1985)

Also in 1985, Felsenstein published "Phylogenies and the comparative method" in
*The American Naturalist*. This paper solved a fundamental statistical problem:
when comparing traits across species, the data points are not independent
because species share evolutionary history. Treating them as independent
overstates statistical significance.

Felsenstein's solution — phylogenetically independent contrasts (PIC) —
provided a mathematically rigorous way to account for shared ancestry in
comparative analyses. The paper has accumulated over 9,000 citations and is the
second most cited article in the history of *The American Naturalist*. It
spawned an entire subfield of phylogenetic comparative methods.

### The Newick Format (1986)

On June 26, 1986, during the Society for the Study of Evolution meetings in
Durham, New Hampshire, Felsenstein convened an informal committee to standardize
a format for representing phylogenetic trees as text strings. The committee —
James Archie, William H.E. Day, Wayne Maddison, Christopher Meacham, F. James
Rohlf, David Swofford, and Felsenstein — held its final session at Newick's
Lobster House in Dover, New Hampshire.

The format they adopted, based on Christopher Meacham's 1984 work for PHYLIP's
tree-drawing programs, uses nested parentheses to represent tree topology — a
correspondence first noted by mathematician Arthur Cayley in 1857. They named it
the Newick format after the restaurant where they finalized it over lobster.

The Newick format was never formally published. Yet it became — and remains —
the de facto universal standard for representing phylogenetic trees. Every
phylogenetics program reads and writes it. It is one of the most successful
informal standards in the history of science.

## PHYLIP Grows

PHYLIP evolved steadily over decades:

- **Version 1.0** (October 1980) — Initial release. Parsimony and distance
  methods. Distributed on magnetic tapes.
- **Version 3.0** (1987) — Major expansion. Increased program count, improved
  portability, initial C implementations alongside Pascal.
- **Version 3.2** (1989) — Added protein sequence support and maximum
  likelihood methods. The version first cited formally in *Cladistics*.
- **Version 3.3** (1993) — Source code rewritten from Pascal to C, greatly
  expanding portability across operating systems.
- **Version 3.5** — Widely used stable release. Precompiled executables for
  DOS, Windows, PowerMac, and Unix.
- **Version 3.6** — Documentation as web pages with hyperlinks. Faster DNA
  and restriction site likelihood programs. Protein likelihood programs.
  Gamma-distributed rate variation.
- **Version 3.696** — First open-source licensed version, allowing
  redistribution alongside other software.
- **Version 3.698** — Latest version as of this writing. 64-bit Windows
  support, consensus tree bug fix.

The package grew to encompass **65 portable programs** written in C, covering:

- **Parsimony methods**: DNA, protein, and discrete character parsimony
- **Distance methods**: Neighbor-joining, UPGMA, Fitch-Margoliash
- **Maximum likelihood**: DNA and protein likelihood with multiple
  substitution models
- **Bootstrap resampling**: Via seqboot
- **Consensus trees**: Strict, majority-rule, and extended consensus
- **Tree drawing**: Drawgram and Drawtree
- **Sequence simulation**: For testing methods
- **Distance computation**: Multiple evolutionary models

By the 2010s, PHYLIP had accumulated over **30,000 registered users** and more
than **25,000 citations** summed across all versions. Even in 2013, well after
newer tools had gained prominence, it was cited over 1,000 times in a single
year.

## The Software Catalog: A Proto-Package Manager

Perhaps as consequential as PHYLIP itself was Felsenstein's decision to maintain
a comprehensive catalog of *all* phylogenetics software — not just his own.

The catalog, originally hosted at `evolution.genetics.washington.edu`, grew to
document **392 phylogeny packages** and **54 free web servers**. Each entry
included the tool's name, author, capabilities, platform, and URL. Felsenstein
maintained it by hand for decades, adding new tools as they appeared and
(where possible) noting when they went dormant.

This was, in effect, a **hand-curated package manager for an entire scientific
discipline**. It existed years before SourceForge (1999), years before
Bioconductor (2001), years before GitHub (2008), years before conda-forge.
For many researchers in the 1990s and 2000s, Felsenstein's catalog page was
the first — and often only — place to discover that software existed for a
particular type of phylogenetic analysis.

The catalog was democratic: it listed tools without quality judgments, trusting
researchers to evaluate software for themselves. It was comprehensive: from
general-purpose packages to highly specialized programs for host-parasite
coevolution or genome rearrangement phylogeny. And it was persistent:
Felsenstein kept it updated for roughly two decades.

In August 2023, the pages were migrated from their original university server
to GitHub Pages at `phylipweb.github.io/phylip/`. The catalog survives, though
many of the tools it points to do not — their URLs lead to defunct servers,
abandoned university pages, and dead FTP archives. Each dead link represents a
tool that once helped researchers answer evolutionary questions, now at risk of
being lost to history.

## The Textbook

In 2004, Felsenstein published *Inferring Phylogenies* with Sinauer Associates
— 664 pages covering the mathematical foundations of every major phylogenetic
method: parsimony, distance, likelihood, Bayesian inference, the bootstrap,
coalescent trees, comparative methods, and more, with approximately 1,000
references.

The reviews were extraordinary. David Penny called it "the book we have been
waiting for — occasionally a book is a classic by the time it is published and
this is it." A.J. Drummond in *Heredity* called it "an instant classic" and
Felsenstein "the father of statistical phylogenetics." Fredrik Ronquist wrote
in *Science* that "the publication of *Inferring Phylogenies* is a milestone
for evolutionary biology."

*Inferring Phylogenies* is not just a textbook. It is a map of an entire
intellectual landscape, written by the person who helped create much of that
landscape. It remains essential reading for anyone entering phylogenetics.

## Why This Matters

Joe Felsenstein's contributions form a throughline from the earliest days of
computational biology to the genomics era. The pruning algorithm makes modern
phylogenomics possible. The bootstrap gives us a language for uncertainty.
Independent contrasts gave comparative biology statistical rigor. The Newick
format lets every program speak the same language. And the software catalog
connected researchers with tools they didn't know existed.

PHYLIP itself — freely distributed since 1980, written in portable C,
documented meticulously — embodied values that the open-source movement would
later codify: that scientific software should be free, accessible, and
reproducible. Felsenstein practiced these values before they had a name.

The catalog embodied a different value: that a senior scientist with standing
in a field has a responsibility to help others navigate it. Maintaining a list
of 392 software packages is not glamorous work. It does not result in Nature
papers or grant funding. But it may be one of the highest-impact contributions
a scientist can make — a quiet act of service that multiplied the effectiveness
of thousands of researchers around the world.

## This Project

The PHYLIP Archaeology project exists because these contributions deserve to be
preserved, studied, and carried forward. The original C code contains
algorithmic wisdom that should not be lost to bit rot. The catalog documents an
ecosystem that should not fade into dead links. And the ideas — the pruning
algorithm, the bootstrap, independent contrasts — deserve clean, modern
implementations that a new generation of researchers can learn from.

This is archaeological work in the best sense: careful excavation, faithful
documentation, and deep respect for what came before.

---

*Thank you, Joe, for building the foundations.*

---

## References

- Felsenstein, J. (1973). Maximum likelihood and minimum-steps methods for
  estimating evolutionary trees from data on discrete characters. *Systematic
  Biology*, 22, 240-249.
- Felsenstein, J. (1978). Cases in which parsimony or compatibility methods
  will be positively misleading. *Systematic Zoology*, 27(4), 401-410.
- Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum
  likelihood approach. *Journal of Molecular Evolution*, 17, 368-376.
- Felsenstein, J. (1985). Confidence limits on phylogenies: an approach using
  the bootstrap. *Evolution*, 39, 783-791.
- Felsenstein, J. (1985). Phylogenies and the comparative method. *The American
  Naturalist*, 125(1), 1-15.
- Felsenstein, J. (1989). PHYLIP — Phylogeny Inference Package (Version 3.2).
  *Cladistics*, 5, 164-166.
- Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates,
  Sunderland, Massachusetts.
- Carvalho, P., Diniz-Filho, J.A.F., & Bini, L.M. (2005). The impact of
  Felsenstein's "Phylogenies and the comparative method" on evolutionary
  biology. *Scientometrics*, 62(1).
- Revell, L.J. & Chamberlain, S.A. (2014). Rphylip: an R interface for PHYLIP.
  *Methods in Ecology and Evolution*, 5(9), 976-981.
- Van Noorden, R., Maher, B., & Nuzzo, R. (2014). The top 100 papers. *Nature*,
  514, 550-553.

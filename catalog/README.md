# Software Catalog Preservation

This directory contains the structured preservation of Joe Felsenstein's
[phylogenetics software catalog](https://phylipweb.github.io/phylip/software.html),
which documents **392 phylogeny packages** and **54 free web servers**.

## Contents

- `schema.json` — JSON Schema defining the structure of each tool entry
- `tools.json` — Structured data for all 407 cataloged tools
- `snapshots/` — Archived HTML pages from the original catalog
- `parse_catalog.py` — Extracts tool names, categories, platforms from index page
- `enrich_catalog.py` — Extracts URLs and descriptions from detail pages, checks liveness
- `fix_matches.py` — Resolves anchor-to-name mismatches for additional coverage

## Methodology

1. **Snapshot**: Archive the original HTML catalog pages
2. **Parse**: Extract structured tool data into `tools.json`
3. **Verify**: Check the status of every URL (alive/dead/archived)
4. **Archive**: Search the Wayback Machine for archived copies of dead links
5. **Enrich**: Add algorithm classifications, influence lineage, and notes

## Tool Status Categories

| Status | Meaning |
|--------|---------|
| `maintained` | Active development or recent updates |
| `archived` | Code available but no recent activity |
| `dormant` | Site exists but appears abandoned |
| `dead` | URL unreachable, no archive found |
| `unknown` | Not yet assessed |

## Current Statistics (2026-02-17)

| Metric | Count |
|--------|-------|
| Total tools | 407 |
| Assessed (URL checked) | 367 (90.2%) |
| Archived (code available) | 196 (48.2%) |
| Dormant (site up, inactive) | 137 (33.7%) |
| Dead (unreachable, no archive) | 34 (8.4%) |
| Unknown (no URL in catalog) | 40 (9.8%) |
| With Wayback Machine archive | 203 |
| With description | 381 (93.6%) |

## Snapshots

The `snapshots/` directory contains HTML pages downloaded from
`phylipweb.github.io/phylip/` on 2026-02-16:

- `index.html` — PHYLIP home page
- `software.html` — The full software catalog index (407 tools)
- `software.pars.html` — Parsimony program entries (43 tools)
- `software.dist.html` — Distance program entries (68 tools)
- `software.etc1.html` — Maximum likelihood and Bayesian entries (133 tools)
- `software.etc2.html` — Visualization, simulation, and other entries (148 tools)
- `software.serv.html` — Web server entries (43 servers)
- `software.xref.html` — Cross-reference table
- `phylip.html` — PHYLIP programs and documentation
- `faq.html` — Frequently asked questions
- `general.html` — General information

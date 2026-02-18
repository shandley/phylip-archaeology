#!/usr/bin/env python3
"""Compute statistics from enriched catalog data for the Nature Methods paper.

Reads tools_enriched.json, computes summary stats, temporal distributions,
survival rates, category evolution, and writes narrative_stats.json.
"""

import json
from collections import Counter
from pathlib import Path

ENRICHED_FILE = Path(__file__).parent / "tools_enriched.json"
STATS_FILE = Path(__file__).parent / "stats" / "narrative_stats.json"


def main():
    tools = json.loads(ENRICHED_FILE.read_text())
    total = len(tools)

    # ── Summary statistics ──────────────────────────────────
    status_counts = Counter(t["status"] for t in tools)
    with_url = sum(1 for t in tools if t.get("url_original"))
    with_archive = sum(1 for t in tools if t.get("url_archive"))
    with_year = sum(1 for t in tools if t.get("first_release"))
    with_author = sum(1 for t in tools if t.get("author"))
    with_lang = sum(1 for t in tools if t.get("languages"))

    dead_no_archive = sum(1 for t in tools
                          if t["status"] == "dead" and not t.get("url_archive"))

    print(f"PHYLIP Software Catalog Analysis")
    print(f"{'='*60}")
    print(f"Total tools cataloged: {total}")
    print(f"  Archived (code available, not active): {status_counts['archived']}")
    print(f"  Dormant (site up, appears abandoned):  {status_counts['dormant']}")
    print(f"  Dead (unreachable):                    {status_counts['dead']}")
    print(f"  Unknown (no URL to check):             {status_counts['unknown']}")
    print()
    print(f"URLs: {with_url} tools have original URLs")
    print(f"  Wayback Machine archived: {with_archive}")
    print(f"  Dead with no archive:     {dead_no_archive}")
    print()
    print(f"Metadata coverage:")
    print(f"  Release year:  {with_year}/{total} ({with_year/total*100:.0f}%)")
    print(f"  Author:        {with_author}/{total} ({with_author/total*100:.0f}%)")
    print(f"  Language:      {with_lang}/{total} ({with_lang/total*100:.0f}%)")
    print()

    # ── Temporal distribution ───────────────────────────────
    bins_5yr: Counter = Counter()
    for t in tools:
        y = t.get("first_release")
        if y:
            bin_start = (y // 5) * 5
            bins_5yr[bin_start] += 1

    print(f"Temporal Distribution (5-year bins):")
    print(f"  Tools with known year: {with_year}")
    peak_bin = max(bins_5yr, key=bins_5yr.get) if bins_5yr else None
    for b in sorted(bins_5yr):
        bar = '#' * bins_5yr[b]
        marker = " <-- peak" if b == peak_bin else ""
        print(f"  {b}-{b+4}: {bar} ({bins_5yr[b]}){marker}")
    print()

    # ── Loss rate by era ────────────────────────────────────
    print(f"Survival Analysis by Era:")
    era_status: dict[int, Counter] = {}
    for t in tools:
        y = t.get("first_release")
        if y:
            era = (y // 5) * 5
            if era not in era_status:
                era_status[era] = Counter()
            era_status[era][t["status"]] += 1

    for era in sorted(era_status):
        counts = era_status[era]
        era_total = sum(counts.values())
        dead = counts.get("dead", 0)
        dormant = counts.get("dormant", 0)
        archived = counts.get("archived", 0)
        alive_pct = (era_total - dead) / era_total * 100
        print(f"  {era}-{era+4}: {era_total:3d} tools | "
              f"dead: {dead} ({dead/era_total*100:.0f}%) | "
              f"dormant: {dormant} | archived: {archived} | "
              f"survival: {alive_pct:.0f}%")
    print()

    # ── Category analysis ───────────────────────────────────
    cat_counts = Counter()
    for t in tools:
        for c in t.get("categories", []):
            cat_counts[c] += 1

    print(f"Category Distribution (top 15):")
    for cat, count in cat_counts.most_common(15):
        bar = '#' * (count // 2)
        print(f"  {cat:25s}: {bar} ({count})")
    print()

    # ── Category evolution by era ───────────────────────────
    top_cats = ["distance", "parsimony", "likelihood", "bayesian",
                "tree-visualization", "bootstrap", "alignment"]
    era_cats: dict[int, Counter] = {}
    for t in tools:
        y = t.get("first_release")
        if y:
            era = (y // 5) * 5
            if era not in era_cats:
                era_cats[era] = Counter()
            for c in t.get("categories", []):
                if c in top_cats:
                    era_cats[era][c] += 1

    print(f"Category Evolution by Era:")
    header = f"  {'Era':>9s}"
    for c in top_cats:
        header += f" {c[:8]:>8s}"
    print(header)
    for era in sorted(era_cats):
        row = f"  {era}-{era+4}"
        for c in top_cats:
            row += f" {era_cats[era].get(c, 0):8d}"
        print(row)
    print()

    # ── Language distribution ───────────────────────────────
    lang_counts = Counter()
    for t in tools:
        for l in t.get("languages", []):
            lang_counts[l] += 1

    print(f"Programming Language Distribution:")
    for lang, count in lang_counts.most_common():
        bar = '#' * (count // 2)
        print(f"  {lang:15s}: {bar} ({count})")
    print()

    # ── Language evolution by era ───────────────────────────
    top_langs = ["C", "Java", "C++", "Perl", "Python", "R"]
    era_langs: dict[int, Counter] = {}
    for t in tools:
        y = t.get("first_release")
        if y:
            era = (y // 5) * 5
            if era not in era_langs:
                era_langs[era] = Counter()
            for l in t.get("languages", []):
                if l in top_langs:
                    era_langs[era][l] += 1

    print(f"Language Evolution by Era:")
    header = f"  {'Era':>9s}"
    for l in top_langs:
        header += f" {l:>6s}"
    print(header)
    for era in sorted(era_langs):
        row = f"  {era}-{era+4}"
        for l in top_langs:
            row += f" {era_langs[era].get(l, 0):6d}"
        print(row)
    print()

    # ── Platform distribution ───────────────────────────────
    plat_counts = Counter()
    for t in tools:
        for p in t.get("platforms", []):
            plat_counts[p] += 1

    print(f"Platform Distribution:")
    for plat, count in plat_counts.most_common():
        print(f"  {plat:15s}: {'#' * (count // 3)} ({count})")
    print()

    # ── Build narrative statistics ──────────────────────────
    stats = {
        "total_tools": total,
        "status_archived": status_counts["archived"],
        "status_dormant": status_counts["dormant"],
        "status_dead": status_counts["dead"],
        "status_unknown": status_counts["unknown"],
        "tools_with_url": with_url,
        "wayback_archived": with_archive,
        "dead_no_archive": dead_no_archive,
        "tools_with_year": with_year,
        "tools_with_author": with_author,
        "tools_with_language": with_lang,
        "peak_era": f"{peak_bin}-{peak_bin+4}" if peak_bin else None,
        "peak_era_count": bins_5yr[peak_bin] if peak_bin else 0,
        "temporal_bins_5yr": {f"{k}-{k+4}": v for k, v in sorted(bins_5yr.items())},
        "category_counts": dict(cat_counts.most_common()),
        "language_counts": dict(lang_counts.most_common()),
        "platform_counts": dict(plat_counts.most_common()),
        "era_status": {
            f"{k}-{k+4}": dict(v) for k, v in sorted(era_status.items())
        },
        "era_categories": {
            f"{k}-{k+4}": dict(v) for k, v in sorted(era_cats.items())
        },
        "era_languages": {
            f"{k}-{k+4}": dict(v) for k, v in sorted(era_langs.items())
        },
        # Pre-formatted sentences for paper text
        "sentences": {
            "total": f"Felsenstein's catalog documents {total} phylogenetics software tools spanning over four decades of development.",
            "status": f"Of {total} tools, {status_counts['dead']} ({status_counts['dead']/total*100:.1f}%) are completely unreachable, {status_counts['dormant']} ({status_counts['dormant']/total*100:.1f}%) appear dormant, and {status_counts['archived']} ({status_counts['archived']/total*100:.1f}%) are archived but no longer actively maintained.",
            "preservation": f"The Wayback Machine has preserved snapshots of {with_archive} tools, but {dead_no_archive} tools have no archived copy and may be permanently lost.",
            "peak": f"Tool publication peaked in the {peak_bin}-{peak_bin+4} period with {bins_5yr[peak_bin]} tools, coinciding with the genomics revolution and the rise of Bayesian phylogenetics." if peak_bin else "",
            "languages": f"C dominates the catalog ({lang_counts.get('C', 0)} tools), followed by Java ({lang_counts.get('Java', 0)}), reflecting the computational demands of phylogenetic inference and the era in which most tools were written.",
        },
    }

    STATS_FILE.parent.mkdir(parents=True, exist_ok=True)
    STATS_FILE.write_text(json.dumps(stats, indent=2))
    print(f"Narrative statistics written to {STATS_FILE}")


if __name__ == "__main__":
    main()

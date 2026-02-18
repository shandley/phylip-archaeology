#!/usr/bin/env python3
"""Extract structured metadata (author, year, language) from PHYLIP catalog.

Reads full HTML detail pages from catalog/snapshots/ to get untruncated text,
then extracts metadata via regex, merges with curated_years.json, and writes
tools_enriched.json. Original tools.json is not modified.
"""

import json
import re
from html import unescape
from pathlib import Path

CATALOG_DIR = Path(__file__).parent.parent
TOOLS_FILE = CATALOG_DIR / "tools.json"
SNAPSHOT_DIR = CATALOG_DIR / "snapshots"
CURATED_FILE = Path(__file__).parent / "curated_years.json"
OUTPUT_FILE = Path(__file__).parent / "tools_enriched.json"

DETAIL_PAGES = [
    "software.pars.html",
    "software.dist.html",
    "software.etc1.html",
    "software.etc2.html",
    "software.serv.html",
]


def strip_html(text: str) -> str:
    """Remove HTML tags and normalize whitespace."""
    text = re.sub(r'<[^>]+>', ' ', text)
    text = unescape(text)
    text = re.sub(r'\s+', ' ', text)
    text = text.replace('(at)', '@').replace(' @ ', '@')
    return text.strip()


def extract_full_entries() -> dict[str, str]:
    """Extract full-length text entries from HTML detail pages.

    Returns dict: anchor_id (lowercased) -> full plain text
    """
    skip_anchors = {'top', 'bottom', 'methods', 'systems', 'datatypes',
                    'recent', 'new', 'news', 'changes', 'otherlists',
                    'ftpservers', 'bayesian', 'phylip-servers'}

    anchor_pattern = re.compile(r'<A\s+NAME="([^"]+)"[^>]*>', re.IGNORECASE)
    entries: dict[str, str] = {}

    for page in DETAIL_PAGES:
        filepath = SNAPSHOT_DIR / page
        if not filepath.exists():
            continue
        html = filepath.read_text(encoding="utf-8", errors="replace")

        all_anchors = [(m.start(), m.group(1)) for m in anchor_pattern.finditer(html)]

        for i, (pos, anchor_id) in enumerate(all_anchors):
            if anchor_id.lower() in skip_anchors:
                continue
            end = all_anchors[i + 1][0] if i + 1 < len(all_anchors) else len(html)
            full_text = strip_html(html[pos:end])

            key = anchor_id.lower()
            # Keep longer version if duplicate
            if key not in entries or len(full_text) > len(entries[key]):
                entries[key] = full_text

    return entries


def build_tool_to_fulltext(tools: list[dict], entries: dict[str, str]) -> dict[str, str]:
    """Map tool names to their full-text entries."""
    mapping: dict[str, str] = {}

    for tool in tools:
        name = tool["name"]
        key = name.lower()

        # Try direct match
        if key in entries:
            mapping[name] = entries[key]
            continue

        # Try variations
        for variant in [key.replace(' ', ''), key.replace('-', ''),
                        key.replace('_', ''), key.replace('*', '')]:
            if variant in entries:
                mapping[name] = entries[variant]
                break

    return mapping


def extract_year(text: str, name: str, curated: dict) -> tuple[int | None, str]:
    """Extract first release year. Returns (year, method)."""
    if name in curated:
        return curated[name]["first_release"], "curated"

    if not text:
        return None, "no_text"

    # "released X in YYYY"
    m = re.search(r'released\b.{0,60}?\b((?:19|20)\d{2})\b', text, re.IGNORECASE)
    if m:
        y = int(m.group(1))
        if 1970 <= y <= 2025:
            return y, "released_in"

    # "version X.X (YYYY)" or "version X.X, YYYY"
    m = re.search(r'version\s+[\d.]+\w*[\s,]*[\(,]\s*((?:19|20)\d{2})\b', text, re.IGNORECASE)
    if m:
        y = int(m.group(1))
        if 1970 <= y <= 2025:
            return y, "version_year"

    # "in YYYY" within first 800 chars
    m = re.search(r'\bin\s+((?:19|20)\d{2})\b', text[:800])
    if m:
        y = int(m.group(1))
        if 1980 <= y <= 2025:
            return y, "in_year"

    # Citation years: "(YYYY)" or ", YYYY."
    m = re.search(r'[,(]\s*((?:19|20)\d{2})\s*[.)]\s', text)
    if m:
        y = int(m.group(1))
        if 1970 <= y <= 2025:
            return y, "citation_year"

    # Earliest year mention in full text
    years = [int(y) for y in re.findall(r'\b((?:19|20)\d{2})\b', text)
             if 1970 <= int(y) <= 2025]
    if years:
        return min(years), "earliest_mention"

    return None, "not_found"


def extract_authors(text: str) -> str | None:
    """Extract primary author(s) from entry text."""
    if not text:
        return None

    # Grab text before first "of the" or email parenthesis
    m = re.match(r'^(.+?)\s+(?:of\s+the|of\s+[A-Z]|\()', text)
    if m:
        raw = m.group(1).strip().rstrip(',')
        if re.search(r'[A-Z][a-z]+', raw) and 5 < len(raw) < 300:
            raw = re.sub(r'\s*,?\s*(?:and|&)\s*$', '', raw)
            return raw

    return None


def extract_institution(text: str) -> str | None:
    """Extract institution from entry text."""
    if not text:
        return None

    # "of the Department/Institute/... University/..."
    m = re.search(
        r'(?:of\s+the|at\s+the)\s+((?:Department|School|Institute?|Laboratory|Center|Centre|'
        r'Faculty|Division|College|Museum|Lehrstuhl|Institut)[\w\s,\-\u00C0-\u024F]+?'
        r'(?:University|Institut\w*|College|Museum|Foundation|Center|Centre|School|Lab\w*)'
        r'[^(]*?)\s*(?:\(|\.|\bhas\b|\bhave\b)',
        text, re.IGNORECASE,
    )
    if m:
        inst = m.group(1).strip().rstrip(',')
        if 10 < len(inst) < 300:
            return inst

    # "of the University of X"
    m = re.search(
        r'(?:of\s+the|at\s+the|at)\s+(University\s+of\s+[\w\s,\-\u00C0-\u024F]+?)(?:\s*\(|\.|,\s*(?:has|have))',
        text, re.IGNORECASE,
    )
    if m:
        inst = m.group(1).strip().rstrip(',')
        if 10 < len(inst) < 200:
            return inst

    return None


LANG_PATTERNS: list[tuple[str, re.Pattern]] = [
    ("C++", re.compile(r'\bC\+\+\b', re.IGNORECASE)),
    ("Java", re.compile(r'\bJava\b(?!\s*Script)', re.IGNORECASE)),
    ("JavaScript", re.compile(r'\bJavaScript\b', re.IGNORECASE)),
    ("Python", re.compile(r'\bPython\b', re.IGNORECASE)),
    ("Perl", re.compile(r'\bPerl\b', re.IGNORECASE)),
    ("R", re.compile(r'\bR\s+package\b|\bCRAN\b|\bwritten\s+in\s+R\b|\bR\s+library\b|\bR\s+functions?\b', re.IGNORECASE)),
    ("Fortran", re.compile(r'\bFORTRAN\b|\bFortran\b')),
    ("Pascal", re.compile(r'\bPascal\b|\bDelphi\b', re.IGNORECASE)),
    ("C", re.compile(r'\bwritten\s+in\s+C\b|\bANSI\s+C\b|\b[Cc]\s+program\b|\bC\s+source\b|\bC\s+code\b|\bportable\s+C\b|\bstandard\s+C\b|\bGNU\s+C\b')),
    ("Ruby", re.compile(r'\bRuby\b', re.IGNORECASE)),
    ("Matlab", re.compile(r'\bMATLAB\b|\bMatlab\b', re.IGNORECASE)),
    ("Visual Basic", re.compile(r'\bVisual\s+Basic\b', re.IGNORECASE)),
]


LANG_NORMALIZE = {
    "java": "Java", "c++": "C++", "c": "C", "python": "Python",
    "perl": "Perl", "r": "R", "fortran": "Fortran", "pascal": "Pascal",
    "ruby": "Ruby", "matlab": "Matlab", "javascript": "JavaScript",
    "visual basic": "Visual Basic",
}


def extract_languages(text: str, existing: list[str] | None) -> list[str]:
    """Extract programming languages mentioned in text."""
    # Normalize existing entries
    langs = set()
    for l in (existing or []):
        langs.add(LANG_NORMALIZE.get(l.lower(), l))
    if not text:
        return sorted(langs)
    for lang_name, pattern in LANG_PATTERNS:
        if pattern.search(text):
            langs.add(lang_name)
    return sorted(langs)


def main():
    tools = json.loads(TOOLS_FILE.read_text())
    curated = json.loads(CURATED_FILE.read_text())

    # Extract full-length text from HTML detail pages
    print("Extracting full entries from HTML snapshots...")
    entries = extract_full_entries()
    print(f"  Found {len(entries)} entries across {len(DETAIL_PAGES)} detail pages")

    # Map tool names to full text
    fulltext = build_tool_to_fulltext(tools, entries)
    print(f"  Matched {len(fulltext)}/{len(tools)} tools to full-text entries")
    print()

    year_methods: dict[str, int] = {}
    year_count = 0
    author_count = 0
    lang_count = 0
    inst_count = 0

    for tool in tools:
        name = tool["name"]
        # Use full text from HTML if available, fall back to description
        text = fulltext.get(name, tool.get("description", ""))

        # Also store the full description if we got a longer one
        if name in fulltext and len(fulltext[name]) > len(tool.get("description", "")):
            tool["description_full"] = fulltext[name]

        # Year
        year, method = extract_year(text, name, curated)
        if year is not None:
            tool["first_release"] = year
            tool["year_extraction_method"] = method
            year_count += 1
            year_methods[method] = year_methods.get(method, 0) + 1

        # Author
        author = extract_authors(text)
        if author:
            tool["author"] = author
            author_count += 1

        # Institution
        inst = extract_institution(text)
        if inst:
            tool["institution"] = inst
            inst_count += 1

        # Languages
        existing_langs = tool.get("languages")
        langs = extract_languages(text, existing_langs)
        if langs:
            tool["languages"] = langs
            lang_count += 1

    # Write enriched output
    OUTPUT_FILE.write_text(json.dumps(tools, indent=2, ensure_ascii=False))

    # Report
    total = len(tools)
    print(f"Metadata Extraction Report")
    print(f"{'='*50}")
    print(f"Total tools: {total}")
    print()
    print(f"Years extracted: {year_count}/{total} ({year_count/total*100:.1f}%)")
    for method, count in sorted(year_methods.items(), key=lambda x: -x[1]):
        print(f"  {method}: {count}")
    print()
    print(f"Authors extracted: {author_count}/{total} ({author_count/total*100:.1f}%)")
    print(f"Institutions extracted: {inst_count}/{total} ({inst_count/total*100:.1f}%)")
    print(f"Languages extracted: {lang_count}/{total} ({lang_count/total*100:.1f}%)")
    print()

    # Show year distribution preview
    if year_count:
        from collections import Counter
        decades = Counter()
        for tool in tools:
            y = tool.get("first_release")
            if y:
                decade = (y // 5) * 5
                decades[decade] += 1
        print("Year distribution (5-year bins):")
        for decade in sorted(decades):
            bar = '#' * decades[decade]
            print(f"  {decade}-{decade+4}: {bar} ({decades[decade]})")
    print()
    print(f"Output: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()

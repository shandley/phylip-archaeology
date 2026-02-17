#!/usr/bin/env python3
"""Enrich tools.json with URLs extracted from catalog detail sub-pages,
then check URL liveness and classify tool status.

Phase 1: Extract URLs from software.pars.html, software.dist.html,
         software.etc1.html, software.etc2.html, software.serv.html
Phase 2: Check URL liveness (HEAD/GET) and classify status
Phase 3: Check Wayback Machine for dead URLs
Phase 4: Update tools.json with enriched data
"""

import json
import re
import os
import sys
import time
import urllib.request
import urllib.error
import ssl
from collections import defaultdict
from datetime import date
from html import unescape

SNAPSHOT_DIR = os.path.join(os.path.dirname(__file__), "snapshots")
TOOLS_FILE = os.path.join(os.path.dirname(__file__), "tools.json")
TODAY = date.today().isoformat()

# Detail sub-pages to parse
DETAIL_PAGES = [
    "software.pars.html",
    "software.dist.html",
    "software.etc1.html",
    "software.etc2.html",
    "software.serv.html",
]

# ============================================================
# PHASE 1: Extract tool entries from detail sub-pages
# ============================================================

def strip_html(text):
    """Remove HTML tags and normalize whitespace."""
    text = re.sub(r'<[^>]+>', ' ', text)
    text = unescape(text)
    # Normalize whitespace
    text = re.sub(r'\s+', ' ', text)
    # Clean up email obfuscation
    text = text.replace('(at)', '@').replace(' @ ', '@')
    return text.strip()


def extract_entries_from_page(filepath):
    """Extract tool entries from a detail sub-page.

    Each tool entry starts with <A NAME="anchor_id"> and runs until
    the next <A NAME=...> anchor. We use anchor-to-anchor boundaries
    because some entries have <HR> immediately after the anchor tag
    (e.g., <A NAME="MrBayes"><HR>).

    Returns dict: anchor_id -> {urls: [...], description: str, raw_html: str}
    """
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        html = f.read()

    # Find all anchor positions
    anchor_pattern = re.compile(
        r'<A\s+NAME="([^"]+)"[^>]*>', re.IGNORECASE
    )

    skip_anchors = {'top', 'bottom', 'methods', 'systems', 'datatypes',
                    'recent', 'new', 'news', 'changes', 'otherlists',
                    'ftpservers', 'Bayesian', 'PHYLIP-servers',
                    'phylip-servers'}

    # Collect all anchor positions (including skipped ones, for boundary detection)
    all_anchors = []
    for m in anchor_pattern.finditer(html):
        all_anchors.append((m.start(), m.group(1)))

    # Extract entries: from one anchor to the next anchor
    entries = {}
    for i, (pos, anchor_id) in enumerate(all_anchors):
        if anchor_id in skip_anchors:
            continue

        # End is at the next anchor, or end of file
        if i + 1 < len(all_anchors):
            end = all_anchors[i + 1][0]
        else:
            end = len(html)

        raw_html = html[pos:end]

        # Extract all external URLs from this entry
        url_pattern = re.compile(
            r'<A\s+HREF="(https?://[^"]+)"[^>]*>',
            re.IGNORECASE
        )
        urls = []
        for um in url_pattern.finditer(raw_html):
            url = um.group(1).strip()
            # Skip common non-tool URLs
            if any(skip in url for skip in [
                'phylipweb.github.io', 'evolution.gs.washington.edu/phylip',
                'software-form', 'wikipedia.org', 'ncbi.nlm.nih.gov/pubmed',
                'doi.org', 'scholar.google', 'amazon.com', 'sinauer.com',
            ]):
                # But keep PHYLIP's own URL
                if anchor_id == 'PHYLIP' and 'phylip.html' in url:
                    urls.append(url)
                continue
            urls.append(url)

        # Extract description (first ~300 chars of plain text)
        desc_text = strip_html(raw_html)
        # Remove the tool name from start if present
        if desc_text:
            # Trim to reasonable length for a description
            if len(desc_text) > 500:
                desc_text = desc_text[:500].rsplit(' ', 1)[0] + '...'

        entries[anchor_id] = {
            'urls': urls,
            'description': desc_text,
        }

    return entries


def pick_best_url(urls):
    """Pick the most likely homepage URL from a list of URLs."""
    if not urls:
        return None
    if len(urls) == 1:
        return urls[0]

    # Prefer URLs that look like homepages
    for url in urls:
        lower = url.lower()
        # Prefer GitHub repos
        if 'github.com' in lower or 'github.io' in lower:
            return url
        # Prefer SourceForge
        if 'sourceforge.net' in lower:
            return url
        # Prefer CRAN
        if 'cran.r-project.org' in lower:
            return url
        # Prefer Bioconductor
        if 'bioconductor.org' in lower:
            return url

    # Prefer "web site" or "web page" links (the first one usually)
    return urls[0]


def extract_all_urls():
    """Extract URLs from all detail sub-pages."""
    all_entries = {}

    for page in DETAIL_PAGES:
        filepath = os.path.join(SNAPSHOT_DIR, page)
        if not os.path.exists(filepath):
            print(f"  WARNING: {page} not found, skipping")
            continue
        entries = extract_entries_from_page(filepath)
        print(f"  {page}: {len(entries)} entries, "
              f"{sum(1 for e in entries.values() if e['urls'])} with URLs")

        for anchor_id, entry in entries.items():
            if anchor_id not in all_entries:
                all_entries[anchor_id] = entry
            else:
                # Merge URLs if tool appears in multiple pages
                existing = all_entries[anchor_id]
                for url in entry['urls']:
                    if url not in existing['urls']:
                        existing['urls'].append(url)
                # Use longer description
                if len(entry['description']) > len(existing['description']):
                    existing['description'] = entry['description']

    return all_entries


# ============================================================
# PHASE 2: Match extracted entries to tools.json
# ============================================================

def build_name_to_anchor_map(tools, entries):
    """Build a mapping from tool names to anchor IDs.

    The challenge: tools.json uses display names, but entries use anchor IDs.
    We need to match them.
    """
    # Direct name matches (case-insensitive)
    anchor_to_name = {}
    name_lower_map = {t['name'].lower(): t['name'] for t in tools}

    for anchor_id in entries:
        # Try exact match first
        if anchor_id in name_lower_map:
            anchor_to_name[anchor_id] = name_lower_map[anchor_id]
        elif anchor_id.lower() in name_lower_map:
            anchor_to_name[anchor_id] = name_lower_map[anchor_id.lower()]
        else:
            # Try with common variations
            variations = [
                anchor_id,
                anchor_id.upper(),
                anchor_id.lower(),
                anchor_id.replace('_', ' '),
                anchor_id.replace('-', ' '),
            ]
            matched = False
            for v in variations:
                if v.lower() in name_lower_map:
                    anchor_to_name[anchor_id] = name_lower_map[v.lower()]
                    matched = True
                    break

    return anchor_to_name


# ============================================================
# PHASE 3: URL liveness checking
# ============================================================

def check_url(url, timeout=15):
    """Check if a URL is reachable. Returns (status_code, final_url, error)."""
    # Create SSL context that doesn't verify (many old academic sites have bad certs)
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) '
                       'AppleWebKit/537.36 (KHTML, like Gecko) '
                       'Chrome/120.0.0.0 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,*/*',
    }

    # Try HEAD first, then GET
    for method in ['HEAD', 'GET']:
        try:
            req = urllib.request.Request(url, headers=headers, method=method)
            resp = urllib.request.urlopen(req, timeout=timeout, context=ctx)
            return (resp.status, resp.url, None)
        except urllib.error.HTTPError as e:
            if e.code in (405, 403) and method == 'HEAD':
                continue  # Try GET
            return (e.code, url, str(e))
        except urllib.error.URLError as e:
            return (None, url, str(e.reason))
        except Exception as e:
            return (None, url, str(e))

    return (None, url, "All methods failed")


def check_wayback(url, timeout=15):
    """Check Wayback Machine for an archived version of a URL."""
    api_url = f"https://archive.org/wayback/available?url={url}"
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    try:
        req = urllib.request.Request(api_url, headers={
            'User-Agent': 'phylip-archaeology/1.0'
        })
        resp = urllib.request.urlopen(req, timeout=timeout, context=ctx)
        data = json.loads(resp.read().decode('utf-8'))
        snapshots = data.get('archived_snapshots', {})
        closest = snapshots.get('closest', {})
        if closest.get('available'):
            return closest.get('url')
    except Exception:
        pass
    return None


def classify_status(status_code, error, url):
    """Classify a tool's status based on URL check results."""
    if status_code is not None and 200 <= status_code < 400:
        # Site is reachable
        url_lower = url.lower()
        # GitHub/SourceForge/CRAN = likely maintained or at least archived
        if any(host in url_lower for host in [
            'github.com', 'github.io', 'gitlab.com', 'bitbucket.org',
            'cran.r-project.org', 'bioconductor.org',
        ]):
            return 'archived'  # Conservative: we can't tell maintained from code alone
        return 'dormant'  # Site up but we don't know if actively maintained
    elif status_code == 404:
        return 'dead'
    elif status_code in (301, 302, 307, 308):
        return 'dormant'  # Redirect, could be alive
    else:
        return 'dead'


# ============================================================
# MAIN
# ============================================================

def main():
    # Load current tools.json
    with open(TOOLS_FILE, "r", encoding="utf-8") as f:
        tools = json.load(f)

    print(f"Loaded {len(tools)} tools from tools.json")
    print()

    # Phase 1: Extract URLs from detail pages
    print("=== Phase 1: Extracting URLs from detail pages ===")
    entries = extract_all_urls()
    print(f"\nTotal entries extracted: {len(entries)}")
    with_urls = sum(1 for e in entries.values() if e['urls'])
    print(f"Entries with URLs: {with_urls}")

    # Phase 2: Match to tools.json
    print("\n=== Phase 2: Matching entries to tools.json ===")
    anchor_to_name = build_name_to_anchor_map(tools, entries)
    print(f"Matched {len(anchor_to_name)} entries to tool names")

    # Build name -> tool index
    name_to_idx = {t['name']: i for i, t in enumerate(tools)}

    # Enrich tools with URLs and descriptions
    enriched_count = 0
    url_count = 0
    desc_count = 0
    for anchor_id, entry in entries.items():
        tool_name = anchor_to_name.get(anchor_id)
        if not tool_name or tool_name not in name_to_idx:
            continue

        idx = name_to_idx[tool_name]
        tool = tools[idx]
        enriched_count += 1

        # Add URL if not already present
        if not tool.get('url_original') and entry['urls']:
            best_url = pick_best_url(entry['urls'])
            if best_url:
                tool['url_original'] = best_url
                url_count += 1

        # Add description if not already present or if new one is better
        if entry['description'] and (
            not tool.get('description') or
            len(entry['description']) > len(tool.get('description', ''))
        ):
            tool['description'] = entry['description']
            desc_count += 1

    print(f"Enriched {enriched_count} tools")
    print(f"  Added URLs to {url_count} tools")
    print(f"  Added/improved descriptions for {desc_count} tools")

    # Count tools with URLs now
    tools_with_urls = [(i, t) for i, t in enumerate(tools) if t.get('url_original')]
    print(f"\nTotal tools with URLs: {len(tools_with_urls)}")

    # Phase 3: Check URL liveness
    print("\n=== Phase 3: Checking URL liveness ===")
    total = len(tools_with_urls)
    alive = 0
    dead = 0
    archive_found = 0

    for count, (idx, tool) in enumerate(tools_with_urls, 1):
        url = tool['url_original']
        name = tool['name']

        sys.stdout.write(f"\r  [{count}/{total}] Checking {name[:30]:30s}")
        sys.stdout.flush()

        status_code, final_url, error = check_url(url)
        status = classify_status(status_code, error, url)

        tool['status'] = status
        tool['status_checked'] = TODAY

        if status in ('archived', 'dormant', 'maintained'):
            alive += 1
        else:
            dead += 1
            # Check Wayback Machine for dead URLs
            wayback_url = check_wayback(url)
            if wayback_url:
                tool['url_archive'] = wayback_url
                tool['status'] = 'archived'
                archive_found += 1
                alive += 1
                dead -= 1

        # Rate limiting
        time.sleep(0.3)

    print(f"\n\nResults:")
    print(f"  Reachable: {alive}")
    print(f"  Dead: {dead}")
    print(f"  Recovered via Wayback: {archive_found}")

    # Mark tools without URLs as 'unknown' (already are, but set checked date)
    no_url_count = 0
    for tool in tools:
        if not tool.get('url_original') and tool['status'] == 'unknown':
            tool['status_checked'] = TODAY
            no_url_count += 1

    print(f"  No URL (remain unknown): {no_url_count}")

    # Phase 4: Summary and save
    print("\n=== Phase 4: Saving enriched tools.json ===")
    status_counts = defaultdict(int)
    for tool in tools:
        status_counts[tool['status']] += 1

    print("Final status distribution:")
    for status, count in sorted(status_counts.items(), key=lambda x: -x[1]):
        print(f"  {status}: {count}")

    with open(TOOLS_FILE, "w", encoding="utf-8") as f:
        json.dump(tools, f, indent=2, ensure_ascii=False)

    print(f"\nSaved {len(tools)} tools to {TOOLS_FILE}")


if __name__ == "__main__":
    main()

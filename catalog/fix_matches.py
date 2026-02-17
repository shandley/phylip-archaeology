#!/usr/bin/env python3
"""Fix unmatched anchor IDs and extract additional URLs.

This script handles:
1. Close name variations (e.g., "PAUP" -> "PAUP*", "Kakusan" -> "Kakusan4")
2. Web servers that should be separate entries or merged with existing tools
3. Adds url_original from detail pages for previously unmatched tools
"""

import json
import os
import re
import ssl
import sys
import time
import urllib.request
import urllib.error
from datetime import date
from html import unescape

SNAPSHOT_DIR = os.path.join(os.path.dirname(__file__), "snapshots")
TOOLS_FILE = os.path.join(os.path.dirname(__file__), "tools.json")
TODAY = date.today().isoformat()

# Manual anchor -> tools.json name mappings
MANUAL_MATCHES = {
    "PAUP": "PAUP*",
    "GelCompar": "GelCompar II",
    "Kakusan": "Kakusan4",
    "Permute": "Permute!",
    "PhySIC": "PhySIC_IST",
    "IM": "IMa2",
    "PI": "Phylogenetic Investigator",
}

DETAIL_PAGES = [
    "software.pars.html",
    "software.dist.html",
    "software.etc1.html",
    "software.etc2.html",
    "software.serv.html",
]

def extract_entry(html, anchor_id):
    """Extract a single entry from HTML by anchor ID."""
    pattern = re.compile(
        rf'<A\s+NAME="{re.escape(anchor_id)}"[^>]*>(.*?)(?=<HR|<A\s+NAME=)',
        re.IGNORECASE | re.DOTALL
    )
    m = pattern.search(html)
    if not m:
        return None, []

    raw = m.group(0)

    # Extract URLs
    url_pattern = re.compile(r'<A\s+HREF="(https?://[^"]+)"', re.IGNORECASE)
    urls = []
    for um in url_pattern.finditer(raw):
        url = um.group(1).strip()
        if not any(skip in url for skip in [
            'phylipweb.github.io', 'software-form', 'wikipedia.org',
            'ncbi.nlm.nih.gov/pubmed', 'doi.org', 'scholar.google',
            'amazon.com', 'sinauer.com',
        ]):
            urls.append(url)

    # Extract description
    desc = re.sub(r'<[^>]+>', ' ', raw)
    desc = unescape(desc)
    desc = re.sub(r'\s+', ' ', desc).strip()
    if len(desc) > 500:
        desc = desc[:500].rsplit(' ', 1)[0] + '...'

    return desc, urls


def check_url(url, timeout=15):
    """Check if a URL is reachable."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) '
                       'AppleWebKit/537.36 Chrome/120.0.0.0 Safari/537.36',
    }
    for method in ['HEAD', 'GET']:
        try:
            req = urllib.request.Request(url, headers=headers, method=method)
            resp = urllib.request.urlopen(req, timeout=timeout, context=ctx)
            return resp.status
        except urllib.error.HTTPError as e:
            if e.code in (405, 403) and method == 'HEAD':
                continue
            return e.code
        except Exception:
            return None
    return None


def check_wayback(url, timeout=15):
    """Check Wayback Machine for an archived version."""
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
        closest = data.get('archived_snapshots', {}).get('closest', {})
        if closest.get('available'):
            return closest.get('url')
    except Exception:
        pass
    return None


def main():
    with open(TOOLS_FILE) as f:
        tools = json.load(f)

    name_to_idx = {t['name']: i for i, t in enumerate(tools)}
    name_lower = {t['name'].lower(): t['name'] for t in tools}

    # Load all detail page HTML
    pages_html = {}
    for page in DETAIL_PAGES:
        path = os.path.join(SNAPSHOT_DIR, page)
        if os.path.exists(path):
            with open(path) as f:
                pages_html[page] = f.read()

    # Find all anchor IDs
    anchor_pattern = re.compile(r'<A\s+NAME="([^"]+)"', re.IGNORECASE)
    skip = {'top','bottom','methods','systems','datatypes','recent','new','news',
            'changes','otherlists','ftpservers','Bayesian','PHYLIP-servers',
            'phylip-servers'}

    new_urls = 0
    new_tools = 0
    urls_to_check = []

    for page, html in pages_html.items():
        for m in anchor_pattern.finditer(html):
            aid = m.group(1)
            if aid in skip:
                continue

            # Determine target tool name
            target = None

            # Check manual matches
            if aid in MANUAL_MATCHES:
                target = MANUAL_MATCHES[aid]
            # Check exact match
            elif aid in name_to_idx:
                target = aid
            # Check case-insensitive
            elif aid.lower() in name_lower:
                target = name_lower[aid.lower()]

            if target and target in name_to_idx:
                idx = name_to_idx[target]
                tool = tools[idx]

                # Only add URL if tool doesn't have one
                if not tool.get('url_original'):
                    desc, urls = extract_entry(html, aid)
                    if urls:
                        tool['url_original'] = urls[0]
                        urls_to_check.append((idx, urls[0]))
                        new_urls += 1
                    if desc and (not tool.get('description') or
                                 len(desc) > len(tool.get('description', ''))):
                        tool['description'] = desc

    print(f"Added {new_urls} new URLs from manual matching")
    print(f"Checking {len(urls_to_check)} new URLs...")

    for count, (idx, url) in enumerate(urls_to_check, 1):
        tool = tools[idx]
        sys.stdout.write(f"\r  [{count}/{len(urls_to_check)}] {tool['name'][:30]:30s}")
        sys.stdout.flush()

        status_code = check_url(url)
        if status_code and 200 <= status_code < 400:
            url_lower = url.lower()
            if any(h in url_lower for h in ['github.com', 'github.io', 'gitlab.com',
                                             'cran.r-project.org', 'bioconductor.org']):
                tool['status'] = 'archived'
            else:
                tool['status'] = 'dormant'
        else:
            wayback = check_wayback(url)
            if wayback:
                tool['url_archive'] = wayback
                tool['status'] = 'archived'
            else:
                tool['status'] = 'dead'

        tool['status_checked'] = TODAY
        time.sleep(0.3)

    print()

    # Save
    with open(TOOLS_FILE, 'w') as f:
        json.dump(tools, f, indent=2, ensure_ascii=False)

    from collections import Counter
    statuses = Counter(t['status'] for t in tools)
    print("\nFinal status distribution:")
    for s, c in statuses.most_common():
        print(f"  {s}: {c}")

    print(f"\nTotal tools with URLs: {sum(1 for t in tools if t.get('url_original'))}")
    print(f"Saved to {TOOLS_FILE}")


if __name__ == "__main__":
    main()

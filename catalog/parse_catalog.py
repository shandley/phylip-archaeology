#!/usr/bin/env python3
"""Parse Felsenstein's software catalog HTML into structured tools.json.

This script parses the main software.html index page from Felsenstein's
Phylogeny Programs catalog. It extracts tool names, categories (from the
"by methods" sections), platforms (from the "by computer systems" sections),
data types, and any available URLs and descriptions.
"""

import json
import re
from collections import defaultdict

HTML_FILE = "/Users/scotthandley/Code/exploration/phylip-archaeology/catalog/snapshots/software.html"
OUTPUT_FILE = "/Users/scotthandley/Code/exploration/phylip-archaeology/catalog/tools.json"

# Read the HTML
with open(HTML_FILE, "r", encoding="utf-8", errors="replace") as f:
    html = f.read()

# ============================================================
# CONFIGURATION: Section-to-category mappings
# ============================================================

# Section anchors in the "by methods" portion -> schema categories
method_section_map = {
    "General": "general-purpose",
    "Parsimony": "parsimony",
    "Distance": "distance",
    "Distances": "distance",
    "Likelihood": "likelihood",
    "Bayesian": "bayesian",
    "Quartets": "quartets",
    "AI": "general-purpose",
    "Invariants": "parsimony",
    "Interactive": "tree-visualization",
    "Recombinant": "recombination",
    "Bootstrapping": "bootstrap",
    "Compatibility": "compatibility",
    "Consensus": "consensus",
    "Alignment": "alignment",
    "Genomics": "gene-duplication",
    "Biogeography": "biogeography",
    "Comparative": "comparative-methods",
    "Simulation": "simulation",
    "Shapes": "tree-comparison",
    "Stratigraphy": "clock-dating",
    "Modelselection": "model-selection",
    "Description": "comparative-methods",
    "Plotting": "tree-visualization",
    "Management": "utility",
    "Teaching": "teaching",
}

platform_section_map = {
    "Unix": "unix",
    "interpreters": "cross-platform",
    "Windows": "windows",
    "DOS": "dos",
    "Macintosh": "macos",
}

# Known anchor ID merges: anchor IDs that refer to the same tool
# (due to case inconsistencies in the HTML)
anchor_merges = {
    "tracer": "Tracer",       # lowercase anchor used in some sections
    "NoNA": "NONA",           # alternate anchor for NONA
    "poptrfdos": "POPTREE2",  # anchor refers to POPTREE2
    "vCEBL": "PEBBLE",        # anchor refers to PEBBLE
}

# Force canonical names for anchors with conflicting display names
# (due to typos in the HTML, e.g., #mesquite -> "PHYLIP" on one line)
forced_canonical = {
    "mesquite": "Mesquite",      # typo: one link says #mesquite -> "PHYLIP"
    "MUST": "MUST",              # also listed as "MUST 2000"
    "COMPONENT": "COMPONENT",   # also listed as "COMPONENT Lite"
    "Permute": "Permute!",      # some links drop the "!"
    "Modeltest": "Modeltest",   # one occurrence says "ModelTest" but most say "Modeltest"
}

# Non-tool anchor IDs to ignore completely
ignore_anchors = {
    "methods", "systems", "datatypes", "recent", "new", "news",
    "changes", "otherlists", "ftpservers",
    "Bayesian",  # "Bayesian inference list" is a section link, not a tool
}

# Tools in the "Mysteries" section and non-tool links to exclude
exclude_names = {
    "here", "submission form", "look here",
    "an extensive table", "List of sequence alignment software",
    "That site", "Bioinformatics software for Mac OS X",
    "bioinformatics.org", "Molecular Linux",
    "List of phylogenetic tree visualization software",
    "Phylogenetics Software Resources", "list of software",
    "a list: Programas para Filogenia",
    "Free Software in Molecular Biology for Macintosh and MS Windows computers",
    "SoftwareSeek searchable index of bioinformatics software",
    "Tree and Tree-software for visualisation and manipulations dealing with phylogenetic trees",
}

# Regex to match tool links in the catalog
link_pattern = re.compile(r'<[Aa]\s+[Hh][Rr][Ee][Ff]="([^"]*?#([^"]+))"[^>]*>([^<]+)</[Aa]>', re.DOTALL)

# Regex for full URL links (waiting list)
full_url_pattern = re.compile(r'<[Aa]\s+[Hh][Rr][Ee][Ff]="(https?://[^"]+)"[^>]*>([^<]+)</[Aa]>', re.DOTALL)

# ============================================================
# STEP 1: Build anchor_id -> canonical name mapping
# ============================================================

anchor_to_names = defaultdict(set)
for match in link_pattern.finditer(html):
    href, anchor_id, display = match.groups()
    display = display.strip()
    display = re.sub(r'\s+', ' ', display)
    if display and anchor_id and display not in exclude_names:
        anchor_to_names[anchor_id].add(display)

def pick_canonical_name(names_set):
    """Pick the best canonical name from a set of display names for a tool."""
    names = list(names_set)
    if len(names) == 1:
        return names[0]
    # Filter out generic-looking names
    filtered = [n for n in names if len(n) > 1 and n not in exclude_names]
    if not filtered:
        return names[0]
    # Prefer names that start with uppercase (proper program names)
    proper = [n for n in filtered if n[0].isupper()]
    if proper:
        filtered = proper
    # Return shortest
    return min(filtered, key=len)

anchor_canonical = {}
for aid, names in anchor_to_names.items():
    if aid in ignore_anchors:
        continue
    # Apply merges
    if aid in anchor_merges:
        anchor_canonical[aid] = anchor_merges[aid]
    elif aid in forced_canonical:
        anchor_canonical[aid] = forced_canonical[aid]
    else:
        anchor_canonical[aid] = pick_canonical_name(names)

# ============================================================
# STEP 2: Find section boundaries
# ============================================================

methods_marker = html.find('<A NAME="methods">')
systems_marker = html.find('<A NAME="systems">')
datatypes_marker = html.find('<A NAME="datatypes">')
recent_marker = html.find('<A NAME="recent">')
new_marker = html.find('<A NAME="new">')
news_marker = html.find('<A NAME="news">')

def find_section_ranges(html_text, section_map, area_start, area_end):
    """Find section ranges within a given area of the HTML."""
    sections = []
    for section_name, category in section_map.items():
        pattern = rf'<A\s+NAME="{re.escape(section_name)}"'
        for match in re.finditer(pattern, html_text, re.IGNORECASE):
            pos = match.start()
            if area_start <= pos < area_end:
                sections.append((pos, section_name, category))
                break

    sections.sort(key=lambda x: x[0])

    result = []
    for i, (start, name, cat) in enumerate(sections):
        if i + 1 < len(sections):
            end = sections[i + 1][0]
        else:
            end = area_end
        result.append((start, end, name, cat))
    return result

def extract_tools_from_range(html_text, start, end):
    """Extract all tool references from a range of HTML."""
    segment = html_text[start:end]
    tools = []
    for match in link_pattern.finditer(segment):
        href, anchor_id, display = match.groups()
        display = display.strip()
        display = re.sub(r'\s+', ' ', display)
        if anchor_id and display and display not in exclude_names and anchor_id not in ignore_anchors:
            # Apply anchor merges and forced canonical names
            canonical = anchor_canonical.get(anchor_id, display)
            tools.append((anchor_id, canonical, href))
    return tools

# ============================================================
# STEP 3: Parse method sections -> categories
# ============================================================

tool_categories = defaultdict(set)
tool_names = {}
tool_hrefs = {}

method_sections = find_section_ranges(html, method_section_map, methods_marker, systems_marker)

for start, end, section_name, category in method_sections:
    tools = extract_tools_from_range(html, start, end)
    for anchor_id, name, href in tools:
        tool_categories[name].add(category)
        tool_names[anchor_id] = name
        tool_hrefs[anchor_id] = href

# ============================================================
# STEP 4: Parse platform sections -> platforms
# ============================================================

tool_platforms = defaultdict(set)

platform_sections = find_section_ranges(html, platform_section_map, systems_marker, datatypes_marker)

for start, end, section_name, platform in platform_sections:
    tools = extract_tools_from_range(html, start, end)
    for anchor_id, name, href in tools:
        tool_platforms[name].add(platform)
        if anchor_id not in tool_names:
            tool_names[anchor_id] = name
        if anchor_id not in tool_hrefs:
            tool_hrefs[anchor_id] = href

# ============================================================
# STEP 5: Parse Java subsection -> languages
# ============================================================

tool_languages = defaultdict(set)

java_label = html.find('<strong>Java</strong>', systems_marker)
if java_label >= 0:
    pcs_section = html.find('<A NAME="PCs">', java_label)
    if pcs_section < 0:
        pcs_section = datatypes_marker
    java_tools = extract_tools_from_range(html, java_label, pcs_section)
    for anchor_id, name, href in java_tools:
        tool_languages[name].add("java")

# ============================================================
# STEP 6: Parse data type sections -> data_types
# ============================================================

tool_datatypes = defaultdict(set)

if datatypes_marker >= 0:
    dt_end = recent_marker if recent_marker >= 0 else len(html)
    dt_segment = html[datatypes_marker:dt_end]

    dt_sections = {
        r'<B>Microsatellite data</B>': "microsatellites",
        r'<B>RAPDs, RFLPs, or AFLPs</B>': "rapd-rflp-aflp",
        r'<B>Continuous quantitative</B>': "continuous-characters",
        r'<B>Gene frequencies</B>': "gene-frequencies",
    }

    for pattern, dtype in dt_sections.items():
        match = re.search(pattern, dt_segment)
        if match:
            ul_start = dt_segment.find('<UL>', match.start())
            ul_end = dt_segment.find('</UL>', ul_start)
            if ul_start >= 0 and ul_end >= 0:
                section = dt_segment[ul_start:ul_end]
                for m in link_pattern.finditer(section):
                    href, anchor_id, display = m.groups()
                    name = anchor_canonical.get(anchor_id, display.strip())
                    tool_datatypes[name].add(dtype)

# ============================================================
# STEP 7: Parse waiting list -> additional tools
# ============================================================

waiting_tools = {}

if new_marker >= 0:
    w_end = news_marker if news_marker >= 0 else len(html)
    waiting_segment = html[new_marker:w_end]

    # Extract each LI entry
    li_pattern = re.compile(
        r'<[Ll][Ii]>\s*<[Aa]\s+[Hh][Rr][Ee][Ff]="(https?://[^"]+)"[^>]*>([^<]+)</[Aa]>[,\s]*([^<]*)',
        re.DOTALL
    )
    for m in li_pattern.finditer(waiting_segment):
        url, name, desc = m.groups()
        name = name.strip()
        desc = desc.strip()
        # Clean description
        desc = re.sub(r'\s+', ' ', desc).strip()
        if desc.endswith('.'):
            desc = desc[:-1].strip()

        # Skip non-tools
        if name.lower() in ('here', 'submission form') or 'software-form' in url:
            continue

        waiting_tools[name] = {
            "url": url,
            "description": desc if desc and len(desc) > 5 else None,
        }

    # Also handle the "also described" pattern for MP-EST
    also_pattern = re.compile(
        r'<[Aa]\s+[Hh][Rr][Ee][Ff]="(https?://[^"]+)"[^>]*>([^<]+)</[Aa]>\s*\(also',
        re.DOTALL
    )

# ============================================================
# STEP 8: Build the consolidated tools dictionary
# ============================================================

tools_dict = {}

# Collect all unique tool names
all_names = set()
for aid, name in tool_names.items():
    all_names.add(name)
for name in tool_categories:
    all_names.add(name)
for name in tool_platforms:
    all_names.add(name)

for name in all_names:
    if name in exclude_names:
        continue

    entry = {
        "name": name,
        "status": "unknown",
        "categories": sorted(tool_categories.get(name, set())),
    }

    # Set default category if none found
    if not entry["categories"]:
        entry["categories"] = ["general-purpose"]

    # Platforms
    platforms = sorted(tool_platforms.get(name, set()))
    if platforms:
        entry["platforms"] = platforms

    # Languages
    languages = sorted(tool_languages.get(name, set()))
    if languages:
        entry["languages"] = languages

    # Data types
    data_types = sorted(tool_datatypes.get(name, set()))
    if data_types:
        entry["data_types"] = data_types

    tools_dict[name] = entry

# Add waiting list tools
for name, info in waiting_tools.items():
    if name in tools_dict:
        # Add URL
        if info["url"]:
            tools_dict[name]["url_original"] = info["url"]
        if info["description"]:
            tools_dict[name]["description"] = info["description"]
    else:
        entry = {
            "name": name,
            "url_original": info["url"],
            "status": "unknown",
            "categories": ["general-purpose"],
        }
        if info["description"]:
            entry["description"] = info["description"]
        tools_dict[name] = entry

# ============================================================
# STEP 9: Post-processing refinements
# ============================================================

# Add "supertrees" category to supertree-specific tools
consensus_start = html.find('<A NAME="Consensus">')
if consensus_start >= 0:
    alignment_start = html.find('<A NAME="Alignment">', consensus_start)
    if alignment_start < 0:
        alignment_start = systems_marker
    consensus_segment = html[consensus_start:alignment_start]
    for m in link_pattern.finditer(consensus_segment):
        href, anchor_id, display = m.groups()
        display = display.strip()
        if 'supertree' in display.lower() or 'SuperTree' in display:
            name = anchor_canonical.get(anchor_id, display)
            if name in tools_dict:
                cats = set(tools_dict[name]["categories"])
                cats.add("supertrees")
                tools_dict[name]["categories"] = sorted(cats)

# Add "host-parasite" for known host-parasite tools in biogeography section
bio_start = html.find('<A NAME="Biogeography">')
if bio_start >= 0:
    comp_start = html.find('<A NAME="Comparative">', bio_start)
    if comp_start < 0:
        comp_start = systems_marker
    bio_segment = html[bio_start:comp_start]
    host_parasite_anchors = {'COMPONENT', 'TREEMAP', 'TreeFitter', 'Tarzan', 'ParaFit',
                             'AxParafit', 'CoRe-PA', 'Jane'}
    for m in link_pattern.finditer(bio_segment):
        href, anchor_id, display = m.groups()
        if anchor_id in host_parasite_anchors:
            name = anchor_canonical.get(anchor_id, display.strip())
            if name in tools_dict:
                cats = set(tools_dict[name]["categories"])
                cats.add("host-parasite")
                tools_dict[name]["categories"] = sorted(cats)

# Add "network" category
network_tools = ['Network', 'SplitsTree', 'Spectronet', 'TCS']
for nt in network_tools:
    if nt in tools_dict:
        cats = set(tools_dict[nt]["categories"])
        cats.add("network")
        tools_dict[nt]["categories"] = sorted(cats)

# Add "population-genetics" category
pop_gen_tools = ['POPGENE', 'Arlequin', 'GenoDive', 'FSTAT', 'Genepop', 'DnaSP',
                 'SPAGeDi', 'TFPGA', 'RSTCALC', 'Genetix', 'IMa2', 'Populations',
                 'GDA', 'DIVAGE']
for pg in pop_gen_tools:
    if pg in tools_dict:
        cats = set(tools_dict[pg]["categories"])
        cats.add("population-genetics")
        tools_dict[pg]["categories"] = sorted(cats)

# Add "sequence-editing" category
seq_edit_tools = ['BioEdit', 'GeneDoc', 'SeqPup', 'Phyledit', 'GeneStudio Pro',
                  'ARB', 'PHYDIT', 'ProSeq']
for se in seq_edit_tools:
    if se in tools_dict:
        cats = set(tools_dict[se]["categories"])
        cats.add("sequence-editing")
        tools_dict[se]["categories"] = sorted(cats)

# Add "genomics" for genomics tools
genomics_start = html.find('<A NAME="Genomics">')
if genomics_start >= 0:
    genomics_end = html.find('<A NAME="Biogeography">', genomics_start)
    if genomics_end < 0:
        genomics_end = systems_marker
    genomics_segment = html[genomics_start:genomics_end]
    for m in link_pattern.finditer(genomics_segment):
        href, anchor_id, display = m.groups()
        name = anchor_canonical.get(anchor_id, display.strip())
        if name in tools_dict:
            cats = set(tools_dict[name]["categories"])
            cats.add("genomics")
            tools_dict[name]["categories"] = sorted(cats)

# Recombination: ensure all tools in recombination section have the category
recomb_start = html.find('<A NAME="Recombinant">')
if recomb_start >= 0:
    recomb_end = html.find('<A NAME="Bootstrapping">', recomb_start)
    if recomb_end < 0:
        recomb_end = systems_marker
    recomb_segment = html[recomb_start:recomb_end]
    for m in link_pattern.finditer(recomb_segment):
        href, anchor_id, display = m.groups()
        name = anchor_canonical.get(anchor_id, display.strip())
        if name in tools_dict:
            cats = set(tools_dict[name]["categories"])
            cats.add("recombination")
            tools_dict[name]["categories"] = sorted(cats)

# ============================================================
# STEP 10: Add descriptions for well-known tools
# These are brief descriptions based on the section names and common knowledge
# ============================================================

# Map catalog section info for tools
for anchor_id, name in tool_names.items():
    if name in tools_dict and "catalog_section" not in tools_dict[name]:
        href = tool_hrefs.get(anchor_id, "")
        if "software.pars.html" in href:
            tools_dict[name]["catalog_section"] = "Parsimony programs"
        elif "software.dist.html" in href:
            tools_dict[name]["catalog_section"] = "Distance matrix and computation programs"
        elif "software.etc1.html" in href:
            tools_dict[name]["catalog_section"] = "Maximum likelihood and Bayesian programs"
        elif "software.etc2.html" in href:
            tools_dict[name]["catalog_section"] = "Tree visualization, simulation, and other programs"
        elif "software.serv.html" in href:
            tools_dict[name]["catalog_section"] = "Web servers"

# ============================================================
# STEP 11: Sort and output
# ============================================================

tools_list = sorted(tools_dict.values(), key=lambda x: x["name"].lower())

# Final cleanup
for tool in tools_list:
    if not tool.get("categories"):
        tool["categories"] = ["general-purpose"]
    for key in ["platforms", "data_types", "languages"]:
        if key in tool and not tool[key]:
            del tool[key]

print(f"Total tools extracted: {len(tools_list)}")

# Summary stats
cat_counts = defaultdict(int)
for tool in tools_list:
    for cat in tool.get("categories", []):
        cat_counts[cat] += 1

print("\nCategory distribution:")
for cat, count in sorted(cat_counts.items(), key=lambda x: -x[1]):
    print(f"  {cat}: {count}")

plat_counts = defaultdict(int)
for tool in tools_list:
    for plat in tool.get("platforms", []):
        plat_counts[plat] += 1

print("\nPlatform distribution:")
for plat, count in sorted(plat_counts.items(), key=lambda x: -x[1]):
    print(f"  {plat}: {count}")

with_cats = sum(1 for t in tools_list if t.get("categories"))
with_plat = sum(1 for t in tools_list if t.get("platforms"))
with_url = sum(1 for t in tools_list if t.get("url_original"))
with_dt = sum(1 for t in tools_list if t.get("data_types"))
with_lang = sum(1 for t in tools_list if t.get("languages"))
with_desc = sum(1 for t in tools_list if t.get("description"))
with_section = sum(1 for t in tools_list if t.get("catalog_section"))

print(f"\nTools with categories:      {with_cats}/{len(tools_list)}")
print(f"Tools with platforms:        {with_plat}/{len(tools_list)}")
print(f"Tools with url_original:     {with_url}/{len(tools_list)}")
print(f"Tools with data_types:       {with_dt}/{len(tools_list)}")
print(f"Tools with languages:        {with_lang}/{len(tools_list)}")
print(f"Tools with description:      {with_desc}/{len(tools_list)}")
print(f"Tools with catalog_section:  {with_section}/{len(tools_list)}")

with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
    json.dump(tools_list, f, indent=2, ensure_ascii=False)

print(f"\nOutput written to {OUTPUT_FILE}")

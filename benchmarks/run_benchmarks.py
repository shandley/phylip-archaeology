#!/usr/bin/env python3
"""
Benchmark orchestrator for phylip-rs vs. modern phylogenetics tools.

Reads datasets/manifest.json, runs each tool on each dataset, collects
timing/accuracy results, and writes results/benchmark_results.csv.

All tools are run single-threaded under JC69 model for fair comparison.
No external Python dependencies required (stdlib only).
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PHYLIP_RS = Path("/Users/scotthandley/.cargo/target/release/phylip-rs")
GTIME = "/opt/homebrew/bin/gtime"

BENCHMARKS_DIR = Path(__file__).resolve().parent
DATASETS_DIR = BENCHMARKS_DIR / "datasets"
RESULTS_DIR = BENCHMARKS_DIR / "results"
MANIFEST_PATH = DATASETS_DIR / "manifest.json"

TIMEOUT_S = 600          # 10 minutes default
ML_LARGE_TIMEOUT_S = 120 # Shorter timeout for phylip-rs ML on large datasets
ML_LARGE_TAXA_THRESH = 200

CSV_PATH = RESULTS_DIR / "benchmark_results.csv"
CSV_COLUMNS = [
    "dataset_id", "ntaxa", "nsites", "replicate", "tool",
    "wall_time_s", "peak_mem_kb", "lnl_jc69",
    "rf_distance", "nrf_distance", "timed_out",
]

# Tool command templates. {input}, {outfile}, {prefix} are substituted per-run.
TOOLS: dict[str, str] = {
    "phylip-rs-ml": (
        f"{PHYLIP_RS} ml --model jc69 -s 42 -o {{outfile}} {{input}}"
    ),
    "phylip-rs-nj": (
        f"{PHYLIP_RS} distance --model jc69 --method nj -o {{outfile}} {{input}}"
    ),
    "iqtree": (
        "conda run -n base iqtree -s {input} -m JC -T 1 --seed 42 "
        "-redo --prefix {prefix}"
    ),
    "raxml-ng": (
        "raxml-ng --msa {input} --model JC --threads 1 --seed 42 "
        "--prefix {prefix} --force perf_threads"
    ),
    "veryfasttree": (
        "OMP_NUM_THREADS=1 VeryFastTree -nt -nocat -nosupport -threads 1 {input}"
    ),
}


# ---------------------------------------------------------------------------
# Helpers: tree parsing
# ---------------------------------------------------------------------------

def strip_internal_labels(newick: str) -> str:
    """Remove internal node labels (e.g. bootstrap values) from a Newick string.

    Internal labels appear after ')' and before ':' or ')' or ',' or ';'.
    We keep the colon and branch length but drop the label text.
    """
    # Match: ) followed by a label (non-special chars) then either : or , or ) or ;
    return re.sub(r"\)([^:;,)\s]+)", ")", newick)


def extract_newick_from_phyliprs(output_path: Path) -> str | None:
    """Parse Newick tree from phylip-rs output file.

    Looks for lines following 'Best tree (Newick):' or 'Tree (Newick):'.
    """
    text = output_path.read_text()
    for marker in ("Best tree (Newick):", "Tree (Newick):"):
        if marker in text:
            idx = text.index(marker) + len(marker)
            rest = text[idx:].strip()
            # The Newick string is the next non-empty line
            for line in rest.splitlines():
                line = line.strip()
                if line:
                    return line
    return None


def extract_newick_from_file(path: Path) -> str | None:
    """Read first non-empty line from a file as a Newick string."""
    if not path.exists():
        return None
    for line in path.read_text().splitlines():
        line = line.strip()
        if line:
            return line
    return None


# ---------------------------------------------------------------------------
# Helpers: phylip-rs evaluation
# ---------------------------------------------------------------------------

def evaluate_likelihood(alignment_path: Path, tree_newick: str,
                        tmpdir: Path) -> float | None:
    """Evaluate the log-likelihood of a tree under JC69 using phylip-rs evaluate."""
    tree_file = tmpdir / "eval_tree.nwk"
    tree_file.write_text(tree_newick + "\n")
    eval_out = tmpdir / "eval_output.txt"

    cmd = [
        str(PHYLIP_RS), "evaluate",
        "--model", "jc69",
        "--tree", str(tree_file),
        "-o", str(eval_out),
        str(alignment_path),
    ]
    try:
        subprocess.run(cmd, capture_output=True, timeout=120, check=True)
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as exc:
        print(f"    [warn] evaluate failed: {exc}", file=sys.stderr)
        return None

    text = eval_out.read_text() if eval_out.exists() else ""
    m = re.search(r"Log-likelihood:\s+([-\d.]+)", text)
    return float(m.group(1)) if m else None


def compute_rf_distance(inferred_newick: str, true_tree_path: Path,
                        tmpdir: Path) -> tuple[int | None, float | None]:
    """Compute RF and normalized RF distance between inferred and true tree."""
    inferred_file = tmpdir / "inferred.nwk"
    inferred_file.write_text(inferred_newick + "\n")

    cmd = [
        str(PHYLIP_RS), "treedist",
        "--tree1", str(inferred_file),
        "--tree2", str(true_tree_path),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, timeout=60, text=True,
                                check=True)
        output = result.stdout
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as exc:
        print(f"    [warn] treedist failed: {exc}", file=sys.stderr)
        return None, None

    rf = nrf = None
    m = re.search(r"Robinson-Foulds distance:\s+(\d+)", output)
    if m:
        rf = int(m.group(1))
    m = re.search(r"Normalized Robinson-Foulds distance:\s+([\d.]+)", output)
    if m:
        nrf = float(m.group(1))
    return rf, nrf


# ---------------------------------------------------------------------------
# Helpers: memory measurement
# ---------------------------------------------------------------------------

def run_with_memory(cmd_str: str, timeout: int,
                    tmpdir: Path) -> tuple[subprocess.CompletedProcess | None,
                                           float, int | None, bool]:
    """Run a command wrapped in gtime to capture peak memory.

    Returns (CompletedProcess | None, wall_time_s, peak_mem_kb | None, timed_out).
    """
    # gtime writes to stderr; we capture both stdout and stderr separately
    # by routing gtime's output to a file.
    gtime_log = tmpdir / "gtime.log"
    # gtime -v writes to stderr; redirect gtime stderr to a file while
    # keeping the child's stderr separate by using shell syntax.
    # We wrap the actual command so gtime measures just the tool.
    full_cmd = f'{GTIME} -v sh -c {_shell_quote(cmd_str)} 2> {gtime_log}'

    timed_out = False
    wall_start = time.perf_counter()
    try:
        proc = subprocess.run(
            full_cmd,
            shell=True,
            capture_output=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        wall_time = time.perf_counter() - wall_start
        return None, wall_time, None, True
    wall_time = time.perf_counter() - wall_start

    peak_mem_kb = None
    if gtime_log.exists():
        gtime_text = gtime_log.read_text()
        m = re.search(r"Maximum resident set size.*?:\s+(\d+)", gtime_text)
        if m:
            peak_mem_kb = int(m.group(1))

    return proc, wall_time, peak_mem_kb, timed_out


def _shell_quote(s: str) -> str:
    """Single-quote a string for safe shell interpolation."""
    return "'" + s.replace("'", "'\"'\"'") + "'"


# ---------------------------------------------------------------------------
# Per-tool runners
# ---------------------------------------------------------------------------

def run_tool(tool_name: str, dataset: dict, tmpdir: Path
             ) -> dict:
    """Run a single tool on a single dataset and return a result dict."""
    input_path = DATASETS_DIR / dataset["fasta"]
    true_tree_path = DATASETS_DIR / dataset["tree"]
    ntaxa = dataset["ntaxa"]
    nsites = dataset["nsites"]
    replicate = dataset.get("replicate", 0)
    dataset_id = dataset["prefix"]

    result = {
        "dataset_id": dataset_id,
        "ntaxa": ntaxa,
        "nsites": nsites,
        "replicate": replicate,
        "tool": tool_name,
        "wall_time_s": "",
        "peak_mem_kb": "",
        "lnl_jc69": "",
        "rf_distance": "",
        "nrf_distance": "",
        "timed_out": False,
    }

    # Determine timeout
    timeout = TIMEOUT_S
    if tool_name == "phylip-rs-ml" and ntaxa >= ML_LARGE_TAXA_THRESH:
        timeout = ML_LARGE_TIMEOUT_S

    # Build the command
    outfile = tmpdir / f"{tool_name}_output.txt"
    prefix = tmpdir / f"{tool_name}_prefix"
    cmd_template = TOOLS[tool_name]
    cmd_str = cmd_template.format(
        input=str(input_path),
        outfile=str(outfile),
        prefix=str(prefix),
    )

    # Run with timing and memory measurement
    proc, wall_time, peak_mem_kb, timed_out = run_with_memory(
        cmd_str, timeout, tmpdir
    )
    result["wall_time_s"] = f"{wall_time:.4f}"
    if peak_mem_kb is not None:
        result["peak_mem_kb"] = peak_mem_kb
    result["timed_out"] = timed_out

    if timed_out:
        print(f"    TIMED OUT after {wall_time:.1f}s")
        return result

    if proc is None or proc.returncode != 0:
        stderr_msg = ""
        if proc is not None:
            stderr_msg = proc.stderr.decode(errors="replace")[:200]
        print(f"    FAILED (rc={proc.returncode if proc else '?'}): {stderr_msg}")
        return result

    # ----- Extract the inferred tree -----
    newick = None

    if tool_name in ("phylip-rs-ml", "phylip-rs-nj"):
        newick = extract_newick_from_phyliprs(outfile)

    elif tool_name == "iqtree":
        treefile = Path(f"{prefix}.treefile")
        newick = extract_newick_from_file(treefile)

    elif tool_name == "raxml-ng":
        treefile = Path(f"{prefix}.raxml.bestTree")
        newick = extract_newick_from_file(treefile)

    elif tool_name == "veryfasttree":
        # stdout IS the Newick tree
        stdout_text = proc.stdout.decode(errors="replace").strip()
        # VeryFastTree may print progress to stderr; stdout should be the tree
        for line in stdout_text.splitlines():
            line = line.strip()
            if line.startswith("("):
                newick = line
                break

    if newick is None:
        print("    [warn] Could not extract Newick tree from output")
        return result

    # Strip internal node labels (bootstrap values, etc.)
    clean_newick = strip_internal_labels(newick)

    # ----- Evaluate log-likelihood under JC69 -----
    lnl = evaluate_likelihood(input_path, clean_newick, tmpdir)
    if lnl is not None:
        result["lnl_jc69"] = f"{lnl:.6f}"

    # ----- Compute RF distance to true tree -----
    rf, nrf = compute_rf_distance(clean_newick, true_tree_path, tmpdir)
    if rf is not None:
        result["rf_distance"] = rf
    if nrf is not None:
        result["nrf_distance"] = f"{nrf:.6f}"

    return result


# ---------------------------------------------------------------------------
# CSV I/O
# ---------------------------------------------------------------------------

def init_csv(path: Path) -> None:
    """Write the CSV header (overwriting any existing file)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()


def append_csv(path: Path, row: dict) -> None:
    """Append a single result row to the CSV."""
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writerow(row)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    # Validate prerequisites
    if not PHYLIP_RS.exists():
        sys.exit(f"Error: phylip-rs binary not found at {PHYLIP_RS}\n"
                 f"Run: cd phylip-rs && cargo build --release")
    if not MANIFEST_PATH.exists():
        sys.exit(f"Error: manifest not found at {MANIFEST_PATH}\n"
                 f"Run the dataset generator first.")

    with open(MANIFEST_PATH) as f:
        manifest = json.load(f)

    datasets: list[dict] = manifest["datasets"]
    tool_names = list(TOOLS.keys())

    total_runs = len(datasets) * len(tool_names)
    print(f"Benchmark: {len(datasets)} datasets x {len(tool_names)} tools "
          f"= {total_runs} runs")
    print(f"Tools: {', '.join(tool_names)}")
    print(f"Output: {CSV_PATH}")
    print()

    init_csv(CSV_PATH)

    run_idx = 0
    for ds in datasets:
        for tool_name in tool_names:
            run_idx += 1
            label = f"[{run_idx}/{total_runs}]"
            ds_id = ds["prefix"]
            print(f"{label} Running {tool_name} on {ds_id}...")

            with tempfile.TemporaryDirectory(prefix="phylip_bench_") as tmpdir:
                try:
                    row = run_tool(tool_name, ds, Path(tmpdir))
                except Exception as exc:
                    print(f"    ERROR: {exc}", file=sys.stderr)
                    row = {
                        "dataset_id": ds["prefix"],
                        "ntaxa": ds["ntaxa"],
                        "nsites": ds["nsites"],
                        "replicate": ds.get("replicate", 0),
                        "tool": tool_name,
                        "wall_time_s": "",
                        "peak_mem_kb": "",
                        "lnl_jc69": "",
                        "rf_distance": "",
                        "nrf_distance": "",
                        "timed_out": False,
                    }

                append_csv(CSV_PATH, row)

                # Print summary for this run
                t = row.get("wall_time_s", "")
                mem = row.get("peak_mem_kb", "")
                lnl = row.get("lnl_jc69", "")
                rf = row.get("rf_distance", "")
                to = row.get("timed_out", False)
                parts = []
                if t:
                    parts.append(f"time={t}s")
                if mem:
                    parts.append(f"mem={mem}KB")
                if lnl:
                    parts.append(f"lnL={lnl}")
                if rf != "":
                    parts.append(f"RF={rf}")
                if to:
                    parts.append("TIMEOUT")
                print(f"    => {', '.join(parts) if parts else 'no data'}")

    print(f"\nDone. Results written to {CSV_PATH}")


if __name__ == "__main__":
    main()

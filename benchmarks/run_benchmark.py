#!/usr/bin/env python3
"""Benchmark protein2genomic at increasing query scales.

Measures wall-clock time and peak RSS for --output {all,isoform} at
N = 100, 1k, 10k, 100k random domain queries against the prebuilt
human GENCODE v49 index. Also times the one-shot --build-index pass.

Outputs:
  benchmarks/results.tsv     - one row per (mode, N, replicate)
  benchmarks/scaling.png     - wall time and peak RSS vs N
"""

import argparse
import json
import os
import random
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BIN = ROOT.parent / "build" / "protein2genomic"
GTF = Path("/data1/peerd/sotougl/protein2genomic_work/gencode.v49.primary_assembly.annotation.gtf")
IDX = Path("/data1/peerd/sotougl/protein2genomic_work/human.idx")
PROTEINS_LIST = ROOT / "proteins_list.txt"
WORK = ROOT / "work"
RESULTS_TSV = ROOT / "results.tsv"
PLOT_PNG = ROOT / "scaling.png"
INDEX_BUILD_JSON = ROOT / "index_build.json"

PEAK_RE = re.compile(r"^BENCH_PEAK_RSS_MB\s+(\d+)")

def parse_peak_mb(stderr: str) -> int:
    peak = 0
    for line in stderr.splitlines():
        m = PEAK_RE.match(line)
        if m:
            peak = max(peak, int(m.group(1)))
    return peak

def run(args, *, capture=True):
    t0 = time.perf_counter()
    p = subprocess.run(args, capture_output=capture, text=True)
    wall_s = time.perf_counter() - t0
    if p.returncode != 0:
        sys.stderr.write(p.stderr or "")
        raise SystemExit(f"binary failed: {' '.join(args)}")
    return wall_s, p.stderr or ""

def generate_bed(proteins, n, seed, path: Path):
    rng = random.Random(seed)
    sampled = rng.choices(proteins, k=n)
    with path.open("w") as f:
        for i, pid in enumerate(sampled):
            # Random aa range, length 10..300. Some will exceed protein length
            # and be reported as unmapped — that's realistic.
            start = rng.randint(1, 400)
            end = start + rng.randint(10, 300)
            f.write(f"{pid}\t{start}\t{end}\tD{i}\n")

def bench_index_build():
    # Build the index into a scratch path so we don't clobber the prebuilt one.
    target = WORK / "human_bench.idx"
    if target.exists():
        target.unlink()
    args = [str(BIN), "--gtf", str(GTF), "--build-index", "--index", str(target)]
    wall, err = run(args)
    peak = parse_peak_mb(err)
    return {"wall_s": wall, "peak_mb": peak, "index_size_mb": target.stat().st_size / 1024 / 1024,
            "gtf_size_mb": GTF.stat().st_size / 1024 / 1024}

def bench_query(mode: str, n: int, bed_path: Path, threads: int, rep: int):
    out_dir = WORK / f"out_{mode}_{n}_{rep}"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    args = [str(BIN), "--index", str(IDX), "--bed", str(bed_path),
            "--out-dir", str(out_dir), "--output", mode,
            "--threads", str(threads)]
    wall, err = run(args)
    peak = parse_peak_mb(err)
    return {"mode": mode, "n": n, "rep": rep, "threads": threads,
            "wall_s": wall, "peak_mb": peak}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=int, nargs="+", default=[100, 1000, 10000, 100000])
    ap.add_argument("--modes", nargs="+", default=["all", "isoform"])
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--skip-build", action="store_true",
                    help="Skip the --build-index benchmark (saves ~10 s).")
    args = ap.parse_args()

    if not BIN.exists():
        raise SystemExit(f"binary not found: {BIN}; run cmake + make first")
    if not IDX.exists():
        raise SystemExit(f"index not found: {IDX}")
    if not GTF.exists() and not args.skip_build:
        raise SystemExit(f"GTF not found: {GTF}; pass --skip-build to skip the build benchmark")
    if not PROTEINS_LIST.exists():
        raise SystemExit(f"missing {PROTEINS_LIST}; extract protein IDs from the GTF first")

    WORK.mkdir(exist_ok=True)
    proteins = PROTEINS_LIST.read_text().splitlines()
    print(f"# {len(proteins):,} unique protein IDs available", file=sys.stderr)

    # Pre-generate BED files (one per size, shared across reps).
    bed_paths = {}
    for n in args.sizes:
        p = WORK / f"queries_{n}.bed"
        if not p.exists():
            generate_bed(proteins, n, seed=42 * n, path=p)
            print(f"# generated {p}", file=sys.stderr)
        bed_paths[n] = p

    # Index build benchmark.
    if not args.skip_build:
        print("# benchmarking --build-index ...", file=sys.stderr)
        build_info = bench_index_build()
        INDEX_BUILD_JSON.write_text(json.dumps(build_info, indent=2))
        print(f"# build: {build_info['wall_s']:.1f}s, peak {build_info['peak_mb']} MB, "
              f"index {build_info['index_size_mb']:.1f} MB", file=sys.stderr)

    # Query benchmarks.
    rows = []
    with RESULTS_TSV.open("w") as f:
        f.write("mode\tn\trep\tthreads\twall_s\tpeak_mb\n")
        for mode in args.modes:
            for n in args.sizes:
                for rep in range(args.reps):
                    print(f"# {mode} n={n:,} rep={rep+1}/{args.reps} ...",
                          file=sys.stderr, flush=True)
                    r = bench_query(mode, n, bed_paths[n], args.threads, rep + 1)
                    rows.append(r)
                    f.write(f"{r['mode']}\t{r['n']}\t{r['rep']}\t{r['threads']}\t"
                            f"{r['wall_s']:.4f}\t{r['peak_mb']}\n")
                    print(f"  wall={r['wall_s']:.2f}s peak={r['peak_mb']}MB",
                          file=sys.stderr)

    print(f"# wrote {RESULTS_TSV}", file=sys.stderr)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

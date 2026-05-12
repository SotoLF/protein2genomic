#!/usr/bin/env python3
"""Plot benchmark scaling curves from benchmarks/results.tsv.

Produces a 1x2 matplotlib figure: wall time vs N and peak RSS vs N,
log-scaled on both axes, with per-mode lines and error bars (min/max).
"""

import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
RESULTS_TSV = ROOT / "results.tsv"
PLOT_PNG = ROOT / "scaling.png"
INDEX_BUILD_JSON = ROOT / "index_build.json"

def load_results(path):
    import csv
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            r["n"] = int(r["n"]); r["rep"] = int(r["rep"])
            r["wall_s"] = float(r["wall_s"]); r["peak_mb"] = float(r["peak_mb"])
            rows.append(r)
    return rows

def aggregate(rows, mode):
    by_n = {}
    for r in rows:
        if r["mode"] != mode:
            continue
        by_n.setdefault(r["n"], {"wall": [], "peak": []})
        by_n[r["n"]]["wall"].append(r["wall_s"])
        by_n[r["n"]]["peak"].append(r["peak_mb"])
    ns = sorted(by_n.keys())
    wall_med = [np.median(by_n[n]["wall"]) for n in ns]
    wall_min = [min(by_n[n]["wall"]) for n in ns]
    wall_max = [max(by_n[n]["wall"]) for n in ns]
    peak_med = [np.median(by_n[n]["peak"]) for n in ns]
    peak_min = [min(by_n[n]["peak"]) for n in ns]
    peak_max = [max(by_n[n]["peak"]) for n in ns]
    return ns, wall_med, wall_min, wall_max, peak_med, peak_min, peak_max

def main():
    rows = load_results(RESULTS_TSV)
    modes = sorted(set(r["mode"] for r in rows))
    colors = {"all": "#1f77b4", "isoform": "#d62728",
              "coding": "#2ca02c", "introns": "#9467bd", "span": "#ff7f0e"}

    fig, (ax_t, ax_m) = plt.subplots(1, 2, figsize=(12, 4.5))

    for mode in modes:
        ns, wm, wlo, whi, pm, plo, phi = aggregate(rows, mode)
        c = colors.get(mode, None)
        ax_t.errorbar(ns, wm,
                      yerr=[np.array(wm) - np.array(wlo), np.array(whi) - np.array(wm)],
                      marker="o", capsize=3, label=f"--output {mode}", color=c)
        ax_m.errorbar(ns, pm,
                      yerr=[np.array(pm) - np.array(plo), np.array(phi) - np.array(pm)],
                      marker="o", capsize=3, label=f"--output {mode}", color=c)

    for ax, ylab, ttl in [(ax_t, "wall time (s)", "Wall time vs query count"),
                          (ax_m, "peak RSS (MB)", "Peak RSS vs query count")]:
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("number of queries")
        ax.set_ylabel(ylab)
        ax.set_title(ttl)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="best", fontsize=9)

    suptitle = "protein2genomic scaling on human GENCODE v49"
    if INDEX_BUILD_JSON.exists():
        info = json.loads(INDEX_BUILD_JSON.read_text())
        suptitle += (f"\n(index build: {info['wall_s']:.1f}s, "
                     f"peak {info['peak_mb']} MB, "
                     f"index {info['index_size_mb']:.0f} MB from {info['gtf_size_mb']:.0f} MB GTF)")
    fig.suptitle(suptitle, fontsize=11)
    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=130)
    print(f"# wrote {PLOT_PNG}")

if __name__ == "__main__":
    main()

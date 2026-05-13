#!/usr/bin/env python3
"""Render `isoform_structure.tsv` into a transcript-architecture plot.

Usage examples
--------------
# Single domain to a PDF
protein2genomic plot --isoform isoform_structure.tsv --input-id RD1 --out RD1.pdf

# Every query in the TSV to a multipage PDF
protein2genomic plot --isoform isoform_structure.tsv --all --out queries.pdf

# Interactive HTML (requires plotly)
protein2genomic plot --isoform isoform_structure.tsv --input-id RD1 \
                    --html RD1.html

Design notes
------------
- The C++ binary writes a tidy "one row per structural segment" table; this
  script just groups by `input_id` and draws boxes. We never re-derive
  coordinates from the genome here.
- Strand is honored: genomic coordinates are always on the X axis (so the plot
  matches IGV), but a small arrow at the top indicates the 5'→3' direction of
  translation.
- Color scheme mirrors `plot_group`:
    five_prime_UTR        -> light grey
    three_prime_UTR       -> light grey
    CDS / CDS_no_domain   -> blue
    CDS_domain            -> red
    intron / intron_domain_span
                          -> thin line, dashed if outside domain span
- `--no-introns` collapses the X axis to a concatenated CDS-only view (the
  "spliced transcript" style) using `feature_order_transcript` ordering.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable

# matplotlib is required; plotly is optional and only imported when --html is set.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages


# --------------------------------------------------------------------------- #
# Row parsing
# --------------------------------------------------------------------------- #

# A single row of isoform_structure.tsv. We only keep the columns the plotter
# actually uses. NA-able numerics become `None`.
@dataclass
class Segment:
    input_id: str
    gene_name: str
    transcript_id: str
    protein_id: str
    domain_id: str
    chrom: str
    strand: str
    feature_type: str
    feature_id: str
    feature_part: int
    feature_genomic_start: int
    feature_genomic_end: int
    feature_order_transcript: int
    overlaps_domain: str
    plot_group: str
    is_mane_select: str = "NA"
    is_ensembl_canonical: str = "NA"


def _i_or_none(x: str) -> int | None:
    if x == "NA" or x == "":
        return None
    try:
        return int(x)
    except ValueError:
        return None


def load_isoform_tsv(path: str) -> dict[str, list[Segment]]:
    """Group rows by input_id, preserving file order within each group."""
    by_id: dict[str, list[Segment]] = defaultdict(list)
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {
            "input_id", "feature_type", "feature_genomic_start",
            "feature_genomic_end", "plot_group", "strand", "chrom",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f"isoform tsv {path!r} is missing columns: {sorted(missing)}"
            )
        for row in reader:
            seg = Segment(
                input_id=row.get("input_id", ""),
                gene_name=row.get("gene_name", "") or "",
                transcript_id=row.get("transcript_id", "") or "",
                protein_id=row.get("protein_id", "") or "",
                domain_id=row.get("domain_id", "") or "",
                chrom=row.get("chrom", "") or "",
                strand=row.get("strand", "+") or "+",
                feature_type=row.get("feature_type", "") or "",
                feature_id=row.get("feature_id", "") or "",
                feature_part=_i_or_none(row.get("feature_part", "1")) or 1,
                feature_genomic_start=int(row["feature_genomic_start"]),
                feature_genomic_end=int(row["feature_genomic_end"]),
                feature_order_transcript=(
                    _i_or_none(row.get("feature_order_transcript", "0")) or 0
                ),
                overlaps_domain=row.get("overlaps_domain", "NA") or "NA",
                plot_group=row.get("plot_group", "") or "",
                is_mane_select=row.get("is_mane_select", "NA") or "NA",
                is_ensembl_canonical=row.get("is_ensembl_canonical", "NA") or "NA",
            )
            by_id[seg.input_id].append(seg)
    return by_id


# --------------------------------------------------------------------------- #
# Colors / heights
# --------------------------------------------------------------------------- #

PLOT_GROUP_COLORS = {
    "five_prime_UTR":     "#BDBDBD",
    "three_prime_UTR":    "#9E9E9E",
    "CDS":                "#1F77B4",
    "CDS_no_domain":      "#1F77B4",
    "CDS_domain":         "#D62728",
    "intron":             "#666666",
    "intron_domain_span": "#D62728",
}

# Height (in y units) for each feature kind. UTRs are slightly thinner than
# CDS, introns are a thin line.
def feature_height(feature_type: str) -> float:
    if feature_type in ("CDS",):
        return 0.6
    if feature_type in ("five_prime_UTR", "three_prime_UTR"):
        return 0.4
    return 0.02  # intron line


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #

def _title_for(segs: list[Segment]) -> str:
    s = segs[0]
    bits = []
    if s.gene_name: bits.append(s.gene_name)
    if s.protein_id: bits.append(s.protein_id)
    if s.transcript_id: bits.append(s.transcript_id)
    if s.domain_id: bits.append(f"domain={s.domain_id}")
    flags = []
    if s.is_mane_select == "true": flags.append("MANE_Select")
    if s.is_ensembl_canonical == "true": flags.append("Ensembl_canonical")
    if flags: bits.append("[" + ",".join(flags) + "]")
    bits.append(f"({s.chrom} {s.strand})")
    return "  ".join(bits)


def _draw_genomic(ax, segs: list[Segment], *, show_introns: bool,
                  show_utr: bool, highlight_domain: bool) -> None:
    """Draw on genomic-coordinate X axis (matches IGV)."""
    y_center = 0.5
    for s in segs:
        if s.feature_type == "intron" and not show_introns:
            continue
        if s.feature_type in ("five_prime_UTR", "three_prime_UTR") and not show_utr:
            continue
        color = PLOT_GROUP_COLORS.get(s.plot_group, "#777777")
        if not highlight_domain:
            # Collapse the domain coloring to the plain feature color.
            if s.plot_group == "CDS_domain":
                color = PLOT_GROUP_COLORS["CDS"]
            elif s.plot_group == "intron_domain_span":
                color = PLOT_GROUP_COLORS["intron"]
        h = feature_height(s.feature_type)
        x = s.feature_genomic_start - 0.5  # center on integer base
        w = s.feature_genomic_end - s.feature_genomic_start + 1
        if s.feature_type == "intron":
            # Single thin horizontal line spanning the intron.
            ax.plot([x, x + w], [y_center, y_center], color=color, lw=1.2,
                    solid_capstyle="butt", zorder=1)
            continue
        rect = mpatches.Rectangle(
            (x, y_center - h / 2.0), w, h,
            facecolor=color, edgecolor="black", linewidth=0.4, zorder=2,
        )
        ax.add_patch(rect)

    # Strand arrow.
    xmin = min(s.feature_genomic_start for s in segs)
    xmax = max(s.feature_genomic_end for s in segs)
    arrow_y = y_center + 0.55
    pad = (xmax - xmin) * 0.04 if xmax > xmin else 1
    if (segs[0].strand or "+") == "+":
        ax.annotate("", xy=(xmax - pad, arrow_y), xytext=(xmin + pad, arrow_y),
                    arrowprops=dict(arrowstyle="->", color="#444444", lw=1))
        ax.text(xmin, arrow_y + 0.05, "5'", color="#444444", va="bottom")
        ax.text(xmax, arrow_y + 0.05, "3'", color="#444444", va="bottom",
                ha="right")
    else:
        ax.annotate("", xy=(xmin + pad, arrow_y), xytext=(xmax - pad, arrow_y),
                    arrowprops=dict(arrowstyle="->", color="#444444", lw=1))
        ax.text(xmax, arrow_y + 0.05, "5'", color="#444444", va="bottom",
                ha="right")
        ax.text(xmin, arrow_y + 0.05, "3'", color="#444444", va="bottom")

    ax.set_xlim(xmin - (xmax - xmin) * 0.02 if xmax > xmin else xmin - 1,
                xmax + (xmax - xmin) * 0.02 if xmax > xmin else xmax + 1)
    ax.set_ylim(0, 1.6)
    ax.set_yticks([])
    ax.set_xlabel(f"genomic position ({segs[0].chrom})")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)


def _draw_spliced(ax, segs: list[Segment], *, show_utr: bool,
                  highlight_domain: bool) -> None:
    """Concatenate non-intron features in translation order. Useful when
    introns dwarf the rest of the figure."""
    kept = [s for s in segs if s.feature_type != "intron"]
    if not show_utr:
        kept = [s for s in kept if s.feature_type not in
                ("five_prime_UTR", "three_prime_UTR")]
    kept = sorted(kept, key=lambda s: s.feature_order_transcript)

    y_center = 0.5
    cursor = 0.0
    for s in kept:
        color = PLOT_GROUP_COLORS.get(s.plot_group, "#777777")
        if not highlight_domain and s.plot_group == "CDS_domain":
            color = PLOT_GROUP_COLORS["CDS"]
        w = s.feature_genomic_end - s.feature_genomic_start + 1
        h = feature_height(s.feature_type)
        rect = mpatches.Rectangle(
            (cursor, y_center - h / 2.0), w, h,
            facecolor=color, edgecolor="black", linewidth=0.4,
        )
        ax.add_patch(rect)
        cursor += w
    ax.set_xlim(0, cursor)
    ax.set_ylim(0, 1.6)
    ax.set_yticks([])
    ax.set_xlabel("spliced position (nt, translation order)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)


def _legend(ax, *, show_introns: bool, show_utr: bool,
            highlight_domain: bool) -> None:
    handles = []
    if show_utr:
        handles.append(mpatches.Patch(facecolor=PLOT_GROUP_COLORS["five_prime_UTR"],
                                      edgecolor="black", label="UTR"))
    handles.append(mpatches.Patch(facecolor=PLOT_GROUP_COLORS["CDS"],
                                  edgecolor="black", label="CDS"))
    if highlight_domain:
        handles.append(mpatches.Patch(facecolor=PLOT_GROUP_COLORS["CDS_domain"],
                                      edgecolor="black", label="CDS (domain)"))
    if show_introns:
        handles.append(mpatches.Patch(facecolor=PLOT_GROUP_COLORS["intron"],
                                      label="intron"))
    ax.legend(handles=handles, loc="upper center",
              bbox_to_anchor=(0.5, -0.25), ncol=len(handles), frameon=False)


def render_one(segs: list[Segment], *, out: str | None = None,
               width: float = 12.0, height: float = 2.2,
               title: str | None = None,
               show_introns: bool = True,
               show_utr: bool = True,
               highlight_domain: bool = True,
               spliced: bool = False,
               fig=None):
    """Render one query. If `fig` is provided we draw into it (PdfPages)."""
    own_fig = fig is None
    if own_fig:
        fig, ax = plt.subplots(figsize=(width, height))
    else:
        ax = fig.add_subplot(1, 1, 1)
    if spliced:
        _draw_spliced(ax, segs, show_utr=show_utr,
                      highlight_domain=highlight_domain)
    else:
        _draw_genomic(ax, segs, show_introns=show_introns, show_utr=show_utr,
                      highlight_domain=highlight_domain)
    _legend(ax, show_introns=show_introns and not spliced, show_utr=show_utr,
            highlight_domain=highlight_domain)
    ax.set_title(title or _title_for(segs))
    fig.tight_layout()
    if own_fig and out:
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out}", file=sys.stderr)
    return fig


# --------------------------------------------------------------------------- #
# Optional plotly HTML export
# --------------------------------------------------------------------------- #

def render_html(segs: list[Segment], out: str, *, highlight_domain: bool = True,
                show_introns: bool = True, show_utr: bool = True) -> None:
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise SystemExit(
            "plotly is required for --html output. Install with: pip install plotly"
        )
    fig = go.Figure()
    y = 0
    for s in segs:
        if s.feature_type == "intron" and not show_introns:
            continue
        if s.feature_type in ("five_prime_UTR", "three_prime_UTR") and not show_utr:
            continue
        color = PLOT_GROUP_COLORS.get(s.plot_group, "#777777")
        if not highlight_domain and s.plot_group == "CDS_domain":
            color = PLOT_GROUP_COLORS["CDS"]
        x0, x1 = s.feature_genomic_start - 0.5, s.feature_genomic_end + 0.5
        h = feature_height(s.feature_type)
        if s.feature_type == "intron":
            fig.add_shape(type="line", x0=x0, x1=x1, y0=y, y1=y,
                          line=dict(color=color, width=2))
        else:
            fig.add_shape(type="rect", x0=x0, x1=x1, y0=y - h / 2, y1=y + h / 2,
                          fillcolor=color, line=dict(color="black", width=0.5))
            fig.add_trace(go.Scatter(
                x=[(x0 + x1) / 2], y=[y],
                mode="markers",
                marker=dict(size=1, color=color),
                hovertemplate=(
                    f"{s.feature_id} ({s.feature_type})<br>"
                    f"{s.chrom}:{s.feature_genomic_start}-{s.feature_genomic_end}"
                    f"<br>plot_group={s.plot_group}<extra></extra>"
                ),
                showlegend=False,
            ))
    fig.update_layout(
        title=_title_for(segs),
        xaxis_title=f"genomic position ({segs[0].chrom})",
        yaxis=dict(visible=False, range=[-0.6, 0.6]),
        height=260, margin=dict(l=40, r=40, t=60, b=40),
    )
    fig.write_html(out, include_plotlyjs="cdn")
    print(f"Wrote {out}", file=sys.stderr)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="protein2genomic plot",
        description="Render isoform_structure.tsv to PDF/PNG/SVG/HTML.",
    )
    p.add_argument("--isoform", required=True,
                   help="Path to isoform_structure.tsv (the plot-ready table).")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--input-id", help="A single input_id to render.")
    g.add_argument("--all", action="store_true",
                   help="Render every input_id (multipage PDF if --out ends in .pdf).")
    p.add_argument("--out", help="Output file (.pdf/.png/.svg).")
    p.add_argument("--html", help="Optional interactive HTML output (plotly).")
    p.add_argument("--title", help="Override the figure title.")
    p.add_argument("--width", type=float, default=12.0, help="Figure width in inches.")
    p.add_argument("--height", type=float, default=2.2, help="Figure height in inches.")
    p.add_argument("--no-highlight", dest="highlight", action="store_false",
                   help="Do not color CDS_domain segments differently.")
    p.add_argument("--no-introns", dest="introns", action="store_false",
                   help="Hide intron lines.")
    p.add_argument("--no-utr", dest="utr", action="store_false",
                   help="Hide UTR boxes.")
    p.add_argument("--spliced", action="store_true",
                   help="Concatenate non-intron features in translation order.")
    p.set_defaults(highlight=True, introns=True, utr=True)
    return p


def main(argv: list[str] | None = None) -> int:
    args = _argparser().parse_args(argv)
    by_id = load_isoform_tsv(args.isoform)
    if not by_id:
        print(f"No rows in {args.isoform}", file=sys.stderr)
        return 1

    if args.input_id:
        if args.input_id not in by_id:
            keys = ", ".join(sorted(by_id)[:10])
            print(f"input_id {args.input_id!r} not found. Sample: {keys}",
                  file=sys.stderr)
            return 2
        ids = [args.input_id]
    else:
        ids = list(by_id.keys())

    if args.out:
        if args.out.lower().endswith(".pdf") and len(ids) > 1:
            with PdfPages(args.out) as pdf:
                for i in ids:
                    fig = plt.figure(figsize=(args.width, args.height))
                    render_one(by_id[i], fig=fig,
                               width=args.width, height=args.height,
                               title=(args.title if len(ids) == 1 else None),
                               show_introns=args.introns, show_utr=args.utr,
                               highlight_domain=args.highlight,
                               spliced=args.spliced)
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)
            print(f"Wrote {args.out} ({len(ids)} pages)", file=sys.stderr)
        else:
            if len(ids) > 1:
                # Non-PDF can only hold one figure; emit one file per id.
                base, ext = os.path.splitext(args.out)
                for i in ids:
                    out_i = f"{base}.{i}{ext}"
                    render_one(by_id[i], out=out_i,
                               width=args.width, height=args.height,
                               title=args.title,
                               show_introns=args.introns, show_utr=args.utr,
                               highlight_domain=args.highlight,
                               spliced=args.spliced)
            else:
                render_one(by_id[ids[0]], out=args.out,
                           width=args.width, height=args.height,
                           title=args.title,
                           show_introns=args.introns, show_utr=args.utr,
                           highlight_domain=args.highlight,
                           spliced=args.spliced)

    if args.html:
        if len(ids) > 1:
            base, ext = os.path.splitext(args.html)
            for i in ids:
                render_html(by_id[i], f"{base}.{i}{ext}",
                            highlight_domain=args.highlight,
                            show_introns=args.introns, show_utr=args.utr)
        else:
            render_html(by_id[ids[0]], args.html,
                        highlight_domain=args.highlight,
                        show_introns=args.introns, show_utr=args.utr)

    if not args.out and not args.html:
        print("Nothing to do: pass --out and/or --html.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Convert UniProt feature annotations into the BED-like format prot2exon eats.

Two input formats are supported, autodetected by extension or via `--format`:

* `.dat` / `.txt` / `.gz` — the classic UniProtKB **flat-file** with `FT`
  lines:

      FT   DOMAIN          93..312
      FT                   /note="Protein kinase"
      FT                   /evidence="ECO:0000255|PROSITE-ProRule:PRU00159"
      FT   ZN_FING         50..73
      FT                   /note="C4-type"

* `.json` — UniProtKB REST output (e.g.
  `https://rest.uniprot.org/uniprotkb/P04637.json`, or a JSONLines dump). We
  read the `features` array. Each feature has `type`, `location`,
  `description` and optionally `featureId`.

The `FT` flat-file format is the one most papers and most pipelines that
predate the REST API still use, so the parser handles it natively.

By default the parser keeps **DOMAIN / REPEAT / ZN_FING / DNA_BIND /
COILED / TRANSMEM / REGION** features (configurable via `--feature-types`).

Mapping UniProt → ENSP is handled by `_mapping.UniProtToEnsp`. UniProt entries
typically carry their own Ensembl cross-references in `DR Ensembl;` lines —
the parser will use those *first* (no external mapping needed) and only fall
back to the Ensembl xref TSV / `--simple-mapping` for entries that don't
carry an Ensembl xref.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
from _mapping import (UniProtToEnsp, looks_like_ensp, looks_like_uniprot,
                      open_text, dedup_and_sort_rows, write_bed_row,
                      strip_version)


DEFAULT_FEATURE_TYPES = {
    "DOMAIN", "REPEAT", "REGION", "ZN_FING", "DNA_BIND",
    "COILED", "TRANSMEM", "MOTIF", "TOPO_DOM",
}


# --------------------------------------------------------------------------- #
# Flat-file (.dat) parser
# --------------------------------------------------------------------------- #

# Match: "FT   DOMAIN          93..312" or "FT   ZN_FING         50..73"
# location can also be "?..312", "<1..312", "312..>500" or "?" — we accept the
# strict integer form and drop the rest with a counter.
FT_HEADER_RE = re.compile(
    r"^FT\s+(\w+)\s+(?P<start><?\??\d+)\.\.(?P<end>>?\??\d+)"
)
FT_NOTE_RE = re.compile(r'/note="([^"]+)"')
FT_ID_RE   = re.compile(r'/id="([^"]+)"')
DR_ENSEMBL_RE = re.compile(r"^DR\s+Ensembl;\s+([^;\s]+);\s+([^;\s]+);\s+([^.\s]+)")
AC_RE   = re.compile(r"^AC\s+([^;]+);")


def parse_dat(path: str, *, feature_types: set[str]
              ) -> tuple[list, dict[str, list[str]], list]:
    """Yield (rows, dr_ensp_by_acc, rejected) for a UniProt .dat file.

    rows: list of (uniprot_acc, aa_start, aa_end, domain_id, description)
    dr_ensp_by_acc: per-accession list of ENSPs found in DR Ensembl lines
    rejected: list of (line, reason) tuples
    """
    rows: list[tuple[str, int, int, str, str]] = []
    dr_by_acc: dict[str, list[str]] = {}
    rejected: list[tuple[str, str]] = []

    cur_acc: str | None = None
    cur_dr: list[str] = []
    pending: dict | None = None  # currently-open FT feature awaiting /note

    def flush_pending():
        nonlocal pending
        if not pending or cur_acc is None:
            pending = None
            return
        rows.append((cur_acc, pending["s"], pending["e"], pending["domain"],
                     pending["desc"]))
        pending = None

    with open_text(path) as fh:
        for line in fh:
            if line.startswith("ID "):
                # New entry — flush state
                flush_pending()
                if cur_acc is not None:
                    dr_by_acc[cur_acc] = cur_dr
                cur_acc = None
                cur_dr = []
                continue
            if line.startswith("AC "):
                m = AC_RE.match(line)
                if m and cur_acc is None:
                    cur_acc = m.group(1).strip()
                continue
            if line.startswith("DR "):
                m = DR_ENSEMBL_RE.match(line)
                if m:
                    enst, ensp, ensg = m.group(1), m.group(2), m.group(3)
                    cur_dr.append(strip_version(ensp))
                continue
            if line.startswith("//"):
                flush_pending()
                if cur_acc is not None:
                    dr_by_acc[cur_acc] = list(cur_dr)
                cur_acc, cur_dr = None, []
                continue
            if not line.startswith("FT"):
                continue
            # A new FT header line resets `pending`.
            m = FT_HEADER_RE.match(line)
            if m:
                flush_pending()
                ftype = re.match(r"^FT\s+(\w+)", line).group(1)
                if ftype not in feature_types:
                    pending = None
                    continue
                s_str = m.group("start").lstrip("<>?")
                e_str = m.group("end").lstrip("<>?")
                if not s_str.isdigit() or not e_str.isdigit():
                    rejected.append((line.rstrip("\n"),
                                     "fuzzy location, skipped"))
                    pending = None
                    continue
                s, e = int(s_str), int(e_str)
                if e < s:
                    rejected.append((line.rstrip("\n"), "end<start"))
                    pending = None
                    continue
                pending = {"s": s, "e": e, "type": ftype,
                           "domain": ftype, "desc": ""}
                continue
            # Continuation lines: pull /note=, /id= if present.
            if pending is not None:
                m_note = FT_NOTE_RE.search(line)
                if m_note:
                    pending["desc"] = re.sub(r"\s+", " ", m_note.group(1)).strip()
                    # Use the note as the domain_id when it's short and ASCII
                    if pending["desc"] and len(pending["desc"]) <= 40:
                        pending["domain"] = (
                            pending["type"] + "_" +
                            re.sub(r"[^A-Za-z0-9_.-]+", "_", pending["desc"]))
                m_id = FT_ID_RE.search(line)
                if m_id:
                    # Some features have a stable /id (e.g. Zn finger
                    # accession) — prefer that as domain_id.
                    pending["domain"] = pending["type"] + "_" + m_id.group(1)

        # EOF
        flush_pending()
        if cur_acc is not None:
            dr_by_acc[cur_acc] = list(cur_dr)

    return rows, dr_by_acc, rejected


# --------------------------------------------------------------------------- #
# JSON parser (REST or JSONLines)
# --------------------------------------------------------------------------- #

def parse_json(path: str, *, feature_types: set[str]
               ) -> tuple[list, dict[str, list[str]], list]:
    rows: list[tuple[str, int, int, str, str]] = []
    dr_by_acc: dict[str, list[str]] = {}
    rejected: list[tuple[str, str]] = []

    def consume_entry(entry: dict):
        acc = entry.get("primaryAccession", "")
        if not acc:
            return
        # DR Ensembl
        ensps: list[str] = []
        for xref in (entry.get("uniProtKBCrossReferences", []) or []):
            if xref.get("database") == "Ensembl":
                for prop in (xref.get("properties", []) or []):
                    if prop.get("key") in ("ProteinId", "EnsemblProteinId"):
                        v = prop.get("value", "")
                        if v: ensps.append(strip_version(v))
        if ensps:
            dr_by_acc[acc] = ensps
        # Features
        for feat in (entry.get("features", []) or []):
            ftype = (feat.get("type") or "").upper().replace(" ", "_")
            if ftype not in feature_types and ftype.replace("-", "_") not in feature_types:
                continue
            loc = feat.get("location", {}) or {}
            s = (loc.get("start", {}) or {}).get("value")
            e = (loc.get("end", {}) or {}).get("value")
            if s is None or e is None or not isinstance(s, int) or not isinstance(e, int):
                rejected.append((json.dumps(feat), "fuzzy/missing location"))
                continue
            if e < s:
                rejected.append((json.dumps(feat), "end<start"))
                continue
            desc = (feat.get("description") or "").strip()
            domain_id = ftype
            fid = feat.get("featureId")
            if fid:
                domain_id = f"{ftype}_{fid}"
            elif desc:
                domain_id = (f"{ftype}_" +
                             re.sub(r"[^A-Za-z0-9_.-]+", "_", desc)[:40])
            rows.append((acc, s, e, domain_id,
                         re.sub(r"\s+", " ", desc).strip()))

    with open_text(path) as fh:
        # JSONLines first (one entry per line); fall back to a single JSON
        # document (single entry OR `{"results": [...]}`).
        first = fh.readline()
        fh.seek(0)
        if first.lstrip().startswith("{") and not first.lstrip().startswith('{"results"'):
            # Try JSONLines.
            ok = True
            try:
                json.loads(first)
            except json.JSONDecodeError:
                ok = False
            if ok:
                for line in fh:
                    line = line.strip()
                    if not line: continue
                    try:
                        consume_entry(json.loads(line))
                    except json.JSONDecodeError as e:
                        rejected.append((line[:120], f"bad json: {e}"))
                return rows, dr_by_acc, rejected
        # Single JSON document.
        doc = json.load(fh)
        if isinstance(doc, dict) and "results" in doc:
            for entry in doc["results"]:
                consume_entry(entry)
        elif isinstance(doc, list):
            for entry in doc:
                consume_entry(entry)
        elif isinstance(doc, dict):
            consume_entry(doc)
        else:
            rejected.append((str(type(doc)), "unrecognized JSON shape"))
    return rows, dr_by_acc, rejected


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="UniProt feature table (.dat or REST .json) → "
                    "prot2exon BED-like.")
    p.add_argument("--in", dest="in_path", required=True,
                   help="UniProt .dat / .dat.gz / .txt / .json file.")
    p.add_argument("--format", choices=("auto", "dat", "json"), default="auto")
    p.add_argument("--out", default="-")
    p.add_argument("--feature-types",
                   help=f"Comma-separated. Default: {','.join(sorted(DEFAULT_FEATURE_TYPES))}")
    p.add_argument("--mapping",
                   help="Ensembl UniProt xref TSV. Used as a fallback when an "
                        "entry has no DR Ensembl cross-reference.")
    p.add_argument("--simple-mapping",
                   help="Two-column TSV mapping (uniprot\\tensp). Alternative.")
    p.add_argument("--min-length", type=int, default=5)
    p.add_argument("--keep-unmapped",
                   help="Write rejected/unmapped rows here.")
    args = p.parse_args(argv)

    feat_types = (set(s.strip().upper() for s in args.feature_types.split(","))
                  if args.feature_types else DEFAULT_FEATURE_TYPES)

    fmt = args.format
    if fmt == "auto":
        low = args.in_path.lower().rstrip(".gz")
        fmt = "json" if low.endswith(".json") else "dat"
        print(f"[uniprot] auto-detected format = {fmt}", file=sys.stderr)

    if fmt == "dat":
        raw, dr_by_acc, rejected = parse_dat(args.in_path,
                                             feature_types=feat_types)
    else:
        raw, dr_by_acc, rejected = parse_json(args.in_path,
                                              feature_types=feat_types)

    # Optional fallback mapping for entries with no DR Ensembl.
    fallback = None
    if args.mapping:
        fallback = UniProtToEnsp.from_ensembl_xref_tsv(args.mapping)
    elif args.simple_mapping:
        fallback = UniProtToEnsp.from_simple_tsv(args.simple_mapping)

    out_rows: list[tuple[str, int, int, str, str]] = []
    unresolved = 0
    for acc, s, e, did, desc in raw:
        if (e - s + 1) < args.min_length:
            continue
        ensps = dr_by_acc.get(acc, [])
        if not ensps and fallback is not None:
            ensps = fallback.lookup(acc)
        if not ensps:
            unresolved += 1
            rejected.append((f"{acc} {s}-{e} {did}", "no ENSP"))
            continue
        for ensp in ensps:
            out_rows.append((ensp, s, e, did, desc))

    out_rows = dedup_and_sort_rows(out_rows)

    out_fh = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        out_fh.write("# prot2exon BED-like, generated from UniProt features\n")
        out_fh.write("# columns: ENSP\\taa_start\\taa_end\\tdomain_id\\tdescription\n")
        for r in out_rows:
            write_bed_row(out_fh, ensp=r[0], aa_start=r[1], aa_end=r[2],
                          domain_id=r[3], source=r[4])
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()
    print(f"[uniprot] wrote {len(out_rows)} rows ({unresolved} unresolved, "
          f"{len(rejected)} total rejected) to {args.out}", file=sys.stderr)
    if args.keep_unmapped and rejected:
        with open(args.keep_unmapped, "w") as f:
            f.write("raw\treason\n")
            for raw, why in rejected:
                f.write(f"{raw}\t{why}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

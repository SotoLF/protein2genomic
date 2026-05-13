"""Shared utilities for the prepare_from_* parsers.

Most of the heavy lifting these scripts do is converting **UniProt accessions**
(P04637, Q9Y6K1, ...) into **Ensembl protein IDs** (ENSP00000269305, ...)
that prot2exon indexes. This module loads a single mapping table and
exposes one lookup function.

Two issues you should know about going in
-----------------------------------------

1. A single UniProt accession can map to **multiple ENSPs** — different
   transcripts of the same gene. Pick a strategy (`--strategy`):

     all          emit one row per (uniprot, ensp) pair (default)
     canonical    keep only the MANE Select / Ensembl_canonical ENSP
                  (requires that information in the mapping file or
                  cross-referenced via the index)
     longest      keep the ENSP with the longest CDS (requires --index)

2. A single ENSP can map to **multiple UniProt accessions** — usually the
   canonical and TrEMBL forms of the same protein. The reverse direction is
   less common in our workflow but the loader keeps both edges.

Mapping file format
-------------------

Ensembl distributes per-release UniProt cross-reference files. The relevant
columns are:

    gene_stable_id  transcript_stable_id  protein_stable_id  xref  db_name

Where `db_name` is one of `Uniprot/SWISSPROT`, `Uniprot/SPTREMBL`,
`Uniprot_gn`, etc. We only consume rows with a SWISSPROT or SPTREMBL `db_name`
by default. Versions on either side are stripped on load.

Download (human, Ensembl 110+, ~20 MB gzipped):

    rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/tsv/homo_sapiens/Homo_sapiens.GRCh38.110.uniprot.tsv.gz .
    gunzip Homo_sapiens.GRCh38.110.uniprot.tsv.gz

Or via HTTP:
    https://ftp.ensembl.org/pub/release-110/tsv/homo_sapiens/

If you don't have / can't get this file, the `--api` flag in each parser will
query UniProt's REST ID-mapping endpoint (slow, rate-limited).
"""

from __future__ import annotations

import csv
import gzip
import re
import sys
from collections import defaultdict
from typing import Iterable


# UniProt accession regex: P04637, Q9Y6K1, A0A0B4J2A2, P04637-2, ...
UNIPROT_ACC_RE = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9](-\d+)?$|^[OPQ][0-9][A-Z0-9]{3}[0-9](-\d+)?$")
ENSP_RE = re.compile(r"^ENSP\d+(\.\d+)?$")
ENST_RE = re.compile(r"^ENST\d+(\.\d+)?$")


def strip_version(x: str) -> str:
    """Strip the .N version suffix Ensembl puts on its IDs."""
    return x.split(".", 1)[0]


def looks_like_uniprot(x: str) -> bool:
    return bool(UNIPROT_ACC_RE.match(x))


def looks_like_ensp(x: str) -> bool:
    return bool(ENSP_RE.match(x))


def looks_like_enst(x: str) -> bool:
    return bool(ENST_RE.match(x))


class UniProtToEnsp:
    """Bidirectional UniProt accession ↔ ENSP lookup."""

    def __init__(self) -> None:
        self.up_to_ensp: dict[str, set[str]] = defaultdict(set)
        self.ensp_to_up: dict[str, set[str]] = defaultdict(set)

    @classmethod
    def from_ensembl_xref_tsv(cls, path: str,
                              keep_dbs: Iterable[str] = ("Uniprot/SWISSPROT",
                                                          "Uniprot/SPTREMBL")) -> "UniProtToEnsp":
        """Load the Ensembl per-release UniProt xref TSV.

        Accepts both the modern multi-column TSV (with a header line) and
        older versions. We only look at columns named ``protein_stable_id``,
        ``xref`` and ``db_name``. Versioned IDs are stripped.
        """
        keep = set(keep_dbs)
        m = cls()
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            # Try header-driven first.
            sample = fh.readline()
            fh.seek(0)
            reader: csv.DictReader | None = None
            if sample.startswith("gene_stable_id") or sample.startswith("#"):
                reader = csv.DictReader(
                    (l for l in fh if not l.startswith("#")), delimiter="\t")
                for row in reader:
                    db = row.get("db_name", "")
                    if db not in keep:
                        continue
                    pid = strip_version(row.get("protein_stable_id", "") or "")
                    up = (row.get("xref", "") or "").strip()
                    if not pid or not up or not looks_like_uniprot(up):
                        continue
                    m.up_to_ensp[up].add(pid)
                    m.ensp_to_up[pid].add(up)
            else:
                # Headerless / unknown layout: best-effort column hunt.
                for line in fh:
                    if line.startswith("#"): continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5: continue
                    db = parts[4]
                    if db not in keep: continue
                    pid = strip_version(parts[2])
                    up = parts[3].strip()
                    if not pid or not up or not looks_like_uniprot(up):
                        continue
                    m.up_to_ensp[up].add(pid)
                    m.ensp_to_up[pid].add(up)
        print(
            f"[mapping] loaded {len(m.up_to_ensp)} UniProt accessions covering "
            f"{len(m.ensp_to_up)} ENSPs from {path}",
            file=sys.stderr,
        )
        return m

    @classmethod
    def from_simple_tsv(cls, path: str, *, up_col: int = 0,
                        ensp_col: int = 1) -> "UniProtToEnsp":
        """Fallback loader for a 2-column TSV (`uniprot\\tensp`). Stable
        across mapping-source choices."""
        m = cls()
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("#") or line.strip() == "":
                    continue
                parts = line.rstrip("\n").split("\t")
                if max(up_col, ensp_col) >= len(parts):
                    continue
                up = parts[up_col].strip()
                pid = strip_version(parts[ensp_col].strip())
                if not looks_like_uniprot(up) or not looks_like_ensp(pid):
                    continue
                m.up_to_ensp[up].add(pid)
                m.ensp_to_up[pid].add(up)
        print(
            f"[mapping] loaded {len(m.up_to_ensp)} UniProt accessions covering "
            f"{len(m.ensp_to_up)} ENSPs from {path}",
            file=sys.stderr,
        )
        return m

    def lookup(self, up: str) -> list[str]:
        """Return all ENSPs for a given UniProt accession. Also tries the
        non-isoform form (`P04637-2` → `P04637`)."""
        candidates = [up]
        if "-" in up:
            candidates.append(up.split("-", 1)[0])
        for c in candidates:
            if c in self.up_to_ensp:
                return sorted(self.up_to_ensp[c])
        return []


def open_text(path: str):
    """Open a possibly-gzipped text file as a context manager."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def write_bed_row(out, *, ensp: str, aa_start: int, aa_end: int,
                  domain_id: str, source: str = "") -> None:
    """Emit one BED-like row in the format prot2exon expects.

    Columns: ENSP, aa_start, aa_end, domain_id, source (free-form 5th col)
    """
    fields = [ensp, str(aa_start), str(aa_end), domain_id]
    if source:
        fields.append(source)
    out.write("\t".join(fields) + "\n")


def dedup_and_sort_rows(rows: list[tuple[str, int, int, str, str]]
                       ) -> list[tuple[str, int, int, str, str]]:
    """Remove exact duplicates (same ENSP, aa range, domain_id) and sort by
    (ENSP, aa_start, aa_end) for stable output."""
    seen = set()
    out = []
    for r in rows:
        key = (r[0], r[1], r[2], r[3])
        if key in seen: continue
        seen.add(key)
        out.append(r)
    return sorted(out, key=lambda r: (r[0], r[1], r[2], r[3]))

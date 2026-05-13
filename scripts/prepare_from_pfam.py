#!/usr/bin/env python3
"""Convert HMMER domtblout (Pfam scan) to the BED-like format prot2exon eats.

HMMER's `--domtblout` is the right file: one row per *domain hit* (not per
sequence hit), so we get per-domain coordinates. Format reference:
http://eddylab.org/software/hmmer/Userguide.pdf §10.

Important columns (space-separated, ragged):

    1   target name           (depends on direction, see below)
    2   target accession
    3   tlen                  (target length)
    4   query name
    5   query accession
    6   qlen                  (query length)
    7   E-value (full)        full-sequence E-value
    8   score (full)          full-sequence bit score
    9   bias (full)
    10  #                     this-domain index
    11  of                    total domains for this protein
    12  c-Evalue              conditional E-value of this domain
    13  i-Evalue              independent E-value of this domain
    14  score (this)          this-domain bit score
    15  bias (this)
    16  hmm from
    17  hmm to
    18  ali from              ← we use these
    19  ali to                ← we use these
    20  env from
    21  env to
    22  acc
    23+ description (rest of the line)

Direction
---------

If you ran `hmmscan` with the Pfam HMM database as the database (typical):
    hmmscan --domtblout out.dom Pfam-A.hmm proteins.fa
then **target = HMM (Pfam family), query = your protein**. Pass `--mode scan`
(this is the default).

If you ran `hmmsearch` with your protein FASTA as the database:
    hmmsearch --domtblout out.dom one_hmm.hmm proteins.fa
then **target = your protein, query = HMM**. Pass `--mode search`.

We compute the protein-coordinate range from `ali from / ali to` of the
**protein** column (not the HMM). hmm_from / hmm_to are HMM-coordinate
positions and would be wrong here.

ID mapping
----------

If your protein FASTA used ENSPs already (the GENCODE-derived FASTA is the
typical case), no mapping is needed — `--id-type ensp`. If you used UniProt
accessions, supply `--mapping` (Ensembl xref TSV) or `--simple-mapping`.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
from _mapping import (UniProtToEnsp, looks_like_ensp, looks_like_uniprot,
                      open_text, dedup_and_sort_rows, write_bed_row,
                      strip_version)


def parse(path: str, *, mode: str, min_score: float, min_length: int,
          id_type: str, mapping: UniProtToEnsp | None
          ) -> tuple[list, list]:
    rows: list[tuple[str, int, int, str, str]] = []
    rejected: list[tuple[str, str]] = []

    if id_type == "auto":
        # Sniff the first non-comment line.
        samples: list[str] = []
        with open_text(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip(): continue
                fields = line.split()
                protein_col = 3 if mode == "scan" else 0
                if len(fields) > protein_col:
                    samples.append(fields[protein_col])
                if len(samples) >= 50: break
        n_ensp = sum(1 for s in samples if looks_like_ensp(s))
        n_up = sum(1 for s in samples if looks_like_uniprot(s))
        id_type = "ensp" if n_ensp >= max(n_up, 1) else (
                  "uniprot" if n_up else "unknown")
        print(f"[pfam] auto-detected id_type = {id_type}", file=sys.stderr)

    if id_type == "uniprot" and mapping is None:
        raise SystemExit(
            "domtblout uses UniProt accessions but no --mapping given.")

    with open_text(path) as fh:
        for ln, line in enumerate(fh, start=1):
            if line.startswith("#") or not line.strip():
                continue
            # domtblout is whitespace-separated with ragged trailing
            # description. The first 22 fields are positional.
            fields = re.split(r"\s+", line.rstrip("\n"), maxsplit=22)
            if len(fields) < 22:
                rejected.append((line.rstrip("\n"), "too few columns"))
                continue
            if mode == "scan":
                hmm_name = fields[0]
                hmm_acc  = fields[1]
                protein  = fields[3]
            else:
                protein  = fields[0]
                hmm_name = fields[3]
                hmm_acc  = fields[4]
            try:
                this_score = float(fields[13])
                ali_from   = int(fields[17])
                ali_to     = int(fields[18])
            except (ValueError, IndexError):
                rejected.append((line.rstrip("\n"), "unparseable numerics"))
                continue
            description = fields[22] if len(fields) > 22 else ""

            if this_score < min_score:
                continue
            if ali_to < ali_from:
                rejected.append((line.rstrip("\n"), "ali_to<ali_from"))
                continue
            if (ali_to - ali_from + 1) < min_length:
                continue

            # Resolve ID.
            ensps: list[str]
            if id_type == "ensp":
                ensps = [strip_version(protein)] if looks_like_ensp(protein) else []
            elif id_type == "uniprot":
                ensps = mapping.lookup(protein) if mapping else []
            else:
                ensps = ([strip_version(protein)] if looks_like_ensp(protein)
                         else (mapping.lookup(protein) if mapping and looks_like_uniprot(protein) else []))
            if not ensps:
                rejected.append((line.rstrip("\n"), f"no ENSP for {protein}"))
                continue

            # Domain id: prefer Pfam accession (with version stripped), then
            # name. We always tag with `_Pfam` so it's clear in the BED.
            pfam_acc = strip_version(hmm_acc) if hmm_acc and hmm_acc != "-" else ""
            tag = pfam_acc or hmm_name or "Pfam_hit"
            domain_id = f"{hmm_name}_{tag}" if (hmm_name and pfam_acc) else f"{tag}_Pfam"
            desc_safe = re.sub(r"\s+", " ", description).strip()

            for ensp in ensps:
                rows.append((ensp, ali_from, ali_to, domain_id, desc_safe))

    return rows, rejected


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="HMMER domtblout (Pfam scan) → prot2exon BED-like.")
    p.add_argument("--in", dest="in_path", required=True,
                   help="HMMER --domtblout file.")
    p.add_argument("--out", default="-")
    p.add_argument("--mode", choices=("scan", "search"), default="scan",
                   help="hmmscan: target=HMM, query=protein (default). "
                        "hmmsearch: target=protein, query=HMM.")
    p.add_argument("--min-score", type=float, default=0.0,
                   help="Drop hits below this per-domain bit score.")
    p.add_argument("--min-length", type=int, default=5)
    p.add_argument("--id-type", choices=("auto", "ensp", "uniprot"), default="auto")
    p.add_argument("--mapping",
                   help="Ensembl UniProt xref TSV (only required when protein IDs are UniProt).")
    p.add_argument("--simple-mapping",
                   help="Two-column TSV mapping (uniprot\\tensp). Alternative.")
    p.add_argument("--keep-unmapped",
                   help="Write rejected rows to this file.")
    args = p.parse_args(argv)

    mapping = None
    if args.mapping:
        mapping = UniProtToEnsp.from_ensembl_xref_tsv(args.mapping)
    elif args.simple_mapping:
        mapping = UniProtToEnsp.from_simple_tsv(args.simple_mapping)

    rows, rejected = parse(args.in_path,
                           mode=args.mode,
                           min_score=args.min_score,
                           min_length=args.min_length,
                           id_type=args.id_type,
                           mapping=mapping)
    rows = dedup_and_sort_rows(rows)

    out_fh = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        out_fh.write("# prot2exon BED-like, generated from HMMER domtblout (Pfam)\n")
        out_fh.write("# columns: ENSP\\taa_start\\taa_end\\tdomain_id\\tdescription\n")
        for r in rows:
            write_bed_row(out_fh, ensp=r[0], aa_start=r[1], aa_end=r[2],
                          domain_id=r[3], source=r[4])
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()
    print(f"[pfam] wrote {len(rows)} rows "
          f"({len(rejected)} rejected) to {args.out}", file=sys.stderr)
    if args.keep_unmapped and rejected:
        with open(args.keep_unmapped, "w") as f:
            f.write("raw_row\treason\n")
            for raw, why in rejected:
                f.write(f"{raw}\t{why}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

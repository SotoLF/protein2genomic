#!/usr/bin/env python3
"""Generate synthetic GTFs and BEDs used by run_tests.py.

We build two GTFs so we can test both "GTF carries tag attributes" and
"GTF carries none" without changing the rest of the schema.

Transcript inventory (with_tags.gtf)
====================================

ENST1 / ENSP1 — plus strand, MANE_Select + Ensembl_canonical, 3 exons, clean CDS
    chrA  +   exon1 100..149   (50 nt, 5'UTR + CDS_1)
                                  5'UTR = 100..119   (20 nt)
                                  CDS_1 = 120..149   (30 nt → aa 1..10)
                exon2 200..250   (51 nt, CDS_2)
                                  CDS_2 = 200..250   (51 nt → aa 11..27)
                exon3 300..400   (101 nt, CDS_3 + 3'UTR)
                                  CDS_3 = 300..350   (51 nt → aa 28..44)
                                  3'UTR = 351..400   (50 nt)
    total CDS = 132 nt = 44 aa   (divisible by 3 — no mismatch)

ENST2 / ENSP2 — minus strand, NO MANE/canonical (but other transcripts in the
    same GTF carry tags, so the GTF has tags overall → is_mane_select = false,
    not NA). 3 exons. Translation order: high genomic → low genomic.
    chrB  -   exon3 (5' in translation = genomic 800..850), 51 nt
                                  5'UTR = 851..??  none (exon ends at 850)
                                  Actually we want some UTR: pad exon to 800..860
    Final layout (genomic):
                exon1_genomic 400..450   = exon_3 in translation = 3'UTR + CDS_3
                exon2_genomic 600..650   = exon_2 in translation = CDS_2
                exon3_genomic 800..860   = exon_1 in translation = CDS_1 + 5'UTR
                CDS_3 (translation, genomic = last)  = 400..429 (30 nt aa 35..44)
                CDS_2 (translation, genomic = middle)= 600..650 (51 nt aa 18..34)
                CDS_1 (translation, genomic = first 5')= 800..850 (51 nt aa 1..17)
                5'UTR genomic 851..860   (10 nt, on '-' strand it's 5' UTR)
                3'UTR genomic 430..450 ? no — 3'UTR genomic 400..429 is CDS,
                so 3'UTR is on the other side: actually the 3'UTR sits 5' of
                CDS in genomic coords for '-' strand transcripts.
    For simplicity:
                exon at 400..450 contains 3'UTR 400..429 (30 nt) + CDS 430..450 (21 nt)
                exon at 600..650 is all CDS (51 nt)
                exon at 800..860 contains CDS 800..850 (51 nt) + 5'UTR 851..860 (10 nt)
                Total CDS = 21+51+51 = 123 → aa = 41, divisible by 3.

ENST3 / ENSP3 — selenoprotein-like, plus strand, CDS%3 != 0. Single exon for
    simplicity. CDS length 25 nt → remainder 1.
    chrC  +   exon1 1000..1049 (50 nt)
                CDS  1000..1024 (25 nt) — aa 1..8 (and an extra nt)
                3'UTR 1025..1049 (25 nt)

ENST4 / ENSP4 — codon split 1+2 (first base of codon 2 sits in CDS_1, last 2
    bases of codon 2 sit in CDS_2). Plus strand.
    chrD  +   exon1 100..103 (4 nt, all CDS)
                CDS_1  = 100..103   nt 1..4 of CDS, encodes aa 1 and base 1 of aa 2
              intron 104..199
              exon2 200..204 (5 nt, all CDS)
                CDS_2  = 200..204   nt 5..9 of CDS, encodes bases 2..3 of aa 2 and aa 3
                Total CDS = 9 nt = 3 aa (no mismatch).

ENST5 / ENSP5 — codon split 2+1 (first 2 bases of codon 2 in CDS_1, last base
    in CDS_2). Plus strand.
    chrE  +   exon1 100..104 (5 nt, all CDS)
                CDS_1 = 100..104   nt 1..5 — aa 1 + bases 1..2 of aa 2
              intron 105..199
              exon2 200..203 (4 nt, all CDS)
                CDS_2 = 200..203   nt 6..9 — base 3 of aa 2 + aa 3
                Total CDS = 9 nt = 3 aa.

ENST6 — non-coding (no protein_id, no CDS records). One exon. Used to test
    that an ENST query against a non-coding transcript yields
    no_CDS_for_protein.
    chrF  +   exon1 500..600 (101 nt)

Transcript inventory (no_tags.gtf)
==================================

ENST7 / ENSP7 — plus strand, deliberately no `tag "..."` attribute on any
    line. Same shape as ENST1. Used to verify that is_mane_select / is_canonical
    are reported as NA in this case.
"""

import os
import sys

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def attrs(**kw):
    """Format GTF attributes column. Empty values are skipped."""
    parts = []
    # Preserve a stable order so the GTF diffs cleanly.
    for k in ("gene_id", "transcript_id", "gene_type", "gene_name",
              "transcript_type", "transcript_name", "exon_number",
              "exon_id", "protein_id"):
        if k in kw and kw[k] is not None:
            parts.append(f'{k} "{kw[k]}"')
    # tags last (we may emit multiple)
    for t in kw.get("tags", []):
        parts.append(f'tag "{t}"')
    return "; ".join(parts) + ";"


def row(chrom, src, feat, start, end, strand, phase, attrs_str):
    return f"{chrom}\t{src}\t{feat}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs_str}"


def emit_simple(rows, *, chrom, strand, gid, tid, pid, gname,
                exons, cds_intervals, tags=None):
    """Write a transcript with the given exon and CDS intervals. `exons` and
    `cds_intervals` are lists of (start, end) in genomic order. exon_number is
    assigned in translation order (matching GENCODE)."""
    tags = tags or []
    n_exons = len(exons)
    # Pseudo transcript line (we don't strictly need it, but real GTFs have it
    # and we want the tag-detection check to cover the case where tags ride on
    # exon/CDS).
    rows.append(row(chrom, "p2g_test", "transcript",
                    min(s for s, _ in exons), max(e for _, e in exons),
                    strand, 0,
                    attrs(gene_id=gid, transcript_id=tid, gene_type="protein_coding",
                          gene_name=gname, transcript_type="protein_coding",
                          protein_id=pid, tags=tags)))
    # Exon numbers in translation order. On '+' strand: 1..N along genomic
    # ascending. On '-' strand: 1..N along genomic descending.
    sorted_exons = sorted(exons, key=lambda se: se[0])
    if strand == "-":
        ordered = list(reversed(sorted_exons))
    else:
        ordered = sorted_exons
    exon_num_by_pair = {pair: i + 1 for i, pair in enumerate(ordered)}

    for (s, e) in sorted_exons:
        rows.append(row(chrom, "p2g_test", "exon", s, e, strand, 0,
                        attrs(gene_id=gid, transcript_id=tid,
                              gene_name=gname, transcript_type="protein_coding",
                              exon_number=exon_num_by_pair[(s, e)],
                              exon_id=f"{tid}_E{exon_num_by_pair[(s, e)]}",
                              protein_id=pid, tags=tags)))
    # CDS rows. exon_number is the enclosing exon's number.
    for (cs, ce) in cds_intervals:
        for (es, ee) in sorted_exons:
            if cs >= es and ce <= ee:
                ex_num = exon_num_by_pair[(es, ee)]
                break
        else:
            raise SystemExit(f"CDS {cs}-{ce} not contained in any exon for {tid}")
        rows.append(row(chrom, "p2g_test", "CDS", cs, ce, strand, 0,
                        attrs(gene_id=gid, transcript_id=tid,
                              gene_name=gname, transcript_type="protein_coding",
                              exon_number=ex_num,
                              protein_id=pid, tags=tags)))


def emit_noncoding(rows, *, chrom, strand, gid, tid, gname, exons):
    rows.append(row(chrom, "p2g_test", "transcript",
                    min(s for s, _ in exons), max(e for _, e in exons),
                    strand, 0,
                    attrs(gene_id=gid, transcript_id=tid, gene_type="lncRNA",
                          gene_name=gname, transcript_type="lncRNA")))
    for i, (s, e) in enumerate(sorted(exons, key=lambda se: se[0]), start=1):
        rows.append(row(chrom, "p2g_test", "exon", s, e, strand, 0,
                        attrs(gene_id=gid, transcript_id=tid, gene_name=gname,
                              transcript_type="lncRNA",
                              exon_number=i, exon_id=f"{tid}_E{i}")))


def build_with_tags(path):
    rows = []

    # ENST1 — plus, MANE+canonical, 3 exons.
    emit_simple(rows,
        chrom="chrA", strand="+",
        gid="ENSG1", tid="ENST1", pid="ENSP1", gname="TEST1",
        exons=[(100, 149), (200, 250), (300, 400)],
        cds_intervals=[(120, 149), (200, 250), (300, 350)],
        tags=["MANE_Select", "Ensembl_canonical"])

    # ENST2 — minus, no tags carried by this transcript but tags exist in GTF.
    emit_simple(rows,
        chrom="chrB", strand="-",
        gid="ENSG2", tid="ENST2", pid="ENSP2", gname="TEST2",
        exons=[(400, 450), (600, 650), (800, 860)],
        cds_intervals=[(430, 450), (600, 650), (800, 850)])

    # ENST3 — selenoprotein-like, CDS%3 != 0 (25 nt remainder 1)
    emit_simple(rows,
        chrom="chrC", strand="+",
        gid="ENSG3", tid="ENST3", pid="ENSP3", gname="TEST3",
        exons=[(1000, 1049)],
        cds_intervals=[(1000, 1024)])

    # ENST4 — codon split 1+2
    emit_simple(rows,
        chrom="chrD", strand="+",
        gid="ENSG4", tid="ENST4", pid="ENSP4", gname="TEST4",
        exons=[(100, 103), (200, 204)],
        cds_intervals=[(100, 103), (200, 204)])

    # ENST5 — codon split 2+1
    emit_simple(rows,
        chrom="chrE", strand="+",
        gid="ENSG5", tid="ENST5", pid="ENSP5", gname="TEST5",
        exons=[(100, 104), (200, 203)],
        cds_intervals=[(100, 104), (200, 203)])

    # ENST6 — non-coding (lncRNA), one exon
    emit_noncoding(rows,
        chrom="chrF", strand="+",
        gid="ENSG6", tid="ENST6", gname="TEST6_LNC",
        exons=[(500, 600)])

    with open(path, "w") as f:
        f.write("#!genome-build synthetic_p2g_test\n")
        for r in rows:
            f.write(r + "\n")
    print(f"wrote {path} ({len(rows)} rows)")


def build_no_tags(path):
    rows = []
    emit_simple(rows,
        chrom="chrA", strand="+",
        gid="ENSG7", tid="ENST7", pid="ENSP7", gname="TEST7",
        exons=[(100, 149), (200, 250), (300, 400)],
        cds_intervals=[(120, 149), (200, 250), (300, 350)],
        tags=None)  # no tags at all
    # Sanity: remove any accidental `tag "..."` (none expected).
    rows = [r for r in rows if 'tag "' not in r]
    with open(path, "w") as f:
        f.write("#!genome-build synthetic_p2g_test_no_tags\n")
        for r in rows:
            f.write(r + "\n")
    print(f"wrote {path} ({len(rows)} rows)")


if __name__ == "__main__":
    build_with_tags(os.path.join(OUT_DIR, "with_tags.gtf"))
    build_no_tags(os.path.join(OUT_DIR, "no_tags.gtf"))

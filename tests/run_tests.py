#!/usr/bin/env python3
"""End-to-end test suite for protein2genomic.

Builds the synthetic GTFs via make_synthetic_gtf.py, builds an index, runs
several BED queries through the C++ binary, and asserts on the produced
outputs.

Run from the repo root after a successful `cmake --build`:
    python3 tests/run_tests.py

Exits non-zero on any failure.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BIN = REPO_ROOT / "build" / "protein2genomic"
WRAPPER = REPO_ROOT / "bin" / "protein2genomic"
TESTS_DIR = Path(__file__).resolve().parent

WITH_TAGS_GTF = TESTS_DIR / "with_tags.gtf"
NO_TAGS_GTF   = TESTS_DIR / "no_tags.gtf"


# --------------------------------------------------------------------------- #
# Tiny test framework
# --------------------------------------------------------------------------- #

PASSED: list[str] = []
FAILED: list[tuple[str, str]] = []


def assert_eq(name: str, expected, actual):
    if expected == actual:
        PASSED.append(name)
    else:
        FAILED.append((name, f"expected {expected!r}, got {actual!r}"))


def assert_in(name: str, needle, haystack):
    if needle in haystack:
        PASSED.append(name)
    else:
        FAILED.append((name, f"expected {needle!r} to be in {haystack!r}"))


def assert_true(name: str, cond, hint: str = ""):
    if cond:
        PASSED.append(name)
    else:
        FAILED.append((name, hint or "condition was false"))


def run(*args, expect_zero: bool = True) -> subprocess.CompletedProcess:
    proc = subprocess.run([str(a) for a in args],
                          capture_output=True, text=True)
    if expect_zero and proc.returncode != 0:
        raise SystemExit(f"command failed ({proc.returncode}): {args}\n"
                         f"stdout: {proc.stdout}\nstderr: {proc.stderr}")
    return proc


def read_tsv(path: Path) -> list[dict]:
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# --------------------------------------------------------------------------- #
# Test cases
# --------------------------------------------------------------------------- #

def main() -> int:
    # Step 0: regenerate the synthetic GTFs.
    run(sys.executable, TESTS_DIR / "make_synthetic_gtf.py")
    if not BIN.exists():
        raise SystemExit(f"binary not found: {BIN}. Build it first.")

    work = Path(tempfile.mkdtemp(prefix="p2g_test_"))
    print(f"work dir: {work}", file=sys.stderr)

    # ---- Build the with-tags index ---------------------------------------
    idx = work / "with_tags.idx"
    run(BIN, "--gtf", WITH_TAGS_GTF, "--build-index", "--index", idx)
    assert_true("with_tags index exists", idx.exists())

    # ---- Build the no-tags index -----------------------------------------
    idx_no = work / "no_tags.idx"
    run(BIN, "--gtf", NO_TAGS_GTF, "--build-index", "--index", idx_no)
    assert_true("no_tags index exists", idx_no.exists())

    # ---- BED 1: cover ENSP vs ENST, MANE flags, structure_only, beyond,
    #            negative strand, codon splits, selenoprotein, non-coding ---
    bed = work / "queries.bed"
    bed.write_text("""\
# 1) ENSP query, MANE+canonical, plus strand, aa 1..10 = CDS_1 fully
ENSP1\t1\t10\tQ1_ENSP
# 2) ENST query for the same transcript — should produce identical mapping
ENST1\t1\t10\tQ2_ENST
# 3) Versioned ENSP — version is stripped
ENSP1.5\t1\t10\tQ3_VER
# 4) Structure-only (no aa). aa_start = aa_end = 0 keeps the domain_id column
#    intact (whitespace-tokenized parser collapses empty fields).
ENSP1\t0\t0\tQ4_STRUCT
ENSP2\t0\t0\tQ4b_STRUCT_M
# 5) Domain beyond protein length (ENSP1 has 44 aa)
ENSP1\t100\t110\tQ5_BEYOND
# 6) Negative strand, full coverage of CDS_1 in translation order (aa 1..17)
ENSP2\t1\t17\tQ6_NEG
# 7) Selenoprotein-like, domain at the end (aa 8..8)
ENSP3\t8\t8\tQ7_SEC
# 8) Codon split 1+2 — domain on the split aa (aa 2..2)
ENSP4\t2\t2\tQ8_SPLIT12
# 9) Codon split 2+1 — domain on the split aa (aa 2..2)
ENSP5\t2\t2\tQ9_SPLIT21
# 10) Non-coding ENST — no CDS
ENST6\t1\t10\tQ10_NONCODING
# 11) Protein not in index
ENSP99\t1\t10\tQ11_NOTFOUND
""")

    out = work / "out_all"
    run(BIN, "--index", idx, "--bed", bed, "--out-dir", out, "--output", "all")

    summary = {r["input_id"]: r for r in read_tsv(out / "domain_mapping_summary.tsv")}
    unmapped = {r["input_id"]: r for r in read_tsv(out / "unmapped_domains.tsv")}

    # Q1: ENSP positive case
    q1 = summary["Q1_ENSP"]
    assert_eq("Q1 status ok", "ok", q1["status"])
    assert_eq("Q1 input_id_type ENSP", "ENSP", q1["input_id_type"])
    assert_eq("Q1 is_mane_select", "true", q1["is_mane_select"])
    assert_eq("Q1 is_ensembl_canonical", "true", q1["is_ensembl_canonical"])
    assert_eq("Q1 cds_length_mismatch", "false", q1["cds_length_mismatch"])
    assert_eq("Q1 cds_nt_remainder", "0", q1["cds_nt_remainder"])
    assert_eq("Q1 protein_length_aa = 44", "44", q1["protein_length_aa"])
    assert_eq("Q1 domain_genomic_start = 120", "120", q1["domain_genomic_start"])
    assert_eq("Q1 domain_genomic_end = 149", "149", q1["domain_genomic_end"])
    assert_eq("Q1 n_coding_segments = 1", "1", q1["n_coding_segments"])

    # Q2: ENST input gives same result, except input_id_type
    q2 = summary["Q2_ENST"]
    assert_eq("Q2 input_id_type ENST", "ENST", q2["input_id_type"])
    assert_eq("Q2 domain_genomic_start same", q1["domain_genomic_start"],
              q2["domain_genomic_start"])
    assert_eq("Q2 domain_genomic_end same", q1["domain_genomic_end"],
              q2["domain_genomic_end"])
    assert_eq("Q2 protein_id resolved to ENSP1", "ENSP1", q2["protein_id"])

    # Q3: versioned id
    q3 = summary["Q3_VER"]
    assert_eq("Q3 versioned stripped, protein_id = ENSP1", "ENSP1", q3["protein_id"])

    # Q4 / Q4b: structure_only
    assert_eq("Q4 status structure_only", "structure_only", summary["Q4_STRUCT"]["status"])
    assert_eq("Q4 no_domain_mode true", "true", summary["Q4_STRUCT"]["no_domain_mode"])
    assert_eq("Q4b status structure_only", "structure_only",
              summary["Q4b_STRUCT_M"]["status"])

    # Q5: beyond protein length → unmapped (still in summary with the reason)
    assert_true("Q5 unmapped row exists", "Q5_BEYOND" in unmapped)
    assert_eq("Q5 reason", "domain_beyond_protein_length",
              unmapped["Q5_BEYOND"]["reason"])

    # Q6: negative strand, aa 1..17 = CDS_1 in translation order (genomic
    # 800..850). protein_length_aa = 41 (123 nt / 3).
    q6 = summary["Q6_NEG"]
    assert_eq("Q6 strand", "-", q6["strand"])
    assert_eq("Q6 domain_genomic_start = 800", "800", q6["domain_genomic_start"])
    assert_eq("Q6 domain_genomic_end = 850", "850", q6["domain_genomic_end"])
    assert_eq("Q6 protein_length_aa = 41", "41", q6["protein_length_aa"])
    # MANE/canonical: this transcript has none, but the GTF has tags overall.
    assert_eq("Q6 is_mane_select false (GTF has tags)", "false", q6["is_mane_select"])
    assert_eq("Q6 is_canonical false (GTF has tags)", "false",
              q6["is_ensembl_canonical"])

    # Q7: selenoprotein-like, status ok_cds_mismatch
    q7 = summary["Q7_SEC"]
    assert_eq("Q7 cds_length_mismatch true", "true", q7["cds_length_mismatch"])
    assert_eq("Q7 cds_nt_remainder", "1", q7["cds_nt_remainder"])
    assert_in("Q7 status carries _cds_mismatch", "_cds_mismatch", q7["status"])

    # Q8: codon split 1+2 — aa 2 covered by both CDS rows. Fractions sum to 1.
    isoform_rows = read_tsv(out / "isoform_structure.tsv")
    q8_cds = [r for r in isoform_rows
              if r["input_id"] == "Q8_SPLIT12" and r["feature_type"] == "CDS"
              and r["overlaps_domain"] == "coding_overlap"]
    assert_eq("Q8 two coding_overlap rows", 2, len(q8_cds))
    q8_total = sum(float(r["domain_overlap_fraction_of_domain"]) for r in q8_cds)
    assert_true("Q8 fractions sum to ~1.0",
                abs(q8_total - 1.0) < 1e-6,
                f"sum = {q8_total}")

    # Q9: codon split 2+1, symmetric check
    q9_cds = [r for r in isoform_rows
              if r["input_id"] == "Q9_SPLIT21" and r["feature_type"] == "CDS"
              and r["overlaps_domain"] == "coding_overlap"]
    assert_eq("Q9 two coding_overlap rows", 2, len(q9_cds))
    q9_total = sum(float(r["domain_overlap_fraction_of_domain"]) for r in q9_cds)
    assert_true("Q9 fractions sum to ~1.0",
                abs(q9_total - 1.0) < 1e-6,
                f"sum = {q9_total}")

    # Q10: non-coding ENST → no_CDS_for_protein
    assert_true("Q10 unmapped row exists", "Q10_NONCODING" in unmapped)
    assert_eq("Q10 reason", "no_CDS_for_protein", unmapped["Q10_NONCODING"]["reason"])

    # Q11: protein not in index
    assert_true("Q11 unmapped row exists", "Q11_NOTFOUND" in unmapped)
    assert_eq("Q11 reason", "protein_not_in_index",
              unmapped["Q11_NOTFOUND"]["reason"])

    # ---- BED 2: same query against the no-tags GTF ------------------------
    bed2 = work / "queries_notags.bed"
    bed2.write_text("ENSP7\t1\t10\tQ_NOTAGS\n")
    out2 = work / "out_notags"
    run(BIN, "--index", idx_no, "--bed", bed2, "--out-dir", out2, "--output", "all")
    no_summary = {r["input_id"]: r for r in
                  read_tsv(out2 / "domain_mapping_summary.tsv")}
    qn = no_summary["Q_NOTAGS"]
    assert_eq("No-tags GTF: is_mane_select = NA", "NA", qn["is_mane_select"])
    assert_eq("No-tags GTF: is_ensembl_canonical = NA", "NA",
              qn["is_ensembl_canonical"])

    # ---- BED12 sanity ----------------------------------------------------
    bed12 = (out / "domain_blocks.bed12").read_text().splitlines()
    # We expect a BED12 row for every successful query that has a domain.
    # Q1, Q2, Q3 (ENSP1 aa1..10), Q6 (ENSP2 aa1..17), Q7 (ENSP3 aa8),
    # Q8 (ENSP4), Q9 (ENSP5). 7 rows.
    assert_eq("BED12 row count", 7, len(bed12))
    # Q1's BED12: one block, size 30, start 0, chromStart 119 (genomic 120-1).
    q1_line = [l for l in bed12 if "Q1_ENSP" in l]
    assert_eq("Q1 BED12 lines == 1", 1, len(q1_line))
    parts = q1_line[0].split("\t")
    assert_eq("Q1 BED12 chrom", "chrA", parts[0])
    assert_eq("Q1 BED12 chromStart", "119", parts[1])
    assert_eq("Q1 BED12 chromEnd", "149", parts[2])
    assert_eq("Q1 BED12 blockCount", "1", parts[9])
    assert_eq("Q1 BED12 blockSizes", "30,", parts[10])
    assert_eq("Q1 BED12 blockStarts", "0,", parts[11])
    # Q8's BED12: split 1+2, two blocks.
    q8_line = [l for l in bed12 if "Q8_SPLIT12" in l][0].split("\t")
    assert_eq("Q8 BED12 blockCount", "2", q8_line[9])

    # ---- Plotter smoke test ----------------------------------------------
    pdf = work / "Q1.pdf"
    proc = subprocess.run(
        [str(WRAPPER), "plot", "--isoform", str(out / "isoform_structure.tsv"),
         "--input-id", "Q1_ENSP", "--out", str(pdf)],
        capture_output=True, text=True,
    )
    assert_eq("plot exit 0", 0, proc.returncode)
    assert_true("PDF created", pdf.exists() and pdf.stat().st_size > 0)

    # ---- Wrap up ---------------------------------------------------------
    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    for n, why in FAILED:
        print(f"  FAIL  {n}: {why}", file=sys.stderr)
    if FAILED:
        return 1
    # Clean up.
    shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

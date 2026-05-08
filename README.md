# Protein2Genomic

Map protein domain coordinates to genomic / transcript structure using a GTF annotation.

For each input domain (a `protein_id` and an aa range), the tool answers two related but distinct questions:

1. **Mapping** — *which exact genomic bases code this domain?*
2. **Structure / visualization** — *how is the whole transcript organized in 5′UTR / CDS / 3′UTR / intron, and where does the domain fall on it?*

The two questions correspond to different output modes (`coding`, `span`, `isoform`, `all`). The `isoform` mode is the plot-ready one: it returns a single tidy table with one row per structural feature of the transcript and an explicit overlap classification per row.

---

## Table of contents

1. [Build](#build)
2. [Quickstart](#quickstart)
3. [Input format](#input-format)
4. [Output modes](#output-modes)
5. [Coordinate conventions](#coordinate-conventions)
6. [File-by-file schemas](#file-by-file-schemas)
   - [domain_mapping_summary.tsv](#domain_mapping_summarytsv)
   - [domain_cds_segments.bed](#domain_cds_segmentsbed)
   - [domain_cds_segments.tsv](#domain_cds_segmentstsv)
   - [domain_span_with_introns.bed](#domain_span_with_intronsbed)
   - [isoform_structure.tsv](#isoform_structuretsv)
   - [unmapped_domains.tsv](#unmapped_domainstsv)
   - [run_metadata.json](#run_metadatajson)
7. [Worked example](#worked-example)
8. [GTF compatibility & index](#gtf-compatibility--index)
9. [CLI reference](#cli-reference)

---

## Build

Requirements: C++17, CMake ≥ 3.16, OpenMP (optional, parallelizes domain processing).

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

The binary is produced at `build/protein2genomic`.

## Quickstart

```bash
# 1. Build a binary index from a GTF (one-time per annotation)
./protein2genomic --gtf gencode.v49.basic.annotation.gtf \
                  --build-index --index human.idx

# 2. Map a BED of domains using that index
./protein2genomic --index human.idx \
                  --bed domains.bed \
                  --out-dir results \
                  --output all
```

`--out-dir` is created if missing.

## Input format

Whitespace-separated, BED-like. Lines beginning with `#` are ignored.

```
ENSP00000269305    10    50    AD1    TF1
ENSP00000306245     5   100    RD1    TF2
```

| Column | Required | Meaning |
|---|---|---|
| 1 | yes | `protein_id` (versioned or unversioned; the tool strips the suffix) |
| 2 | yes | `aa_start` — 1-based inclusive |
| 3 | yes | `aa_end` — 1-based inclusive |
| 4 | no | `domain_id` (used as `input_id` for tracking) |
| 5+ | no | ignored |

If column 4 is missing, `input_id` falls back to `protein_id:aa_start-aa_end`, so every row remains identifiable in the outputs.

## Output modes

Pass `--output KIND` together with `--out-dir DIR`. All modes additionally write `domain_mapping_summary.tsv` (one row per input domain, ok or not) and `unmapped_domains.tsv` (only if at least one row failed).

| `--output` | Question it answers | Files written (in addition to summary / unmapped) |
|---|---|---|
| `coding` | Which exact genomic bases code each domain? | `domain_cds_segments.bed`, `domain_cds_segments.tsv` |
| `span` | What is the genomic envelope of the domain (incl. introns between coding CDS)? | `domain_span_with_introns.bed` |
| `isoform` | How is the whole transcript organized, and where does the domain fall on it? | `isoform_structure.tsv` |
| `all` (default) | Everything above. | All of the above + `run_metadata.json` |

### `coding`
For each domain, the **exact CDS genomic segments** that encode it. A domain that spans 3 CDS exons produces 3 rows. Rows are emitted in **translation order** (so `segment_index_in_domain = 1` is the most 5′ CDS slice of the protein).

### `span`
For each domain, the **single genomic envelope** from the first to the last domain-coding base, including any introns between them. Useful for locus-level views (IGV, browser tracks).

### `isoform`
The **plot-ready transcript structure**. One row per structural feature (`five_prime_UTR`, `CDS`, `three_prime_UTR`, `intron`). CDS exons that are partially covered by the domain are split into separate rows so that `domain` vs `no-domain` portions are independent rows. This is the table to drive an isoform-architecture figure.

### `all`
Everything above plus `run_metadata.json` (tool version, output kind, annotation source, coordinate conventions, mapped / unmapped counts, timestamp, full CLI invocation).

## Coordinate conventions

| Output | System |
|---|---|
| `*.bed` | 0-based half-open (BED standard). `start = 0-based inclusive`, `end = 0-based exclusive`. |
| `*.tsv` | 1-based inclusive for genomic, CDS-nt, and aa coordinates. Same as GTF. |

`NA` means *not applicable to this row* (e.g. CDS-nt fields on a UTR row).

> **Why does the BED differ by 1 from the TSV?** They describe the same interval in different conventions. A CDS at GTF positions `7676219..7676272` (1-based inclusive, length 54) is BED `7676218..7676272` (0-based half-open, length still 54). `end - start` matches the length in both systems; only `start` shifts.

---

## File-by-file schemas

### `domain_mapping_summary.tsv`

One row per input domain, written for every `--output` mode.

| Column | Type | Meaning |
|---|---|---|
| `input_id` | string | User-supplied identifier (BED column 4 if present, else `protein_id:aa_start-aa_end`) |
| `protein_id` | string | Normalized (version suffix stripped) |
| `transcript_id` | string | Transcript that this protein belongs to |
| `gene_id` | string | Ensembl gene id |
| `gene_name` | string | HGNC-style gene symbol (if present in the GTF) |
| `domain_id` | string | BED column 4 (if any) |
| `chrom` | string | Chromosome of the transcript |
| `strand` | char | `+` or `−` |
| `aa_start` | int | Input domain start (1-based inclusive aa) |
| `aa_end` | int | Input domain end |
| `domain_length_aa` | int | `aa_end − aa_start + 1` |
| `domain_length_nt` | int | `domain_length_aa × 3` |
| `protein_length_aa` | int | Total CDS length / 3 for this protein |
| `domain_genomic_start` | int | Min genomic coord of any domain-coding base (1-based) |
| `domain_genomic_end` | int | Max genomic coord of any domain-coding base (1-based) |
| `n_coding_segments` | int | How many CDS slices the domain spans |
| `fully_mapped` | bool | `true` if the entire aa range fits inside the CDS, else `false` (the range was clipped to fit) |
| `status` | string | `ok` / `partial` / unmapped reason |

### `domain_cds_segments.bed`

6-column BED, one row per CDS slice that codes the domain.

| Column | Meaning |
|---|---|
| 1 | `chrom` |
| 2 | `start` — 0-based |
| 3 | `end` — 0-based exclusive |
| 4 | `name` — `protein_id[_domain_id]_aa_start-aa_end` |
| 5 | `score` — always `0` |
| 6 | `strand` |

Rows are emitted in **translation order** for a given domain; on `−` strand this means high genomic coord first. If you need genomic-ascending order, run `sort -k1,1 -k2,2n` after.

### `domain_cds_segments.tsv`

Same segments as the BED, with full annotation.

| Column | Type | Meaning |
|---|---|---|
| `input_id`, `protein_id`, `transcript_id`, `gene_id`, `gene_name`, `domain_id` | string | Identity columns |
| `chrom`, `strand` | | |
| `genomic_start` | int | 1-based inclusive |
| `genomic_end` | int | 1-based inclusive |
| `segment_index_in_domain` | int | 1..N in translation order |
| `cds_nt_start`, `cds_nt_end` | int | CDS-relative nt offsets (1-based) of this slice |
| `aa_start_encoded`, `aa_end_encoded` | int | Range of aa that this slice encodes |
| `aa_start`, `aa_end` | int | The input domain bounds (repeated on every row of the same domain, for joining) |

### `domain_span_with_introns.bed`

6-column BED, one row per input domain. `start..end` covers from the first to the last domain-coding base, so introns between coding CDS are inside the interval.

| Column | Meaning |
|---|---|
| 1 | `chrom` |
| 2 | `start` — 0-based |
| 3 | `end` — 0-based exclusive |
| 4 | `name` — `protein_id[_domain_id]_aa_start-aa_end` |
| 5 | `score` — always `0` |
| 6 | `strand` |

### `isoform_structure.tsv`

The plot-ready table. One row per structural feature of the transcript, in genomic order. CDS exons are **split** when the domain only overlaps part of them, so a single original CDS exon can produce up to 3 rows: the no-domain part before, the coding-overlap part, and the no-domain part after.

Columns are grouped by purpose. `NA` means the column does not apply to that row type.

#### Identity columns

| Column | Meaning |
|---|---|
| `input_id` | User-supplied identifier; trace each row back to its BED row |
| `gene_id` | Ensembl gene id |
| `gene_name` | Gene symbol |
| `transcript_id` | Transcript id |
| `protein_id` | Normalized protein id |
| `domain_id` | Domain id from BED column 4 |

#### Location columns

| Column | Meaning |
|---|---|
| `chrom` | Chromosome |
| `strand` | `+` or `−` |
| `feature_genomic_start` | 1-based inclusive |
| `feature_genomic_end` | 1-based inclusive |
| `feature_length_nt` | `end − start + 1` |

#### Feature-type columns

| Column | Meaning |
|---|---|
| `feature_type` | One of `five_prime_UTR`, `CDS`, `three_prime_UTR`, `intron` |
| `feature_id` | `<feature_type>_<n>` numbered in **translation order**. `CDS_1` is the most 5′ CDS in the protein, `CDS_2` the next, etc. After CDS-splitting, two adjacent rows may share the same `feature_id` because they are slices of the same original CDS exon |
| `exon_number` | Source GTF `exon_number` for UTR / CDS rows. `NA` for introns |

#### Ordering columns

| Column | Meaning |
|---|---|
| `feature_order_genomic` | 1..N along the chromosome (low → high coord). Always equals row position |
| `feature_order_transcript` | 1..N in translation direction. Equals `feature_order_genomic` on `+` strand; reversed on `−` strand |

For `−` strand genes, plot using `feature_genomic_start/end` on the X axis (genomic coords), but use `feature_order_transcript` to interpret biological order (5′ → 3′ of the protein).

#### CDS-coordinate columns (NA on UTR / intron rows)

| Column | Meaning |
|---|---|
| `cds_nt_start` | CDS-relative nt offset (1-based) of this slice's first base |
| `cds_nt_end` | CDS-relative nt offset of the last base |
| `aa_start_encoded` | First aa that this slice encodes (1-based) |
| `aa_end_encoded` | Last aa that this slice encodes |

The mapping is `aa = ⌈cds_nt / 3⌉`. A 1-nt CDS slice that contains only the third base of an aa still reports that aa in `aa_start_encoded` / `aa_end_encoded`.

#### Domain-overlap columns

`overlaps_domain` is **not a yes/no flag** — it discriminates *coding* overlap from *intronic* overlap inside the domain envelope:

| Value | Meaning |
|---|---|
| `no` | The row is outside the domain entirely. UTR rows always carry `no` (UTRs cannot encode the domain even if they fall inside the genomic envelope) |
| `coding_overlap` | CDS row whose genomic interval overlaps a domain-coding range. The row encodes part of the domain |
| `inside_domain_genomic_span` | Intron located between two `coding_overlap` CDS rows. The intron itself does not encode the domain, but it lies inside the domain's genomic envelope and should usually be drawn the same colour as the surrounding domain CDS |

The companion columns are filled only for `coding_overlap` rows:

| Column | Meaning |
|---|---|
| `domain_overlap_genomic_start` / `_end` | The sub-interval of this row that codes the domain (1-based inclusive). For a fully coding row this equals `feature_genomic_start..feature_genomic_end` |
| `domain_overlap_cds_nt_start` / `_end` | Same overlap projected to CDS-nt (1-based) |
| `domain_overlap_aa_start` / `_end` | Same overlap projected to aa |
| `domain_overlap_fraction_of_feature` | overlap_length / `feature_length_nt`. `1.0000` if the row is fully inside the domain |
| `domain_overlap_fraction_of_domain` | overlap_length / `domain_length_nt`. Sums to `1.0` across all `coding_overlap` rows of the same domain |

#### Plotting column

`plot_group` is a single string suitable for direct mapping to a colour scale:

| Value | Maps to |
|---|---|
| `five_prime_UTR` | 5′ UTR exon segment |
| `three_prime_UTR` | 3′ UTR exon segment |
| `CDS_no_domain` | CDS segment outside the domain |
| `CDS_domain` | CDS segment that encodes the domain |
| `intron` | Intron outside the domain genomic span |
| `intron_domain_span` | Intron between two `CDS_domain` rows |

In `ggplot`/`ggtranscript`-style code:
```r
ggplot(rows, aes(xmin = feature_genomic_start, xmax = feature_genomic_end,
                 y = transcript_id, fill = plot_group)) +
  geom_rect()
```

### `unmapped_domains.tsv`

Written only when at least one row failed.

| Column | Meaning |
|---|---|
| `input_id`, `protein_id`, `aa_start`, `aa_end`, `domain_id` | identity |
| `reason` | one of: `protein_not_in_index`, `no_CDS_for_protein`, `domain_beyond_protein_length`, `no_overlap` |

### `run_metadata.json`

Written only with `--output all`. Records:

- `tool`, `version`, `timestamp_utc`
- `output_kind` — the value of `--output`
- `annotation_source` — path to GTF or index used
- `index_format_version`
- `coordinate_conventions` — restated explicitly
- `domain_counts` — `{ total, mapped, unmapped }`
- `cli` — full argv of the invocation

---

## Worked example

Input BED:
```
ENSP00000306245    5    100    RD1    TF2
```

`ENSP00000306245` is a `+` strand transcript with this layout:
```
[CDS slice 1]   intron 1   [CDS slice 2]   intron 2   [CDS slice 3]   ...
75278983..94    75279124..  75279877..035   ...
```

Domain `RD1` (aa 5..100) spans the end of slice 1 and most of slice 2. The relevant rows of `isoform_structure.tsv` look like (selected columns shown):

| feature_type | feature_id | feature_genomic_start | feature_genomic_end | aa_start_encoded | aa_end_encoded | overlaps_domain | plot_group |
|---|---|---|---|---|---|---|---|
| CDS | CDS_1 | 75278983 | 75278994 | 1 | 4 | no | CDS_no_domain |
| CDS | CDS_2 | 75278995 | 75279123 | 5 | 47 | coding_overlap | CDS_domain |
| intron | intron_1 | 75279124 | 75279876 | NA | NA | inside_domain_genomic_span | intron_domain_span |
| CDS | CDS_3 | 75279877 | 75280035 | 48 | 100 | coding_overlap | CDS_domain |
| CDS | CDS_4 | 75280036 | 75280128 | 101 | 131 | no | CDS_no_domain |
| intron | intron_2 | 75280129 | 75280559 | NA | NA | no | intron |
| ... | ... | ... | ... | ... | ... | ... | ... |
| three_prime_UTR | three_prime_UTR_1 | 75281422 | 75282230 | NA | NA | no | three_prime_UTR |

Reading this row by row:
- `CDS_1` exists in the transcript but encodes aa 1..4 — outside the domain (which starts at aa 5).
- `CDS_2` is fully inside the domain; it encodes aa 5..47.
- `intron_1` is between two `CDS_domain` rows, so it's marked `inside_domain_genomic_span` even though introns can't code anything.
- `CDS_3` is fully inside the domain; it encodes aa 48..100.
- `CDS_4` and beyond are outside the domain.

To plot, fill by `plot_group`. To highlight the domain, the rows you want are `plot_group ∈ {CDS_domain, intron_domain_span}`.

---

## GTF compatibility & index

Both major GTF flavours work without configuration:

- **GENCODE** — versioned IDs (`protein_id "ENSP00000306245.4"`)
- **Ensembl** — unversioned IDs (`protein_id "ENSP00000306245"`)

The tool normalizes IDs by stripping the `.<version>` suffix on both sides (GTF and BED), so a BED with versioned protein IDs works against an Ensembl-derived index, and vice versa.

The binary index stores per-protein exon and CDS intervals plus `gene_id`, `gene_name`, and `exon_number`. The format is versioned (`INDEX_FORMAT_VERSION = 2`); loading an older index returns an explicit error asking you to rebuild it. Rebuild is fast (~10 s for human GENCODE).

---

## CLI reference

```
USAGE
  protein2genomic --gtf FILE --build-index --index FILE
  protein2genomic (--gtf FILE | --index FILE) --bed FILE --out-dir DIR [--output KIND]

OPTIONS
  --bed FILE       BED with rows: protein_id  aa_start  aa_end  [domain_id]
                   aa coordinates are 1-based inclusive.
  --out-dir DIR    Output directory (created if missing).
  --gtf FILE       GTF annotation, parsed on the fly.
  --index FILE     Pre-built binary index (faster, recommended).
  --output KIND    {coding, span, isoform, all}. Default: all.
  --build-index    Build a binary index from --gtf into --index.
  --threads NUM    Process domains in parallel via OpenMP. Default: 1.
  --verbose        Log progress to stderr.
  --version        Print version and exit.
  --help           Full help including all output schemas.
```

Run `protein2genomic --help` for the full schema reference at the terminal.

## License

MIT.

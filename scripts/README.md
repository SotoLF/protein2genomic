# Input preparation scripts

`prot2exon` expects a BED-like file with `ENSP   aa_start   aa_end   domain_id`. Most users don't have that file lying around — they have an InterProScan run, a UniProt entry, or an HMMER Pfam scan. These scripts convert each of those into the BED-like format the mapper takes.

| Script | Input format | Source |
|---|---|---|
| [prepare_from_interpro.py](prepare_from_interpro.py) | InterProScan TSV (`-f TSV`) | https://github.com/ebi-pf-team/interproscan |
| [prepare_from_uniprot_features.py](prepare_from_uniprot_features.py) | UniProtKB flat-file (`.dat`) or REST JSON | https://www.uniprot.org/ |
| [prepare_from_pfam.py](prepare_from_pfam.py) | HMMER `--domtblout` (against Pfam-A.hmm) | http://eddylab.org/software/hmmer/ |

All three emit the same 5-column BED-like (`ENSP   aa_start   aa_end   domain_id   description`). Column 5 is free-form context (signature description / Pfam description / UniProt /note=). The mapper ignores it but it makes `head` outputs much easier to read.

## The UniProt → ENSP mapping problem

UniProt and Ensembl don't share identifiers. A UniProt accession (`P04637`) can map to several Ensembl proteins (`ENSP00000269305`, `ENSP00000455263`, … — one per transcript of *TP53*), and a few ENSPs don't map to any UniProt accession at all. The scripts handle this in two layers:

1. **Inline xrefs (preferred).** UniProt flat-files and REST JSON carry `DR Ensembl;` lines with the matching ENSP. `prepare_from_uniprot_features.py` reads those *first* — no extra file needed for any entry that has them.
2. **Mapping file fallback.** For inputs that don't carry inline xrefs (notably InterProScan output), pass `--mapping` pointing at Ensembl's per-release UniProt cross-reference TSV. Same column order across releases.

Download once per Ensembl release (~20 MB compressed):

```bash
# Pick the release that matches your GTF. For GENCODE v49 / Ensembl 115:
RELEASE=115
curl -O https://ftp.ensembl.org/pub/release-${RELEASE}/tsv/homo_sapiens/Homo_sapiens.GRCh38.${RELEASE}.uniprot.tsv.gz
```

If your fasta used a custom mapping already (e.g. you ran HMMER on the GENCODE pc_translations FASTA which is ENSP-headered), pass `--id-type ensp` and skip the mapping entirely.

If you have a tiny project-specific mapping (a 2-column `uniprot \t ensp`) you can use `--simple-mapping FILE` instead of `--mapping`.

## Multiple ENSPs per UniProt

Default behavior: emit one BED row per (UniProt, ENSP) pair. That gives you a domain mapped against every transcript that codes for the protein. Pair the output with `--output isoform` and you'll see how a single domain looks across all of the protein's isoforms.

If you only want the canonical mapping, post-filter on `is_mane_select` or `is_ensembl_canonical` in `domain_mapping_summary.tsv` after running `prot2exon`.

## Examples

### InterProScan

```bash
# Common case: InterProScan was run on a UniProt proteome (column 1 = UniProt acc)
python3 scripts/prepare_from_interpro.py \
    --in proteome.interpro.tsv \
    --mapping Homo_sapiens.GRCh38.115.uniprot.tsv.gz \
    --analyses Pfam,SMART,PROSITE_PROFILES \
    --min-length 20 \
    --out domains.bed

# InterProScan was run on the GENCODE pc_translations FASTA (column 1 = ENSP)
python3 scripts/prepare_from_interpro.py \
    --in proteome.interpro.tsv \
    --id-type ensp \
    --out domains.bed

# Map the result
./prot2exon --index human.idx --bed domains.bed --out-dir results --output all
```

### UniProt features (flat-file)

```bash
# Download the TP53 entry as a flat file
curl -o P04637.dat https://rest.uniprot.org/uniprotkb/P04637.txt

# Convert all DOMAIN/REPEAT/REGION/ZN_FING features to BED-like
python3 scripts/prepare_from_uniprot_features.py \
    --in P04637.dat \
    --out p53_domains.bed

# REST JSON works too
curl -o P04637.json https://rest.uniprot.org/uniprotkb/P04637.json
python3 scripts/prepare_from_uniprot_features.py \
    --in P04637.json \
    --feature-types DOMAIN,REGION,ZN_FING \
    --out p53_domains.bed
```

For bulk runs, the REST API also serves JSONLines:

```bash
curl 'https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=jsonl' \
    -o human_reviewed.jsonl
python3 scripts/prepare_from_uniprot_features.py \
    --in human_reviewed.jsonl --out all_domains.bed
```

### Pfam via HMMER

```bash
# 1. Run hmmscan against Pfam-A.hmm (Pfam HMMs as the database, your
#    proteins as the query). Output domtblout, not the per-sequence tblout.
hmmscan --cut_ga --domtblout pfam_hits.dom Pfam-A.hmm proteins.fa > /dev/null

# 2. Convert to BED-like
python3 scripts/prepare_from_pfam.py \
    --in pfam_hits.dom \
    --mode scan \
    --id-type ensp \
    --min-score 25 \
    --out pfam.bed

# If you ran hmmsearch instead (one HMM at a time, your protein FASTA as DB):
python3 scripts/prepare_from_pfam.py --in pfam_hits.dom --mode search --out pfam.bed
```

## Output sanity check

Every script writes a 2-line comment header so you can eyeball what it produced:

```
# prot2exon BED-like, generated from InterProScan
# columns: ENSP   aa_start   aa_end   domain_id   signature_description
ENSP00000269305   95   288   IPR011615_Pfam   P53 transactivation motif
ENSP00000269305   323  356   IPR010991_Pfam   p53, tetramerisation
...
```

The mapper ignores `#` lines and the 5th column, so the file is immediately consumable:

```bash
./prot2exon --index human.idx --bed pfam.bed --out-dir results
```

## Common failure modes

| Symptom | Likely cause | Fix |
|---|---|---|
| Most rows say `protein_not_in_index` | Your BED contains UniProt accessions but the parser was run with `--id-type ensp` | Re-run with `--id-type uniprot` and a mapping file |
| BED is empty after parsing | All hits filtered by `--min-length` / `--min-score` or wrong `--mode` | Try lowering thresholds, or check `--keep-unmapped FILE` for the rejection reasons |
| Several ENSPs per query in the summary | One UniProt accession maps to multiple Ensembl proteins (alt. isoforms) | Expected. Filter on `is_mane_select=true` post-mapping if you want only the canonical |
| `no ENSP for P0DTC2` | UniProt accession isn't in the Ensembl mapping (newer than your TSV, or non-human) | Update the mapping file, or use a different mapping for non-human species |

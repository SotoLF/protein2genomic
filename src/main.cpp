#include "common.hpp"
#include "gtf_parser.hpp"
#include "domain_mapper.hpp"
#include "output_writer.hpp"
#include "utils.hpp"
#include <iostream>
#include <getopt.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

struct Config {
    std::string gtf_file;
    std::string index_file;
    std::string bed_file;
    std::string out_dir;
    OutputKind output_kind = OutputKind::ALL;
    bool build_index = false;
    bool use_index = false;
    bool verbose = false;
    int threads = 1;
};

void print_usage(const char* program_name) {
    std::cout
<< "Protein2Genomic - map protein domain coordinates to genomic / transcript\n"
<< "structure using a GTF annotation.\n"
<< "\n"
<< "USAGE\n"
<< "  " << program_name << " --gtf FILE --build-index --index FILE\n"
<< "  " << program_name << " (--gtf FILE | --index FILE) --bed FILE --out-dir DIR [--output KIND]\n"
<< "\n"
<< "OPTIONS\n"
<< "  Mapping inputs (required):\n"
<< "    --bed FILE       Query file (TSV/BED-like). See BED INPUT FORMAT.\n"
<< "                     Rows can be:\n"
<< "                       id  aa_start  aa_end  [domain_id]\n"
<< "                     or, for whole-transcript structure with no domain:\n"
<< "                       id\n"
<< "                     `id` can be an ENSP (protein) or ENST (transcript).\n"
<< "    --out-dir DIR    Output directory (created if missing).\n"
<< "    --gtf FILE       GTF annotation, parsed on the fly.\n"
<< "    --index FILE     Pre-built binary index (faster, recommended).\n"
<< "                     Use --gtf or --index, not both.\n"
<< "\n"
<< "  Output selection:\n"
<< "    --output KIND    One of {coding, introns, span, isoform, bed12, all}.\n"
<< "                     Default: all. See OUTPUT MODES below.\n"
<< "\n"
<< "  Index management:\n"
<< "    --build-index    Build a binary index from --gtf into --index, then exit\n"
<< "                     (or proceed to mapping if --bed is also given).\n"
<< "\n"
<< "  Other:\n"
<< "    --threads NUM    Process queries in parallel via OpenMP. Default: 1.\n"
<< "    --verbose        Log progress to stderr.\n"
<< "    --version        Print version and exit.\n"
<< "    --help           This message.\n"
<< "\n"
<< "OUTPUT MODES\n"
<< "  Every mode emits all features of the requested type (e.g. every CDS or\n"
<< "  every intron of the transcript) and adds an overlap classification column\n"
<< "  that says which features intersect the domain. The companion BED is the\n"
<< "  subset that actually overlaps with the domain.\n"
<< "\n"
<< "  --output coding\n"
<< "    All CDS rows of the transcript, with domain overlap classification.\n"
<< "    Files written:\n"
<< "      domain_cds_segments.tsv   every CDS of the transcript + overlaps_domain\n"
<< "      domain_cds_segments.bed   subset: CDS rows with coding_overlap only\n"
<< "    A CDS exon that the domain only partially covers is split into multiple\n"
<< "    rows that share the same feature_id and differ by feature_part (1..K).\n"
<< "\n"
<< "  --output introns\n"
<< "    All introns of the transcript, with domain overlap classification.\n"
<< "    Files written:\n"
<< "      domain_introns.tsv        every intron of the transcript + overlaps_domain\n"
<< "      domain_introns.bed        subset: introns inside the domain genomic span\n"
<< "    An intron lying between two domain-coding CDS rows is labeled\n"
<< "    inside_domain_genomic_span (it does not encode the domain but is inside\n"
<< "    its envelope).\n"
<< "\n"
<< "  --output span\n"
<< "    Genomic envelope of the domain: first to last domain-coding base,\n"
<< "    introns included.\n"
<< "    Files written:\n"
<< "      domain_span_with_introns.bed   one BED row per domain (chrom, min..max)\n"
<< "\n"
<< "  --output isoform\n"
<< "    Full transcript structure: 5'UTR / CDS / 3'UTR / intron, in one tidy table.\n"
<< "    Files written:\n"
<< "      isoform_structure.tsv   one row per structural segment of the transcript\n"
<< "    Plot-ready. CDS rows that the domain partially covers are split, sharing\n"
<< "    the same feature_id with different feature_part.\n"
<< "\n"
<< "  --output bed12\n"
<< "    One BED12 row per domain, ready to drop into IGV / UCSC.\n"
<< "    Files written:\n"
<< "      domain_blocks.bed12   chromStart..chromEnd is the genomic envelope\n"
<< "                            of the domain (introns included). Blocks are the\n"
<< "                            CDS slices that code the domain. thickStart and\n"
<< "                            thickEnd equal chromStart/chromEnd so the entire\n"
<< "                            feature is drawn thick. itemRgb defaults to red.\n"
<< "\n"
<< "  --output all\n"
<< "    Writes coding + introns + span + isoform + bed12 + run_metadata.json.\n"
<< "\n"
<< "  Every mode also writes:\n"
<< "    domain_mapping_summary.tsv   one row per input query (success or partial)\n"
<< "    unmapped_domains.tsv         only present when at least one row failed\n"
<< "\n"
<< "  No-domain mode: when a BED row provides only protein_id (or aa_start =\n"
<< "  aa_end = 0), the row is processed as 'just give me the structure of this\n"
<< "  transcript'. All overlap-related columns are NA, plot_group becomes the\n"
<< "  feature_type (e.g. CDS, intron), and the *.bed subset files are empty for\n"
<< "  those rows (there is no domain to subset to).\n"
<< "\n"
<< "COORDINATE CONVENTIONS\n"
<< "  *.bed files : 0-based half-open (BED standard)\n"
<< "  *.tsv files : 1-based inclusive for genomic, CDS-nt, and aa coordinates\n"
<< "  Values reported as NA mean 'does not apply to this row'.\n"
<< "\n"
<< "SCHEMA: domain_mapping_summary.tsv\n"
<< "  One row per input query. status is 'ok' if the entire aa range was\n"
<< "  mapped to CDS, 'partial' if it was clipped, 'structure_only' for no-domain\n"
<< "  rows, or the unmapped reason.\n"
<< "  Columns:\n"
<< "    input_id              user-supplied row identifier (column 4 of BED, else\n"
<< "                          protein_id:aa_start-aa_end, or protein_id alone)\n"
<< "    protein_id            normalized (version stripped)\n"
<< "    transcript_id         transcript that the protein belongs to\n"
<< "    gene_id               Ensembl gene id\n"
<< "    gene_name             HGNC-style gene symbol (if present in GTF)\n"
<< "    domain_id             column 4 of the BED, if any\n"
<< "    chrom, strand         chromosome and strand of the transcript\n"
<< "    aa_start, aa_end      input domain bounds (1-based inclusive aa)\n"
<< "                          NA in no-domain mode\n"
<< "    domain_length_aa      aa_end - aa_start + 1 (NA in no-domain mode)\n"
<< "    domain_length_nt      domain_length_aa * 3 (NA in no-domain mode)\n"
<< "    protein_length_aa     total CDS length / 3 for this protein\n"
<< "    domain_genomic_start  min genomic coord of any domain-coding base\n"
<< "    domain_genomic_end    max genomic coord of any domain-coding base\n"
<< "    n_coding_segments     number of CDS exon slices the domain spans\n"
<< "    fully_mapped          true if the full aa range fits inside the CDS\n"
<< "                          (always true in no-domain mode)\n"
<< "    no_domain_mode        true when the BED row had no aa range\n"
<< "    status                ok | partial | structure_only | <unmapped reason>\n"
<< "\n"
<< "SCHEMA: feature TSV files (domain_cds_segments.tsv, domain_introns.tsv,\n"
<< "                          isoform_structure.tsv)\n"
<< "  All three share the same column layout, differing only in which feature\n"
<< "  types they include:\n"
<< "    domain_cds_segments.tsv  -> CDS rows only\n"
<< "    domain_introns.tsv       -> intron rows only\n"
<< "    isoform_structure.tsv    -> 5'UTR + CDS + 3'UTR + intron rows\n"
<< "\n"
<< "  Identity:\n"
<< "    input_id, gene_id, gene_name, transcript_id, protein_id, domain_id\n"
<< "  Location:\n"
<< "    chrom, strand\n"
<< "    feature_genomic_start, feature_genomic_end   1-based inclusive\n"
<< "    feature_length_nt                            end - start + 1\n"
<< "  Feature type:\n"
<< "    feature_type    five_prime_UTR | CDS | three_prime_UTR | intron\n"
<< "    feature_id      stable across CDS splits, e.g. CDS_3 for all sub-parts of\n"
<< "                    GTF CDS #3 (translation-order numbering)\n"
<< "    feature_part    1..K within rows sharing the same feature_id (CDS splits),\n"
<< "                    in translation order. Always 1 for UTR and intron rows\n"
<< "    exon_number     source GTF exon_number; NA for introns\n"
<< "  Ordering:\n"
<< "    feature_order_genomic     1..N along the chromosome\n"
<< "    feature_order_transcript  1..N in translation direction\n"
<< "                              (= feature_order_genomic on +, reversed on -)\n"
<< "  CDS coordinates (NA for UTR / intron rows):\n"
<< "    cds_nt_start, cds_nt_end       CDS-relative nt offsets (1-based)\n"
<< "    aa_start_encoded, aa_end_encoded\n"
<< "                                   aa range that this row encodes\n"
<< "  Domain overlap (NA when no_domain_mode):\n"
<< "    overlaps_domain   one of:\n"
<< "      no                          - the row does not intersect the domain\n"
<< "      coding_overlap              - CDS row that encodes part of the domain\n"
<< "      inside_domain_genomic_span  - intron between two coding_overlap rows\n"
<< "    domain_overlap_genomic_start, domain_overlap_genomic_end\n"
<< "    domain_overlap_cds_nt_start,  domain_overlap_cds_nt_end\n"
<< "    domain_overlap_aa_start,      domain_overlap_aa_end\n"
<< "    domain_overlap_fraction_of_feature   overlap_length / feature_length_nt\n"
<< "    domain_overlap_fraction_of_domain    overlap_length / total_domain_length_nt\n"
<< "  Plotting:\n"
<< "    plot_group   color-coding shortcut. Values:\n"
<< "      With a domain:\n"
<< "        five_prime_UTR | three_prime_UTR\n"
<< "        CDS_no_domain | CDS_domain\n"
<< "        intron | intron_domain_span\n"
<< "      No-domain mode (just the feature type, no overlap info):\n"
<< "        five_prime_UTR | three_prime_UTR | CDS | intron\n"
<< "\n"
<< "SCHEMA: companion BED files\n"
<< "  domain_cds_segments.bed     CDS rows with overlaps_domain == coding_overlap\n"
<< "  domain_introns.bed          intron rows inside_domain_genomic_span\n"
<< "  domain_span_with_introns.bed   single row per domain: first to last\n"
<< "                              domain-coding base, introns included\n"
<< "  All are 6-column BED: chrom, start_0based, end, name, score=0, strand.\n"
<< "  name = protein_id[_domain_id]_<aa_start>-<aa_end> (or '_no_domain' suffix\n"
<< "  in no-domain mode). All three BEDs are empty for no-domain queries.\n"
<< "\n"
<< "SCHEMA: unmapped_domains.tsv\n"
<< "  Only written when at least one BED row could not be mapped.\n"
<< "  Columns:\n"
<< "    input_id, protein_id, aa_start, aa_end, domain_id\n"
<< "    reason   one of:\n"
<< "      protein_not_in_index        protein_id is absent from the GTF/index\n"
<< "      no_CDS_for_protein          protein has no CDS records (rare)\n"
<< "      domain_beyond_protein_length\n"
<< "                                  aa_start exceeds protein length\n"
<< "      no_overlap                  aa range does not intersect the CDS\n"
<< "\n"
<< "SCHEMA: run_metadata.json (only with --output all)\n"
<< "  Records tool version, output_kind, annotation source, index format\n"
<< "  version, coordinate conventions, total/mapped/unmapped/no_domain counts,\n"
<< "  timestamp, and the full CLI invocation.\n"
<< "\n"
<< "BED INPUT FORMAT\n"
<< "  Whitespace-separated. Lines starting with '#' are comments.\n"
<< "    column 1   id: ENSP (protein) or ENST (transcript), with or without\n"
<< "               version suffix. ENST and ENSP that point to the same\n"
<< "               transcript produce identical mapping output.\n"
<< "    column 2   aa_start (1-based inclusive). Optional.\n"
<< "    column 3   aa_end   (1-based inclusive). Optional.\n"
<< "    column 4   domain_id (optional, used as input_id)\n"
<< "    column 5+  ignored\n"
<< "  No-domain mode: rows with only column 1, or with aa_start = aa_end = 0,\n"
<< "  are processed for whole-transcript structure (overlap columns are NA).\n"
<< "  Versions are stripped on both sides, so a BED with\n"
<< "  ENSP00000269305.4 (or ENST00000269305.9) will resolve against an index\n"
<< "  built from a GTF with unversioned ids (and vice versa).\n"
<< "  Non-coding ENST queries surface as 'no_CDS_for_protein' in\n"
<< "  unmapped_domains.tsv.\n"
<< "\n"
<< "EXAMPLES\n"
<< "  Build the index once:\n"
<< "    " << program_name << " --gtf gencode.v49.basic.annotation.gtf \\\n"
<< "                          --build-index --index human.idx\n"
<< "\n"
<< "  Map domains, all outputs:\n"
<< "    " << program_name << " --index human.idx --bed domains.bed \\\n"
<< "                          --out-dir results --output all\n"
<< "\n"
<< "  Plot-ready table only:\n"
<< "    " << program_name << " --index human.idx --bed domains.bed \\\n"
<< "                          --out-dir results --output isoform\n"
<< "\n"
<< "  Whole-transcript structure (no domain):\n"
<< "    printf 'ENSP00000269305\\nENSP00000306245\\n' > proteins.txt\n"
<< "    " << program_name << " --index human.idx --bed proteins.txt \\\n"
<< "                          --out-dir results --output isoform\n"
        ;
}

void print_version() {
    std::cout << "Protein2Genomic 2.2.0 (index format v" << INDEX_FORMAT_VERSION << ")\n";
}

bool parse_output_kind(const std::string& s, OutputKind& out) {
    if (s == "coding")  { out = OutputKind::CODING;  return true; }
    if (s == "introns") { out = OutputKind::INTRONS; return true; }
    if (s == "span")    { out = OutputKind::SPAN;    return true; }
    if (s == "isoform") { out = OutputKind::ISOFORM; return true; }
    if (s == "bed12")   { out = OutputKind::BED12;   return true; }
    if (s == "all")     { out = OutputKind::ALL;     return true; }
    return false;
}

int parse_arguments(int argc, char* argv[], Config& config) {
    static struct option long_options[] = {
        {"gtf",         required_argument, 0, 'g'},
        {"bed",         required_argument, 0, 'b'},
        {"out-dir",     required_argument, 0, 'O'},
        {"output",      required_argument, 0, 'o'},
        {"build-index", no_argument,       0, 'i'},
        {"index",       required_argument, 0, 'x'},
        {"threads",     required_argument, 0, 't'},
        {"verbose",     no_argument,       0, 'v'},
        {"help",        no_argument,       0, 'h'},
        {"version",     no_argument,       0, 'V'},
        {0, 0, 0, 0}
    };
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "g:b:O:o:ix:t:vhV", long_options, &option_index)) != -1) {
        switch (c) {
            case 'g': config.gtf_file = optarg; break;
            case 'b': config.bed_file = optarg; break;
            case 'O': config.out_dir = optarg; break;
            case 'o': {
                if (!parse_output_kind(optarg, config.output_kind)) {
                    std::cerr << "Error: --output must be one of "
                                 "{coding,introns,span,isoform,bed12,all}\n";
                    return 1;
                }
                break;
            }
            case 'i': config.build_index = true; break;
            case 'x': config.index_file = optarg; config.use_index = true; break;
            case 't':
                config.threads = std::stoi(optarg);
                if (config.threads < 1) { std::cerr << "Error: --threads must be >= 1\n"; return 1; }
                break;
            case 'v': config.verbose = true; break;
            case 'h': print_usage(argv[0]); return 0;
            case 'V': print_version(); return 0;
            case '?': return 1;
            default:  return 1;
        }
    }
    return -1;
}

int validate_config(const Config& config) {
    if (config.build_index) {
        if (config.gtf_file.empty()) {
            std::cerr << "Error: --build-index requires --gtf\n"; return 1;
        }
        if (config.index_file.empty()) {
            std::cerr << "Error: --build-index requires --index\n"; return 1;
        }
        if (!utils::file_exists(config.gtf_file)) {
            std::cerr << "Error: GTF file not found: " << config.gtf_file << "\n"; return 1;
        }
        return 0;
    }

    // Mapping mode.
    if (config.gtf_file.empty() && !config.use_index) {
        std::cerr << "Error: provide --gtf or --index\n"; return 1;
    }
    if (config.bed_file.empty()) {
        std::cerr << "Error: --bed is required\n"; return 1;
    }
    if (config.out_dir.empty()) {
        std::cerr << "Error: --out-dir is required\n"; return 1;
    }
    if (!config.gtf_file.empty() && !utils::file_exists(config.gtf_file)) {
        std::cerr << "Error: GTF file not found: " << config.gtf_file << "\n"; return 1;
    }
    if (config.use_index && !utils::file_exists(config.index_file)) {
        std::cerr << "Error: index file not found: " << config.index_file << "\n"; return 1;
    }
    if (!utils::file_exists(config.bed_file)) {
        std::cerr << "Error: BED file not found: " << config.bed_file << "\n"; return 1;
    }
    return 0;
}

int main(int argc, char* argv[]) {
    Config config;
    int parse_result = parse_arguments(argc, argv, config);
    if (parse_result >= 0) return parse_result;
    if (validate_config(config) != 0) return 1;

#ifdef USE_OPENMP
    omp_set_num_threads(config.threads);
    if (config.verbose) std::cerr << "Threads: " << config.threads << "\n";
#endif

    try {
        GTFParser parser;

        if (config.build_index) {
            Timer timer("Index building");
            ErrorCode rc = parser.load_gtf(config.gtf_file);
            if (rc != ErrorCode::SUCCESS) {
                std::cerr << "Error building index: " << utils::error_code_to_string(rc) << "\n";
                return 1;
            }
            rc = parser.save_index(config.index_file);
            if (rc != ErrorCode::SUCCESS) {
                std::cerr << "Error saving index: " << utils::error_code_to_string(rc) << "\n";
                return 1;
            }
            std::cerr << "Index saved to: " << config.index_file << "\n";
            MemoryTracker::print_memory_stats();
            // If no --bed was given, we're done.
            if (config.bed_file.empty()) return 0;
        } else if (config.use_index) {
            Timer timer("Index loading");
            ErrorCode rc = parser.load_index(config.index_file);
            if (rc != ErrorCode::SUCCESS) {
                std::cerr << "Error loading index: " << utils::error_code_to_string(rc) << "\n";
                return 1;
            }
        } else {
            Timer timer("GTF parsing");
            ErrorCode rc = parser.load_gtf(config.gtf_file);
            if (rc != ErrorCode::SUCCESS) {
                std::cerr << "Error parsing GTF: " << utils::error_code_to_string(rc) << "\n";
                return 1;
            }
        }

        DomainMapper mapper(parser, config.output_kind);
        {
            Timer t("Domain loading");
            ErrorCode rc = mapper.load_domains(config.bed_file);
            if (rc != ErrorCode::SUCCESS) return 1;
        }
        {
            Timer t("Domain mapping");
            ErrorCode rc = mapper.process_domains();
            if (rc != ErrorCode::SUCCESS) return 1;
        }

        std::vector<std::string> cli_args;
        for (int i = 0; i < argc; ++i) cli_args.emplace_back(argv[i]);
        std::string source = config.use_index ? config.index_file : config.gtf_file;

        ErrorCode rc = output::write_all(config.out_dir, config.output_kind,
                                         mapper.results(), source, cli_args);
        if (rc != ErrorCode::SUCCESS) return 1;

        if (config.verbose) {
            std::cerr << "\n=== Final ===\n"
                      << "Proteins in index: " << parser.get_protein_count() << "\n"
                      << "Queries processed: " << mapper.domain_count() << "\n"
                      << "Mapped: " << mapper.mapped_count()
                      << " / Unmapped: " << mapper.unmapped_count() << "\n";
        }
        // Always emit peak RSS so benchmarks can parse it without --verbose.
        MemoryTracker::print_memory_stats();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}

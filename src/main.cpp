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
<< "    --bed FILE       BED with rows: protein_id  aa_start  aa_end  [domain_id]\n"
<< "                     aa coordinates are 1-based inclusive.\n"
<< "    --out-dir DIR    Output directory (created if missing).\n"
<< "    --gtf FILE       GTF annotation, parsed on the fly.\n"
<< "    --index FILE     Pre-built binary index (faster, recommended).\n"
<< "                     Use --gtf or --index, not both.\n"
<< "\n"
<< "  Output selection:\n"
<< "    --output KIND    One of {coding, span, isoform, all}. Default: all.\n"
<< "                     See OUTPUT MODES below.\n"
<< "\n"
<< "  Index management:\n"
<< "    --build-index    Build a binary index from --gtf into --index, then exit\n"
<< "                     (or proceed to mapping if --bed is also given).\n"
<< "\n"
<< "  Other:\n"
<< "    --threads NUM    Process domains in parallel via OpenMP. Default: 1.\n"
<< "    --verbose        Log progress to stderr.\n"
<< "    --version        Print version and exit.\n"
<< "    --help           This message.\n"
<< "\n"
<< "OUTPUT MODES\n"
<< "  Each --output value selects a different question to answer.\n"
<< "\n"
<< "  --output coding\n"
<< "    \"Which exact genomic bases code this protein domain?\"\n"
<< "    Files written:\n"
<< "      domain_cds_segments.bed   one BED row per coding genomic segment\n"
<< "      domain_cds_segments.tsv   same segments with full identity / aa-nt columns\n"
<< "    Each domain may produce 1..N segments (one per CDS exon it spans).\n"
<< "\n"
<< "  --output span\n"
<< "    \"What is the genomic envelope of the domain, including introns between\n"
<< "     domain-coding CDS segments?\"\n"
<< "    Files written:\n"
<< "      domain_span_with_introns.bed   one BED row per domain (chrom, min..max)\n"
<< "    Use this for locus-level views (IGV, browser tracks).\n"
<< "\n"
<< "  --output isoform\n"
<< "    \"How is the whole transcript organized (UTR / CDS / intron) and where\n"
<< "     does the domain fall on it?\"\n"
<< "    Files written:\n"
<< "      isoform_structure.tsv   one row per structural segment of the transcript\n"
<< "    This is the plot-ready table. CDS rows are split so that domain-overlapping\n"
<< "    and non-overlapping portions are separate rows. Schema below.\n"
<< "\n"
<< "  --output all\n"
<< "    Writes coding + span + isoform + run_metadata.json.\n"
<< "\n"
<< "  Every mode also writes:\n"
<< "    domain_mapping_summary.tsv   one row per input domain (success or partial)\n"
<< "    unmapped_domains.tsv         only present when at least one row failed\n"
<< "\n"
<< "COORDINATE CONVENTIONS\n"
<< "  *.bed files : 0-based half-open (BED standard)\n"
<< "  *.tsv files : 1-based inclusive for genomic, CDS-nt, and aa coordinates\n"
<< "  Values reported as NA mean \"does not apply to this row\".\n"
<< "\n"
<< "SCHEMA: domain_mapping_summary.tsv\n"
<< "  One row per input domain. status is 'ok' if the entire aa range was\n"
<< "  mapped to CDS, 'partial' if it was clipped, or the unmapped reason.\n"
<< "  Columns:\n"
<< "    input_id              user-supplied row identifier (column 4 of BED, else\n"
<< "                          protein_id:aa_start-aa_end)\n"
<< "    protein_id            normalized (version stripped)\n"
<< "    transcript_id         transcript that the protein belongs to\n"
<< "    gene_id               Ensembl gene id\n"
<< "    gene_name             HGNC-style gene symbol (if present in GTF)\n"
<< "    domain_id             column 4 of the BED, if any\n"
<< "    chrom, strand         chromosome and strand of the transcript\n"
<< "    aa_start, aa_end      input domain bounds (1-based inclusive aa)\n"
<< "    domain_length_aa      aa_end - aa_start + 1\n"
<< "    domain_length_nt      domain_length_aa * 3\n"
<< "    protein_length_aa     total CDS length / 3 for this protein\n"
<< "    domain_genomic_start  min genomic coord of any domain-coding base\n"
<< "    domain_genomic_end    max genomic coord of any domain-coding base\n"
<< "    n_coding_segments     number of CDS exon slices the domain spans\n"
<< "    fully_mapped          true if the full aa range fits inside the CDS\n"
<< "    status                ok | partial | <reason> (e.g. protein_not_in_index)\n"
<< "\n"
<< "SCHEMA: domain_cds_segments.tsv\n"
<< "  One row per CDS slice that codes part of the domain.\n"
<< "  Columns:\n"
<< "    input_id, protein_id, transcript_id, gene_id, gene_name, domain_id\n"
<< "    chrom, strand            same as summary\n"
<< "    genomic_start, genomic_end\n"
<< "                             1-based inclusive coords of this segment\n"
<< "    segment_index_in_domain  1..N in translation order (5' -> 3' of protein)\n"
<< "    cds_nt_start, cds_nt_end CDS-relative nt offsets (1-based) for this slice\n"
<< "    aa_start_encoded         first aa encoded by this slice (1-based)\n"
<< "    aa_end_encoded           last aa encoded by this slice\n"
<< "    aa_start, aa_end         the input domain bounds (repeated for each row)\n"
<< "  domain_cds_segments.bed contains the same segments as a 6-column BED:\n"
<< "    chrom  start_0based  end  name  score  strand\n"
<< "  where name = protein_id[_domain_id]_<aa_start>-<aa_end>.\n"
<< "\n"
<< "SCHEMA: domain_span_with_introns.bed\n"
<< "  One BED row per input domain. Columns: chrom, start_0based, end, name,\n"
<< "  score (always 0), strand. start..end covers from the first to the last\n"
<< "  domain-coding base, so introns between coding CDS are included.\n"
<< "\n"
<< "SCHEMA: isoform_structure.tsv\n"
<< "  One row per structural feature of the transcript that is relevant to\n"
<< "  the domain. CDS exons are split when the domain only covers part of them.\n"
<< "  Rows are sorted by feature_order_genomic (low to high genomic coord).\n"
<< "  Columns:\n"
<< "    Identity:\n"
<< "      input_id, gene_id, gene_name, transcript_id, protein_id, domain_id\n"
<< "      Trace each row back to its input.\n"
<< "    Location:\n"
<< "      chrom, strand\n"
<< "      feature_genomic_start, feature_genomic_end   1-based inclusive\n"
<< "      feature_length_nt                            end - start + 1\n"
<< "    Feature type:\n"
<< "      feature_type   one of {five_prime_UTR, CDS, three_prime_UTR, intron}\n"
<< "      feature_id     <feature_type>_<n> numbered in translation order\n"
<< "                     (CDS_1 is the most 5' CDS in the protein)\n"
<< "      exon_number    source GTF exon_number for UTR/CDS rows; NA for introns\n"
<< "    Ordering:\n"
<< "      feature_order_genomic    1..N along the chromosome (low -> high)\n"
<< "      feature_order_transcript 1..N in translation direction\n"
<< "                               (= feature_order_genomic on +, reversed on -)\n"
<< "    CDS coordinates (NA for UTR / intron rows):\n"
<< "      cds_nt_start, cds_nt_end       CDS-relative nt offsets (1-based)\n"
<< "      aa_start_encoded, aa_end_encoded\n"
<< "                                     aa range that this row encodes\n"
<< "    Domain overlap:\n"
<< "      overlaps_domain   one of:\n"
<< "        no                          - this row does not encode the domain\n"
<< "        coding_overlap              - CDS row that encodes part of the domain\n"
<< "        inside_domain_genomic_span  - intron located between two\n"
<< "                                      domain-coding CDS rows; the intron\n"
<< "                                      itself does not encode the domain but\n"
<< "                                      lies inside its genomic envelope\n"
<< "      domain_overlap_genomic_start, domain_overlap_genomic_end\n"
<< "                                     sub-interval that codes the domain\n"
<< "                                     (set only for coding_overlap rows)\n"
<< "      domain_overlap_cds_nt_start, domain_overlap_cds_nt_end\n"
<< "      domain_overlap_aa_start, domain_overlap_aa_end\n"
<< "                                     same overlap projected to CDS-nt and aa\n"
<< "      domain_overlap_fraction_of_feature\n"
<< "                                     overlap_length / feature_length_nt\n"
<< "      domain_overlap_fraction_of_domain\n"
<< "                                     overlap_length / total_domain_length_nt\n"
<< "                                     (sums to 1.0 across all coding_overlap\n"
<< "                                     rows of the same domain)\n"
<< "    Plotting:\n"
<< "      plot_group   shorthand for direct color mapping. Values:\n"
<< "        five_prime_UTR        - 5' UTR exon segment\n"
<< "        three_prime_UTR       - 3' UTR exon segment\n"
<< "        CDS_no_domain         - CDS segment outside the domain\n"
<< "        CDS_domain            - CDS segment that encodes the domain\n"
<< "        intron                - intron outside the domain genomic span\n"
<< "        intron_domain_span    - intron between two CDS_domain rows\n"
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
<< "  Records the tool version, output_kind, annotation source, index format\n"
<< "  version, coordinate conventions, mapped/unmapped counts, timestamp,\n"
<< "  and the full CLI invocation.\n"
<< "\n"
<< "BED INPUT FORMAT\n"
<< "  Whitespace-separated. Lines starting with '#' are comments.\n"
<< "    column 1   protein_id (with or without version suffix)\n"
<< "    column 2   aa_start (1-based inclusive)\n"
<< "    column 3   aa_end (1-based inclusive)\n"
<< "    column 4   domain_id (optional)\n"
<< "    column 5+  ignored\n"
<< "  protein_id versions are stripped on both sides, so a BED with\n"
<< "  ENSP00000269305.4 will resolve against an index built from a GTF\n"
<< "  with ENSP00000269305 (and vice versa).\n"
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
        ;
}

void print_version() {
    std::cout << "Protein2Genomic 2.0.0 (index format v" << INDEX_FORMAT_VERSION << ")\n";
}

bool parse_output_kind(const std::string& s, OutputKind& out) {
    if (s == "coding")  { out = OutputKind::CODING;  return true; }
    if (s == "span")    { out = OutputKind::SPAN;    return true; }
    if (s == "isoform") { out = OutputKind::ISOFORM; return true; }
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
                    std::cerr << "Error: --output must be one of {coding,span,isoform,all}\n";
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
                      << "Domains processed: " << mapper.domain_count() << "\n"
                      << "Mapped: " << mapper.mapped_count()
                      << " / Unmapped: " << mapper.unmapped_count() << "\n";
            MemoryTracker::print_memory_stats();
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}

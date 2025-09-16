#include "common.hpp"
#include "gtf_parser.hpp"
#include "domain_mapper.hpp"
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
    std::string output_file;
    std::string output_format = "simple";  // simple|detailed
    std::string mode = "basic";             // basic|full
    bool build_index = false;
    bool use_index = false;
    bool verbose = false;
    int threads = 1;
};


void print_usage(const char* program_name) {
    std::cout << "Protein2Genomic - High-Performance Domain Mapping Tool\n"
              << "Usage: " << program_name << " [OPTIONS]\n\n"
              << "Options:\n"
              << "  --gtf FILE         GTF annotation file\n"
              << "  --bed FILE         BED file with protein domains (protein_id start end [domain_id])\n"
              << "  --output FILE      Output file (default: stdout)\n"
              << "  --mode MODE        Mapping mode: basic|full (default: basic)\n"
              << "  --format FORMAT    Output format: simple|detailed (default: simple)\n"
              << "  --build-index      Build binary index from GTF file\n"
              << "  --index FILE       Use pre-built binary index\n"
              << "  --threads NUM      Number of threads (default: 1)\n"
              << "  --verbose          Verbose output\n"
              << "  --help             Show this help message\n\n"
              << "Mapping Modes:\n"
              << "  basic  - CDS regions only within domain boundaries\n"
              << "  full   - CDS and introns within domain boundaries\n\n"
              << "Output Formats:\n"
              << "  simple   - Only regions within domain boundaries\n"
              << "  detailed - All transcript CDS with domain overlap annotation (Yes/No)\n\n"
              << "Examples:\n"
              << "  # Build index (one-time setup)\n"
              << "  " << program_name << " --gtf annotations.gtf --build-index --index annotations.idx\n\n"
              << "  # Basic mapping (CDS only within domain)\n"
              << "  " << program_name << " --index annotations.idx --bed domains.bed --mode basic --format simple\n\n"
              << "  # Full mapping with introns\n"
              << "  " << program_name << " --index annotations.idx --bed domains.bed --mode full --format simple\n\n"
              << "  # Detailed view of all transcript CDS\n"
              << "  " << program_name << " --index annotations.idx --bed domains.bed --format detailed\n";
}

void print_version() {
    std::cout << "Protein2Genomic version 1.0.0\n"
              << "Built with C++17, optimized for performance\n";
}

int parse_arguments(int argc, char* argv[], Config& config) {
    static struct option long_options[] = {
        {"gtf",         required_argument, 0, 'g'},
        {"bed",         required_argument, 0, 'b'},
        {"output",      required_argument, 0, 'o'},
        {"format",      required_argument, 0, 'f'},
        {"mode",        required_argument, 0, 'm'},
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
    
    while ((c = getopt_long(argc, argv, "g:b:o:f:m:ix:t:vhV", long_options, &option_index)) != -1) {
        switch (c) {
            case 'g':
                config.gtf_file = optarg;
                break;
            case 'b':
                config.bed_file = optarg;
                break;
            case 'o':
                config.output_file = optarg;
                break;
            case 'f':
                config.output_format = optarg;
                if (config.output_format != "simple" && config.output_format != "detailed") {
                    std::cerr << "Error: Invalid output format. Use 'simple' or 'detailed'\n";
                    return 1;
                }
                break;
            case 'm':
                config.mode = optarg;
                if (config.mode != "basic" && config.mode != "full") {
                    std::cerr << "Error: Invalid mode. Use 'basic' or 'full'\n";
                    return 1;
                }
                break;
            case 'i':
                config.build_index = true;
                break;
            case 'x':
                config.index_file = optarg;
                config.use_index = true;
                break;
            case 't':
                config.threads = std::stoi(optarg);
                if (config.threads < 1) {
                    std::cerr << "Error: Number of threads must be positive\n";
                    return 1;
                }
                break;
            case 'v':
                config.verbose = true;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            case 'V':
                print_version();
                return 0;
            case '?':
                return 1;
            default:
                return 1;
        }
    }
    
    return -1; // Continue execution
}

int validate_config(const Config& config) {
    // Check required files
    if (!config.build_index && !config.use_index && config.gtf_file.empty()) {
        std::cerr << "Error: GTF file is required (use --gtf or --index)\n";
        return 1;
    }
    
    if (!config.build_index && config.bed_file.empty()) {
        std::cerr << "Error: BED file is required (use --bed)\n";
        return 1;
    }
    
    if (config.build_index && config.index_file.empty()) {
        std::cerr << "Error: Index file is required when building index (use --index)\n";
        return 1;
    }
    
    // Check file existence
    if (!config.gtf_file.empty() && !utils::file_exists(config.gtf_file)) {
        std::cerr << "Error: GTF file not found: " << config.gtf_file << std::endl;
        return 1;
    }
    
    if (!config.bed_file.empty() && !utils::file_exists(config.bed_file)) {
        std::cerr << "Error: BED file not found: " << config.bed_file << std::endl;
        return 1;
    }
    
    // CORREGIR ESTA LÍNEA: Solo verificar si existe cuando NO estamos construyendo el índice
    if (config.use_index && !config.build_index && !utils::file_exists(config.index_file)) {
        std::cerr << "Error: Index file not found: " << config.index_file << std::endl;
        return 1;
    }
    
    return 0;
}


int main(int argc, char* argv[]) {
    Config config;
    
    // Parse command line arguments
    int parse_result = parse_arguments(argc, argv, config);
    if (parse_result >= 0) {
        return parse_result;
    }
    
    // Validate configuration
    if (validate_config(config) != 0) {
        return 1;
    }
    
    try {
        GTFParser parser;
        
        // Configure feature types based on mode
        FeatureType requested_features;
        if (config.mode == "basic") {
            requested_features = FeatureType::CDS;
        } else { // full
            requested_features = FeatureType::CDS | FeatureType::INTRON;
        }
        
        if (config.verbose) {
            std::cerr << "Mode: " << config.mode << std::endl;
            std::cerr << "Format: " << config.output_format << std::endl;
        }
        
        // Set number of threads for OpenMP
#ifdef USE_OPENMP
        omp_set_num_threads(config.threads);
        if (config.verbose) {
            std::cerr << "Using " << config.threads << " threads" << std::endl;
        }
#endif
        
        // Build or load index
        if (config.build_index) {
            if (config.verbose) {
                std::cerr << "Building universal index from GTF file..." << std::endl;
            }
            
            Timer timer("Index building");
            // Always build complete index with all features for maximum compatibility
            FeatureType index_features = FeatureType::CDS | FeatureType::INTRON;
            ErrorCode result = parser.load_gtf(config.gtf_file, index_features, true);
            if (result != ErrorCode::SUCCESS) {
                std::cerr << "Error building index: " << utils::error_code_to_string(result) << std::endl;
                return 1;
            }
            
            result = parser.save_index(config.index_file);
            if (result != ErrorCode::SUCCESS) {
                std::cerr << "Error saving index: " << utils::error_code_to_string(result) << std::endl;
                return 1;
            }
            
            std::cerr << "Index saved to: " << config.index_file << std::endl;
            MemoryTracker::print_memory_stats();
            
            if (config.bed_file.empty()) {
                return 0; // Only building index
            }
        } else if (config.use_index) {
            if (config.verbose) {
                std::cerr << "Loading pre-built index..." << std::endl;
            }
            
            Timer timer("Index loading");
            ErrorCode result = parser.load_index(config.index_file);
            if (result != ErrorCode::SUCCESS) {
                std::cerr << "Error loading index: " << utils::error_code_to_string(result) << std::endl;
                return 1;
            }
        } else {
            if (config.verbose) {
                std::cerr << "Parsing GTF file directly..." << std::endl;
            }
            
            Timer timer("GTF parsing");
            ErrorCode result = parser.load_gtf(config.gtf_file, requested_features, true);
            if (result != ErrorCode::SUCCESS) {
                std::cerr << "Error parsing GTF: " << utils::error_code_to_string(result) << std::endl;
                return 1;
            }
        }
        
        // Process domains
        DomainMapper mapper(parser, config.mode, config.output_format);
        
        {
            Timer timer("Domain loading");
            ErrorCode result = mapper.load_domains(config.bed_file);
            if (result != ErrorCode::SUCCESS) {
                std::cerr << "Error loading domains: " << utils::error_code_to_string(result) << std::endl;
                return 1;
            }
        }
        
        {
            Timer timer("Domain mapping");
            ErrorCode result = mapper.process_domains();
            if (result != ErrorCode::SUCCESS) {
                std::cerr << "Error processing domains: " << utils::error_code_to_string(result) << std::endl;
                return 1;
            }
        }
        
        // Write results
        std::string output_file = config.output_file.empty() ? "/dev/stdout" : config.output_file;
        
        {
            Timer timer("Results writing");
            ErrorCode result = mapper.write_results(output_file, config.output_format);
            if (result != ErrorCode::SUCCESS) {
                std::cerr << "Error writing results: " << utils::error_code_to_string(result) << std::endl;
                return 1;
            }
        }
        
        // Print final statistics
        if (config.verbose) {
            std::cerr << "\n=== Final Statistics ===" << std::endl;
            std::cerr << "Proteins in index: " << parser.get_protein_count() << std::endl;
            std::cerr << "Domains processed: " << mapper.get_domain_count() << std::endl;
            std::cerr << "Genomic intervals generated: " << mapper.get_mapped_count() << std::endl;
            std::cerr << "Unmapped domains: " << mapper.get_unmapped_count() << std::endl;
            MemoryTracker::print_memory_stats();
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

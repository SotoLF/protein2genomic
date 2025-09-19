#include "gtf_parser.hpp"
#include "utils.hpp"
#include <fstream>
#include <sstream>
#include <regex>
#include <algorithm>
#include <set>

GTFParser::GTFParser() : lines_processed_(0), requested_features_(FeatureType::CDS), resolve_overlaps_(true) {
    // Pre-allocate containers for better performance
    protein_index_.reserve(INITIAL_PROTEIN_CAPACITY);
    protein_to_transcript_.reserve(INITIAL_PROTEIN_CAPACITY);
    temp_fields_.reserve(20); // GTF has 9 fields + attributes
}

ErrorCode GTFParser::load_gtf(const std::string& filename, 
                      FeatureType features,
                      bool resolve_overlaps) {
    requested_features_ = features;
    resolve_overlaps_ = resolve_overlaps;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open GTF file: " << filename << std::endl;
        return ErrorCode::FILE_NOT_FOUND;
    }
    
    current_file_ = filename;
    lines_processed_ = 0;
    
    std::string line;
    line.reserve(MAX_LINE_LENGTH);
    
    GTFRecord record;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    while (std::getline(file, line)) {
        ++lines_processed_;
        
        // Skip header lines and empty lines
        if (line.empty() || line[0] == '#') continue;
        
        if (parse_line(line, record)) {
            if (is_relevant_feature(record.feature)) {
                process_record(record);
            }
        }
        
        // Progress indicator every 100k lines
        if (lines_processed_ % 100000 == 0) {
            std::cerr << "Processed " << lines_processed_ << " lines...\r" << std::flush;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cerr << "\nCompleted parsing GTF file in " << duration.count() << "ms" << std::endl;
    std::cerr << "Total lines processed: " << lines_processed_ << std::endl;
    std::cerr << "Proteins indexed: " << protein_index_.size() << std::endl;
    std::cerr << "Estimated memory usage: " << estimate_memory_usage() / (1024*1024) << " MB" << std::endl;
    
    // If features beyond CDS are requested, build gene structures
    if (requested_features_ != FeatureType::CDS) {
        std::cerr << "Building gene structures for extended features..." << std::endl;
        auto structure_start = std::chrono::high_resolution_clock::now();
        
        build_gene_structures();
        build_protein_index_from_structures();
        
        auto structure_end = std::chrono::high_resolution_clock::now();
        auto structure_duration = std::chrono::duration_cast<std::chrono::milliseconds>(structure_end - structure_start);
        std::cerr << "Gene structure building completed in " << structure_duration.count() << "ms" << std::endl;
        std::cerr << "Gene structures: " << gene_structures_.size() << std::endl;
    }
    
    // Sort intervals for each protein for better cache locality
    // For negative strand proteins, sort in reverse order to maintain translation order
    for (auto& [protein_id, intervals] : protein_index_) {
        std::sort(intervals.begin(), intervals.end());
        
        // For negative strand, reverse the order to match translation direction
        if (!intervals.empty() && intervals[0].strand == '-') {
            std::reverse(intervals.begin(), intervals.end());
        }
    }
    
    return ErrorCode::SUCCESS;
}

bool GTFParser::parse_line(const std::string& line, GTFRecord& record) const {
    temp_fields_.clear();
    temp_fields_ = split_fast(line, '\t');
    
    if (temp_fields_.size() < 9) {
        return false;
    }
    
    try {
        record.chromosome = temp_fields_[0];
        record.source = temp_fields_[1];
        record.feature = temp_fields_[2];
        record.start = std::stoul(temp_fields_[3]);
        record.end = std::stoul(temp_fields_[4]);
        record.score = temp_fields_[5];
        record.strand = temp_fields_[6][0];
        record.phase = (temp_fields_[7] == ".") ? 0 : std::stoi(temp_fields_[7]);
        record.attributes = temp_fields_[8];
        
        return true;
    } catch (const std::exception& e) {
        return false;
    }
}

std::vector<std::string> GTFParser::split_fast(const std::string& str, char delimiter) const {
    std::vector<std::string> result;
    result.reserve(12); // Most GTF lines have ~9-12 fields
    
    size_t start = 0;
    size_t end = str.find(delimiter);
    
    while (end != std::string::npos) {
        result.emplace_back(str.substr(start, end - start));
        start = end + 1;
        end = str.find(delimiter, start);
    }
    result.emplace_back(str.substr(start));
    
    return result;
}

std::string GTFParser::extract_attribute(const std::string& attributes, const std::string& key) const {
    // Fast attribute extraction using string search
    std::string search_key = key + " \"";
    size_t start = attributes.find(search_key);
    
    if (start == std::string::npos) {
        return "";
    }
    
    start += search_key.length();
    size_t end = attributes.find("\"", start);
    
    if (end == std::string::npos) {
        return "";
    }
    
    return attributes.substr(start, end - start);
}

void GTFParser::process_record(const GTFRecord& record) {
    // Process ALL relevant features, not just those with protein_id
    if (!is_relevant_feature(record.feature)) {
        return;
    }
    
    // Convert string feature to FeatureType
    FeatureType feature_type = string_to_feature_type(record.feature);
    
    // Extract gene_id and transcript_id
    std::string gene_id = extract_attribute(record.attributes, "gene_id");
    std::string transcript_id = extract_attribute(record.attributes, "transcript_id");
    std::string protein_id = extract_attribute(record.attributes, "protein_id");
    
    // Skip if no gene_id or transcript_id
    if (gene_id.empty() || transcript_id.empty()) {
        return;
    }
    
    // Create genomic interval
    GenomicInterval interval(
        record.chromosome,
        record.start,
        record.end,
        record.strand,
        record.phase,
        feature_type,
        transcript_id,
        gene_id
    );
    
    // Store in temp_transcript_data_ for later processing
    temp_transcript_data_[transcript_id].push_back(interval);
    
    // If this record has protein_id, store the protein -> transcript mapping
    if (!protein_id.empty()) {
        protein_to_transcript_[protein_id] = transcript_id;
        
        // For CDS features with protein_id, also add directly to protein_index_ (fallback)
        if (feature_type == FeatureType::CDS) {
            protein_index_[protein_id].push_back(interval);
        }
    }
}

bool GTFParser::is_relevant_feature(const std::string& feature) const {
    // Always include exon and CDS for basic functionality
    if (feature == "exon" || feature == "CDS") {
        return true;
    }
    
    // Include additional features if requested
    if (requested_features_ != FeatureType::CDS) {
        return (feature == "start_codon" || feature == "stop_codon");
    }
    
    return false;
}

const std::vector<GenomicInterval>* GTFParser::get_protein_intervals(const std::string& protein_id) const {
    auto it = protein_index_.find(protein_id);
    return (it != protein_index_.end()) ? &it->second : nullptr;
}

const std::string* GTFParser::get_transcript_id(const std::string& protein_id) const {
    auto it = protein_to_transcript_.find(protein_id);
    return (it != protein_to_transcript_.end()) ? &it->second : nullptr;
}

size_t GTFParser::estimate_memory_usage() const {
    size_t total = 0;
    
    // Protein index
    for (const auto& [protein_id, intervals] : protein_index_) {
        total += protein_id.size();
        total += intervals.size() * sizeof(GenomicInterval);
        for (const auto& interval : intervals) {
            total += interval.chromosome.size();
            // feature_type is an enum, not a string, so no need to add its size
        }
    }
    
    // Protein to transcript mapping
    for (const auto& [protein_id, transcript_id] : protein_to_transcript_) {
        total += protein_id.size() + transcript_id.size();
    }
    
    return total;
}

ErrorCode GTFParser::save_index(const std::string& filename) const {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return ErrorCode::FILE_NOT_FOUND;
    }
    
    // Write header
    uint32_t num_proteins = protein_index_.size();
    file.write(reinterpret_cast<const char*>(&num_proteins), sizeof(num_proteins));
    
    // Write protein data
    for (const auto& [protein_id, intervals] : protein_index_) {
        // Write protein_id length and data
        uint16_t id_len = protein_id.length();
        file.write(reinterpret_cast<const char*>(&id_len), sizeof(id_len));
        file.write(protein_id.c_str(), id_len);
        
        // Write number of intervals
        uint32_t num_intervals = intervals.size();
        file.write(reinterpret_cast<const char*>(&num_intervals), sizeof(num_intervals));
        
        // Write intervals
        for (const auto& interval : intervals) {
            uint16_t chr_len = interval.chromosome.length();
            file.write(reinterpret_cast<const char*>(&chr_len), sizeof(chr_len));
            file.write(interval.chromosome.c_str(), chr_len);
            
            file.write(reinterpret_cast<const char*>(&interval.start), sizeof(interval.start));
            file.write(reinterpret_cast<const char*>(&interval.end), sizeof(interval.end));
            file.write(reinterpret_cast<const char*>(&interval.strand), sizeof(interval.strand));
            file.write(reinterpret_cast<const char*>(&interval.phase), sizeof(interval.phase));
            
            // Write feature_type as enum value
            file.write(reinterpret_cast<const char*>(&interval.feature_type), sizeof(interval.feature_type));
            
            // Write transcript_id and gene_id lengths and data
            uint16_t transcript_len = interval.transcript_id.length();
            file.write(reinterpret_cast<const char*>(&transcript_len), sizeof(transcript_len));
            file.write(interval.transcript_id.c_str(), transcript_len);
            
            uint16_t gene_len = interval.gene_id.length();
            file.write(reinterpret_cast<const char*>(&gene_len), sizeof(gene_len));
            file.write(interval.gene_id.c_str(), gene_len);
        }
    }
    
    // Save gene structures for full mode support
    uint32_t num_genes = gene_structures_.size();
    file.write(reinterpret_cast<const char*>(&num_genes), sizeof(num_genes));
    
    for (const auto& [gene_id, gene_structure] : gene_structures_) {
        // Write gene_id
        uint16_t gene_id_len = gene_id.size();
        file.write(reinterpret_cast<const char*>(&gene_id_len), sizeof(gene_id_len));
        file.write(gene_id.c_str(), gene_id_len);
        
        // Write gene structure data
        uint16_t chr_len = gene_structure.chromosome.size();
        file.write(reinterpret_cast<const char*>(&chr_len), sizeof(chr_len));
        file.write(gene_structure.chromosome.c_str(), chr_len);
        
        file.write(reinterpret_cast<const char*>(&gene_structure.gene_start), sizeof(gene_structure.gene_start));
        file.write(reinterpret_cast<const char*>(&gene_structure.gene_end), sizeof(gene_structure.gene_end));
        file.write(reinterpret_cast<const char*>(&gene_structure.strand), sizeof(gene_structure.strand));
        file.write(reinterpret_cast<const char*>(&gene_structure.cds_start), sizeof(gene_structure.cds_start));
        file.write(reinterpret_cast<const char*>(&gene_structure.cds_end), sizeof(gene_structure.cds_end));
        
        // Write protein_id
        uint16_t protein_id_len = gene_structure.protein_id.size();
        file.write(reinterpret_cast<const char*>(&protein_id_len), sizeof(protein_id_len));
        file.write(gene_structure.protein_id.c_str(), protein_id_len);
        
        // Write exons
        uint32_t num_exons = gene_structure.exons.size();
        file.write(reinterpret_cast<const char*>(&num_exons), sizeof(num_exons));
        for (const auto& exon : gene_structure.exons) {
            file.write(reinterpret_cast<const char*>(&exon.start), sizeof(exon.start));
            file.write(reinterpret_cast<const char*>(&exon.end), sizeof(exon.end));
        }
        
        // Write CDS regions
        uint32_t num_cds = gene_structure.cds_regions.size();
        file.write(reinterpret_cast<const char*>(&num_cds), sizeof(num_cds));
        for (const auto& cds : gene_structure.cds_regions) {
            file.write(reinterpret_cast<const char*>(&cds.start), sizeof(cds.start));
            file.write(reinterpret_cast<const char*>(&cds.end), sizeof(cds.end));
        }
    }    return ErrorCode::SUCCESS;
}

ErrorCode GTFParser::load_index(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return ErrorCode::FILE_NOT_FOUND;
    }
    
    clear();
    
    uint32_t num_proteins;
    file.read(reinterpret_cast<char*>(&num_proteins), sizeof(num_proteins));
    
    protein_index_.reserve(num_proteins);
    
    for (uint32_t i = 0; i < num_proteins; ++i) {
        // Read protein_id
        uint16_t id_len;
        file.read(reinterpret_cast<char*>(&id_len), sizeof(id_len));
        std::string protein_id(id_len, '\0');
        file.read(&protein_id[0], id_len);
        
        // Read intervals
        uint32_t num_intervals;
        file.read(reinterpret_cast<char*>(&num_intervals), sizeof(num_intervals));
        
        std::vector<GenomicInterval> intervals;
        intervals.reserve(num_intervals);
        
        for (uint32_t j = 0; j < num_intervals; ++j) {
            GenomicInterval interval;
            
            uint16_t chr_len;
            file.read(reinterpret_cast<char*>(&chr_len), sizeof(chr_len));
            interval.chromosome.resize(chr_len);
            file.read(&interval.chromosome[0], chr_len);
            
            file.read(reinterpret_cast<char*>(&interval.start), sizeof(interval.start));
            file.read(reinterpret_cast<char*>(&interval.end), sizeof(interval.end));
            file.read(reinterpret_cast<char*>(&interval.strand), sizeof(interval.strand));
            file.read(reinterpret_cast<char*>(&interval.phase), sizeof(interval.phase));
            
            // Read feature_type as enum value
            file.read(reinterpret_cast<char*>(&interval.feature_type), sizeof(interval.feature_type));
            
            // Read transcript_id and gene_id
            uint16_t transcript_len;
            file.read(reinterpret_cast<char*>(&transcript_len), sizeof(transcript_len));
            interval.transcript_id.resize(transcript_len);
            file.read(&interval.transcript_id[0], transcript_len);
            
            uint16_t gene_len;
            file.read(reinterpret_cast<char*>(&gene_len), sizeof(gene_len));
            interval.gene_id.resize(gene_len);
            file.read(&interval.gene_id[0], gene_len);
            
            intervals.emplace_back(std::move(interval));
        }
        
        protein_index_[protein_id] = std::move(intervals);
    }
    
    // Read gene structures for full mode support
    uint32_t num_genes;
    file.read(reinterpret_cast<char*>(&num_genes), sizeof(num_genes));
    
    gene_structures_.reserve(num_genes);
    
    for (uint32_t i = 0; i < num_genes; ++i) {
        // Read gene_id
        uint16_t gene_id_len;
        file.read(reinterpret_cast<char*>(&gene_id_len), sizeof(gene_id_len));
        std::string gene_id(gene_id_len, '\0');
        file.read(&gene_id[0], gene_id_len);
        
        GeneStructure gene_structure;
        gene_structure.gene_id = gene_id;
        
        // Read gene structure data
        uint16_t chr_len;
        file.read(reinterpret_cast<char*>(&chr_len), sizeof(chr_len));
        gene_structure.chromosome.resize(chr_len);
        file.read(&gene_structure.chromosome[0], chr_len);
        
        file.read(reinterpret_cast<char*>(&gene_structure.gene_start), sizeof(gene_structure.gene_start));
        file.read(reinterpret_cast<char*>(&gene_structure.gene_end), sizeof(gene_structure.gene_end));
        file.read(reinterpret_cast<char*>(&gene_structure.strand), sizeof(gene_structure.strand));
        file.read(reinterpret_cast<char*>(&gene_structure.cds_start), sizeof(gene_structure.cds_start));
        file.read(reinterpret_cast<char*>(&gene_structure.cds_end), sizeof(gene_structure.cds_end));
        
        // Read protein_id
        uint16_t protein_id_len;
        file.read(reinterpret_cast<char*>(&protein_id_len), sizeof(protein_id_len));
        gene_structure.protein_id.resize(protein_id_len);
        file.read(&gene_structure.protein_id[0], protein_id_len);
        
        // Read exons
        uint32_t num_exons;
        file.read(reinterpret_cast<char*>(&num_exons), sizeof(num_exons));
        gene_structure.exons.reserve(num_exons);
        
        for (uint32_t j = 0; j < num_exons; ++j) {
            GenomicInterval exon;
            file.read(reinterpret_cast<char*>(&exon.start), sizeof(exon.start));
            file.read(reinterpret_cast<char*>(&exon.end), sizeof(exon.end));
            gene_structure.exons.emplace_back(exon);
        }
        
        // Read CDS regions
        uint32_t num_cds;
        file.read(reinterpret_cast<char*>(&num_cds), sizeof(num_cds));
        gene_structure.cds_regions.reserve(num_cds);
        
        for (uint32_t j = 0; j < num_cds; ++j) {
            GenomicInterval cds;
            file.read(reinterpret_cast<char*>(&cds.start), sizeof(cds.start));
            file.read(reinterpret_cast<char*>(&cds.end), sizeof(cds.end));
            gene_structure.cds_regions.emplace_back(cds);
        }
        
        gene_structures_[gene_id] = std::move(gene_structure);
    }
    
    return ErrorCode::SUCCESS;
}
std::string feature_type_to_string(FeatureType type) {
    switch (type) {
        case FeatureType::EXON: return "exon";
        case FeatureType::CDS: return "CDS";
        case FeatureType::INTRON: return "intron";
        case FeatureType::START_CODON: return "start_codon";
        case FeatureType::STOP_CODON: return "stop_codon";
        default: return "unknown";
    }
}

FeatureType string_to_feature_type(const std::string& str) {
    if (str == "exon") return FeatureType::EXON;
    if (str == "CDS" || str == "cds") return FeatureType::CDS;
    if (str == "intron") return FeatureType::INTRON;
    if (str == "start_codon") return FeatureType::START_CODON;
    if (str == "stop_codon") return FeatureType::STOP_CODON;
    if (str == "all") return FeatureType::ALL;
    return FeatureType::EXON; // default
}

FeatureType parse_feature_types(const std::string& types_str) {
    if (types_str == "all") {
        return FeatureType::ALL;
    }
    
    FeatureType result = static_cast<FeatureType>(0);
    std::istringstream iss(types_str);
    std::string type;
    
    while (std::getline(iss, type, ',')) {
        // Trim whitespace
        type.erase(0, type.find_first_not_of(" \t"));
        type.erase(type.find_last_not_of(" \t") + 1);
        
        if (type == "exon") {
            result = result | FeatureType::EXON;
        } else if (type == "cds" || type == "CDS") {
            result = result | FeatureType::CDS;
        } else if (type == "intron") {
            result = result | FeatureType::INTRON;
        } else if (type == "start_codon") {
            result = result | FeatureType::START_CODON;
        } else if (type == "stop_codon") {
            result = result | FeatureType::STOP_CODON;
        }
    }
    
    // Default to CDS if nothing specified
    if (static_cast<int>(result) == 0) {
        result = FeatureType::CDS;
    }
    
    return result;
}

std::string get_feature_type_name(FeatureType type) {
    switch (type) {
        case FeatureType::EXON: return "exon";
        case FeatureType::CDS: return "CDS";
        case FeatureType::INTRON: return "intron";
        case FeatureType::START_CODON: return "start_codon";
        case FeatureType::STOP_CODON: return "stop_codon";
        default: return "unknown";
    }
}

// Implementar métodos de GeneStructure

std::vector<GenomicInterval> GeneStructure::get_introns() const {
    std::vector<GenomicInterval> introns;
    
    if (cds_regions.size() < 2) {
        return introns; // Need at least 2 CDS for introns
    }
    
    // Sort CDS by position
    auto sorted_cds = cds_regions;
    std::sort(sorted_cds.begin(), sorted_cds.end());
    
    // Generate introns between consecutive CDS
    for (size_t i = 0; i < sorted_cds.size() - 1; ++i) {
        uint32_t intron_start = sorted_cds[i].end + 1;
        uint32_t intron_end = sorted_cds[i + 1].start - 1;
        
        if (intron_start <= intron_end) { // Valid intron
            introns.emplace_back(chromosome, intron_start, intron_end, 
                               strand, 0, FeatureType::INTRON,
                               sorted_cds[i].transcript_id, gene_id);
        }
    }
    
    return introns;
}


// GTFParser methods for building gene structures
void GTFParser::build_gene_structures() {
    // Simplified: we don't need complex gene structures anymore
    // The heavy lifting is done in build_protein_index_from_structures()
    gene_structures_.clear();
    
    // We still need some basic gene structures for compatibility
    // Group by gene_id for basic statistics
    std::unordered_map<std::string, std::string> gene_to_first_transcript;
    
    for (const auto& [transcript_id, intervals] : temp_transcript_data_) {
        if (!intervals.empty()) {
            const std::string& gene_id = intervals[0].gene_id;
            if (gene_to_first_transcript.find(gene_id) == gene_to_first_transcript.end()) {
                gene_to_first_transcript[gene_id] = transcript_id;
                
                // Create minimal gene structure for stats
                GeneStructure gene_struct;
                gene_struct.gene_id = gene_id;
                gene_struct.chromosome = intervals[0].chromosome;
                gene_struct.strand = intervals[0].strand;
                
                // Calculate bounds
                auto minmax = std::minmax_element(intervals.begin(), intervals.end(),
                    [](const GenomicInterval& a, const GenomicInterval& b) {
                        return a.start < b.start;
                    });
                gene_struct.gene_start = minmax.first->start;
                gene_struct.gene_end = std::max_element(intervals.begin(), intervals.end(),
                    [](const GenomicInterval& a, const GenomicInterval& b) {
                        return a.end < b.end;
                    })->end;
                
                gene_structures_[gene_id] = std::move(gene_struct);
            }
        }
    }
}

void GTFParser::build_protein_index_from_structures() {
    // For each protein, rebuild its intervals from the specific transcript's features
    for (const auto& [protein_id, transcript_id] : protein_to_transcript_) {
        // Find the transcript's intervals
        auto transcript_it = temp_transcript_data_.find(transcript_id);
        if (transcript_it == temp_transcript_data_.end()) {
            continue; // No data for this transcript
        }
        
        const auto& transcript_intervals = transcript_it->second;
        
        // Clear existing protein intervals and rebuild
        auto& protein_intervals = protein_index_[protein_id];
        protein_intervals.clear();
        
        // Collect features by type from THIS transcript only
        std::vector<GenomicInterval> exons, cdss, start_codons, stop_codons;
        
        for (const auto& interval : transcript_intervals) {
            switch (interval.feature_type) {
                case FeatureType::EXON:
                    exons.push_back(interval);
                    break;
                case FeatureType::CDS:
                    cdss.push_back(interval);
                    break;
                case FeatureType::START_CODON:
                    start_codons.push_back(interval);
                    break;
                case FeatureType::STOP_CODON:
                    stop_codons.push_back(interval);
                    break;
                default:
                    break;
            }
        }
        
        // Sort exons by position for proper intron calculation
        std::sort(exons.begin(), exons.end(), 
            [](const GenomicInterval& a, const GenomicInterval& b) {
                return a.start < b.start;
            });
        
        // Add requested features to protein intervals
        if (has_feature(requested_features_, FeatureType::EXON)) {
            protein_intervals.insert(protein_intervals.end(), exons.begin(), exons.end());
        }
        if (has_feature(requested_features_, FeatureType::CDS)) {
            protein_intervals.insert(protein_intervals.end(), cdss.begin(), cdss.end());
        }
        if (has_feature(requested_features_, FeatureType::START_CODON)) {
            protein_intervals.insert(protein_intervals.end(), start_codons.begin(), start_codons.end());
        }
        if (has_feature(requested_features_, FeatureType::STOP_CODON)) {
            protein_intervals.insert(protein_intervals.end(), stop_codons.begin(), stop_codons.end());
        }
        
        // Generate introns if requested and we have multiple CDS
        if (has_feature(requested_features_, FeatureType::INTRON) && cdss.size() > 1) {
            // Sort CDS by position for proper intron calculation
            auto sorted_cds = cdss;
            std::sort(sorted_cds.begin(), sorted_cds.end(), 
                [](const GenomicInterval& a, const GenomicInterval& b) {
                    return a.start < b.start;
                });
            
            for (size_t i = 0; i < sorted_cds.size() - 1; ++i) {
                uint32_t intron_start = sorted_cds[i].end + 1;
                uint32_t intron_end = sorted_cds[i + 1].start - 1;
                
                // Only create intron if coordinates are valid
                if (intron_start <= intron_end) {
                    GenomicInterval intron(
                        sorted_cds[i].chromosome,
                        intron_start,
                        intron_end,
                        sorted_cds[i].strand,
                        -1, // No phase for introns
                        FeatureType::INTRON,
                        transcript_id,
                        sorted_cds[i].gene_id
                    );
                    protein_intervals.push_back(intron);
                }
            }
        }
        
        // Apply overlap resolution if requested
        if (resolve_overlaps_) {
            protein_intervals = resolve_overlapping_features(protein_intervals);
        }
        
        // Sort final intervals
        std::sort(protein_intervals.begin(), protein_intervals.end());
        
        // If no intervals were found, don't remove the protein from the index
        // This can happen if the requested features don't exist for this protein
        if (protein_intervals.empty()) {
            // Keep the protein in the index but with empty intervals
            // This allows the mapping to still report it as unmapped rather than not found
        }
    }
    
    // Clear temporary data to save memory
    temp_transcript_data_.clear();
}

std::vector<GenomicInterval> GTFParser::resolve_overlapping_features(
    const std::vector<GenomicInterval>& intervals) const {
    
    if (intervals.empty()) return intervals;
    
    std::vector<GenomicInterval> resolved;
    auto sorted_intervals = intervals;
    std::sort(sorted_intervals.begin(), sorted_intervals.end());
    
    for (const auto& interval : sorted_intervals) {
        bool overlaps = false;
        
        // Check for overlaps with existing intervals
        for (auto& existing : resolved) {
            if (existing.chromosome == interval.chromosome &&
                existing.start <= interval.end && existing.end >= interval.start) {
                
                // Overlap detected - keep higher priority feature
                if (interval.priority > existing.priority) {
                    // Replace existing with higher priority interval
                    existing = interval;
                }
                overlaps = true;
                break;
            }
        }
        
        if (!overlaps) {
            resolved.push_back(interval);
        }
    }
    
    return resolved;
}

const GeneStructure* GTFParser::get_gene_structure(const std::string& gene_id) const {
    auto it = gene_structures_.find(gene_id);
    return (it != gene_structures_.end()) ? &it->second : nullptr;
}


#ifndef GTF_PARSER_HPP
#define GTF_PARSER_HPP

#include "common.hpp"

struct GenomicInterval {
    std::string chromosome;
    uint32_t start;
    uint32_t end;
    char strand;
    uint8_t phase;
    FeatureType feature_type;
    FeaturePriority priority;
    std::string transcript_id;
    std::string gene_id;
    
    GenomicInterval() : start(0), end(0), strand('+'), phase(0), 
                       feature_type(FeatureType::EXON), priority(FeaturePriority::EXON) {}
    
    GenomicInterval(const std::string& chr, uint32_t s, uint32_t e, 
                   char str, uint8_t ph, FeatureType type, 
                   const std::string& tx_id = "", const std::string& g_id = "")
        : chromosome(chr), start(s), end(e), strand(str), phase(ph), 
          feature_type(type), transcript_id(tx_id), gene_id(g_id) {
        set_priority_from_type();
    }
    
    void set_priority_from_type() {
        switch (feature_type) {
            case FeatureType::INTRON: priority = FeaturePriority::INTRON; break;
            case FeatureType::EXON: priority = FeaturePriority::EXON; break;
            case FeatureType::CDS: priority = FeaturePriority::CDS; break;
            case FeatureType::START_CODON: priority = FeaturePriority::START_CODON; break;
            case FeatureType::STOP_CODON: priority = FeaturePriority::STOP_CODON; break;
            default: priority = FeaturePriority::EXON; break;
        }
    }
    
    bool operator<(const GenomicInterval& other) const {
        if (chromosome != other.chromosome) return chromosome < other.chromosome;
        if (start != other.start) return start < other.start;
        if (end != other.end) return end < other.end;
        return static_cast<int>(priority) < static_cast<int>(other.priority);
    }
};

struct GeneStructure {
    std::string gene_id;
    std::string chromosome;
    char strand;
    std::vector<GenomicInterval> cds_regions;  // Only store CDS for intron generation
    std::string protein_id;
    
    GeneStructure() : strand('+') {}
    
    // Generate derived features
    std::vector<GenomicInterval> get_introns() const;
};


struct GTFRecord {
    std::string chromosome;
    std::string source;
    std::string feature;
    uint32_t start;
    uint32_t end;
    std::string score;
    char strand;
    uint8_t phase;
    std::string attributes;
    
    GTFRecord() : start(0), end(0), strand('+'), phase(0) {}
};

class GTFParser {
private:
    ProteinMap protein_index_;
    TranscriptMap protein_to_transcript_;
    GeneMap gene_structures_; // NEW: Store complete gene structures
    FeatureType requested_features_; // NEW: What features to include
    bool resolve_overlaps_; // NEW: Whether to resolve overlapping features
    std::string current_file_;
    size_t lines_processed_;
    
    // Memory optimization: reuse string objects
    mutable std::string temp_string_;
    mutable std::vector<std::string> temp_fields_;
    
    // Temporary storage for building gene structures
    std::unordered_map<std::string, std::vector<GenomicInterval>> temp_transcript_data_;
    
    // Parsing methods
    bool parse_line(const std::string& line, GTFRecord& record) const;
    std::string extract_attribute(const std::string& attributes, const std::string& key) const;
    std::vector<std::string> split_fast(const std::string& str, char delimiter) const;
    
    // Index building methods
    void process_record(const GTFRecord& record);
    bool is_relevant_feature(const std::string& feature) const;
    
    // NEW: Gene structure methods
    void build_gene_structures();
    void build_protein_index_from_structures();
    std::vector<GenomicInterval> resolve_overlapping_features(
        const std::vector<GenomicInterval>& intervals) const;
    
public:
    GTFParser();
    ~GTFParser() = default;
    
    // Main interface
    ErrorCode load_gtf(const std::string& filename, 
                      FeatureType features = FeatureType::CDS,
                      bool resolve_overlaps = true);
    ErrorCode save_index(const std::string& filename) const;
    ErrorCode load_index(const std::string& filename);
    
    // Query interface
    const std::vector<GenomicInterval>* get_protein_intervals(const std::string& protein_id) const;
    const std::string* get_transcript_id(const std::string& protein_id) const;
    const GeneStructure* get_gene_structure(const std::string& gene_id) const; // NEW
    
    // Configuration
    void set_feature_types(FeatureType features) { requested_features_ = features; }
    void set_resolve_overlaps(bool resolve) { resolve_overlaps_ = resolve; }
    
    // Statistics
    size_t get_protein_count() const { return protein_index_.size(); }
    size_t get_gene_count() const { return gene_structures_.size(); }
    size_t get_lines_processed() const { return lines_processed_; }
    
    // Memory management
    void clear() { 
        protein_index_.clear(); 
        protein_to_transcript_.clear(); 
        gene_structures_.clear();
    }
    size_t estimate_memory_usage() const;
};

#endif // GTF_PARSER_HPP

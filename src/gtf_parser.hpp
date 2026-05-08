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
    std::string gene_name;
    uint32_t exon_number; // 0 if not set; populated for exon/CDS via GTF attribute when present

    GenomicInterval() : start(0), end(0), strand('+'), phase(0),
                       feature_type(FeatureType::EXON), priority(FeaturePriority::EXON),
                       exon_number(0) {}

    GenomicInterval(const std::string& chr, uint32_t s, uint32_t e,
                   char str, uint8_t ph, FeatureType type,
                   const std::string& tx_id = "", const std::string& g_id = "",
                   const std::string& g_name = "", uint32_t ex_num = 0)
        : chromosome(chr), start(s), end(e), strand(str), phase(ph),
          feature_type(type), transcript_id(tx_id), gene_id(g_id),
          gene_name(g_name), exon_number(ex_num) {
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
    std::vector<GenomicInterval> cds_regions;
    std::string protein_id;

    GeneStructure() : strand('+') {}
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
    std::unordered_map<std::string, std::string> protein_to_gene_name_;
    std::unordered_map<std::string, std::string> protein_to_gene_id_;
    std::string current_file_;
    size_t lines_processed_;

    mutable std::string temp_string_;
    mutable std::vector<std::string> temp_fields_;

    std::unordered_map<std::string, std::vector<GenomicInterval>> temp_transcript_data_;

    bool parse_line(const std::string& line, GTFRecord& record) const;
    std::string extract_attribute(const std::string& attributes, const std::string& key) const;
    std::vector<std::string> split_fast(const std::string& str, char delimiter) const;

    void process_record(const GTFRecord& record);
    bool is_relevant_feature(const std::string& feature) const;

    void build_protein_index_from_structures();

public:
    GTFParser();
    ~GTFParser() = default;

    ErrorCode load_gtf(const std::string& filename);
    ErrorCode save_index(const std::string& filename) const;
    ErrorCode load_index(const std::string& filename);

    const std::vector<GenomicInterval>* get_protein_intervals(const std::string& protein_id) const;
    const std::string* get_transcript_id(const std::string& protein_id) const;
    std::string get_gene_name(const std::string& protein_id) const;
    std::string get_gene_id(const std::string& protein_id) const;

    size_t get_protein_count() const { return protein_index_.size(); }
    size_t get_lines_processed() const { return lines_processed_; }

    void clear() {
        protein_index_.clear();
        protein_to_transcript_.clear();
        protein_to_gene_name_.clear();
        protein_to_gene_id_.clear();
    }
    size_t estimate_memory_usage() const;
};

#endif // GTF_PARSER_HPP

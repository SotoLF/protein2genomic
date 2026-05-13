#ifndef COMMON_HPP
#define COMMON_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <thread>
#include <mutex>

using ProteinMap = std::unordered_map<std::string, std::vector<class GenomicInterval>>;
using TranscriptMap = std::unordered_map<std::string, std::string>; // protein_id -> transcript_id
using GeneMap = std::unordered_map<std::string, class GeneStructure>;

struct GenomicInterval;
struct ProteinDomain;
struct GTFRecord;
struct GeneStructure;

// GTF/derived feature types stored in the index.
enum class FeatureType {
    EXON = 1,
    CDS = 2,
    INTRON = 4,
    START_CODON = 8,
    STOP_CODON = 16,
    ALL = 31
};

enum class FeaturePriority {
    INTRON = 1,
    EXON = 2,
    CDS = 3,
    START_CODON = 4,
    STOP_CODON = 5
};

// What the user wants the package to produce.
// Selected via --output {coding,introns,span,isoform,bed12,all}.
enum class OutputKind {
    CODING,    // domain_cds_segments.{bed,tsv} — ALL CDS rows + overlap column
    INTRONS,   // domain_introns.{bed,tsv}     — ALL intron rows + overlap column
    SPAN,      // domain_span_with_introns.bed — single envelope row per domain
    ISOFORM,   // isoform_structure.tsv         — UTR/CDS/intron, plot-ready
    BED12,     // domain_blocks.bed12          — one BED12 per domain, IGV-ready
    ALL        // everything plus summary + unmapped + run_metadata
};

// Tri-state for transcript-level flags that may be unknown if the GTF doesn't
// carry tag attributes at all (e.g. plain Ensembl base GTF). We emit NA in
// that case rather than silently writing `false`.
enum class TriBool : uint8_t {
    NA   = 0,
    TRUE_ = 1,
    FALSE_ = 2
};
inline const char* tribool_to_string(TriBool b) {
    switch (b) {
        case TriBool::TRUE_:  return "true";
        case TriBool::FALSE_: return "false";
        case TriBool::NA:
        default:              return "NA";
    }
}

// Plot-ready feature kind for isoform_structure.tsv.
// This is *derived* from GTF EXON/CDS intersections, not raw GTF features.
enum class PlotFeatureType {
    FIVE_PRIME_UTR,
    CDS,
    THREE_PRIME_UTR,
    INTRON
};

// Domain-overlap classification for a feature in isoform_structure.tsv.
enum class DomainOverlapKind {
    NO,                          // feature is outside the domain entirely
    CODING_OVERLAP,              // CDS feature that codes part of the domain
    INSIDE_DOMAIN_GENOMIC_SPAN   // intron between two domain-coding CDS segments
};

constexpr size_t INITIAL_PROTEIN_CAPACITY = 30000;
constexpr size_t INITIAL_INTERVAL_CAPACITY = 200000;
constexpr size_t READ_BUFFER_SIZE = 1024 * 1024;
constexpr size_t MAX_LINE_LENGTH = 8192;

// v3: added transcript->protein reverse map, MANE_Select / Ensembl_canonical
// tags per transcript, gtf_has_tags flag, cds_total_nt per protein. Old (v1/v2)
// indices will be rejected with an explicit "rebuild" message.
// Magic is "P2E8" (project renamed prot2exon, May 2026). Indices built under
// the old "P2G8" magic must be rebuilt.
constexpr uint32_t INDEX_FORMAT_MAGIC = 0x50324538; // "P2E8"
constexpr uint32_t INDEX_FORMAT_VERSION = 3;

enum class ErrorCode {
    SUCCESS = 0,
    FILE_NOT_FOUND,
    PARSE_ERROR,
    MEMORY_ERROR,
    INVALID_FORMAT,
    INDEX_VERSION_MISMATCH
};

constexpr inline FeatureType operator|(FeatureType a, FeatureType b) {
    return static_cast<FeatureType>(static_cast<int>(a) | static_cast<int>(b));
}

constexpr inline FeatureType operator&(FeatureType a, FeatureType b) {
    return static_cast<FeatureType>(static_cast<int>(a) & static_cast<int>(b));
}

constexpr inline bool has_feature(FeatureType mask, FeatureType feature) {
    return (mask & feature) == feature;
}

std::string feature_type_to_string(FeatureType type);
FeatureType string_to_feature_type(const std::string& str);
FeatureType parse_feature_types(const std::string& types_str);
std::string get_feature_type_name(FeatureType type);

std::string plot_feature_to_string(PlotFeatureType t);
std::string overlap_kind_to_string(DomainOverlapKind k);
std::string plot_group_for(PlotFeatureType t, DomainOverlapKind k);

#endif // COMMON_HPP

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

// Performance optimized types
using ProteinMap = std::unordered_map<std::string, std::vector<class GenomicInterval>>;
using TranscriptMap = std::unordered_map<std::string, std::string>; // protein_id -> transcript_id
using GeneMap = std::unordered_map<std::string, class GeneStructure>; // gene_id -> structure

// Forward declarations
struct GenomicInterval;
struct ProteinDomain;
struct GTFRecord;
struct GeneStructure;

// Feature types (simplified, removed UTR references)
enum class FeatureType {
    EXON = 1,
    CDS = 2,
    INTRON = 4,
    START_CODON = 8,
    STOP_CODON = 16,
    ALL = 31  // All bits set
};

// Feature priority (higher number = higher priority for overlap resolution)
enum class FeaturePriority {
    INTRON = 1,
    EXON = 2,
    CDS = 3,
    START_CODON = 4,
    STOP_CODON = 5
};

// Constants
constexpr size_t INITIAL_PROTEIN_CAPACITY = 30000;
constexpr size_t INITIAL_INTERVAL_CAPACITY = 200000;
constexpr size_t READ_BUFFER_SIZE = 1024 * 1024;
constexpr size_t MAX_LINE_LENGTH = 8192;

// Error codes
enum class ErrorCode {
    SUCCESS = 0,
    FILE_NOT_FOUND,
    PARSE_ERROR,
    MEMORY_ERROR,
    INVALID_FORMAT
};

// Utility functions for feature types
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

#endif // COMMON_HPP

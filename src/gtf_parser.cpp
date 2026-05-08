#include "gtf_parser.hpp"
#include "utils.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>

GTFParser::GTFParser() : lines_processed_(0) {
    protein_index_.reserve(INITIAL_PROTEIN_CAPACITY);
    protein_to_transcript_.reserve(INITIAL_PROTEIN_CAPACITY);
    protein_to_gene_name_.reserve(INITIAL_PROTEIN_CAPACITY);
    protein_to_gene_id_.reserve(INITIAL_PROTEIN_CAPACITY);
    temp_fields_.reserve(20);
}

ErrorCode GTFParser::load_gtf(const std::string& filename) {
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
        if (line.empty() || line[0] == '#') continue;
        if (parse_line(line, record)) {
            if (is_relevant_feature(record.feature)) {
                process_record(record);
            }
        }
        if (lines_processed_ % 100000 == 0) {
            std::cerr << "Processed " << lines_processed_ << " lines...\r" << std::flush;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cerr << "\nCompleted parsing GTF file in " << duration.count() << "ms" << std::endl;
    std::cerr << "Total lines processed: " << lines_processed_ << std::endl;

    build_protein_index_from_structures();
    std::cerr << "Proteins indexed: " << protein_index_.size() << std::endl;
    std::cerr << "Estimated memory usage: " << estimate_memory_usage() / (1024*1024) << " MB" << std::endl;

    // Sort intervals per protein. Negative-strand proteins keep genomic order
    // here; the mapper handles translation order explicitly via the strand field.
    for (auto& [protein_id, intervals] : protein_index_) {
        std::sort(intervals.begin(), intervals.end());
    }

    return ErrorCode::SUCCESS;
}

bool GTFParser::parse_line(const std::string& line, GTFRecord& record) const {
    temp_fields_.clear();
    temp_fields_ = split_fast(line, '\t');
    if (temp_fields_.size() < 9) return false;
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
    } catch (const std::exception&) {
        return false;
    }
}

std::vector<std::string> GTFParser::split_fast(const std::string& str, char delimiter) const {
    std::vector<std::string> result;
    result.reserve(12);
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
    std::string search_key = key + " \"";
    size_t start = attributes.find(search_key);
    if (start == std::string::npos) return "";
    start += search_key.length();
    size_t end = attributes.find("\"", start);
    if (end == std::string::npos) return "";
    return attributes.substr(start, end - start);
}

void GTFParser::process_record(const GTFRecord& record) {
    FeatureType feature_type = string_to_feature_type(record.feature);

    std::string gene_id = extract_attribute(record.attributes, "gene_id");
    std::string transcript_id = extract_attribute(record.attributes, "transcript_id");
    std::string protein_id = extract_attribute(record.attributes, "protein_id");
    std::string gene_name = extract_attribute(record.attributes, "gene_name");
    std::string exon_number_str = extract_attribute(record.attributes, "exon_number");

    if (gene_id.empty() || transcript_id.empty()) return;

    auto remove_version = [](const std::string& id) -> std::string {
        size_t dot_pos = id.find('.');
        return (dot_pos != std::string::npos) ? id.substr(0, dot_pos) : id;
    };

    gene_id = remove_version(gene_id);
    transcript_id = remove_version(transcript_id);
    if (!protein_id.empty()) protein_id = remove_version(protein_id);

    uint32_t exon_number = 0;
    if (!exon_number_str.empty()) {
        try { exon_number = static_cast<uint32_t>(std::stoul(exon_number_str)); }
        catch (...) { exon_number = 0; }
    }

    GenomicInterval interval(record.chromosome, record.start, record.end,
                             record.strand, record.phase, feature_type,
                             transcript_id, gene_id, gene_name, exon_number);

    temp_transcript_data_[transcript_id].push_back(interval);

    if (!protein_id.empty()) {
        protein_to_transcript_[protein_id] = transcript_id;
        if (!gene_name.empty()) protein_to_gene_name_[protein_id] = gene_name;
        protein_to_gene_id_[protein_id] = gene_id;
    }
}

bool GTFParser::is_relevant_feature(const std::string& feature) const {
    return feature == "exon" || feature == "CDS";
}

const std::vector<GenomicInterval>* GTFParser::get_protein_intervals(const std::string& protein_id) const {
    std::string normalized = protein_id;
    size_t dot_pos = normalized.find('.');
    if (dot_pos != std::string::npos) normalized = normalized.substr(0, dot_pos);
    auto it = protein_index_.find(normalized);
    return (it != protein_index_.end()) ? &it->second : nullptr;
}

const std::string* GTFParser::get_transcript_id(const std::string& protein_id) const {
    std::string normalized = protein_id;
    size_t dot_pos = normalized.find('.');
    if (dot_pos != std::string::npos) normalized = normalized.substr(0, dot_pos);
    auto it = protein_to_transcript_.find(normalized);
    return (it != protein_to_transcript_.end()) ? &it->second : nullptr;
}

std::string GTFParser::get_gene_name(const std::string& protein_id) const {
    std::string normalized = protein_id;
    size_t dot_pos = normalized.find('.');
    if (dot_pos != std::string::npos) normalized = normalized.substr(0, dot_pos);
    auto it = protein_to_gene_name_.find(normalized);
    return (it != protein_to_gene_name_.end()) ? it->second : std::string();
}

std::string GTFParser::get_gene_id(const std::string& protein_id) const {
    std::string normalized = protein_id;
    size_t dot_pos = normalized.find('.');
    if (dot_pos != std::string::npos) normalized = normalized.substr(0, dot_pos);
    auto it = protein_to_gene_id_.find(normalized);
    return (it != protein_to_gene_id_.end()) ? it->second : std::string();
}

size_t GTFParser::estimate_memory_usage() const {
    size_t total = 0;
    for (const auto& [protein_id, intervals] : protein_index_) {
        total += protein_id.size();
        total += intervals.size() * sizeof(GenomicInterval);
        for (const auto& interval : intervals) {
            total += interval.chromosome.size() + interval.transcript_id.size()
                   + interval.gene_id.size() + interval.gene_name.size();
        }
    }
    return total;
}

void GTFParser::build_protein_index_from_structures() {
    // Build per-protein interval lists: keep both EXON and CDS from the
    // transcript that this protein belongs to.
    for (const auto& [protein_id, transcript_id] : protein_to_transcript_) {
        auto it = temp_transcript_data_.find(transcript_id);
        if (it == temp_transcript_data_.end()) continue;

        auto& protein_intervals = protein_index_[protein_id];
        protein_intervals.clear();
        protein_intervals.reserve(it->second.size());
        for (const auto& iv : it->second) {
            if (iv.feature_type == FeatureType::EXON || iv.feature_type == FeatureType::CDS) {
                protein_intervals.push_back(iv);
            }
        }
    }
    temp_transcript_data_.clear();
}

ErrorCode GTFParser::save_index(const std::string& filename) const {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) return ErrorCode::FILE_NOT_FOUND;

    uint32_t magic = INDEX_FORMAT_MAGIC;
    uint32_t version = INDEX_FORMAT_VERSION;
    file.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    file.write(reinterpret_cast<const char*>(&version), sizeof(version));

    auto write_str = [&](const std::string& s) {
        uint16_t len = static_cast<uint16_t>(s.size());
        file.write(reinterpret_cast<const char*>(&len), sizeof(len));
        if (len) file.write(s.data(), len);
    };

    uint32_t num_proteins = protein_index_.size();
    file.write(reinterpret_cast<const char*>(&num_proteins), sizeof(num_proteins));

    for (const auto& [protein_id, intervals] : protein_index_) {
        write_str(protein_id);

        auto tx_it = protein_to_transcript_.find(protein_id);
        write_str(tx_it != protein_to_transcript_.end() ? tx_it->second : "");
        auto gn_it = protein_to_gene_name_.find(protein_id);
        write_str(gn_it != protein_to_gene_name_.end() ? gn_it->second : "");
        auto gi_it = protein_to_gene_id_.find(protein_id);
        write_str(gi_it != protein_to_gene_id_.end() ? gi_it->second : "");

        uint32_t num_intervals = intervals.size();
        file.write(reinterpret_cast<const char*>(&num_intervals), sizeof(num_intervals));
        for (const auto& iv : intervals) {
            write_str(iv.chromosome);
            file.write(reinterpret_cast<const char*>(&iv.start), sizeof(iv.start));
            file.write(reinterpret_cast<const char*>(&iv.end), sizeof(iv.end));
            file.write(reinterpret_cast<const char*>(&iv.strand), sizeof(iv.strand));
            file.write(reinterpret_cast<const char*>(&iv.phase), sizeof(iv.phase));
            file.write(reinterpret_cast<const char*>(&iv.feature_type), sizeof(iv.feature_type));
            file.write(reinterpret_cast<const char*>(&iv.exon_number), sizeof(iv.exon_number));
        }
    }
    return ErrorCode::SUCCESS;
}

ErrorCode GTFParser::load_index(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) return ErrorCode::FILE_NOT_FOUND;

    clear();

    uint32_t magic = 0, version = 0;
    file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (magic != INDEX_FORMAT_MAGIC || version != INDEX_FORMAT_VERSION) {
        std::cerr << "Error: index file '" << filename
                  << "' is not in the expected format (magic=" << std::hex << magic
                  << " version=" << std::dec << version
                  << "). Rebuild with --build-index." << std::endl;
        return ErrorCode::INDEX_VERSION_MISMATCH;
    }

    auto read_str = [&]() {
        uint16_t len = 0;
        file.read(reinterpret_cast<char*>(&len), sizeof(len));
        std::string s(len, '\0');
        if (len) file.read(&s[0], len);
        return s;
    };

    uint32_t num_proteins = 0;
    file.read(reinterpret_cast<char*>(&num_proteins), sizeof(num_proteins));
    protein_index_.reserve(num_proteins);

    for (uint32_t i = 0; i < num_proteins; ++i) {
        std::string protein_id = read_str();
        std::string transcript_id = read_str();
        std::string gene_name = read_str();
        std::string gene_id = read_str();

        protein_to_transcript_[protein_id] = transcript_id;
        if (!gene_name.empty()) protein_to_gene_name_[protein_id] = gene_name;
        if (!gene_id.empty()) protein_to_gene_id_[protein_id] = gene_id;

        uint32_t num_intervals = 0;
        file.read(reinterpret_cast<char*>(&num_intervals), sizeof(num_intervals));
        std::vector<GenomicInterval> intervals;
        intervals.reserve(num_intervals);
        for (uint32_t j = 0; j < num_intervals; ++j) {
            GenomicInterval iv;
            uint16_t chr_len = 0;
            file.read(reinterpret_cast<char*>(&chr_len), sizeof(chr_len));
            iv.chromosome.resize(chr_len);
            if (chr_len) file.read(&iv.chromosome[0], chr_len);
            file.read(reinterpret_cast<char*>(&iv.start), sizeof(iv.start));
            file.read(reinterpret_cast<char*>(&iv.end), sizeof(iv.end));
            file.read(reinterpret_cast<char*>(&iv.strand), sizeof(iv.strand));
            file.read(reinterpret_cast<char*>(&iv.phase), sizeof(iv.phase));
            file.read(reinterpret_cast<char*>(&iv.feature_type), sizeof(iv.feature_type));
            file.read(reinterpret_cast<char*>(&iv.exon_number), sizeof(iv.exon_number));
            iv.transcript_id = transcript_id;
            iv.gene_id = gene_id;
            iv.gene_name = gene_name;
            iv.set_priority_from_type();
            intervals.emplace_back(std::move(iv));
        }
        protein_index_[protein_id] = std::move(intervals);
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
    return FeatureType::EXON;
}

FeatureType parse_feature_types(const std::string& types_str) {
    if (types_str == "all") return FeatureType::ALL;
    FeatureType result = static_cast<FeatureType>(0);
    std::istringstream iss(types_str);
    std::string type;
    while (std::getline(iss, type, ',')) {
        type.erase(0, type.find_first_not_of(" \t"));
        type.erase(type.find_last_not_of(" \t") + 1);
        if (type == "exon") result = result | FeatureType::EXON;
        else if (type == "cds" || type == "CDS") result = result | FeatureType::CDS;
        else if (type == "intron") result = result | FeatureType::INTRON;
        else if (type == "start_codon") result = result | FeatureType::START_CODON;
        else if (type == "stop_codon") result = result | FeatureType::STOP_CODON;
    }
    if (static_cast<int>(result) == 0) result = FeatureType::CDS;
    return result;
}

std::string get_feature_type_name(FeatureType type) {
    return feature_type_to_string(type);
}

std::string plot_feature_to_string(PlotFeatureType t) {
    switch (t) {
        case PlotFeatureType::FIVE_PRIME_UTR:  return "five_prime_UTR";
        case PlotFeatureType::CDS:             return "CDS";
        case PlotFeatureType::THREE_PRIME_UTR: return "three_prime_UTR";
        case PlotFeatureType::INTRON:          return "intron";
    }
    return "unknown";
}

std::string overlap_kind_to_string(DomainOverlapKind k) {
    switch (k) {
        case DomainOverlapKind::NO:                         return "no";
        case DomainOverlapKind::CODING_OVERLAP:             return "coding_overlap";
        case DomainOverlapKind::INSIDE_DOMAIN_GENOMIC_SPAN: return "inside_domain_genomic_span";
    }
    return "no";
}

std::string plot_group_for(PlotFeatureType t, DomainOverlapKind k) {
    switch (t) {
        case PlotFeatureType::FIVE_PRIME_UTR:  return "five_prime_UTR";
        case PlotFeatureType::THREE_PRIME_UTR: return "three_prime_UTR";
        case PlotFeatureType::CDS:
            return (k == DomainOverlapKind::CODING_OVERLAP) ? "CDS_domain" : "CDS_no_domain";
        case PlotFeatureType::INTRON:
            return (k == DomainOverlapKind::INSIDE_DOMAIN_GENOMIC_SPAN)
                ? "intron_domain_span" : "intron";
    }
    return "unknown";
}

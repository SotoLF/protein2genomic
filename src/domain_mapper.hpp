#ifndef DOMAIN_MAPPER_HPP
#define DOMAIN_MAPPER_HPP

#include "common.hpp"
#include "gtf_parser.hpp"

struct ProteinDomain {
    std::string protein_id;
    uint32_t start;
    uint32_t end;
    std::string domain_id;
    
    ProteinDomain() : start(0), end(0) {}
    
    ProteinDomain(const std::string& pid, uint32_t s, uint32_t e, const std::string& did = "")
        : protein_id(pid), start(s), end(e), domain_id(did) {}
};

struct GenomicDomain {
    std::string chromosome;
    uint32_t start;
    uint32_t end;
    std::string protein_id;
    std::string domain_id;
    std::string feature_type;
    char strand;
    uint32_t protein_start;
    uint32_t protein_end;
    std::string cds_id;           // For detailed format: CDS_1, CDS_2, etc.
    std::string domain_overlap;   // For detailed format: Yes/No
    uint32_t region_index;        // Index of the region (1, 2, 3, ...)
    std::string feature_id;       // Feature identifier: CDS_1, intron_1, etc.
    
    GenomicDomain() : start(0), end(0), strand('+'), protein_start(0), protein_end(0), domain_overlap("No"), region_index(0) {}
};

class DomainMapper {
private:
    const GTFParser& gtf_parser_;
    std::vector<ProteinDomain> domains_;
    std::vector<GenomicDomain> results_;
    std::string mode_;           // "basic" or "full"
    std::string format_;         // "simple" or "detailed"
    
    // Coordinate conversion methods
    std::vector<GenomicDomain> map_domain_to_genomic(const ProteinDomain& domain) const;
    std::vector<GenomicDomain> map_domain_detailed(const ProteinDomain& domain) const;
    uint32_t protein_to_genomic_position(uint32_t protein_pos, 
                                        const std::vector<GenomicInterval>& intervals,
                                        char strand) const;
    
    // Optimization: batch processing
    void process_domains_batch(const std::vector<ProteinDomain>& batch);
    
public:
    explicit DomainMapper(const GTFParser& parser, const std::string& mode = "basic", const std::string& format = "simple");
    ~DomainMapper() = default;
    
    // Main interface
    ErrorCode load_domains(const std::string& bed_filename);
    ErrorCode process_domains();
    ErrorCode write_results(const std::string& output_filename, const std::string& format = "bed") const;
    
    // Statistics
    size_t get_domain_count() const { return domains_.size(); }
    size_t get_mapped_count() const;
    size_t get_unmapped_count() const;
    
    void clear() { domains_.clear(); results_.clear(); }
};

#endif // DOMAIN_MAPPER_HPP

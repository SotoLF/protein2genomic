#ifndef DOMAIN_MAPPER_HPP
#define DOMAIN_MAPPER_HPP

#include "common.hpp"
#include "gtf_parser.hpp"

struct ProteinDomain {
    std::string input_id;        // user-supplied row identifier (BED column 4 if present, else auto)
    std::string protein_id;
    uint32_t start;              // 1-based aa
    uint32_t end;                // 1-based aa, inclusive
    std::string domain_id;       // semantic name (BED column 5 if present)

    ProteinDomain() : start(0), end(0) {}
};

// One row per CDS segment that *codes* part of the domain.
// Backs domain_cds_segments.bed and .tsv.
struct CodingSegmentRow {
    std::string chrom;
    uint32_t genomic_start;      // 1-based, inclusive
    uint32_t genomic_end;        // 1-based, inclusive
    char strand;
    std::string input_id;
    std::string protein_id;
    std::string transcript_id;
    std::string gene_id;
    std::string gene_name;
    std::string domain_id;
    uint32_t aa_start;           // domain aa range
    uint32_t aa_end;
    uint32_t segment_index_in_domain; // 1..N across the domain (in translation order)
    uint32_t cds_nt_start;       // domain-relative nt offsets, 1-based
    uint32_t cds_nt_end;
    uint32_t aa_start_encoded;
    uint32_t aa_end_encoded;
};

// One row per input domain capturing the full genomic envelope (incl. introns).
// Backs domain_span_with_introns.bed.
struct SpanRow {
    std::string chrom;
    uint32_t genomic_start;
    uint32_t genomic_end;
    char strand;
    std::string input_id;
    std::string protein_id;
    std::string domain_id;
    std::string gene_name;
};

// One row per structural feature of the transcript (UTR/CDS/intron),
// split so that domain-overlapping vs non-overlapping CDS are separate rows.
// Backs isoform_structure.tsv.
struct IsoformSegmentRow {
    std::string input_id;
    std::string gene_id;
    std::string gene_name;
    std::string transcript_id;
    std::string protein_id;
    std::string domain_id;

    std::string chrom;
    char strand;

    PlotFeatureType feature_type;
    std::string feature_id;       // five_prime_UTR_1, CDS_1, intron_1, ...
    uint32_t exon_number;         // 0 for introns
    uint32_t feature_genomic_start;
    uint32_t feature_genomic_end;
    uint32_t feature_length_nt;
    uint32_t feature_order_genomic;
    uint32_t feature_order_transcript;

    // CDS-only nucleotide/aa ranges; UTR/intron rows leave these as 0 (printed as NA).
    bool has_cds_coords;
    uint32_t cds_nt_start;
    uint32_t cds_nt_end;
    uint32_t aa_start_encoded;
    uint32_t aa_end_encoded;

    DomainOverlapKind overlap;
    bool has_overlap_coords;
    uint32_t domain_overlap_genomic_start;
    uint32_t domain_overlap_genomic_end;
    uint32_t domain_overlap_cds_nt_start;
    uint32_t domain_overlap_cds_nt_end;
    uint32_t domain_overlap_aa_start;
    uint32_t domain_overlap_aa_end;
    double domain_overlap_fraction_of_feature;
    double domain_overlap_fraction_of_domain;

    std::string plot_group;

    IsoformSegmentRow()
        : strand('+'), feature_type(PlotFeatureType::CDS), exon_number(0),
          feature_genomic_start(0), feature_genomic_end(0), feature_length_nt(0),
          feature_order_genomic(0), feature_order_transcript(0),
          has_cds_coords(false), cds_nt_start(0), cds_nt_end(0),
          aa_start_encoded(0), aa_end_encoded(0),
          overlap(DomainOverlapKind::NO), has_overlap_coords(false),
          domain_overlap_genomic_start(0), domain_overlap_genomic_end(0),
          domain_overlap_cds_nt_start(0), domain_overlap_cds_nt_end(0),
          domain_overlap_aa_start(0), domain_overlap_aa_end(0),
          domain_overlap_fraction_of_feature(0.0),
          domain_overlap_fraction_of_domain(0.0) {}
};

// One row per input domain — successful or partial.
struct MappingSummaryRow {
    std::string input_id;
    std::string protein_id;
    std::string transcript_id;
    std::string gene_id;
    std::string gene_name;
    std::string domain_id;
    std::string chrom;
    char strand;
    uint32_t aa_start;
    uint32_t aa_end;
    uint32_t domain_length_aa;
    uint32_t domain_length_nt;
    uint32_t protein_length_aa;
    uint32_t domain_genomic_start;
    uint32_t domain_genomic_end;
    uint32_t n_coding_segments;
    bool fully_mapped;
};

struct UnmappedDomainRow {
    std::string input_id;
    std::string protein_id;
    uint32_t aa_start;
    uint32_t aa_end;
    std::string domain_id;
    std::string reason;
};

// Per-domain result bundle, computed once and consumed by the writers.
struct DomainResult {
    ProteinDomain domain;
    bool mapped;
    UnmappedDomainRow unmapped;

    MappingSummaryRow summary;
    std::vector<CodingSegmentRow> coding_segments;
    SpanRow span;
    std::vector<IsoformSegmentRow> isoform_segments;

    DomainResult() : mapped(false) {}
};

class DomainMapper {
public:
    DomainMapper(const GTFParser& parser, OutputKind output_kind);

    ErrorCode load_domains(const std::string& bed_filename);
    ErrorCode process_domains();

    const std::vector<DomainResult>& results() const { return results_; }
    size_t domain_count() const { return domains_.size(); }
    size_t mapped_count() const;
    size_t unmapped_count() const { return domain_count() - mapped_count(); }

private:
    const GTFParser& gtf_;
    OutputKind output_kind_;
    std::vector<ProteinDomain> domains_;
    std::vector<DomainResult> results_;

    DomainResult process_one(const ProteinDomain& d) const;
};

#endif // DOMAIN_MAPPER_HPP

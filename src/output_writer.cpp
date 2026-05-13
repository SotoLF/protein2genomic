#include "output_writer.hpp"
#include "utils.hpp"
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <ctime>
#include <iomanip>

namespace {

bool ensure_dir(const std::string& dir) {
    struct stat st;
    if (stat(dir.c_str(), &st) == 0) return S_ISDIR(st.st_mode);
    return mkdir(dir.c_str(), 0755) == 0;
}

std::string join(const std::string& dir, const std::string& name) {
    if (dir.empty()) return name;
    if (dir.back() == '/') return dir + name;
    return dir + "/" + name;
}

void w_int_or_na(std::ofstream& f, uint32_t v) {
    if (v == 0) f << "NA"; else f << v;
}

void w_str_or_na(std::ofstream& f, const std::string& s) {
    if (s.empty()) f << "NA"; else f << s;
}

void w_double(std::ofstream& f, double v) {
    f << std::fixed << std::setprecision(4) << v;
}

std::string iso8601_now() {
    std::time_t t = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return buf;
}

std::string output_kind_to_string(OutputKind k) {
    switch (k) {
        case OutputKind::CODING:  return "coding";
        case OutputKind::INTRONS: return "introns";
        case OutputKind::SPAN:    return "span";
        case OutputKind::ISOFORM: return "isoform";
        case OutputKind::BED12:   return "bed12";
        case OutputKind::ALL:     return "all";
    }
    return "unknown";
}

std::string bed_name(const DomainResult& r) {
    std::ostringstream oss;
    // Prefer the resolved protein id (set in summary) so ENST queries still
    // get a meaningful name; fall back to the raw input.
    const std::string& pid = r.summary.protein_id.empty()
        ? r.domain.protein_id : r.summary.protein_id;
    oss << pid;
    if (!r.domain.domain_id.empty()) oss << "_" << r.domain.domain_id;
    if (r.domain.has_domain()) oss << "_" << r.domain.start << "-" << r.domain.end;
    else oss << "_no_domain";
    return oss.str();
}

ErrorCode write_summary(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    // Existing columns are unchanged in order; new columns are appended at the
    // end so consumers that read by index still work.
    f << "input_id\tprotein_id\ttranscript_id\tgene_id\tgene_name\tdomain_id\tchrom\tstrand\t"
         "aa_start\taa_end\tdomain_length_aa\tdomain_length_nt\tprotein_length_aa\t"
         "domain_genomic_start\tdomain_genomic_end\tn_coding_segments\tfully_mapped\tno_domain_mode\t"
         "input_id_type\tis_mane_select\tis_ensembl_canonical\t"
         "cds_length_mismatch\tcds_nt_remainder\tstatus\n";
    for (const auto& r : results) {
        if (!r.mapped) {
            f << r.domain.input_id << "\t";
            w_str_or_na(f, r.domain.protein_id);
            f << "\tNA\tNA\tNA\t";
            w_str_or_na(f, r.domain.domain_id);
            f << "\tNA\tNA\t";
            if (r.domain.has_domain()) {
                f << r.domain.start << "\t" << r.domain.end
                  << "\t" << (r.domain.end - r.domain.start + 1)
                  << "\t" << ((r.domain.end - r.domain.start + 1) * 3);
            } else {
                f << "NA\tNA\tNA\tNA";
            }
            f << "\tNA\tNA\tNA\t0\tfalse\t"
              << (r.no_domain_mode ? "true" : "false")
              << "\tNA\tNA\tNA\tNA\tNA"
              << "\t" << r.unmapped.reason << "\n";
            continue;
        }
        const auto& s = r.summary;
        f << s.input_id << "\t";
        w_str_or_na(f, s.protein_id);
        f << "\t" << s.transcript_id << "\t";
        w_str_or_na(f, s.gene_id);
        f << "\t"; w_str_or_na(f, s.gene_name);
        f << "\t"; w_str_or_na(f, s.domain_id);
        f << "\t" << s.chrom << "\t" << s.strand << "\t";
        if (s.no_domain_mode) {
            f << "NA\tNA\tNA\tNA";
        } else {
            f << s.aa_start << "\t" << s.aa_end
              << "\t" << s.domain_length_aa << "\t" << s.domain_length_nt;
        }
        f << "\t" << s.protein_length_aa
          << "\t";
        if (s.no_domain_mode) {
            f << "NA\tNA\tNA";
        } else {
            f << s.domain_genomic_start << "\t" << s.domain_genomic_end
              << "\t" << s.n_coding_segments;
        }
        f << "\t" << (s.fully_mapped ? "true" : "false")
          << "\t" << (s.no_domain_mode ? "true" : "false")
          << "\t"; w_str_or_na(f, s.input_id_type);
        f << "\t" << tribool_to_string(s.is_mane_select)
          << "\t" << tribool_to_string(s.is_ensembl_canonical)
          << "\t" << (s.cds_length_mismatch ? "true" : "false")
          << "\t" << static_cast<int>(s.cds_nt_remainder)
          << "\t";
        // Status string. ok / partial get a _cds_mismatch suffix when relevant.
        const char* base = s.no_domain_mode ? "structure_only"
                                             : (s.fully_mapped ? "ok" : "partial");
        if (s.cds_length_mismatch && !s.no_domain_mode) {
            f << base << "_cds_mismatch";
        } else {
            f << base;
        }
        f << "\n";
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

// Header + columns shared by isoform_structure.tsv, domain_cds_segments.tsv, domain_introns.tsv.
const char* kFeatureTsvHeader =
"input_id\tgene_id\tgene_name\ttranscript_id\tprotein_id\tdomain_id\t"
"is_mane_select\tis_ensembl_canonical\tcds_length_mismatch\tcds_nt_remainder\t"
"chrom\tstrand\tfeature_type\tfeature_id\tfeature_part\texon_number\t"
"feature_genomic_start\tfeature_genomic_end\tfeature_length_nt\t"
"feature_order_genomic\tfeature_order_transcript\t"
"cds_nt_start\tcds_nt_end\taa_start_encoded\taa_end_encoded\t"
"overlaps_domain\t"
"domain_overlap_genomic_start\tdomain_overlap_genomic_end\t"
"domain_overlap_cds_nt_start\tdomain_overlap_cds_nt_end\t"
"domain_overlap_aa_start\tdomain_overlap_aa_end\t"
"domain_overlap_fraction_of_feature\tdomain_overlap_fraction_of_domain\t"
"plot_group\n";

void write_feature_tsv_row(std::ofstream& f, const IsoformSegmentRow& s) {
    f << s.input_id << "\t";
    w_str_or_na(f, s.gene_id); f << "\t";
    w_str_or_na(f, s.gene_name); f << "\t";
    w_str_or_na(f, s.transcript_id); f << "\t";
    w_str_or_na(f, s.protein_id); f << "\t";
    w_str_or_na(f, s.domain_id); f << "\t"
      << tribool_to_string(s.is_mane_select) << "\t"
      << tribool_to_string(s.is_ensembl_canonical) << "\t"
      << (s.cds_length_mismatch ? "true" : "false") << "\t"
      << static_cast<int>(s.cds_nt_remainder) << "\t";
    f << s.chrom << "\t" << s.strand << "\t"
      << plot_feature_to_string(s.feature_type) << "\t"
      << s.feature_id << "\t"
      << s.feature_part << "\t";
    w_int_or_na(f, s.exon_number); f << "\t"
      << s.feature_genomic_start << "\t" << s.feature_genomic_end << "\t"
      << s.feature_length_nt << "\t"
      << s.feature_order_genomic << "\t" << s.feature_order_transcript << "\t";
    if (s.has_cds_coords) {
        f << s.cds_nt_start << "\t" << s.cds_nt_end << "\t"
          << s.aa_start_encoded << "\t" << s.aa_end_encoded;
    } else {
        f << "NA\tNA\tNA\tNA";
    }
    f << "\t";
    if (s.no_domain_mode) f << "NA";
    else f << overlap_kind_to_string(s.overlap);
    f << "\t";
    if (!s.no_domain_mode && s.has_overlap_coords) {
        f << s.domain_overlap_genomic_start << "\t" << s.domain_overlap_genomic_end << "\t"
          << s.domain_overlap_cds_nt_start << "\t" << s.domain_overlap_cds_nt_end << "\t"
          << s.domain_overlap_aa_start << "\t" << s.domain_overlap_aa_end << "\t";
        w_double(f, s.domain_overlap_fraction_of_feature); f << "\t";
        w_double(f, s.domain_overlap_fraction_of_domain);
    } else {
        f << "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
    }
    f << "\t" << s.plot_group << "\n";
}

ErrorCode write_isoform_tsv(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    f << kFeatureTsvHeader;
    for (const auto& r : results) {
        if (!r.mapped) continue;
        for (const auto& s : r.isoform_segments) write_feature_tsv_row(f, s);
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_coding_tsv(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    f << kFeatureTsvHeader;
    for (const auto& r : results) {
        if (!r.mapped) continue;
        for (const auto& s : r.isoform_segments) {
            if (s.feature_type == PlotFeatureType::CDS) write_feature_tsv_row(f, s);
        }
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_introns_tsv(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    f << kFeatureTsvHeader;
    for (const auto& r : results) {
        if (!r.mapped) continue;
        for (const auto& s : r.isoform_segments) {
            if (s.feature_type == PlotFeatureType::INTRON) write_feature_tsv_row(f, s);
        }
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

// BED of CDS rows that code part of the domain. Empty in no-domain mode.
ErrorCode write_coding_bed(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    for (const auto& r : results) {
        if (!r.mapped || r.no_domain_mode) continue;
        std::string name = bed_name(r);
        for (const auto& s : r.isoform_segments) {
            if (s.feature_type != PlotFeatureType::CDS) continue;
            if (s.overlap != DomainOverlapKind::CODING_OVERLAP) continue;
            f << s.chrom << "\t" << (s.feature_genomic_start - 1) << "\t" << s.feature_genomic_end
              << "\t" << name << "\t0\t" << s.strand << "\n";
        }
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

// BED of introns located inside the domain genomic span. Empty in no-domain mode.
ErrorCode write_introns_bed(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    for (const auto& r : results) {
        if (!r.mapped || r.no_domain_mode) continue;
        std::string name = bed_name(r);
        for (const auto& s : r.isoform_segments) {
            if (s.feature_type != PlotFeatureType::INTRON) continue;
            if (s.overlap != DomainOverlapKind::INSIDE_DOMAIN_GENOMIC_SPAN) continue;
            f << s.chrom << "\t" << (s.feature_genomic_start - 1) << "\t" << s.feature_genomic_end
              << "\t" << name << "\t0\t" << s.strand << "\n";
        }
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

// BED12 with one row per domain. chromStart..chromEnd is the genomic envelope
// of the domain (introns included). The blocks are the CDS slices that code
// the domain (overlap_kind == CODING_OVERLAP). thickStart=chromStart,
// thickEnd=chromEnd so IGV draws the whole feature thick. Empty in no-domain
// mode. Block lists are ordered by genomic start as BED12 requires.
ErrorCode write_bed12(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    for (const auto& r : results) {
        if (!r.mapped || r.no_domain_mode || !r.has_span) continue;
        // Collect coding-overlap CDS slices in genomic order.
        std::vector<std::pair<uint32_t, uint32_t>> blocks; // 1-based [s,e]
        for (const auto& s : r.isoform_segments) {
            if (s.feature_type != PlotFeatureType::CDS) continue;
            if (s.overlap != DomainOverlapKind::CODING_OVERLAP) continue;
            // Only the domain-overlapping sub-slice of this CDS row.
            if (!s.has_overlap_coords) continue;
            blocks.emplace_back(s.domain_overlap_genomic_start,
                                s.domain_overlap_genomic_end);
        }
        if (blocks.empty()) continue;
        std::sort(blocks.begin(), blocks.end());

        uint32_t chrom_start = r.span.genomic_start - 1; // 0-based
        uint32_t chrom_end   = r.span.genomic_end;       // half-open
        std::string name = bed_name(r);
        f << r.span.chrom << "\t" << chrom_start << "\t" << chrom_end << "\t"
          << name << "\t0\t" << r.span.strand << "\t"
          << chrom_start << "\t" << chrom_end << "\t"
          << "255,0,0\t" << blocks.size() << "\t";
        // blockSizes
        for (size_t i = 0; i < blocks.size(); ++i) {
            if (i) f << ",";
            f << (blocks[i].second - blocks[i].first + 1);
        }
        f << ",\t";
        // blockStarts (offset from chromStart)
        for (size_t i = 0; i < blocks.size(); ++i) {
            if (i) f << ",";
            f << (blocks[i].first - 1 - chrom_start);
        }
        f << ",\n";
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_span_bed(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    for (const auto& r : results) {
        if (!r.mapped || !r.has_span) continue;
        std::string name = bed_name(r);
        const auto& s = r.span;
        f << s.chrom << "\t" << (s.genomic_start - 1) << "\t" << s.genomic_end
          << "\t" << name << "\t0\t" << s.strand << "\n";
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_unmapped(const std::string& path, const std::vector<DomainResult>& results) {
    bool any = false;
    for (const auto& r : results) if (!r.mapped) { any = true; break; }
    if (!any) return ErrorCode::SUCCESS;
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    f << "input_id\tprotein_id\taa_start\taa_end\tdomain_id\treason\n";
    for (const auto& r : results) {
        if (r.mapped) continue;
        const auto& u = r.unmapped;
        f << u.input_id << "\t" << u.protein_id << "\t";
        if (r.domain.has_domain()) f << u.aa_start << "\t" << u.aa_end;
        else f << "NA\tNA";
        f << "\t";
        w_str_or_na(f, u.domain_id);
        f << "\t" << u.reason << "\n";
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_metadata(const std::string& path,
                         OutputKind kind,
                         const std::vector<DomainResult>& results,
                         const std::string& source,
                         const std::vector<std::string>& cli_args) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    size_t mapped = 0, unmapped = 0, no_domain = 0;
    for (const auto& r : results) {
        if (r.mapped) ++mapped; else ++unmapped;
        if (r.no_domain_mode) ++no_domain;
    }
    f << "{\n";
    f << "  \"tool\": \"protein2genomic\",\n";
    f << "  \"version\": \"2.2.0\",\n";
    f << "  \"timestamp_utc\": \"" << iso8601_now() << "\",\n";
    f << "  \"output_kind\": \"" << output_kind_to_string(kind) << "\",\n";
    f << "  \"annotation_source\": \"" << source << "\",\n";
    f << "  \"index_format_version\": " << INDEX_FORMAT_VERSION << ",\n";
    f << "  \"coordinate_conventions\": {\n";
    f << "    \"bed\": \"0-based half-open\",\n";
    f << "    \"tsv\": \"1-based inclusive (genomic and CDS nt)\",\n";
    f << "    \"aa\": \"1-based inclusive\"\n";
    f << "  },\n";
    f << "  \"query_counts\": { \"total\": " << results.size()
      << ", \"mapped\": " << mapped
      << ", \"unmapped\": " << unmapped
      << ", \"no_domain_mode\": " << no_domain << " },\n";
    f << "  \"cli\": [";
    for (size_t i = 0; i < cli_args.size(); ++i) {
        if (i) f << ", ";
        f << "\"";
        for (char c : cli_args[i]) { if (c == '\\' || c == '"') f << '\\'; f << c; }
        f << "\"";
    }
    f << "]\n";
    f << "}\n";
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

} // namespace

namespace output {

ErrorCode write_all(const std::string& out_dir,
                    OutputKind kind,
                    const std::vector<DomainResult>& results,
                    const std::string& gtf_or_index_path,
                    const std::vector<std::string>& cli_args) {
    if (!ensure_dir(out_dir)) {
        std::cerr << "Error: cannot create output directory: " << out_dir << std::endl;
        return ErrorCode::FILE_NOT_FOUND;
    }

    ErrorCode rc = write_summary(join(out_dir, "domain_mapping_summary.tsv"), results);
    if (rc != ErrorCode::SUCCESS) return rc;

    if (kind == OutputKind::CODING || kind == OutputKind::ALL) {
        rc = write_coding_tsv(join(out_dir, "domain_cds_segments.tsv"), results);
        if (rc != ErrorCode::SUCCESS) return rc;
        rc = write_coding_bed(join(out_dir, "domain_cds_segments.bed"), results);
        if (rc != ErrorCode::SUCCESS) return rc;
    }
    if (kind == OutputKind::INTRONS || kind == OutputKind::ALL) {
        rc = write_introns_tsv(join(out_dir, "domain_introns.tsv"), results);
        if (rc != ErrorCode::SUCCESS) return rc;
        rc = write_introns_bed(join(out_dir, "domain_introns.bed"), results);
        if (rc != ErrorCode::SUCCESS) return rc;
    }
    if (kind == OutputKind::SPAN || kind == OutputKind::ALL) {
        rc = write_span_bed(join(out_dir, "domain_span_with_introns.bed"), results);
        if (rc != ErrorCode::SUCCESS) return rc;
    }
    if (kind == OutputKind::ISOFORM || kind == OutputKind::ALL) {
        rc = write_isoform_tsv(join(out_dir, "isoform_structure.tsv"), results);
        if (rc != ErrorCode::SUCCESS) return rc;
    }
    if (kind == OutputKind::BED12 || kind == OutputKind::ALL) {
        rc = write_bed12(join(out_dir, "domain_blocks.bed12"), results);
        if (rc != ErrorCode::SUCCESS) return rc;
    }

    write_unmapped(join(out_dir, "unmapped_domains.tsv"), results);

    if (kind == OutputKind::ALL) {
        rc = write_metadata(join(out_dir, "run_metadata.json"), kind, results,
                            gtf_or_index_path, cli_args);
        if (rc != ErrorCode::SUCCESS) return rc;
    }
    return ErrorCode::SUCCESS;
}

} // namespace output

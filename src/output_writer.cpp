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

// Print a uint that is "0 means missing" as NA, else as the value.
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
        case OutputKind::CODING: return "coding";
        case OutputKind::SPAN: return "span";
        case OutputKind::ISOFORM: return "isoform";
        case OutputKind::ALL: return "all";
    }
    return "unknown";
}

ErrorCode write_summary(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    f << "input_id\tprotein_id\ttranscript_id\tgene_id\tgene_name\tdomain_id\tchrom\tstrand\t"
         "aa_start\taa_end\tdomain_length_aa\tdomain_length_nt\tprotein_length_aa\t"
         "domain_genomic_start\tdomain_genomic_end\tn_coding_segments\tfully_mapped\tstatus\n";
    for (const auto& r : results) {
        if (!r.mapped) {
            f << r.domain.input_id << "\t" << r.domain.protein_id
              << "\tNA\tNA\tNA\t";
            w_str_or_na(f, r.domain.domain_id);
            f << "\tNA\tNA\t" << r.domain.start << "\t" << r.domain.end
              << "\t" << (r.domain.end - r.domain.start + 1)
              << "\t" << ((r.domain.end - r.domain.start + 1) * 3)
              << "\tNA\tNA\tNA\t0\tfalse\t" << r.unmapped.reason << "\n";
            continue;
        }
        const auto& s = r.summary;
        f << s.input_id << "\t" << s.protein_id << "\t" << s.transcript_id << "\t";
        w_str_or_na(f, s.gene_id);
        f << "\t"; w_str_or_na(f, s.gene_name);
        f << "\t"; w_str_or_na(f, s.domain_id);
        f << "\t" << s.chrom << "\t" << s.strand << "\t" << s.aa_start << "\t" << s.aa_end
          << "\t" << s.domain_length_aa << "\t" << s.domain_length_nt
          << "\t" << s.protein_length_aa
          << "\t" << s.domain_genomic_start << "\t" << s.domain_genomic_end
          << "\t" << s.n_coding_segments
          << "\t" << (s.fully_mapped ? "true" : "false")
          << "\t" << (s.fully_mapped ? "ok" : "partial") << "\n";
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

std::string bed_name(const DomainResult& r) {
    std::ostringstream oss;
    oss << r.domain.protein_id;
    if (!r.domain.domain_id.empty()) oss << "_" << r.domain.domain_id;
    oss << "_" << r.domain.start << "-" << r.domain.end;
    return oss.str();
}

ErrorCode write_coding_bed(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    for (const auto& r : results) {
        if (!r.mapped) continue;
        std::string name = bed_name(r);
        for (const auto& c : r.coding_segments) {
            // BED is 0-based half-open.
            f << c.chrom << "\t" << (c.genomic_start - 1) << "\t" << c.genomic_end
              << "\t" << name << "\t0\t" << c.strand << "\n";
        }
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_coding_tsv(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    f << "input_id\tprotein_id\ttranscript_id\tgene_id\tgene_name\tdomain_id\t"
         "chrom\tstrand\tgenomic_start\tgenomic_end\t"
         "segment_index_in_domain\tcds_nt_start\tcds_nt_end\t"
         "aa_start_encoded\taa_end_encoded\taa_start\taa_end\n";
    for (const auto& r : results) {
        if (!r.mapped) continue;
        for (const auto& c : r.coding_segments) {
            f << c.input_id << "\t" << c.protein_id << "\t" << c.transcript_id << "\t";
            w_str_or_na(f, c.gene_id);
            f << "\t"; w_str_or_na(f, c.gene_name);
            f << "\t"; w_str_or_na(f, c.domain_id);
            f << "\t" << c.chrom << "\t" << c.strand
              << "\t" << c.genomic_start << "\t" << c.genomic_end
              << "\t" << c.segment_index_in_domain
              << "\t" << c.cds_nt_start << "\t" << c.cds_nt_end
              << "\t" << c.aa_start_encoded << "\t" << c.aa_end_encoded
              << "\t" << c.aa_start << "\t" << c.aa_end << "\n";
        }
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_span_bed(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    for (const auto& r : results) {
        if (!r.mapped) continue;
        std::string name = bed_name(r);
        const auto& s = r.span;
        f << s.chrom << "\t" << (s.genomic_start - 1) << "\t" << s.genomic_end
          << "\t" << name << "\t0\t" << s.strand << "\n";
    }
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_isoform_tsv(const std::string& path, const std::vector<DomainResult>& results) {
    std::ofstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    f << "input_id\tgene_id\tgene_name\ttranscript_id\tprotein_id\tdomain_id\t"
         "chrom\tstrand\tfeature_type\tfeature_id\texon_number\t"
         "feature_genomic_start\tfeature_genomic_end\tfeature_length_nt\t"
         "feature_order_genomic\tfeature_order_transcript\t"
         "cds_nt_start\tcds_nt_end\taa_start_encoded\taa_end_encoded\t"
         "overlaps_domain\t"
         "domain_overlap_genomic_start\tdomain_overlap_genomic_end\t"
         "domain_overlap_cds_nt_start\tdomain_overlap_cds_nt_end\t"
         "domain_overlap_aa_start\tdomain_overlap_aa_end\t"
         "domain_overlap_fraction_of_feature\tdomain_overlap_fraction_of_domain\t"
         "plot_group\n";
    for (const auto& r : results) {
        if (!r.mapped) continue;
        for (const auto& s : r.isoform_segments) {
            f << s.input_id << "\t";
            w_str_or_na(f, s.gene_id); f << "\t";
            w_str_or_na(f, s.gene_name); f << "\t";
            w_str_or_na(f, s.transcript_id); f << "\t";
            f << s.protein_id << "\t";
            w_str_or_na(f, s.domain_id); f << "\t";
            f << s.chrom << "\t" << s.strand << "\t"
              << plot_feature_to_string(s.feature_type) << "\t"
              << s.feature_id << "\t";
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
            f << "\t" << overlap_kind_to_string(s.overlap) << "\t";
            if (s.has_overlap_coords) {
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
        f << u.input_id << "\t" << u.protein_id << "\t"
          << u.aa_start << "\t" << u.aa_end << "\t";
        w_str_or_na(f, u.domain_id); f << "\t" << u.reason << "\n";
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
    size_t mapped = 0, unmapped = 0;
    for (const auto& r : results) (r.mapped ? mapped : unmapped)++;

    f << "{\n";
    f << "  \"tool\": \"protein2genomic\",\n";
    f << "  \"version\": \"2.0.0\",\n";
    f << "  \"timestamp_utc\": \"" << iso8601_now() << "\",\n";
    f << "  \"output_kind\": \"" << output_kind_to_string(kind) << "\",\n";
    f << "  \"annotation_source\": \"" << source << "\",\n";
    f << "  \"index_format_version\": " << INDEX_FORMAT_VERSION << ",\n";
    f << "  \"coordinate_conventions\": {\n";
    f << "    \"bed\": \"0-based half-open\",\n";
    f << "    \"tsv\": \"1-based inclusive (genomic and CDS nt)\",\n";
    f << "    \"aa\": \"1-based inclusive\"\n";
    f << "  },\n";
    f << "  \"domain_counts\": { \"total\": " << results.size()
      << ", \"mapped\": " << mapped << ", \"unmapped\": " << unmapped << " },\n";
    f << "  \"cli\": [";
    for (size_t i = 0; i < cli_args.size(); ++i) {
        if (i) f << ", ";
        f << "\"";
        for (char c : cli_args[i]) {
            if (c == '\\' || c == '"') f << '\\';
            f << c;
        }
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

    // Always write the summary so users have a single-row-per-domain index.
    ErrorCode rc = write_summary(join(out_dir, "domain_mapping_summary.tsv"), results);
    if (rc != ErrorCode::SUCCESS) return rc;

    if (kind == OutputKind::CODING || kind == OutputKind::ALL) {
        rc = write_coding_bed(join(out_dir, "domain_cds_segments.bed"), results);
        if (rc != ErrorCode::SUCCESS) return rc;
        rc = write_coding_tsv(join(out_dir, "domain_cds_segments.tsv"), results);
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

    write_unmapped(join(out_dir, "unmapped_domains.tsv"), results);

    if (kind == OutputKind::ALL) {
        rc = write_metadata(join(out_dir, "run_metadata.json"), kind, results,
                            gtf_or_index_path, cli_args);
        if (rc != ErrorCode::SUCCESS) return rc;
    }
    return ErrorCode::SUCCESS;
}

} // namespace output

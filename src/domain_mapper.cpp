#include "domain_mapper.hpp"
#include "utils.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <unordered_map>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace {

inline bool intervals_overlap(uint32_t a_s, uint32_t a_e, uint32_t b_s, uint32_t b_e) {
    return a_s <= b_e && b_s <= a_e;
}

std::vector<GenomicInterval> cds_in_translation_order(const std::vector<GenomicInterval>& intervals) {
    std::vector<GenomicInterval> cds;
    for (const auto& iv : intervals) if (iv.feature_type == FeatureType::CDS) cds.push_back(iv);
    std::sort(cds.begin(), cds.end(),
              [](const GenomicInterval& a, const GenomicInterval& b) { return a.start < b.start; });
    if (!cds.empty() && cds[0].strand == '-') std::reverse(cds.begin(), cds.end());
    return cds;
}

std::vector<GenomicInterval> exons_in_genomic_order(const std::vector<GenomicInterval>& intervals) {
    std::vector<GenomicInterval> exons;
    for (const auto& iv : intervals) if (iv.feature_type == FeatureType::EXON) exons.push_back(iv);
    std::sort(exons.begin(), exons.end(),
              [](const GenomicInterval& a, const GenomicInterval& b) { return a.start < b.start; });
    return exons;
}

struct DomainGenomic {
    struct Range {
        uint32_t genomic_start;
        uint32_t genomic_end;
        uint32_t cds_index_in_translation;
        uint32_t cumulative_nt_start;
        uint32_t cumulative_nt_end;
        uint32_t aa_overlap_start;
        uint32_t aa_overlap_end;
    };
    std::vector<Range> ranges;
    uint32_t total_protein_length_aa = 0;
    bool fully_mapped = false;
};

DomainGenomic map_domain_to_genomic(const ProteinDomain& d,
                                    const std::vector<GenomicInterval>& cds_tx_order) {
    DomainGenomic dg;
    if (cds_tx_order.empty()) return dg;

    std::vector<uint32_t> cumulative;
    cumulative.reserve(cds_tx_order.size() + 1);
    cumulative.push_back(0);
    for (const auto& c : cds_tx_order) cumulative.push_back(cumulative.back() + (c.end - c.start + 1));
    uint32_t total_nt = cumulative.back();
    dg.total_protein_length_aa = total_nt / 3;

    if (!d.has_domain()) return dg;

    uint32_t domain_nt_start = (d.start - 1) * 3;
    uint32_t domain_nt_end = d.end * 3 - 1;
    if (domain_nt_end >= total_nt) {
        if (domain_nt_start >= total_nt) return dg;
        domain_nt_end = total_nt - 1;
    } else {
        dg.fully_mapped = true;
    }

    for (size_t i = 0; i < cds_tx_order.size(); ++i) {
        uint32_t cds_nt_start = cumulative[i];
        uint32_t cds_nt_end = cumulative[i + 1] - 1;
        if (!intervals_overlap(domain_nt_start, domain_nt_end, cds_nt_start, cds_nt_end)) continue;

        uint32_t ovl_nt_start = std::max(domain_nt_start, cds_nt_start);
        uint32_t ovl_nt_end = std::min(domain_nt_end, cds_nt_end);
        uint32_t off_s = ovl_nt_start - cds_nt_start;
        uint32_t off_e = ovl_nt_end - cds_nt_start;

        const auto& c = cds_tx_order[i];
        uint32_t g_s, g_e;
        if (c.strand == '-') { g_s = c.end - off_e; g_e = c.end - off_s; }
        else                 { g_s = c.start + off_s; g_e = c.start + off_e; }

        DomainGenomic::Range r;
        r.genomic_start = g_s; r.genomic_end = g_e;
        r.cds_index_in_translation = static_cast<uint32_t>(i);
        r.cumulative_nt_start = ovl_nt_start;
        r.cumulative_nt_end = ovl_nt_end;
        r.aa_overlap_start = ovl_nt_start / 3 + 1;
        r.aa_overlap_end = ovl_nt_end / 3 + 1;
        dg.ranges.push_back(r);
    }
    return dg;
}

bool genomic_overlaps_any_domain_range(uint32_t s, uint32_t e,
                                       const std::vector<DomainGenomic::Range>& ranges,
                                       uint32_t& ovl_s, uint32_t& ovl_e) {
    bool found = false;
    for (const auto& r : ranges) {
        if (intervals_overlap(s, e, r.genomic_start, r.genomic_end)) {
            uint32_t a = std::max(s, r.genomic_start);
            uint32_t b = std::min(e, r.genomic_end);
            if (!found) { ovl_s = a; ovl_e = b; found = true; }
            else { ovl_s = std::min(ovl_s, a); ovl_e = std::max(ovl_e, b); }
        }
    }
    return found;
}

struct CdsCoord { uint32_t cds_nt_start; uint32_t cds_nt_end; uint32_t aa_start; uint32_t aa_end; };
CdsCoord cds_genomic_to_aa(const GenomicInterval& c, uint32_t cum_nt_start,
                           uint32_t g_s, uint32_t g_e) {
    uint32_t off_s, off_e;
    if (c.strand == '-') { off_s = c.end - g_e; off_e = c.end - g_s; }
    else                 { off_s = g_s - c.start; off_e = g_e - c.start; }
    uint32_t nt_s = cum_nt_start + off_s;
    uint32_t nt_e = cum_nt_start + off_e;
    CdsCoord cc;
    cc.cds_nt_start = nt_s + 1;
    cc.cds_nt_end = nt_e + 1;
    cc.aa_start = nt_s / 3 + 1;
    cc.aa_end = nt_e / 3 + 1;
    return cc;
}

} // namespace

DomainMapper::DomainMapper(const GTFParser& parser, OutputKind output_kind)
    : gtf_(parser), output_kind_(output_kind) {
    domains_.reserve(100000);
}

ErrorCode DomainMapper::load_domains(const std::string& bed_filename) {
    std::ifstream file(bed_filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open BED file: " << bed_filename << std::endl;
        return ErrorCode::FILE_NOT_FOUND;
    }
    domains_.clear();
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::vector<std::string> cols;
        std::string tok;
        while (iss >> tok) cols.push_back(tok);
        if (cols.empty()) continue;

        ProteinDomain d;
        d.protein_id = cols[0];

        // Domain coords are optional. Missing -> no-domain mode.
        if (cols.size() >= 3) {
            try {
                d.start = static_cast<uint32_t>(std::stoul(cols[1]));
                d.end   = static_cast<uint32_t>(std::stoul(cols[2]));
            } catch (...) { d.start = 0; d.end = 0; }
        }
        if (cols.size() >= 4) d.domain_id = cols[3];

        if (!d.domain_id.empty()) {
            d.input_id = d.domain_id;
        } else if (d.has_domain()) {
            std::ostringstream oss;
            oss << d.protein_id << ":" << d.start << "-" << d.end;
            d.input_id = oss.str();
        } else {
            d.input_id = d.protein_id;
        }
        domains_.push_back(std::move(d));
    }
    std::cerr << "Loaded " << domains_.size() << " queries from " << bed_filename << std::endl;
    return ErrorCode::SUCCESS;
}

DomainResult DomainMapper::process_one(const ProteinDomain& d) const {
    DomainResult R;
    R.domain = d;
    R.no_domain_mode = !d.has_domain();

    const auto* intervals = gtf_.get_protein_intervals(d.protein_id);
    if (!intervals || intervals->empty()) {
        R.unmapped.input_id = d.input_id;
        R.unmapped.protein_id = d.protein_id;
        R.unmapped.aa_start = d.start;
        R.unmapped.aa_end = d.end;
        R.unmapped.domain_id = d.domain_id;
        R.unmapped.reason = "protein_not_in_index";
        return R;
    }

    auto cds_tx = cds_in_translation_order(*intervals);
    auto exons = exons_in_genomic_order(*intervals);
    if (cds_tx.empty()) {
        R.unmapped.input_id = d.input_id;
        R.unmapped.protein_id = d.protein_id;
        R.unmapped.aa_start = d.start;
        R.unmapped.aa_end = d.end;
        R.unmapped.domain_id = d.domain_id;
        R.unmapped.reason = "no_CDS_for_protein";
        return R;
    }

    char strand = cds_tx.front().strand;
    std::string chrom = cds_tx.front().chromosome;
    std::string transcript_id = cds_tx.front().transcript_id;
    std::string gene_id = gtf_.get_gene_id(d.protein_id);
    if (gene_id.empty()) gene_id = cds_tx.front().gene_id;
    std::string gene_name = gtf_.get_gene_name(d.protein_id);
    if (gene_name.empty()) gene_name = cds_tx.front().gene_name;

    DomainGenomic dg = map_domain_to_genomic(d, cds_tx);

    // Unmapped if a domain was requested but didn't intersect any CDS.
    if (d.has_domain() && dg.ranges.empty()) {
        R.unmapped.input_id = d.input_id;
        R.unmapped.protein_id = d.protein_id;
        R.unmapped.aa_start = d.start;
        R.unmapped.aa_end = d.end;
        R.unmapped.domain_id = d.domain_id;
        R.unmapped.reason = dg.total_protein_length_aa == 0
            ? "no_CDS_for_protein"
            : (d.start > dg.total_protein_length_aa
                ? "domain_beyond_protein_length"
                : "no_overlap");
        return R;
    }

    R.mapped = true;
    uint32_t domain_length_nt = d.has_domain() ? (d.end - d.start + 1) * 3 : 0;
    uint32_t domain_g_min = 0, domain_g_max = 0;
    if (!dg.ranges.empty()) {
        domain_g_min = dg.ranges.front().genomic_start;
        domain_g_max = dg.ranges.front().genomic_end;
        for (const auto& r : dg.ranges) {
            domain_g_min = std::min(domain_g_min, r.genomic_start);
            domain_g_max = std::max(domain_g_max, r.genomic_end);
        }
    }

    // ----- span (only when a domain exists and was mapped) -----
    if (d.has_domain() && !dg.ranges.empty()) {
        SpanRow s;
        s.chrom = chrom;
        s.genomic_start = domain_g_min;
        s.genomic_end = domain_g_max;
        s.strand = strand;
        s.input_id = d.input_id;
        s.protein_id = d.protein_id;
        s.domain_id = d.domain_id;
        s.gene_name = gene_name;
        R.span = s;
        R.has_span = true;
    }

    // ----- summary -----
    {
        MappingSummaryRow s;
        s.input_id = d.input_id;
        s.protein_id = d.protein_id;
        s.transcript_id = transcript_id;
        s.gene_id = gene_id;
        s.gene_name = gene_name;
        s.domain_id = d.domain_id;
        s.chrom = chrom;
        s.strand = strand;
        s.aa_start = d.start;
        s.aa_end = d.end;
        s.domain_length_aa = d.has_domain() ? (d.end - d.start + 1) : 0;
        s.domain_length_nt = domain_length_nt;
        s.protein_length_aa = dg.total_protein_length_aa;
        s.domain_genomic_start = domain_g_min;
        s.domain_genomic_end = domain_g_max;
        s.n_coding_segments = static_cast<uint32_t>(dg.ranges.size());
        s.fully_mapped = d.has_domain() ? dg.fully_mapped : true;
        s.no_domain_mode = !d.has_domain();
        R.summary = s;
    }

    // ----- isoform_structure rows -----
    // We always build them: they back coding, introns, isoform, and (via filtering) all.

    // Cumulative nt offsets per CDS in translation order.
    std::vector<uint32_t> cum_nt(cds_tx.size() + 1, 0);
    for (size_t i = 0; i < cds_tx.size(); ++i) cum_nt[i + 1] = cum_nt[i] + (cds_tx[i].end - cds_tx[i].start + 1);
    // Map CDS genomic start -> translation index.
    std::unordered_map<uint32_t, size_t> cds_start_to_tx_idx;
    for (size_t i = 0; i < cds_tx.size(); ++i) cds_start_to_tx_idx[cds_tx[i].start] = i;

    std::vector<GenomicInterval> cds_gen = cds_tx;
    std::sort(cds_gen.begin(), cds_gen.end(),
              [](const GenomicInterval& a, const GenomicInterval& b) { return a.start < b.start; });
    uint32_t cds_min_g = cds_gen.front().start;
    uint32_t cds_max_g = cds_gen.back().end;

    auto classify_utr = [&](uint32_t s, uint32_t e) -> PlotFeatureType {
        if (strand == '+') {
            return (e < cds_min_g) ? PlotFeatureType::FIVE_PRIME_UTR
                                   : PlotFeatureType::THREE_PRIME_UTR;
        } else {
            return (e < cds_min_g) ? PlotFeatureType::THREE_PRIME_UTR
                                   : PlotFeatureType::FIVE_PRIME_UTR;
        }
    };

    // Initial rows in genomic order: introns between exons, then UTR/CDS within each exon.
    std::vector<IsoformSegmentRow> rows;
    rows.reserve(exons.size() * 3);

    for (size_t ei = 0; ei < exons.size(); ++ei) {
        if (ei > 0) {
            const auto& prev = exons[ei - 1];
            const auto& cur = exons[ei];
            if (prev.end + 1 <= cur.start - 1) {
                IsoformSegmentRow ir;
                ir.chrom = chrom; ir.strand = strand;
                ir.feature_type = PlotFeatureType::INTRON;
                ir.feature_genomic_start = prev.end + 1;
                ir.feature_genomic_end = cur.start - 1;
                ir.feature_length_nt = ir.feature_genomic_end - ir.feature_genomic_start + 1;
                rows.push_back(std::move(ir));
            }
        }
        const auto& exon = exons[ei];
        std::vector<const GenomicInterval*> covering;
        for (const auto& c : cds_gen) {
            if (intervals_overlap(c.start, c.end, exon.start, exon.end)) covering.push_back(&c);
        }
        std::sort(covering.begin(), covering.end(),
                  [](const GenomicInterval* a, const GenomicInterval* b) { return a->start < b->start; });

        uint32_t cursor = exon.start;
        for (const auto* cp : covering) {
            uint32_t cs = std::max(cp->start, exon.start);
            uint32_t ce = std::min(cp->end, exon.end);
            if (cursor < cs) {
                IsoformSegmentRow ir;
                ir.chrom = chrom; ir.strand = strand;
                ir.exon_number = exon.exon_number;
                ir.feature_genomic_start = cursor;
                ir.feature_genomic_end = cs - 1;
                ir.feature_length_nt = ir.feature_genomic_end - ir.feature_genomic_start + 1;
                ir.feature_type = classify_utr(ir.feature_genomic_start, ir.feature_genomic_end);
                rows.push_back(std::move(ir));
            }
            {
                IsoformSegmentRow ir;
                ir.chrom = chrom; ir.strand = strand;
                ir.exon_number = exon.exon_number;
                ir.feature_type = PlotFeatureType::CDS;
                ir.feature_genomic_start = cs;
                ir.feature_genomic_end = ce;
                ir.feature_length_nt = ce - cs + 1;
                ir.has_cds_coords = true;
                auto it = cds_start_to_tx_idx.find(cp->start);
                if (it != cds_start_to_tx_idx.end()) {
                    size_t tx_i = it->second;
                    ir.source_cds_index = static_cast<uint32_t>(tx_i);
                    CdsCoord cc = cds_genomic_to_aa(cds_tx[tx_i], cum_nt[tx_i], cs, ce);
                    ir.cds_nt_start = cc.cds_nt_start;
                    ir.cds_nt_end = cc.cds_nt_end;
                    ir.aa_start_encoded = cc.aa_start;
                    ir.aa_end_encoded = cc.aa_end;
                }
                rows.push_back(std::move(ir));
            }
            cursor = ce + 1;
        }
        if (cursor <= exon.end) {
            IsoformSegmentRow ir;
            ir.chrom = chrom; ir.strand = strand;
            ir.exon_number = exon.exon_number;
            ir.feature_genomic_start = cursor;
            ir.feature_genomic_end = exon.end;
            ir.feature_length_nt = ir.feature_genomic_end - ir.feature_genomic_start + 1;
            ir.feature_type = classify_utr(ir.feature_genomic_start, ir.feature_genomic_end);
            rows.push_back(std::move(ir));
        }
    }

    // Split CDS rows where the domain partially overlaps.
    std::vector<IsoformSegmentRow> rows2;
    rows2.reserve(rows.size() * 2);
    for (auto& ir : rows) {
        if (ir.feature_type != PlotFeatureType::CDS || !d.has_domain()) {
            rows2.push_back(std::move(ir));
            continue;
        }
        uint32_t ovl_s = 0, ovl_e = 0;
        bool any = genomic_overlaps_any_domain_range(ir.feature_genomic_start, ir.feature_genomic_end,
                                                    dg.ranges, ovl_s, ovl_e);
        if (!any) { rows2.push_back(std::move(ir)); continue; }
        auto make_cds_row = [&](uint32_t s, uint32_t e) {
            IsoformSegmentRow nr = ir;
            nr.feature_genomic_start = s;
            nr.feature_genomic_end = e;
            nr.feature_length_nt = e - s + 1;
            size_t tx_i = ir.source_cds_index;
            CdsCoord cc = cds_genomic_to_aa(cds_tx[tx_i], cum_nt[tx_i], s, e);
            nr.cds_nt_start = cc.cds_nt_start;
            nr.cds_nt_end = cc.cds_nt_end;
            nr.aa_start_encoded = cc.aa_start;
            nr.aa_end_encoded = cc.aa_end;
            nr.has_cds_coords = true;
            return nr;
        };
        if (ir.feature_genomic_start < ovl_s) rows2.push_back(make_cds_row(ir.feature_genomic_start, ovl_s - 1));
        rows2.push_back(make_cds_row(ovl_s, ovl_e));
        if (ovl_e < ir.feature_genomic_end) rows2.push_back(make_cds_row(ovl_e + 1, ir.feature_genomic_end));
    }

    // Genomic ordering.
    uint32_t ord_g = 0;
    for (auto& ir : rows2) ir.feature_order_genomic = ++ord_g;
    if (strand == '-') {
        uint32_t n = static_cast<uint32_t>(rows2.size());
        for (auto& ir : rows2) ir.feature_order_transcript = n - ir.feature_order_genomic + 1;
    } else {
        for (auto& ir : rows2) ir.feature_order_transcript = ir.feature_order_genomic;
    }

    // Assign feature_id and feature_part.
    //   CDS:    feature_id = "CDS_<source_cds_index+1>" (stable across splits)
    //           feature_part = 1..K within rows sharing the same source_cds_index,
    //                           in translation order (= order_transcript ascending).
    //   UTR:    numbered separately per UTR kind, in translation order
    //   intron: numbered in translation order
    {
        // Index of rows by translation order.
        std::vector<size_t> idx(rows2.size());
        for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
        std::sort(idx.begin(), idx.end(),
                  [&](size_t a, size_t b) {
                      return rows2[a].feature_order_transcript < rows2[b].feature_order_transcript;
                  });

        std::unordered_map<int, uint32_t> per_utr_counter;       // by PlotFeatureType
        uint32_t intron_counter = 0;
        // For CDS: source_cds_index -> running part counter
        std::unordered_map<uint32_t, uint32_t> cds_part_counter;

        for (size_t k : idx) {
            auto& r = rows2[k];
            switch (r.feature_type) {
                case PlotFeatureType::CDS: {
                    r.feature_id = "CDS_" + std::to_string(r.source_cds_index + 1);
                    r.feature_part = ++cds_part_counter[r.source_cds_index];
                    break;
                }
                case PlotFeatureType::INTRON: {
                    ++intron_counter;
                    r.feature_id = "intron_" + std::to_string(intron_counter);
                    r.feature_part = 1;
                    break;
                }
                case PlotFeatureType::FIVE_PRIME_UTR:
                case PlotFeatureType::THREE_PRIME_UTR: {
                    int key = static_cast<int>(r.feature_type);
                    uint32_t n = ++per_utr_counter[key];
                    r.feature_id = plot_feature_to_string(r.feature_type) + "_" + std::to_string(n);
                    r.feature_part = 1;
                    break;
                }
            }
        }
    }

    // Identity + overlap classification + plot_group.
    for (auto& r : rows2) {
        r.input_id = d.input_id;
        r.gene_id = gene_id;
        r.gene_name = gene_name;
        r.transcript_id = transcript_id;
        r.protein_id = d.protein_id;
        r.domain_id = d.domain_id;
        r.no_domain_mode = !d.has_domain();

        if (!d.has_domain()) {
            r.overlap = DomainOverlapKind::NO;
            r.has_overlap_coords = false;
            // plot_group: just the feature type, no domain coloring.
            r.plot_group = plot_feature_to_string(r.feature_type);
            continue;
        }

        if (r.feature_type == PlotFeatureType::CDS) {
            uint32_t ovl_s = 0, ovl_e = 0;
            bool any = genomic_overlaps_any_domain_range(r.feature_genomic_start, r.feature_genomic_end,
                                                        dg.ranges, ovl_s, ovl_e);
            if (any) {
                r.overlap = DomainOverlapKind::CODING_OVERLAP;
                r.has_overlap_coords = true;
                r.domain_overlap_genomic_start = ovl_s;
                r.domain_overlap_genomic_end = ovl_e;
                size_t tx_i = r.source_cds_index;
                CdsCoord cc = cds_genomic_to_aa(cds_tx[tx_i], cum_nt[tx_i], ovl_s, ovl_e);
                r.domain_overlap_cds_nt_start = cc.cds_nt_start;
                r.domain_overlap_cds_nt_end = cc.cds_nt_end;
                r.domain_overlap_aa_start = cc.aa_start;
                r.domain_overlap_aa_end = cc.aa_end;
                double feat_len = static_cast<double>(r.feature_length_nt);
                double ovl_len = static_cast<double>(ovl_e - ovl_s + 1);
                r.domain_overlap_fraction_of_feature = feat_len > 0 ? ovl_len / feat_len : 0.0;
                r.domain_overlap_fraction_of_domain = domain_length_nt > 0
                    ? ovl_len / static_cast<double>(domain_length_nt) : 0.0;
            } else {
                r.overlap = DomainOverlapKind::NO;
            }
        } else if (r.feature_type == PlotFeatureType::INTRON) {
            r.overlap = (r.feature_genomic_start >= domain_g_min && r.feature_genomic_end <= domain_g_max)
                ? DomainOverlapKind::INSIDE_DOMAIN_GENOMIC_SPAN
                : DomainOverlapKind::NO;
        } else {
            r.overlap = DomainOverlapKind::NO;
        }
        r.plot_group = plot_group_for(r.feature_type, r.overlap);
    }

    R.isoform_segments = std::move(rows2);
    return R;
}

ErrorCode DomainMapper::process_domains() {
    results_.clear();
    results_.resize(domains_.size());
    auto t0 = std::chrono::high_resolution_clock::now();
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic, 64)
#endif
    for (size_t i = 0; i < domains_.size(); ++i) {
        results_[i] = process_one(domains_[i]);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    size_t mapped = 0, unmapped = 0;
    for (const auto& r : results_) (r.mapped ? mapped : unmapped)++;
    std::cerr << "Mapped " << mapped << "/" << domains_.size()
              << " queries (" << unmapped << " unmapped) in " << ms << "ms" << std::endl;
    return ErrorCode::SUCCESS;
}

size_t DomainMapper::mapped_count() const {
    size_t n = 0;
    for (const auto& r : results_) if (r.mapped) ++n;
    return n;
}

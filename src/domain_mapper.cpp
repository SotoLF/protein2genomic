#include "domain_mapper.hpp"
#include "utils.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace {

inline bool intervals_overlap(uint32_t a_s, uint32_t a_e, uint32_t b_s, uint32_t b_e) {
    return a_s <= b_e && b_s <= a_e;
}

// CDS intervals in translation order: + strand keeps genomic order, - strand reversed.
std::vector<GenomicInterval> cds_in_translation_order(const std::vector<GenomicInterval>& intervals) {
    std::vector<GenomicInterval> cds;
    for (const auto& iv : intervals) {
        if (iv.feature_type == FeatureType::CDS) cds.push_back(iv);
    }
    std::sort(cds.begin(), cds.end(),
              [](const GenomicInterval& a, const GenomicInterval& b) { return a.start < b.start; });
    if (!cds.empty() && cds[0].strand == '-') std::reverse(cds.begin(), cds.end());
    return cds;
}

std::vector<GenomicInterval> exons_in_genomic_order(const std::vector<GenomicInterval>& intervals) {
    std::vector<GenomicInterval> exons;
    for (const auto& iv : intervals) {
        if (iv.feature_type == FeatureType::EXON) exons.push_back(iv);
    }
    std::sort(exons.begin(), exons.end(),
              [](const GenomicInterval& a, const GenomicInterval& b) { return a.start < b.start; });
    return exons;
}

struct DomainGenomic {
    // One sub-range per CDS-in-translation-order that the domain overlaps.
    // Each entry stores the genomic [start,end] of the overlap and the
    // translation-frame nt offsets [nt_start_in_cds, nt_end_in_cds] (1-based,
    // domain-relative? -> we use cds-relative cumulative).
    struct Range {
        uint32_t genomic_start;
        uint32_t genomic_end;
        uint32_t cds_index_in_translation;    // 0-based
        uint32_t nt_offset_in_cds_start;      // 0-based offset within that CDS interval
        uint32_t nt_offset_in_cds_end;
        uint32_t aa_overlap_start;            // 1-based
        uint32_t aa_overlap_end;              // 1-based
        uint32_t cumulative_nt_start;         // 0-based across the whole CDS
        uint32_t cumulative_nt_end;
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
    for (const auto& c : cds_tx_order) {
        cumulative.push_back(cumulative.back() + (c.end - c.start + 1));
    }
    uint32_t total_nt = cumulative.back();
    dg.total_protein_length_aa = total_nt / 3;

    // 0-based half-open in nt
    uint32_t domain_nt_start = (d.start - 1) * 3;
    uint32_t domain_nt_end = d.end * 3 - 1;
    if (domain_nt_end >= total_nt) {
        // Domain extends beyond CDS — try to clip and still report what we can.
        if (domain_nt_start >= total_nt) return dg; // entirely beyond
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
        if (c.strand == '-') {
            g_s = c.end - off_e;
            g_e = c.end - off_s;
        } else {
            g_s = c.start + off_s;
            g_e = c.start + off_e;
        }

        DomainGenomic::Range r;
        r.genomic_start = g_s;
        r.genomic_end = g_e;
        r.cds_index_in_translation = static_cast<uint32_t>(i);
        r.nt_offset_in_cds_start = off_s;
        r.nt_offset_in_cds_end = off_e;
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

// Map a genomic [g_s,g_e] within a known CDS interval c to (cds_nt, aa) coordinates,
// assuming CDS i in translation order with cumulative_nt_start = cum_nt_start.
struct CdsCoord { uint32_t cds_nt_start; uint32_t cds_nt_end; uint32_t aa_start; uint32_t aa_end; };
CdsCoord cds_genomic_to_aa(const GenomicInterval& c, uint32_t cum_nt_start,
                           uint32_t g_s, uint32_t g_e) {
    uint32_t off_s, off_e;
    if (c.strand == '-') {
        off_s = c.end - g_e;
        off_e = c.end - g_s;
    } else {
        off_s = g_s - c.start;
        off_e = g_e - c.start;
    }
    uint32_t nt_s = cum_nt_start + off_s; // 0-based
    uint32_t nt_e = cum_nt_start + off_e;
    CdsCoord cc;
    cc.cds_nt_start = nt_s + 1; // 1-based for output
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
    size_t line_no = 0;
    while (std::getline(file, line)) {
        ++line_no;
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        ProteinDomain d;
        std::string col1, col2, col3, col4, col5;
        iss >> col1 >> col2 >> col3;
        if (col1.empty() || col2.empty() || col3.empty()) continue;
        d.protein_id = col1;
        try {
            d.start = static_cast<uint32_t>(std::stoul(col2));
            d.end = static_cast<uint32_t>(std::stoul(col3));
        } catch (...) { continue; }
        if (iss >> col4) d.domain_id = col4;
        // input_id: prefer domain_id; fall back to protein:start-end so each row is identifiable.
        if (!d.domain_id.empty()) d.input_id = d.domain_id;
        else {
            std::ostringstream oss;
            oss << d.protein_id << ":" << d.start << "-" << d.end;
            d.input_id = oss.str();
        }
        domains_.push_back(std::move(d));
    }
    std::cerr << "Loaded " << domains_.size() << " domains from " << bed_filename << std::endl;
    return ErrorCode::SUCCESS;
}

DomainResult DomainMapper::process_one(const ProteinDomain& d) const {
    DomainResult R;
    R.domain = d;

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
    if (dg.ranges.empty()) {
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

    // Genomic envelope (min start / max end across all coding ranges).
    uint32_t domain_g_min = dg.ranges.front().genomic_start;
    uint32_t domain_g_max = dg.ranges.front().genomic_end;
    for (const auto& r : dg.ranges) {
        domain_g_min = std::min(domain_g_min, r.genomic_start);
        domain_g_max = std::max(domain_g_max, r.genomic_end);
    }
    uint32_t domain_length_nt = (d.end - d.start + 1) * 3;

    // ----- coding segments -----
    // dg.ranges is in translation order; each entry corresponds to one CDS slice.
    {
        uint32_t segment_idx = 0;
        for (const auto& r : dg.ranges) {
            ++segment_idx;
            CodingSegmentRow row;
            row.chrom = chrom;
            row.genomic_start = r.genomic_start;
            row.genomic_end = r.genomic_end;
            row.strand = strand;
            row.input_id = d.input_id;
            row.protein_id = d.protein_id;
            row.transcript_id = transcript_id;
            row.gene_id = gene_id;
            row.gene_name = gene_name;
            row.domain_id = d.domain_id;
            row.aa_start = d.start;
            row.aa_end = d.end;
            row.segment_index_in_domain = segment_idx;
            row.cds_nt_start = r.cumulative_nt_start + 1;
            row.cds_nt_end = r.cumulative_nt_end + 1;
            row.aa_start_encoded = r.aa_overlap_start;
            row.aa_end_encoded = r.aa_overlap_end;
            R.coding_segments.push_back(std::move(row));
        }
    }

    // ----- span -----
    {
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
        s.domain_length_aa = d.end - d.start + 1;
        s.domain_length_nt = domain_length_nt;
        s.protein_length_aa = dg.total_protein_length_aa;
        s.domain_genomic_start = domain_g_min;
        s.domain_genomic_end = domain_g_max;
        s.n_coding_segments = static_cast<uint32_t>(dg.ranges.size());
        s.fully_mapped = dg.fully_mapped;
        R.summary = s;
    }

    // ----- isoform structure -----
    if (output_kind_ == OutputKind::ISOFORM || output_kind_ == OutputKind::ALL) {
        // Pre-compute cumulative nt offsets per CDS interval (translation order).
        std::vector<uint32_t> cum_nt(cds_tx.size() + 1, 0);
        for (size_t i = 0; i < cds_tx.size(); ++i) {
            cum_nt[i + 1] = cum_nt[i] + (cds_tx[i].end - cds_tx[i].start + 1);
        }
        // Lookup: genomic CDS interval -> index in translation order.
        // Key by start coord (CDS regions don't overlap within a transcript).
        std::unordered_map<uint32_t, size_t> cds_start_to_tx_idx;
        for (size_t i = 0; i < cds_tx.size(); ++i) cds_start_to_tx_idx[cds_tx[i].start] = i;

        // CDS in genomic order.
        std::vector<GenomicInterval> cds_gen = cds_tx;
        std::sort(cds_gen.begin(), cds_gen.end(),
                  [](const GenomicInterval& a, const GenomicInterval& b) { return a.start < b.start; });
        uint32_t cds_min_g = cds_gen.front().start;
        uint32_t cds_max_g = cds_gen.back().end;

        // Build raw segments in genomic order: introns between exons, then UTR/CDS within each exon.
        std::vector<IsoformSegmentRow> rows;
        rows.reserve(exons.size() * 3);

        auto classify_utr = [&](uint32_t s, uint32_t e) -> PlotFeatureType {
            // Region s..e is wholly outside CDS span by construction (caller ensures this).
            if (strand == '+') {
                return (e < cds_min_g) ? PlotFeatureType::FIVE_PRIME_UTR
                                       : PlotFeatureType::THREE_PRIME_UTR;
            } else {
                return (e < cds_min_g) ? PlotFeatureType::THREE_PRIME_UTR
                                       : PlotFeatureType::FIVE_PRIME_UTR;
            }
        };

        for (size_t ei = 0; ei < exons.size(); ++ei) {
            // Intron before this exon (between previous exon and this one).
            if (ei > 0) {
                const auto& prev = exons[ei - 1];
                const auto& cur = exons[ei];
                if (prev.end + 1 <= cur.start - 1) {
                    IsoformSegmentRow ir;
                    ir.chrom = chrom;
                    ir.strand = strand;
                    ir.feature_type = PlotFeatureType::INTRON;
                    ir.feature_genomic_start = prev.end + 1;
                    ir.feature_genomic_end = cur.start - 1;
                    ir.feature_length_nt = ir.feature_genomic_end - ir.feature_genomic_start + 1;
                    ir.has_cds_coords = false;
                    rows.push_back(std::move(ir));
                }
            }

            const auto& exon = exons[ei];
            // CDS slices that overlap this exon, sorted by genomic start.
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
                    // UTR slice [cursor, cs-1]
                    IsoformSegmentRow ir;
                    ir.chrom = chrom;
                    ir.strand = strand;
                    ir.exon_number = exon.exon_number;
                    ir.feature_genomic_start = cursor;
                    ir.feature_genomic_end = cs - 1;
                    ir.feature_length_nt = ir.feature_genomic_end - ir.feature_genomic_start + 1;
                    ir.feature_type = classify_utr(ir.feature_genomic_start, ir.feature_genomic_end);
                    ir.has_cds_coords = false;
                    rows.push_back(std::move(ir));
                }
                {
                    IsoformSegmentRow ir;
                    ir.chrom = chrom;
                    ir.strand = strand;
                    ir.exon_number = exon.exon_number;
                    ir.feature_type = PlotFeatureType::CDS;
                    ir.feature_genomic_start = cs;
                    ir.feature_genomic_end = ce;
                    ir.feature_length_nt = ce - cs + 1;
                    ir.has_cds_coords = true;
                    auto it = cds_start_to_tx_idx.find(cp->start);
                    if (it != cds_start_to_tx_idx.end()) {
                        size_t tx_i = it->second;
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
                ir.chrom = chrom;
                ir.strand = strand;
                ir.exon_number = exon.exon_number;
                ir.feature_genomic_start = cursor;
                ir.feature_genomic_end = exon.end;
                ir.feature_length_nt = ir.feature_genomic_end - ir.feature_genomic_start + 1;
                ir.feature_type = classify_utr(ir.feature_genomic_start, ir.feature_genomic_end);
                ir.has_cds_coords = false;
                rows.push_back(std::move(ir));
            }
        }

        // Annotate domain overlap, fill orderings + plot_group + identity columns.
        // Need split-CDS rows for partial domain overlap. Walk and, for any CDS row
        // that overlaps the domain partially, split into up to three rows.
        std::vector<IsoformSegmentRow> rows2;
        rows2.reserve(rows.size() * 2);
        for (auto& ir : rows) {
            if (ir.feature_type != PlotFeatureType::CDS) {
                rows2.push_back(std::move(ir));
                continue;
            }
            // Split CDS by overlap with any dg.ranges.
            uint32_t ovl_s = 0, ovl_e = 0;
            bool any = genomic_overlaps_any_domain_range(ir.feature_genomic_start, ir.feature_genomic_end,
                                                        dg.ranges, ovl_s, ovl_e);
            if (!any) { rows2.push_back(std::move(ir)); continue; }
            // Build up to 3 rows: [start, ovl_s-1] no-overlap, [ovl_s, ovl_e] overlap, [ovl_e+1, end] no-overlap.
            auto make_cds_row = [&](uint32_t s, uint32_t e) {
                IsoformSegmentRow nr = ir;
                nr.feature_genomic_start = s;
                nr.feature_genomic_end = e;
                nr.feature_length_nt = e - s + 1;
                // Recompute CDS coords for this slice.
                size_t tx_i = 0;
                for (size_t i = 0; i < cds_tx.size(); ++i) {
                    if (cds_tx[i].start <= s && e <= cds_tx[i].end) { tx_i = i; break; }
                }
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

        // Annotate ordering and overlap classification.
        uint32_t ord_g = 0;
        for (auto& ir : rows2) ir.feature_order_genomic = ++ord_g;

        // Translation order: + strand same as genomic, - strand reversed.
        if (strand == '-') {
            uint32_t n = static_cast<uint32_t>(rows2.size());
            for (auto& ir : rows2) ir.feature_order_transcript = n - ir.feature_order_genomic + 1;
        } else {
            for (auto& ir : rows2) ir.feature_order_transcript = ir.feature_order_genomic;
        }

        // Per-feature counters (CDS_1, intron_1, ...) numbered in translation order.
        // Build lookup of rows sorted by transcript order, assign numeric ids per feature_type.
        std::vector<size_t> idx(rows2.size());
        for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
        std::sort(idx.begin(), idx.end(),
                  [&](size_t a, size_t b) {
                      return rows2[a].feature_order_transcript < rows2[b].feature_order_transcript;
                  });
        std::unordered_map<int, uint32_t> per_type_counter;
        for (size_t k : idx) {
            int key = static_cast<int>(rows2[k].feature_type);
            uint32_t n = ++per_type_counter[key];
            rows2[k].feature_id = plot_feature_to_string(rows2[k].feature_type) + "_" + std::to_string(n);
        }

        // Domain overlap classification + coordinates.
        for (auto& ir : rows2) {
            ir.input_id = d.input_id;
            ir.gene_id = gene_id;
            ir.gene_name = gene_name;
            ir.transcript_id = transcript_id;
            ir.protein_id = d.protein_id;
            ir.domain_id = d.domain_id;

            if (ir.feature_type == PlotFeatureType::CDS) {
                uint32_t ovl_s = 0, ovl_e = 0;
                bool any = genomic_overlaps_any_domain_range(ir.feature_genomic_start, ir.feature_genomic_end,
                                                            dg.ranges, ovl_s, ovl_e);
                if (any) {
                    ir.overlap = DomainOverlapKind::CODING_OVERLAP;
                    ir.has_overlap_coords = true;
                    ir.domain_overlap_genomic_start = ovl_s;
                    ir.domain_overlap_genomic_end = ovl_e;
                    // CDS-relative coords for the overlap slice.
                    size_t tx_i = 0;
                    for (size_t i = 0; i < cds_tx.size(); ++i) {
                        if (cds_tx[i].start <= ovl_s && ovl_e <= cds_tx[i].end) { tx_i = i; break; }
                    }
                    CdsCoord cc = cds_genomic_to_aa(cds_tx[tx_i], cum_nt[tx_i], ovl_s, ovl_e);
                    ir.domain_overlap_cds_nt_start = cc.cds_nt_start;
                    ir.domain_overlap_cds_nt_end = cc.cds_nt_end;
                    ir.domain_overlap_aa_start = cc.aa_start;
                    ir.domain_overlap_aa_end = cc.aa_end;
                    double feat_len = static_cast<double>(ir.feature_length_nt);
                    double ovl_len = static_cast<double>(ovl_e - ovl_s + 1);
                    ir.domain_overlap_fraction_of_feature = feat_len > 0 ? ovl_len / feat_len : 0.0;
                    ir.domain_overlap_fraction_of_domain = domain_length_nt > 0
                        ? ovl_len / static_cast<double>(domain_length_nt) : 0.0;
                } else {
                    ir.overlap = DomainOverlapKind::NO;
                }
            } else if (ir.feature_type == PlotFeatureType::INTRON) {
                if (ir.feature_genomic_start >= domain_g_min && ir.feature_genomic_end <= domain_g_max) {
                    ir.overlap = DomainOverlapKind::INSIDE_DOMAIN_GENOMIC_SPAN;
                } else {
                    ir.overlap = DomainOverlapKind::NO;
                }
            } else {
                ir.overlap = DomainOverlapKind::NO;
            }
            ir.plot_group = plot_group_for(ir.feature_type, ir.overlap);
        }

        R.isoform_segments = std::move(rows2);
    }

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
              << " domains (" << unmapped << " unmapped) in " << ms << "ms" << std::endl;
    return ErrorCode::SUCCESS;
}

size_t DomainMapper::mapped_count() const {
    size_t n = 0;
    for (const auto& r : results_) if (r.mapped) ++n;
    return n;
}

#include "domain_mapper.hpp"
#include "utils.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>

#ifdef USE_OPENMP
#include <omp.h>
#endif

DomainMapper::DomainMapper(const GTFParser& parser, const std::string& mode, const std::string& format) 
    : gtf_parser_(parser), mode_(mode), format_(format) {
    domains_.reserve(100000); // Reserve space for typical domain counts
    results_.reserve(500000); // Estimate 5 genomic intervals per domain
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
        std::string protein_id, domain_id;
        uint32_t start, end;
        
        if (iss >> protein_id >> start >> end) {
            // Optional domain_id
            iss >> domain_id;
            
            domains_.emplace_back(protein_id, start, end, domain_id);
        }
    }
    
    std::cerr << "Loaded " << domains_.size() << " domains from " << bed_filename << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode DomainMapper::process_domains() {
    results_.clear();
    results_.reserve(domains_.size() * 5); // Estimate 5 intervals per domain
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    size_t processed = 0;
    size_t unmapped = 0;
    
#ifdef USE_OPENMP
    // Parallel processing
    #pragma omp parallel for reduction(+:unmapped) schedule(dynamic, 100)
    for (size_t i = 0; i < domains_.size(); ++i) {
        auto genomic_domains = map_domain_to_genomic(domains_[i]);
        
        // Assign region index (1-based)
        for (auto& gd : genomic_domains) {
            gd.region_index = static_cast<uint32_t>(i + 1);
        }
        
        if (genomic_domains.empty()) {
            ++unmapped;
        } else {
            #pragma omp critical
            {
                results_.insert(results_.end(), genomic_domains.begin(), genomic_domains.end());
            }
        }
        
        #pragma omp atomic
        ++processed;
        
        if (processed % 10000 == 0) {
            #pragma omp critical
            {
                std::cerr << "Processed " << processed << "/" << domains_.size() 
                         << " domains...\r" << std::flush;
            }
        }
    }
#else
    // Sequential processing
    for (size_t i = 0; i < domains_.size(); ++i) {
        auto genomic_domains = map_domain_to_genomic(domains_[i]);
        
        // Assign region index (1-based)
        for (auto& gd : genomic_domains) {
            gd.region_index = static_cast<uint32_t>(i + 1);
        }
        
        if (genomic_domains.empty()) {
            ++unmapped;
        } else {
            results_.insert(results_.end(), genomic_domains.begin(), genomic_domains.end());
        }
        
        ++processed;
        if (processed % 10000 == 0) {
            std::cerr << "Processed " << processed << "/" << domains_.size() 
                     << " domains...\r" << std::flush;
        }
    }
#endif
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cerr << "\nCompleted domain mapping in " << duration.count() << "ms" << std::endl;
    std::cerr << "Successfully mapped: " << (domains_.size() - unmapped) << "/" << domains_.size() << std::endl;
    std::cerr << "Total genomic intervals: " << results_.size() << std::endl;
    
    if (unmapped > 0) {
        std::cerr << "Warning: " << unmapped << " domains could not be mapped" << std::endl;
    }
    
    return ErrorCode::SUCCESS;
}

std::vector<GenomicDomain> DomainMapper::map_domain_to_genomic(const ProteinDomain& domain) const {
    if (format_ == "detailed") {
        return map_domain_detailed(domain);
    }
    
    // Simple format logic
    std::vector<GenomicDomain> result;
    
    const auto* intervals = gtf_parser_.get_protein_intervals(domain.protein_id);
    if (!intervals || intervals->empty()) {
        return result; // Protein not found
    }
    
    // Find CDS intervals for coordinate mapping
    std::vector<GenomicInterval> cds_intervals;
    for (const auto& interval : *intervals) {
        if (interval.feature_type == FeatureType::CDS) {
            cds_intervals.push_back(interval);
        }
    }
    
    if (cds_intervals.empty()) {
        return result; // No CDS found for protein coordinate mapping
    }
    
    // CDS intervals should already be in correct order from index
    // For positive strand: genomic order (start to end)
    // For negative strand: translation order (downstream to upstream)
    
    // Calculate cumulative protein positions for CDS
    std::vector<uint32_t> cumulative_lengths;
    cumulative_lengths.reserve(cds_intervals.size() + 1);
    cumulative_lengths.push_back(0);
    
    uint32_t cumulative = 0;
    for (const auto& interval : cds_intervals) {
        uint32_t length = interval.end - interval.start + 1;
        cumulative += length;
        cumulative_lengths.push_back(cumulative);
    }
    
    uint32_t total_mapping_length = cumulative;
    
    // Convert domain coordinates to nucleotide positions
    uint32_t domain_start_0based = domain.start - 1;
    uint32_t domain_end_0based = domain.end - 1;
    uint32_t nt_start = domain_start_0based * 3;
    uint32_t nt_end = (domain_end_0based + 1) * 3 - 1;
    
    if (nt_end >= total_mapping_length) {
        return result; // Domain extends beyond CDS length
    }
    
    // Find genomic regions covered by the domain
    std::vector<std::pair<uint32_t, uint32_t>> domain_genomic_ranges;
    
    for (size_t i = 0; i < cds_intervals.size(); ++i) {
        uint32_t interval_start_nt = cumulative_lengths[i];
        uint32_t interval_end_nt = cumulative_lengths[i + 1] - 1;
        
        if (nt_start <= interval_end_nt && nt_end >= interval_start_nt) {
            uint32_t overlap_start_nt = std::max(nt_start, interval_start_nt);
            uint32_t overlap_end_nt = std::min(nt_end, interval_end_nt);
            
            uint32_t offset_start = overlap_start_nt - interval_start_nt;
            uint32_t offset_end = overlap_end_nt - interval_start_nt;
            
            uint32_t genomic_start, genomic_end;
            
            // For negative strand, calculate coordinates differently
            if (!cds_intervals.empty() && cds_intervals[0].strand == '-') {
                // For negative strand, translation starts from the 'end' of the CDS
                // and moves towards the 'start', so coordinates are reversed within each CDS
                genomic_start = cds_intervals[i].end - offset_end;
                genomic_end = cds_intervals[i].end - offset_start;
            } else {
                // For positive strand, normal calculation
                genomic_start = cds_intervals[i].start + offset_start;
                genomic_end = cds_intervals[i].start + offset_end;
            }
            
            domain_genomic_ranges.emplace_back(genomic_start, genomic_end);
        }
    }
    
    if (domain_genomic_ranges.empty()) {
        return result;
    }
    
    // Get overall genomic span
    uint32_t domain_genomic_start = domain_genomic_ranges.front().first;
    uint32_t domain_genomic_end = domain_genomic_ranges.back().second;
    
    for (const auto& range : domain_genomic_ranges) {
        domain_genomic_start = std::min(domain_genomic_start, range.first);
        domain_genomic_end = std::max(domain_genomic_end, range.second);
    }
    
    // Collect relevant intervals first
    std::vector<GenomicInterval> relevant_intervals;
    for (const auto& interval : *intervals) {
        // Skip features based on mode
        if (mode_ == "basic" && interval.feature_type != FeatureType::CDS) {
            continue;
        }
        if (mode_ == "full" && interval.feature_type != FeatureType::CDS && interval.feature_type != FeatureType::INTRON) {
            continue;
        }
        
        // Check overlap with domain's genomic span
        if (interval.start <= domain_genomic_end && interval.end >= domain_genomic_start) {
            relevant_intervals.push_back(interval);
        }
    }
    
    // Sort intervals by genomic position for consistent output
    std::sort(relevant_intervals.begin(), relevant_intervals.end());
    
    // Process intervals in sorted order
    for (const auto& interval : relevant_intervals) {
        GenomicDomain genomic_domain;
        genomic_domain.chromosome = interval.chromosome;
        
        // Check if this CDS overlaps with the domain precisely
        bool has_precise_overlap = false;
        uint32_t precise_start = interval.start;
        uint32_t precise_end = interval.end;
        
        // For CDS: check precise mapping
        if (interval.feature_type == FeatureType::CDS) {
            for (const auto& [range_start, range_end] : domain_genomic_ranges) {
                if (interval.start <= range_end && interval.end >= range_start) {
                    has_precise_overlap = true;
                    precise_start = std::max(interval.start, range_start);
                    precise_end = std::min(interval.end, range_end);
                    break;
                }
            }
            
            if (has_precise_overlap) {
                genomic_domain.start = precise_start;
                genomic_domain.end = precise_end;
                genomic_domain.domain_overlap = "Yes";
            } else {
                // This CDS is adjacent but not overlapping with domain
                genomic_domain.start = interval.start;
                genomic_domain.end = interval.end;
                genomic_domain.domain_overlap = "No";
            }
        } else {
            // For introns: include if they fall within domain span
            genomic_domain.start = std::max(interval.start, domain_genomic_start);
            genomic_domain.end = std::min(interval.end, domain_genomic_end);
        }
        
        genomic_domain.protein_id = domain.protein_id;
        genomic_domain.domain_id = domain.domain_id;
        genomic_domain.feature_type = feature_type_to_string(interval.feature_type);
        genomic_domain.strand = interval.strand;
        genomic_domain.protein_start = domain.start;
        genomic_domain.protein_end = domain.end;
        
        result.push_back(genomic_domain);
    }
    
    return result;
}

std::vector<GenomicDomain> DomainMapper::map_domain_detailed(const ProteinDomain& domain) const {
    std::vector<GenomicDomain> result;
    
    const auto* intervals = gtf_parser_.get_protein_intervals(domain.protein_id);
    if (!intervals || intervals->empty()) {
        return result; // Protein not found
    }
    
    // Find CDS intervals in translation order (already correct from GTF parser)
    std::vector<GenomicInterval> translation_order_cds;
    for (const auto& interval : *intervals) {
        if (interval.feature_type == FeatureType::CDS) {
            translation_order_cds.push_back(interval);
        }
    }
    
    if (translation_order_cds.empty()) {
        return result; // No CDS found
    }
    
    // Calculate cumulative protein coordinates for translation-ordered CDS
    std::vector<uint32_t> cumulative_lengths;
    cumulative_lengths.reserve(translation_order_cds.size() + 1);
    cumulative_lengths.push_back(0);
    
    uint32_t cumulative = 0;
    for (const auto& interval : translation_order_cds) {
        uint32_t length = interval.end - interval.start + 1;
        cumulative += length;
        cumulative_lengths.push_back(cumulative);
    }
    
    // Domain boundaries in nucleotides
    uint32_t domain_nt_start = (domain.start - 1) * 3;
    uint32_t domain_nt_end = domain.end * 3 - 1;
    
    // Find which genomic regions the domain covers
    std::vector<std::pair<uint32_t, uint32_t>> domain_genomic_ranges;
    
    for (size_t i = 0; i < translation_order_cds.size(); ++i) {
        uint32_t cds_start_nt = cumulative_lengths[i];
        uint32_t cds_end_nt = cumulative_lengths[i + 1] - 1;
        
        // Check if domain overlaps with this CDS
        if (domain_nt_start <= cds_end_nt && domain_nt_end >= cds_start_nt) {
            uint32_t overlap_start_nt = std::max(domain_nt_start, cds_start_nt);
            uint32_t overlap_end_nt = std::min(domain_nt_end, cds_end_nt);
            
            uint32_t offset_start = overlap_start_nt - cds_start_nt;
            uint32_t offset_end = overlap_end_nt - cds_start_nt;
            
            const auto& cds = translation_order_cds[i];
            uint32_t genomic_start, genomic_end;
            
            // Handle strand-specific coordinate calculation
            if (cds.strand == '-') {
                // For negative strand, coordinates are reversed within each CDS
                genomic_start = cds.end - offset_end;
                genomic_end = cds.end - offset_start;
            } else {
                // For positive strand, normal calculation
                genomic_start = cds.start + offset_start;
                genomic_end = cds.start + offset_end;
            }
            
            domain_genomic_ranges.emplace_back(genomic_start, genomic_end);
        }
    }
    
    if (domain_genomic_ranges.empty()) {
        return result;
    }
    
    // Calculate overall genomic span
    uint32_t domain_genomic_start = domain_genomic_ranges[0].first;
    uint32_t domain_genomic_end = domain_genomic_ranges[0].second;
    
    for (const auto& range : domain_genomic_ranges) {
        domain_genomic_start = std::min(domain_genomic_start, range.first);
        domain_genomic_end = std::max(domain_genomic_end, range.second);
    }
    
    // Now process all intervals in GENOMIC order for output
    std::vector<GenomicInterval> output_intervals;
    
    // Add all relevant intervals
    for (const auto& interval : *intervals) {
        if (mode_ == "basic" && interval.feature_type != FeatureType::CDS) {
            continue;
        }
        if (mode_ == "full" && interval.feature_type != FeatureType::CDS && interval.feature_type != FeatureType::INTRON) {
            continue;
        }
        
        // Include if overlaps with domain span
        if (interval.start <= domain_genomic_end && interval.end >= domain_genomic_start) {
            output_intervals.push_back(interval);
        }
    }
    
    // Sort by genomic position for output
    std::sort(output_intervals.begin(), output_intervals.end());
    
    // Count total CDS and introns for proper numbering in negative strand
    int total_cds_count = 0;
    int total_intron_count = 0;
    for (const auto& interval : output_intervals) {
        if (interval.feature_type == FeatureType::CDS) {
            total_cds_count++;
        } else if (interval.feature_type == FeatureType::INTRON) {
            total_intron_count++;
        }
    }
    
    int region_counter = 1;
    int cds_counter = 1;
    int intron_counter = 1;
    
    // For negative strand, we need to reverse the feature numbering
    bool is_negative_strand = !output_intervals.empty() && output_intervals[0].strand == '-';
    
    for (const auto& interval : output_intervals) {
        // Calculate the correct feature counter for this interval
        int current_cds_id = cds_counter;
        int current_intron_id = intron_counter;
        
        if (is_negative_strand) {
            if (interval.feature_type == FeatureType::CDS) {
                current_cds_id = total_cds_count - cds_counter + 1;
            } else if (interval.feature_type == FeatureType::INTRON) {
                current_intron_id = total_intron_count - intron_counter + 1;
            }
        }
        
        if (interval.feature_type == FeatureType::CDS) {
            // For CDS, we need to split into overlapping and non-overlapping parts
            std::vector<std::pair<uint32_t, uint32_t>> cds_parts;
            
            // Find overlapping regions with the domain
            std::vector<std::pair<uint32_t, uint32_t>> overlapping_parts;
            for (const auto& [range_start, range_end] : domain_genomic_ranges) {
                if (interval.start <= range_end && interval.end >= range_start) {
                    uint32_t overlap_start = std::max(interval.start, range_start);
                    uint32_t overlap_end = std::min(interval.end, range_end);
                    overlapping_parts.emplace_back(overlap_start, overlap_end);
                }
            }
            
            // Merge overlapping parts if they touch or overlap
            if (!overlapping_parts.empty()) {
                std::sort(overlapping_parts.begin(), overlapping_parts.end());
                std::vector<std::pair<uint32_t, uint32_t>> merged_overlaps;
                merged_overlaps.push_back(overlapping_parts[0]);
                
                for (size_t i = 1; i < overlapping_parts.size(); ++i) {
                    if (overlapping_parts[i].first <= merged_overlaps.back().second + 1) {
                        // Merge with previous
                        merged_overlaps.back().second = std::max(merged_overlaps.back().second, overlapping_parts[i].second);
                    } else {
                        merged_overlaps.push_back(overlapping_parts[i]);
                    }
                }
                
                // Now create CDS parts: non-overlapping, overlapping, non-overlapping
                uint32_t current_pos = interval.start;
                
                for (const auto& [overlap_start, overlap_end] : merged_overlaps) {
                    // Add non-overlapping part before this overlap (if any)
                    if (current_pos < overlap_start) {
                        GenomicDomain non_overlap;
                        non_overlap.chromosome = interval.chromosome;
                        non_overlap.start = current_pos;
                        non_overlap.end = overlap_start - 1;
                        non_overlap.protein_id = domain.protein_id;
                        non_overlap.domain_id = domain.domain_id;
                        non_overlap.feature_type = "CDS";
                        non_overlap.strand = interval.strand;
                        non_overlap.protein_start = domain.start;
                        non_overlap.protein_end = domain.end;
                        non_overlap.region_index = region_counter++;
                        non_overlap.domain_overlap = "No";
                        non_overlap.feature_id = "CDS_" + std::to_string(current_cds_id);
                        result.push_back(non_overlap);
                    }
                    
                    // Add overlapping part
                    GenomicDomain overlap;
                    overlap.chromosome = interval.chromosome;
                    overlap.start = overlap_start;
                    overlap.end = overlap_end;
                    overlap.protein_id = domain.protein_id;
                    overlap.domain_id = domain.domain_id;
                    overlap.feature_type = "CDS";
                    overlap.strand = interval.strand;
                    overlap.protein_start = domain.start;
                    overlap.protein_end = domain.end;
                    overlap.region_index = region_counter++;
                    overlap.domain_overlap = "Yes";
                    overlap.feature_id = "CDS_" + std::to_string(current_cds_id);
                    result.push_back(overlap);
                    
                    current_pos = overlap_end + 1;
                }
                
                // Add final non-overlapping part (if any)
                if (current_pos <= interval.end) {
                    GenomicDomain non_overlap;
                    non_overlap.chromosome = interval.chromosome;
                    non_overlap.start = current_pos;
                    non_overlap.end = interval.end;
                    non_overlap.protein_id = domain.protein_id;
                    non_overlap.domain_id = domain.domain_id;
                    non_overlap.feature_type = "CDS";
                    non_overlap.strand = interval.strand;
                    non_overlap.protein_start = domain.start;
                    non_overlap.protein_end = domain.end;
                    non_overlap.region_index = region_counter++;
                    non_overlap.domain_overlap = "No";
                    non_overlap.feature_id = "CDS_" + std::to_string(current_cds_id);
                    result.push_back(non_overlap);
                }
            } else {
                // No overlap - entire CDS is non-overlapping
                GenomicDomain non_overlap;
                non_overlap.chromosome = interval.chromosome;
                non_overlap.start = interval.start;
                non_overlap.end = interval.end;
                non_overlap.protein_id = domain.protein_id;
                non_overlap.domain_id = domain.domain_id;
                non_overlap.feature_type = "CDS";
                non_overlap.strand = interval.strand;
                non_overlap.protein_start = domain.start;
                non_overlap.protein_end = domain.end;
                non_overlap.region_index = region_counter++;
                non_overlap.domain_overlap = "No";
                non_overlap.feature_id = "CDS_" + std::to_string(current_cds_id);
                result.push_back(non_overlap);
            }
            
            // Increment CDS counter after processing this CDS
            cds_counter++;
        } else {
            // For introns, handle as before
            GenomicDomain genomic_domain;
            genomic_domain.chromosome = interval.chromosome;
            genomic_domain.protein_id = domain.protein_id;
            genomic_domain.domain_id = domain.domain_id;
            genomic_domain.feature_type = feature_type_to_string(interval.feature_type);
            genomic_domain.strand = interval.strand;
            genomic_domain.protein_start = domain.start;
            genomic_domain.protein_end = domain.end;
            genomic_domain.region_index = region_counter++;
            genomic_domain.feature_id = "intron_" + std::to_string(current_intron_id);
            
            // Check if intron overlaps with domain span
            if (interval.start <= domain_genomic_end && interval.end >= domain_genomic_start) {
                genomic_domain.start = std::max(interval.start, domain_genomic_start);
                genomic_domain.end = std::min(interval.end, domain_genomic_end);
                genomic_domain.domain_overlap = "Yes";
            } else {
                genomic_domain.start = interval.start;
                genomic_domain.end = interval.end;
                genomic_domain.domain_overlap = "No";
            }
            
            result.push_back(genomic_domain);
            
            // Increment intron counter after processing this intron
            intron_counter++;
        }
    }
    
    return result;
}

ErrorCode DomainMapper::write_results(const std::string& output_filename, const std::string& format) const {
    std::ofstream file(output_filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot create output file: " << output_filename << std::endl;
        return ErrorCode::FILE_NOT_FOUND;
    }
    
    if (format == "simple") {
        // Simple BED format output with feature_type and region_index as additional columns
        for (const auto& domain : results_) {
            file << domain.chromosome << "\t"
                 << domain.start << "\t"
                 << domain.end << "\t"
                 << domain.protein_id;
            
            if (!domain.domain_id.empty()) {
                file << "_" << domain.domain_id;
            }
            
            file << "_" << domain.protein_start << "-" << domain.protein_end
                 << "\t0\t" << domain.strand << "\t"
                 << domain.feature_type << "\t"
                 << domain.region_index << "\n";
        }
    } else if (format == "detailed") {
        // Detailed BED format output (no header) with overlap information and region_index
        for (const auto& domain : results_) {
            file << domain.chromosome << "\t"
                 << domain.start << "\t"
                 << domain.end << "\t"
                 << domain.protein_id;
            
            if (!domain.domain_id.empty()) {
                file << "_" << domain.domain_id;
            }
            
            file << "_" << domain.protein_start << "-" << domain.protein_end
                 << "\t0\t" << domain.strand << "\t"
                 << domain.feature_type << "\t"
                 << domain.region_index << "\t"
                 << domain.domain_overlap << "\t"
                 << domain.feature_id << "\n";
        }
    }
    
    std::cerr << "Results written to: " << output_filename << std::endl;
    return ErrorCode::SUCCESS;
}

size_t DomainMapper::get_unmapped_count() const {
    return domains_.size() - get_mapped_count();
}

size_t DomainMapper::get_mapped_count() const {
    std::set<std::string> mapped_proteins;
    for (const auto& result : results_) {
        mapped_proteins.insert(result.protein_id);
    }
    return mapped_proteins.size();
}

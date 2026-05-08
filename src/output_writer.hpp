#ifndef OUTPUT_WRITER_HPP
#define OUTPUT_WRITER_HPP

#include "domain_mapper.hpp"

namespace output {

// Writes the requested set of files into out_dir based on output_kind.
// Always writes domain_mapping_summary.tsv and unmapped_domains.tsv (the latter
// only if there are unmapped rows). For ALL, also writes run_metadata.json.
ErrorCode write_all(const std::string& out_dir,
                    OutputKind kind,
                    const std::vector<DomainResult>& results,
                    const std::string& gtf_or_index_path,
                    const std::vector<std::string>& cli_args);

} // namespace output

#endif // OUTPUT_WRITER_HPP

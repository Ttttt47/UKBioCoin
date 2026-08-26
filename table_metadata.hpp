#ifndef UKBIOCOIN_TABLE_METADATA_HPP
#define UKBIOCOIN_TABLE_METADATA_HPP

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace ukc {

class TableError : public std::runtime_error {
public:
    explicit TableError(const std::string &message) : std::runtime_error(message) {}
};

struct MetadataProjection {
    std::vector<std::string> input_names;
    std::vector<std::string> output_names;
    std::vector<std::size_t> output_indices;
    bool canonical = false;
};

std::vector<std::string> split_table_fields(const std::string &line,
                                            const std::string &path,
                                            std::size_t line_number);

MetadataProjection make_metadata_projection(
    const std::vector<std::string> &input_names,
    const std::string &path);

MetadataProjection parse_metadata_header(const std::string &line,
                                         const std::string &path);

std::string metadata_header_tsv(const MetadataProjection &projection);

std::string project_metadata_row(const std::string &line,
                                 const MetadataProjection &projection,
                                 const std::string &path,
                                 std::size_t line_number);

std::string join_table_fields(const std::vector<std::string> &fields,
                              const std::vector<std::size_t> &indices);

}  // namespace ukc

#endif

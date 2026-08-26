#include "table_metadata.hpp"

#include <cctype>
#include <unordered_map>
#include <unordered_set>

namespace ukc {
namespace {

struct Token {
    const char *begin = nullptr;
    const char *end = nullptr;
};

bool next_token(const char *&cursor, const char *end, Token &token) {
    while (cursor < end && std::isspace(static_cast<unsigned char>(*cursor))) {
        ++cursor;
    }
    if (cursor == end) return false;
    if (*cursor == '"') {
        token.begin = ++cursor;
        while (cursor < end && *cursor != '"') ++cursor;
        if (cursor == end) throw TableError("unterminated quoted field");
        token.end = cursor++;
        if (cursor < end && !std::isspace(static_cast<unsigned char>(*cursor))) {
            throw TableError("unexpected character after quoted field");
        }
        return true;
    }
    token.begin = cursor;
    while (cursor < end && !std::isspace(static_cast<unsigned char>(*cursor))) {
        ++cursor;
    }
    token.end = cursor;
    return true;
}

std::string canonical_name(const std::string &name) {
    if (name == "#CHROM" || name == "CHROM") return "#CHROM";
    if (name == "POS") return "POS";
    if (name == "ID") return "ID";
    if (name == "REF_Allele" || name == "REF") return "REF";
    if (name == "ALT_Allele" || name == "ALT1") return "ALT1";
    if (name == "Effect_Allele" || name == "A1") return "A1";
    if (name == "REF_FREQ") return "REF_FREQ";
    return "";
}

const std::vector<std::string> &canonical_order() {
    static const std::vector<std::string> names{
        "#CHROM", "POS", "ID", "REF", "ALT1", "A1", "REF_FREQ"};
    return names;
}

}  // namespace

std::vector<std::string> split_table_fields(const std::string &line,
                                            const std::string &path,
                                            std::size_t line_number) {
    std::vector<std::string> fields;
    if (line.find('\t') != std::string::npos) {
        std::size_t start = 0;
        while (start <= line.size()) {
            const std::size_t tab = line.find('\t', start);
            const std::size_t end = tab == std::string::npos ? line.size() : tab;
            fields.push_back(line.substr(start, end - start));
            if (tab == std::string::npos) break;
            start = tab + 1;
        }
        if (!fields.empty() && !fields.back().empty() &&
            fields.back().back() == '\r') {
            fields.back().pop_back();
        }
        return fields;
    }
    const char *cursor = line.data();
    const char *end = cursor + line.size();
    Token token;
    try {
        while (next_token(cursor, end, token)) {
            fields.emplace_back(token.begin, token.end);
        }
    } catch (const TableError &error) {
        throw TableError(path + ":" + std::to_string(line_number) + ": " +
                         error.what());
    }
    if (fields.empty()) {
        throw TableError(path + ":" + std::to_string(line_number) +
                         ": empty table row");
    }
    return fields;
}

MetadataProjection make_metadata_projection(
    const std::vector<std::string> &input_names,
    const std::string &path) {
    if (input_names.empty()) throw TableError(path + ":1: empty metadata header");
    MetadataProjection projection;
    projection.input_names = input_names;
    std::unordered_set<std::string> exact_names;
    std::unordered_map<std::string, std::size_t> recognized;
    for (std::size_t i = 0; i < input_names.size(); ++i) {
        const std::string &name = input_names[i];
        if (name.empty() || !exact_names.insert(name).second) {
            throw TableError(path + ":1: duplicate metadata column '" + name + "'");
        }
        const std::string canonical = canonical_name(name);
        if (!canonical.empty()) {
            const auto inserted = recognized.emplace(canonical, i);
            if (!inserted.second) {
                throw TableError(path + ":1: metadata aliases '" +
                                 input_names[inserted.first->second] + "' and '" +
                                 name + "' both map to '" + canonical + "'");
            }
        }
    }

    projection.canonical = recognized.size() == canonical_order().size();
    if (projection.canonical) {
        for (const std::string &name : canonical_order()) {
            projection.output_names.push_back(name);
            projection.output_indices.push_back(recognized.at(name));
        }
        for (std::size_t i = 0; i < input_names.size(); ++i) {
            if (canonical_name(input_names[i]).empty()) {
                projection.output_names.push_back(input_names[i]);
                projection.output_indices.push_back(i);
            }
        }
    } else {
        projection.output_names = input_names;
        if (projection.output_names[0].empty() ||
            projection.output_names[0][0] != '#') {
            projection.output_names[0].insert(projection.output_names[0].begin(), '#');
        }
        for (std::size_t i = 0; i < input_names.size(); ++i) {
            projection.output_indices.push_back(i);
        }
    }
    return projection;
}

MetadataProjection parse_metadata_header(const std::string &line,
                                         const std::string &path) {
    return make_metadata_projection(split_table_fields(line, path, 1), path);
}

std::string join_table_fields(const std::vector<std::string> &fields,
                              const std::vector<std::size_t> &indices) {
    std::string result;
    std::size_t reserve = 0;
    for (std::size_t index : indices) reserve += fields[index].size() + 1;
    result.reserve(reserve);
    for (std::size_t i = 0; i < indices.size(); ++i) {
        if (i != 0) result.push_back('\t');
        result.append(fields[indices[i]]);
    }
    return result;
}

std::string metadata_header_tsv(const MetadataProjection &projection) {
    std::vector<std::size_t> indices(projection.output_names.size());
    for (std::size_t i = 0; i < indices.size(); ++i) indices[i] = i;
    return join_table_fields(projection.output_names, indices);
}

std::string project_metadata_row(const std::string &line,
                                 const MetadataProjection &projection,
                                 const std::string &path,
                                 std::size_t line_number) {
    const std::vector<std::string> fields =
        split_table_fields(line, path, line_number);
    if (fields.size() != projection.input_names.size()) {
        throw TableError(path + ":" + std::to_string(line_number) +
                         ": metadata field count " +
                         std::to_string(fields.size()) +
                         " does not match header field count " +
                         std::to_string(projection.input_names.size()));
    }
    return join_table_fields(fields, projection.output_indices);
}

}  // namespace ukc

#include "nss_v2.hpp"

#include <algorithm>
#include <array>
#include <boost/filesystem.hpp>
#include <boost/json.hpp>
#include <cctype>
#include <cstring>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <sys/resource.h>

namespace fs = boost::filesystem;
namespace json = boost::json;

namespace ukc {
namespace {

bool host_is_little_endian() {
    const std::uint16_t value = 1;
    return *reinterpret_cast<const unsigned char *>(&value) == 1;
}

template <typename T>
void reverse_bytes(T &value) {
    unsigned char *begin = reinterpret_cast<unsigned char *>(&value);
    std::reverse(begin, begin + sizeof(T));
}

std::uint64_t checked_product(const std::vector<std::uint64_t> &shape,
                              const std::string &context) {
    if (shape.empty()) {
        throw NssError(context + ": scalar arrays are not supported");
    }
    std::uint64_t product = 1;
    for (std::uint64_t dimension : shape) {
        if (dimension == 0) {
            throw NssError(context + ": dimensions must be positive");
        }
        if (product > std::numeric_limits<std::uint64_t>::max() / dimension) {
            throw NssError(context + ": array shape overflows uint64");
        }
        product *= dimension;
    }
    return product;
}

std::size_t dtype_size(const std::string &dtype, const std::string &context) {
    if (dtype == "<f4" || dtype == ">f4" || dtype == "=f4") {
        return 4;
    }
    if (dtype == "<f8" || dtype == ">f8" || dtype == "=f8") {
        return 8;
    }
    throw NssError(context + ": unsupported dtype '" + dtype + "'");
}

std::string read_file(const std::string &path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw NssError("Error opening " + path + " for input");
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    if (!input.eof() && input.fail()) {
        throw NssError("Error reading " + path);
    }
    return contents.str();
}

const json::object &require_object(const json::value &value,
                                   const std::string &context) {
    if (!value.is_object()) {
        throw NssError(context + " must be an object");
    }
    return value.as_object();
}

const json::value &require_member(const json::object &object,
                                  const char *key,
                                  const std::string &context) {
    const json::value *value = object.if_contains(key);
    if (value == nullptr) {
        throw NssError(context + ": missing required field '" + key + "'");
    }
    return *value;
}

std::string require_string(const json::object &object,
                           const char *key,
                           const std::string &context) {
    const json::value &value = require_member(object, key, context);
    if (!value.is_string()) {
        throw NssError(context + "." + key + " must be a string");
    }
    return std::string(value.as_string().c_str());
}

std::uint64_t require_positive_integer(const json::object &object,
                                       const char *key,
                                       const std::string &context) {
    const json::value &value = require_member(object, key, context);
    std::uint64_t result = 0;
    if (value.is_uint64()) {
        result = value.as_uint64();
    } else if (value.is_int64() && value.as_int64() > 0) {
        result = static_cast<std::uint64_t>(value.as_int64());
    } else {
        throw NssError(context + "." + key + " must be a positive integer");
    }
    if (result == 0) {
        throw NssError(context + "." + key + " must be a positive integer");
    }
    return result;
}

void validate_relative_path(const std::string &path, const std::string &context) {
    const fs::path parsed(path);
    if (path.empty() || parsed.is_absolute()) {
        throw NssError(context + " must be a non-empty relative path");
    }
    for (const auto &part : parsed) {
        if (part == "..") {
            throw NssError(context + " must not contain '..'");
        }
    }
}

std::vector<std::uint64_t> parse_shape(const json::object &object,
                                       const std::string &context) {
    const json::value &value = require_member(object, "shape", context);
    if (!value.is_array() || value.as_array().empty()) {
        throw NssError(context + ".shape must be a non-empty integer array");
    }
    std::vector<std::uint64_t> shape;
    for (const auto &entry : value.as_array()) {
        std::uint64_t dimension = 0;
        if (entry.is_uint64()) {
            dimension = entry.as_uint64();
        } else if (entry.is_int64() && entry.as_int64() > 0) {
            dimension = static_cast<std::uint64_t>(entry.as_int64());
        }
        if (dimension == 0) {
            throw NssError(context + ".shape values must be positive integers");
        }
        shape.push_back(dimension);
    }
    checked_product(shape, context);
    return shape;
}

ArrayDescriptor parse_array(const json::value &value, const std::string &context) {
    const json::object &object = require_object(value, context);
    ArrayDescriptor descriptor;
    descriptor.path = require_string(object, "path", context);
    descriptor.dtype = require_string(object, "dtype", context);
    descriptor.shape = parse_shape(object, context);
    descriptor.order = require_string(object, "order", context);
    validate_relative_path(descriptor.path, context + ".path");
    if (descriptor.order != "C") {
        throw NssError(context + ".order must be 'C'");
    }
    dtype_size(descriptor.dtype, context);
    return descriptor;
}

json::object array_json(const ArrayDescriptor &array) {
    json::array shape;
    for (std::uint64_t dimension : array.shape) {
        shape.push_back(dimension);
    }
    return {{"path", array.path}, {"dtype", array.dtype},
            {"shape", std::move(shape)}, {"order", array.order}};
}

NpyHeader parse_npy_header(std::ifstream &input, const std::string &path) {
    std::array<unsigned char, 8> prefix{};
    input.read(reinterpret_cast<char *>(prefix.data()), prefix.size());
    const std::array<unsigned char, 6> magic{{0x93, 'N', 'U', 'M', 'P', 'Y'}};
    if (!input || !std::equal(magic.begin(), magic.end(), prefix.begin())) {
        throw NssError(path + ": invalid NPY magic/header");
    }
    const unsigned major = prefix[6];
    const unsigned minor = prefix[7];
    (void) minor;
    std::uint64_t header_length = 0;
    if (major == 1) {
        unsigned char bytes[2]{};
        input.read(reinterpret_cast<char *>(bytes), 2);
        header_length = static_cast<std::uint64_t>(bytes[0]) |
                        (static_cast<std::uint64_t>(bytes[1]) << 8U);
    } else if (major == 2 || major == 3) {
        unsigned char bytes[4]{};
        input.read(reinterpret_cast<char *>(bytes), 4);
        header_length = static_cast<std::uint64_t>(bytes[0]) |
                        (static_cast<std::uint64_t>(bytes[1]) << 8U) |
                        (static_cast<std::uint64_t>(bytes[2]) << 16U) |
                        (static_cast<std::uint64_t>(bytes[3]) << 24U);
    } else {
        throw NssError(path + ": unsupported NPY version " + std::to_string(major));
    }
    if (!input || header_length == 0 || header_length > 1024 * 1024) {
        throw NssError(path + ": invalid NPY header length");
    }
    std::string header(static_cast<std::size_t>(header_length), '\0');
    input.read(&header[0], static_cast<std::streamsize>(header.size()));
    if (!input) {
        throw NssError(path + ": truncated NPY header");
    }

    const auto extract_quoted = [&](const std::string &key) {
        const std::size_t key_pos = header.find(key);
        if (key_pos == std::string::npos) {
            throw NssError(path + ": NPY header lacks " + key);
        }
        const std::size_t colon = header.find(':', key_pos + key.size());
        const std::size_t quote = header.find_first_of("'\"", colon);
        if (colon == std::string::npos || quote == std::string::npos) {
            throw NssError(path + ": malformed NPY " + key);
        }
        const std::size_t end_quote = header.find(header[quote], quote + 1);
        if (end_quote == std::string::npos) {
            throw NssError(path + ": malformed NPY " + key);
        }
        return header.substr(quote + 1, end_quote - quote - 1);
    };

    NpyHeader parsed;
    parsed.dtype = extract_quoted("descr");
    parsed.item_size = dtype_size(parsed.dtype, path);
    const std::size_t fortran_pos = header.find("fortran_order");
    const std::size_t fortran_colon = header.find(':', fortran_pos);
    if (fortran_pos == std::string::npos || fortran_colon == std::string::npos) {
        throw NssError(path + ": malformed NPY fortran_order");
    }
    std::size_t bool_start = fortran_colon + 1;
    while (bool_start < header.size() && std::isspace(static_cast<unsigned char>(header[bool_start]))) {
        ++bool_start;
    }
    if (header.compare(bool_start, 4, "True") == 0) {
        parsed.fortran_order = true;
    } else if (header.compare(bool_start, 5, "False") != 0) {
        throw NssError(path + ": malformed NPY fortran_order");
    }

    const std::size_t shape_pos = header.find("shape");
    const std::size_t left = header.find('(', shape_pos);
    const std::size_t right = header.find(')', left);
    if (shape_pos == std::string::npos || left == std::string::npos || right == std::string::npos) {
        throw NssError(path + ": malformed NPY shape");
    }
    std::string shape_text = header.substr(left + 1, right - left - 1);
    std::size_t cursor = 0;
    while (cursor < shape_text.size()) {
        while (cursor < shape_text.size() &&
               (std::isspace(static_cast<unsigned char>(shape_text[cursor])) || shape_text[cursor] == ',')) {
            ++cursor;
        }
        if (cursor == shape_text.size()) {
            break;
        }
        std::size_t end = cursor;
        while (end < shape_text.size() && std::isdigit(static_cast<unsigned char>(shape_text[end]))) {
            ++end;
        }
        if (end == cursor) {
            throw NssError(path + ": malformed NPY shape");
        }
        const std::uint64_t value = std::stoull(shape_text.substr(cursor, end - cursor));
        if (value == 0) {
            throw NssError(path + ": NPY dimensions must be positive");
        }
        parsed.shape.push_back(value);
        cursor = end;
    }
    parsed.element_count = checked_product(parsed.shape, path);
    parsed.data_offset = input.tellg();
    return parsed;
}

std::string npy_header(const std::string &dtype,
                       const std::vector<std::uint64_t> &shape) {
    std::ostringstream shape_stream;
    shape_stream << '(';
    for (std::size_t i = 0; i < shape.size(); ++i) {
        if (i != 0) {
            shape_stream << ", ";
        }
        shape_stream << shape[i];
    }
    if (shape.size() == 1) {
        shape_stream << ',';
    }
    shape_stream << ')';
    std::string header = "{'descr': '" + dtype +
        "', 'fortran_order': False, 'shape': " + shape_stream.str() + ", }";
    const std::size_t base = 10;
    const std::size_t padding = (16 - ((base + header.size() + 1) % 16)) % 16;
    header.append(padding, ' ');
    header.push_back('\n');
    if (header.size() > 65535) {
        throw NssError("NPY v1 header is too large");
    }
    return header;
}

void write_npy_prefix(std::ofstream &output,
                      const std::string &dtype,
                      const std::vector<std::uint64_t> &shape,
                      const std::string &path) {
    const std::string header = npy_header(dtype, shape);
    const unsigned char prefix[8] = {0x93, 'N', 'U', 'M', 'P', 'Y', 1, 0};
    output.write(reinterpret_cast<const char *>(prefix), sizeof(prefix));
    const std::uint16_t length = static_cast<std::uint16_t>(header.size());
    const unsigned char length_bytes[2] = {
        static_cast<unsigned char>(length & 0xffU),
        static_cast<unsigned char>((length >> 8U) & 0xffU)};
    output.write(reinterpret_cast<const char *>(length_bytes), sizeof(length_bytes));
    output.write(header.data(), static_cast<std::streamsize>(header.size()));
    if (!output) {
        throw NssError("Error writing NPY header to " + path);
    }
}

template <typename T>
void write_little_endian(std::ofstream &output,
                         const std::vector<T> &values,
                         const std::string &path) {
    if (host_is_little_endian()) {
        output.write(reinterpret_cast<const char *>(values.data()),
                     static_cast<std::streamsize>(values.size() * sizeof(T)));
    } else {
        for (T value : values) {
            reverse_bytes(value);
            output.write(reinterpret_cast<const char *>(&value), sizeof(T));
        }
    }
    if (!output) {
        throw NssError("Error writing NPY payload to " + path);
    }
}

void validate_descriptor(const std::string &directory,
                         const ArrayDescriptor &descriptor,
                         const std::string &dtype,
                         const std::vector<std::uint64_t> &shape,
                         const std::string &context) {
    if (descriptor.dtype != dtype) {
        throw NssError(context + ": manifest dtype must be " + dtype);
    }
    if (descriptor.shape != shape) {
        throw NssError(context + ": manifest shape does not match dimensions");
    }
    NpyReader reader(join_path(directory, descriptor.path), dtype, shape);
    (void) reader;
}

}  // namespace

bool file_exists(const std::string &path) {
    return fs::is_regular_file(fs::path(path));
}

bool directory_exists(const std::string &path) {
    return fs::is_directory(fs::path(path));
}

std::string join_path(const std::string &base, const std::string &relative) {
    return (fs::path(base) / fs::path(relative)).string();
}

Manifest read_manifest(const std::string &directory) {
    const std::string path = join_path(directory, "manifest.json");
    boost::system::error_code error;
    const json::value document = json::parse(read_file(path), error);
    if (error) {
        throw NssError(path + ": invalid JSON: " + error.message());
    }
    const json::object &root = require_object(document, path);
    if (require_string(root, "schema", path) != "org.ukbiocoin.nss") {
        throw NssError(path + ": unsupported schema");
    }
    const std::uint64_t version = require_positive_integer(root, "schema_version", path);
    if (version != 2) {
        throw NssError(path + ": unsupported schema_version " + std::to_string(version));
    }
    if (require_string(root, "ukc_version", path) != "2.0") {
        throw NssError(path + ": ukc_version must be '2.0'");
    }
    if (require_string(root, "byte_order", path) != "little") {
        throw NssError(path + ": byte_order must be 'little'");
    }

    Manifest manifest;
    const std::string dimensions_context = path + ".dimensions";
    const json::object &dimensions = require_object(
        require_member(root, "dimensions", path), dimensions_context);
    manifest.sample_count = require_positive_integer(dimensions, "sample_count", dimensions_context);
    manifest.variant_count = require_positive_integer(dimensions, "variant_count", dimensions_context);
    manifest.phenotype_count = require_positive_integer(dimensions, "phenotype_count", dimensions_context);

    const std::string metadata_context = path + ".metadata";
    const json::object &metadata = require_object(
        require_member(root, "metadata", path), metadata_context);
    manifest.meta_path = require_string(metadata, "path", metadata_context);
    validate_relative_path(manifest.meta_path, metadata_context + ".path");
    const std::uint64_t metadata_rows = require_positive_integer(metadata, "rows", metadata_context);
    if (metadata_rows != manifest.variant_count) {
        throw NssError(path + ": metadata.rows does not match dimensions.variant_count");
    }
    if (require_string(metadata, "format", metadata_context) != "tsv") {
        throw NssError(path + ": metadata.format must be 'tsv'");
    }

    const std::string arrays_context = path + ".arrays";
    const json::object &arrays = require_object(
        require_member(root, "arrays", path), arrays_context);
    manifest.cov_yy = parse_array(require_member(arrays, "cov_yy", arrays_context), arrays_context + ".cov_yy");
    manifest.var_x = parse_array(require_member(arrays, "var_x", arrays_context), arrays_context + ".var_x");
    const json::value &x_missing = require_member(arrays, "x_missing", arrays_context);
    if (!x_missing.is_null()) {
        manifest.has_x_missing = true;
        manifest.x_missing = parse_array(x_missing, path + ".arrays.x_missing");
    }
    const json::value &y_missing = require_member(arrays, "y_missing", arrays_context);
    if (!y_missing.is_null()) {
        manifest.has_y_missing = true;
        manifest.y_missing = parse_array(y_missing, path + ".arrays.y_missing");
    }

    const json::value &phenotypes_value = require_member(root, "phenotypes", path);
    if (!phenotypes_value.is_array() || phenotypes_value.as_array().size() != manifest.phenotype_count) {
        throw NssError(path + ": phenotypes length does not match dimensions.phenotype_count");
    }
    std::set<std::string> names;
    std::set<std::string> paths;
    for (std::size_t i = 0; i < phenotypes_value.as_array().size(); ++i) {
        const std::string context = path + ".phenotypes[" + std::to_string(i) + "]";
        const json::object &object = require_object(phenotypes_value.as_array()[i], context);
        const json::value &index_value = require_member(object, "index", context);
        std::uint64_t index = std::numeric_limits<std::uint64_t>::max();
        if (index_value.is_uint64()) index = index_value.as_uint64();
        if (index_value.is_int64() && index_value.as_int64() >= 0) index = static_cast<std::uint64_t>(index_value.as_int64());
        if (index != i) {
            throw NssError(context + ".index must equal its zero-based array position");
        }
        PhenotypeDescriptor phenotype;
        phenotype.index = i;
        phenotype.name = require_string(object, "name", context);
        if (phenotype.name.empty() || !names.insert(phenotype.name).second) {
            throw NssError(context + ".name must be non-empty and unique");
        }
        phenotype.cov_xy = parse_array(require_member(object, "cov_xy", context), context + ".cov_xy");
        if (!paths.insert(phenotype.cov_xy.path).second) {
            throw NssError(context + ".cov_xy.path must be unique");
        }
        manifest.phenotypes.push_back(std::move(phenotype));
    }
    return manifest;
}

void write_manifest(const std::string &directory, const Manifest &manifest) {
    json::object dimensions{{"sample_count", manifest.sample_count},
                            {"variant_count", manifest.variant_count},
                            {"phenotype_count", manifest.phenotype_count}};
    json::object metadata{{"path", manifest.meta_path}, {"format", "tsv"},
                          {"rows", manifest.variant_count}};
    json::object arrays;
    arrays["cov_yy"] = array_json(manifest.cov_yy);
    arrays["var_x"] = array_json(manifest.var_x);
    arrays["x_missing"] = manifest.has_x_missing
        ? json::value(array_json(manifest.x_missing)) : json::value(nullptr);
    arrays["y_missing"] = manifest.has_y_missing
        ? json::value(array_json(manifest.y_missing)) : json::value(nullptr);
    json::array phenotypes;
    for (const auto &phenotype : manifest.phenotypes) {
        phenotypes.emplace_back(json::object{{"index", phenotype.index},
                                             {"name", phenotype.name},
                                             {"cov_xy", array_json(phenotype.cov_xy)}});
    }
    json::object root{{"schema", "org.ukbiocoin.nss"},
                      {"schema_version", 2},
                      {"ukc_version", "2.0"},
                      {"byte_order", "little"},
                      {"dimensions", std::move(dimensions)},
                      {"metadata", std::move(metadata)},
                      {"arrays", std::move(arrays)},
                      {"phenotypes", std::move(phenotypes)}};
    const std::string path = join_path(directory, "manifest.json");
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        throw NssError("Error opening " + path + " for output");
    }
    output << json::serialize(root) << '\n';
    if (!output) {
        throw NssError("Error writing " + path);
    }
}

void validate_manifest(const std::string &directory,
                       const Manifest &manifest,
                       bool validate_all_cov_xy) {
    if (!directory_exists(directory)) {
        throw NssError(directory + ": NSS v2 directory does not exist");
    }
    if (manifest.sample_count == 0 || manifest.variant_count == 0 ||
        manifest.phenotype_count == 0 ||
        manifest.phenotypes.size() != manifest.phenotype_count) {
        throw NssError(directory + ": invalid manifest dimensions");
    }
    const std::vector<std::uint64_t> variants{manifest.variant_count};
    const std::vector<std::uint64_t> phenotypes{manifest.phenotype_count};
    const std::vector<std::uint64_t> covariance{
        manifest.phenotype_count, manifest.phenotype_count};
    validate_descriptor(directory, manifest.cov_yy, "<f8", covariance, "cov_yy");
    validate_descriptor(directory, manifest.var_x, "<f8", variants, "var_x");
    if (manifest.has_x_missing) {
        validate_descriptor(directory, manifest.x_missing, "<f8", variants, "x_missing");
    }
    if (manifest.has_y_missing) {
        validate_descriptor(directory, manifest.y_missing, "<f8", phenotypes, "y_missing");
    }
    if (validate_all_cov_xy) {
        for (const auto &phenotype : manifest.phenotypes) {
            validate_descriptor(directory, phenotype.cov_xy, "<f4", variants,
                                "cov_xy for '" + phenotype.name + "'");
        }
    }
    const std::string meta = join_path(directory, manifest.meta_path);
    const std::uint64_t rows = count_tsv_rows(meta);
    if (rows != manifest.variant_count) {
        throw NssError(meta + ": metadata row count " + std::to_string(rows) +
                       " does not match variant_count " +
                       std::to_string(manifest.variant_count));
    }
}

NpyReader::NpyReader(const std::string &path,
                     const std::string &expected_dtype,
                     const std::vector<std::uint64_t> &expected_shape)
    : path_(path), stream_(path, std::ios::binary) {
    if (!stream_) {
        throw NssError("Error opening " + path + " for input");
    }
    header_ = parse_npy_header(stream_, path);
    if (header_.fortran_order) {
        throw NssError(path + ": Fortran-order arrays are not supported");
    }
    if (header_.dtype != expected_dtype) {
        throw NssError(path + ": dtype " + header_.dtype +
                       " does not match manifest dtype " + expected_dtype);
    }
    if (header_.shape != expected_shape) {
        throw NssError(path + ": NPY shape does not match manifest shape");
    }
    if (header_.element_count >
        static_cast<std::uint64_t>(std::numeric_limits<std::streamoff>::max() /
                                   static_cast<std::streamoff>(header_.item_size))) {
        throw NssError(path + ": NPY payload is too large");
    }
    const std::streamoff expected_size = header_.data_offset +
        static_cast<std::streamoff>(header_.element_count * header_.item_size);
    stream_.seekg(0, std::ios::end);
    const std::streamoff actual_size = stream_.tellg();
    if (actual_size != expected_size) {
        throw NssError(path + ": NPY file size does not match its declared shape (expected " +
                       std::to_string(static_cast<long long>(expected_size)) +
                       " bytes, found " + std::to_string(static_cast<long long>(actual_size)) + ")");
    }
    stream_.seekg(header_.data_offset);
}

void NpyReader::read(double *destination, std::size_t count) {
    if (count > remaining()) {
        throw NssError(path_ + ": attempted to read beyond NPY payload");
    }
    const bool file_little = header_.dtype[0] != '>';
    const bool swap = file_little != host_is_little_endian();
    if (header_.item_size == 4) {
        std::vector<float> values(count);
        stream_.read(reinterpret_cast<char *>(values.data()),
                     static_cast<std::streamsize>(count * sizeof(float)));
        if (!stream_) throw NssError(path_ + ": truncated NPY payload");
        for (std::size_t i = 0; i < count; ++i) {
            if (swap) reverse_bytes(values[i]);
            destination[i] = static_cast<double>(values[i]);
        }
    } else {
        stream_.read(reinterpret_cast<char *>(destination),
                     static_cast<std::streamsize>(count * sizeof(double)));
        if (!stream_) throw NssError(path_ + ": truncated NPY payload");
        if (swap) {
            for (std::size_t i = 0; i < count; ++i) reverse_bytes(destination[i]);
        }
    }
    elements_read_ += count;
}

std::vector<double> NpyReader::read_all() {
    std::vector<double> values(static_cast<std::size_t>(remaining()));
    read(values.data(), values.size());
    return values;
}

void write_npy_f32(const std::string &path, const std::vector<float> &values) {
    std::ofstream output(path, std::ios::binary);
    if (!output) throw NssError("Error opening " + path + " for output");
    write_npy_prefix(output, "<f4", {static_cast<std::uint64_t>(values.size())}, path);
    write_little_endian(output, values, path);
}

void write_npy_f64(const std::string &path,
                   const std::vector<double> &values,
                   const std::vector<std::uint64_t> &shape) {
    if (checked_product(shape, path) != values.size()) {
        throw NssError(path + ": value count does not match requested NPY shape");
    }
    std::ofstream output(path, std::ios::binary);
    if (!output) throw NssError("Error opening " + path + " for output");
    write_npy_prefix(output, "<f8", shape, path);
    write_little_endian(output, values, path);
}

void create_npy_f32(const std::string &path, std::uint64_t count) {
    std::ofstream output(path, std::ios::binary);
    if (!output) throw NssError("Error opening " + path + " for output");
    write_npy_prefix(output, "<f4", {count}, path);
}

void append_f32_payload(const std::string &path, const std::vector<float> &values) {
    std::ofstream output(path, std::ios::binary | std::ios::app);
    if (!output) throw NssError("Error opening " + path + " for append");
    write_little_endian(output, values, path);
}

std::uint64_t count_tsv_rows(const std::string &path) {
    std::ifstream input(path);
    if (!input) throw NssError("Error opening " + path + " for input");
    std::string line;
    if (!std::getline(input, line) || line.empty()) {
        throw NssError(path + ": missing metadata header");
    }
    std::uint64_t rows = 0;
    while (std::getline(input, line)) {
        if (line.empty()) throw NssError(path + ": empty metadata row");
        ++rows;
    }
    return rows;
}

std::uint64_t peak_rss_kib() {
    struct rusage usage {};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#if defined(__APPLE__)
    return static_cast<std::uint64_t>(usage.ru_maxrss / 1024);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#endif
}

}  // namespace ukc

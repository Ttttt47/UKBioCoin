#include "adjustment.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <vector>

#include <boost/charconv.hpp>
#include <unistd.h>

#include "nss_v2.hpp"
#include "statistics.hpp"
#include "table_metadata.hpp"

namespace ukc {
namespace {

constexpr std::size_t kInputResultColumns = 7;
constexpr std::size_t kAdjustmentBatchRows = 32768;

class AtomicAdjustmentOutput {
public:
    explicit AtomicAdjustmentOutput(const std::string &final_path)
        : final_path_(final_path),
          temporary_path_(final_path + ".tmp." +
                          std::to_string(static_cast<long long>(::getpid()))),
          stream_(temporary_path_, std::ios::binary) {
        if (!stream_) {
            throw std::runtime_error("Error opening " + temporary_path_ +
                                     " for output");
        }
    }

    ~AtomicAdjustmentOutput() {
        if (!committed_) {
            stream_.close();
            std::remove(temporary_path_.c_str());
        }
    }

    std::ofstream &stream() { return stream_; }

    void commit() {
        stream_.close();
        if (!stream_) throw std::runtime_error("Error writing " + temporary_path_);
        if (std::rename(temporary_path_.c_str(), final_path_.c_str()) != 0) {
            throw std::runtime_error("Error replacing " + final_path_ + ": " +
                                     std::strerror(errno));
        }
        committed_ = true;
    }

private:
    std::string final_path_;
    std::string temporary_path_;
    std::ofstream stream_;
    bool committed_ = false;
};

struct TableLayout {
    std::vector<std::string> names;
    MetadataProjection metadata_projection;
    std::vector<std::size_t> metadata_output_indices;
    std::size_t t_index = 0;
    std::size_t p_index = 0;
};

std::string statistic_name(const std::string &name) {
    if (name == "BETA") return "BETA";
    if (name == "SE") return "SE";
    if (name == "T_STAT" || name == "T-STAT") return "T_STAT";
    if (name == "NEG_LOG10_P" || name == "-log10_P") return "NEG_LOG10_P";
    if (name == "VIF") return "VIF";
    if (name == "OBS_CT" || name == "nobs") return "OBS_CT";
    if (name == "QUALITY_SCORE" || name == "Quality-Score") {
        return "QUALITY_SCORE";
    }
    return "";
}

TableLayout table_layout(const std::string &header,
                         const std::string &path) {
    TableLayout layout;
    layout.names = split_table_fields(header, path, 1);
    std::unordered_map<std::string, std::size_t> statistics;
    std::vector<std::string> metadata_names;
    std::vector<std::size_t> metadata_source_indices;
    for (std::size_t i = 0; i < layout.names.size(); ++i) {
        if (layout.names[i] == "NEG_LOG10_UNADJ" ||
            layout.names[i] == "NEG_LOG10_GC") {
            throw TableError(path + ":1: input is already an adjusted result table");
        }
        const std::string canonical = statistic_name(layout.names[i]);
        if (canonical.empty()) {
            metadata_names.push_back(layout.names[i]);
            metadata_source_indices.push_back(i);
            continue;
        }
        const auto inserted = statistics.emplace(canonical, i);
        if (!inserted.second) {
            throw TableError(path + ":1: result columns '" +
                             layout.names[inserted.first->second] + "' and '" +
                             layout.names[i] + "' both map to '" + canonical + "'");
        }
    }
    const auto t = statistics.find("T_STAT");
    const auto p = statistics.find("NEG_LOG10_P");
    if (t == statistics.end() || p == statistics.end()) {
        throw TableError(path +
                         ":1: adjustment input requires T_STAT/NEG_LOG10_P "
                         "(or T-STAT/-log10_P) columns");
    }
    layout.t_index = t->second;
    layout.p_index = p->second;
    if (!metadata_names.empty()) {
        layout.metadata_projection =
            make_metadata_projection(metadata_names, path);
        layout.metadata_output_indices.reserve(
            layout.metadata_projection.output_indices.size());
        for (std::size_t relative_index :
             layout.metadata_projection.output_indices) {
            layout.metadata_output_indices.push_back(
                metadata_source_indices[relative_index]);
        }
    }
    return layout;
}

double parse_result_number(const std::string &text,
                           const std::string &path,
                           std::size_t line_number,
                           const std::string &column) {
    if (text == "NA" || text == "NaN" || text == "nan") {
        return std::numeric_limits<double>::quiet_NaN();
    }
    errno = 0;
    char *end = nullptr;
    const double value = std::strtod(text.c_str(), &end);
    if (end == text.c_str() || *end != '\0' || errno == ERANGE) {
        throw TableError(path + ":" + std::to_string(line_number) +
                         ": invalid numeric value for " + column + " ('" +
                         text + "')");
    }
    return value;
}

bool valid_adjustment_values(double t_statistic, double neg_log10_p) {
    if (!std::isfinite(t_statistic) || !std::isfinite(neg_log10_p) ||
        neg_log10_p < 0.0) {
        return false;
    }
    return true;
}

void check_row_width(const std::vector<std::string> &fields,
                     const TableLayout &layout, const std::string &path,
                     std::size_t line_number) {
    if (fields.size() != layout.names.size()) {
        throw TableError(path + ":" + std::to_string(line_number) +
                         ": field count " + std::to_string(fields.size()) +
                         " does not match header field count " +
                         std::to_string(layout.names.size()));
    }
}

std::size_t scan_table(const std::string &path, TableLayout &layout,
                       std::vector<double> &chi_squares) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Error opening " + path + " for input");
    std::string line;
    if (!std::getline(input, line)) throw TableError(path + ": missing header");
    layout = table_layout(line, path);
    std::size_t rows = 0;
    while (std::getline(input, line)) {
        ++rows;
        const std::size_t line_number = rows + 1;
        const std::vector<std::string> fields =
            split_table_fields(line, path, line_number);
        check_row_width(fields, layout, path, line_number);
        const double t = parse_result_number(fields[layout.t_index], path,
                                             line_number, layout.names[layout.t_index]);
        const double p = parse_result_number(fields[layout.p_index], path,
                                             line_number, layout.names[layout.p_index]);
        if (valid_adjustment_values(t, p)) {
            const double chi_square = t * t;
            if (!std::isfinite(chi_square)) {
                throw TableError(path + ":" + std::to_string(line_number) +
                                 ": T_STAT is too large to estimate lambda");
            }
            chi_squares.push_back(chi_square);
        }
    }
    if (rows == 0) throw TableError(path + ": adjustment input has no data rows");
    return rows;
}

void validate_adjustment_npy(const NpyHeader &header,
                             const std::string &path) {
    if (header.dtype != "<f8") {
        throw NssError(path + ": adjustment NPY dtype must be little-endian float64 (<f8)");
    }
    if (header.fortran_order) {
        throw NssError(path + ": adjustment NPY must use C order");
    }
    if (header.shape.size() != 2) {
        throw NssError(path + ": adjustment NPY must be two-dimensional");
    }
    if (header.shape[1] != kInputResultColumns) {
        throw NssError(path + ": adjustment NPY must have exactly 7 columns");
    }
}

std::size_t scan_npy(const std::string &path,
                     std::vector<double> &chi_squares) {
    NpyReader reader(path);
    validate_adjustment_npy(reader.header(), path);
    if (reader.header().shape[0] >
        static_cast<std::uint64_t>((std::numeric_limits<std::size_t>::max)())) {
        throw NssError(path +
                       ": adjustment NPY row count exceeds this platform's limit");
    }
    const std::size_t rows = static_cast<std::size_t>(reader.header().shape[0]);
    std::vector<double> values(kAdjustmentBatchRows * kInputResultColumns);
    std::size_t completed = 0;
    while (completed < rows) {
        const std::size_t count = std::min(kAdjustmentBatchRows, rows - completed);
        reader.read(values.data(), count * kInputResultColumns);
        for (std::size_t i = 0; i < count; ++i) {
            const double t = values[i * kInputResultColumns + 2];
            const double p = values[i * kInputResultColumns + 3];
            if (valid_adjustment_values(t, p)) {
                const double chi_square = t * t;
                if (!std::isfinite(chi_square)) {
                    throw NssError(path + ": T_STAT at row " +
                                   std::to_string(completed + i) +
                                   " is too large to estimate lambda");
                }
                chi_squares.push_back(chi_square);
            }
        }
        completed += count;
    }
    return rows;
}

double median(std::vector<double> &values) {
    const std::size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    if (values.size() % 2 != 0) return values[middle];
    const double lower = *std::max_element(values.begin(), values.begin() + middle);
    return lower + (values[middle] - lower) / 2.0;
}

bool append_double(std::string &output, double value) {
    if (!std::isfinite(value)) {
        output.append("NA");
        return true;
    }
    char buffer[128];
    const auto converted = boost::charconv::to_chars(
        buffer, buffer + sizeof(buffer), value,
        boost::charconv::chars_format::general,
        std::numeric_limits<double>::max_digits10);
    if (converted.ec != std::errc()) return false;
    output.append(buffer, converted.ptr);
    return true;
}

double adjusted_neg_log10(double t_statistic, double lambda_used) {
    return chi_square_one_df_neg_log10_from_z(
        std::abs(t_statistic) / std::sqrt(lambda_used));
}

void write_table_from_table(const AdjustmentOptions &options,
                            const std::string &output_path,
                            const TableLayout &layout,
                            std::size_t expected_rows,
                            double lambda_used) {
    std::ifstream input(options.input_path);
    if (!input) {
        throw std::runtime_error("Error opening " + options.input_path + " for input");
    }
    std::string line;
    std::getline(input, line);
    AtomicAdjustmentOutput file(output_path);
    std::ofstream &output = file.stream();
    if (layout.metadata_output_indices.empty()) {
        output << "#ROW_INDEX";
    } else {
        output << metadata_header_tsv(layout.metadata_projection);
    }
    output << "\tNEG_LOG10_UNADJ\tNEG_LOG10_GC\n";

    std::size_t row = 0;
    while (row < expected_rows && std::getline(input, line)) {
        ++row;
        const std::vector<std::string> fields =
            split_table_fields(line, options.input_path, row + 1);
        check_row_width(fields, layout, options.input_path, row + 1);
        if (layout.metadata_output_indices.empty()) {
            output << row - 1;
        } else {
            output << join_table_fields(fields, layout.metadata_output_indices);
        }
        const double t = parse_result_number(
            fields[layout.t_index], options.input_path, row + 1,
            layout.names[layout.t_index]);
        const double p = parse_result_number(
            fields[layout.p_index], options.input_path, row + 1,
            layout.names[layout.p_index]);
        output.put('\t');
        if (valid_adjustment_values(t, p)) {
            std::string numbers;
            if (!append_double(numbers, p)) {
                throw std::runtime_error("Error formatting unadjusted P value");
            }
            numbers.push_back('\t');
            if (!append_double(numbers, adjusted_neg_log10(t, lambda_used))) {
                throw std::runtime_error("Error formatting GC P value");
            }
            output << numbers;
        } else {
            output << "NA\tNA";
        }
        output.put('\n');
        if (!output) throw std::runtime_error("Error writing " + output_path);
    }
    if (row != expected_rows || std::getline(input, line)) {
        throw TableError(options.input_path +
                         ": input row count changed during adjustment output");
    }
    file.commit();
}

void write_table_from_npy(const AdjustmentOptions &options,
                          const std::string &output_path,
                          std::size_t expected_rows,
                          double lambda_used) {
    NpyReader reader(options.input_path);
    validate_adjustment_npy(reader.header(), options.input_path);
    const std::size_t rows = static_cast<std::size_t>(reader.header().shape[0]);
    if (rows != expected_rows) {
        throw NssError(options.input_path +
                       ": NPY shape changed during adjustment output");
    }
    AtomicAdjustmentOutput file(output_path);
    std::ofstream &output = file.stream();
    output << "#ROW_INDEX\tNEG_LOG10_UNADJ\tNEG_LOG10_GC\n";
    std::vector<double> values(kAdjustmentBatchRows * kInputResultColumns);
    std::size_t completed = 0;
    while (completed < rows) {
        const std::size_t count = std::min(kAdjustmentBatchRows, rows - completed);
        reader.read(values.data(), count * kInputResultColumns);
        std::string buffer;
        buffer.reserve(count * 64);
        for (std::size_t i = 0; i < count; ++i) {
            const double t = values[i * kInputResultColumns + 2];
            const double p = values[i * kInputResultColumns + 3];
            buffer.append(std::to_string(completed + i));
            buffer.push_back('\t');
            if (valid_adjustment_values(t, p)) {
                if (!append_double(buffer, p)) {
                    throw std::runtime_error("Error formatting unadjusted P value");
                }
                buffer.push_back('\t');
                if (!append_double(buffer, adjusted_neg_log10(t, lambda_used))) {
                    throw std::runtime_error("Error formatting GC P value");
                }
            } else {
                buffer.append("NA\tNA");
            }
            buffer.push_back('\n');
        }
        output.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        if (!output) throw std::runtime_error("Error writing " + output_path);
        completed += count;
    }
    file.commit();
}

template <typename RowSource>
void write_npy_adjustment(const std::string &output_path, std::size_t rows,
                          double lambda_used, RowSource source) {
    AtomicAdjustmentOutput file(output_path);
    std::ofstream &output = file.stream();
    write_npy_f64_header(output, {static_cast<std::uint64_t>(rows), 2}, output_path);
    std::vector<double> adjusted;
    adjusted.reserve(kAdjustmentBatchRows * 2);
    for (std::size_t row = 0; row < rows; ++row) {
        double t = std::numeric_limits<double>::quiet_NaN();
        double p = std::numeric_limits<double>::quiet_NaN();
        source(row, t, p);
        if (valid_adjustment_values(t, p)) {
            adjusted.push_back(p);
            adjusted.push_back(adjusted_neg_log10(t, lambda_used));
        } else {
            adjusted.push_back(std::numeric_limits<double>::quiet_NaN());
            adjusted.push_back(std::numeric_limits<double>::quiet_NaN());
        }
        if (adjusted.size() == kAdjustmentBatchRows * 2 || row + 1 == rows) {
            write_npy_f64_payload(output, adjusted.data(), adjusted.size(), output_path);
            adjusted.clear();
        }
    }
    file.commit();
}

void write_npy_from_table(const AdjustmentOptions &options,
                          const std::string &output_path,
                          const TableLayout &layout, std::size_t rows,
                          double lambda_used) {
    std::ifstream input(options.input_path);
    if (!input) {
        throw std::runtime_error("Error opening " + options.input_path + " for input");
    }
    std::string line;
    std::getline(input, line);
    write_npy_adjustment(output_path, rows, lambda_used,
        [&](std::size_t row, double &t, double &p) {
            if (!std::getline(input, line)) {
                throw TableError(options.input_path +
                                 ": input ended during adjustment output");
            }
            const std::vector<std::string> fields =
                split_table_fields(line, options.input_path, row + 2);
            check_row_width(fields, layout, options.input_path, row + 2);
            t = parse_result_number(fields[layout.t_index], options.input_path,
                                    row + 2, layout.names[layout.t_index]);
            p = parse_result_number(fields[layout.p_index], options.input_path,
                                    row + 2, layout.names[layout.p_index]);
            if (row + 1 == rows && std::getline(input, line)) {
                throw TableError(options.input_path +
                                 ": input row count changed during adjustment output");
            }
        });
}

void write_npy_from_npy(const AdjustmentOptions &options,
                        const std::string &output_path, std::size_t rows,
                        double lambda_used) {
    NpyReader reader(options.input_path);
    validate_adjustment_npy(reader.header(), options.input_path);
    if (reader.header().shape[0] != rows) {
        throw NssError(options.input_path +
                       ": NPY shape changed during adjustment output");
    }
    std::vector<double> row(kInputResultColumns);
    write_npy_adjustment(output_path, rows, lambda_used,
        [&](std::size_t, double &t, double &p) {
            reader.read(row.data(), row.size());
            t = row[2];
            p = row[3];
        });
}

bool has_npy_suffix(const std::string &path) {
    return path.size() >= 4 && path.compare(path.size() - 4, 4, ".npy") == 0;
}

}  // namespace

void run_adjustment(const AdjustmentOptions &options) {
    if (options.input_path.empty()) {
        throw std::runtime_error("--adjust-file must not be empty");
    }
    if (options.output_prefix.empty()) {
        throw std::runtime_error("--out must not be empty");
    }
    if (options.output_format != "table" && options.output_format != "npy") {
        throw std::runtime_error("--output-format must be 'table' or 'npy'");
    }
    if (options.lambda_supplied &&
        (!std::isfinite(options.lambda) || options.lambda <= 0.0)) {
        throw std::runtime_error("--lambda must be finite and greater than zero");
    }

    const bool npy_input = has_npy_suffix(options.input_path);
    std::vector<double> chi_squares;
    TableLayout layout;
    const std::size_t rows = npy_input
        ? scan_npy(options.input_path, chi_squares)
        : scan_table(options.input_path, layout, chi_squares);
    if (chi_squares.empty()) {
        throw std::runtime_error(
            options.input_path + ": no valid tests are available for GC adjustment");
    }
    const std::size_t valid_tests = chi_squares.size();
    const double lambda_raw = median(chi_squares) / 0.456;
    const double lambda_used = options.lambda_supplied
        ? options.lambda : std::max(1.0, lambda_raw);
    if (!std::isfinite(lambda_raw)) {
        throw std::runtime_error("Unable to estimate a finite genomic-control lambda");
    }

    const std::string output_path = options.output_prefix +
        (options.output_format == "npy" ? ".adjusted.npy" : ".adjusted");
    if (options.output_format == "table") {
        if (npy_input) {
            write_table_from_npy(options, output_path, rows, lambda_used);
        } else {
            write_table_from_table(options, output_path, layout, rows,
                                   lambda_used);
        }
    } else if (npy_input) {
        write_npy_from_npy(options, output_path, rows, lambda_used);
    } else {
        write_npy_from_table(options, output_path, layout, rows, lambda_used);
    }

    std::cout << "Adjustment mode: genomic control\n"
              << "Input format: " << (npy_input ? "npy" : "table") << '\n'
              << "Output format: " << options.output_format << '\n'
              << "Rows: " << rows << '\n'
              << "Valid tests: " << valid_tests << '\n'
              << std::setprecision(std::numeric_limits<double>::max_digits10)
              << "lambda_raw: " << lambda_raw << '\n'
              << "lambda_used: " << lambda_used << '\n'
              << "Output: " << output_path << '\n'
              << std::defaultfloat;
}

}  // namespace ukc

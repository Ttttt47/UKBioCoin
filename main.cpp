#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <unistd.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/program_options.hpp>

#include "nss_v2.hpp"
#include "tools/cpp_cdfs-master/cdf_chisqt/cdf_base.h"

namespace po = boost::program_options;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::size_t;
using std::string;
using std::vector;

namespace {

constexpr size_t kBatchSize = 4096;
using Clock = std::chrono::steady_clock;

class UkcError : public std::runtime_error {
public:
    explicit UkcError(const string &message) : std::runtime_error(message) {}
};

double elapsed_seconds(const Clock::time_point &start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

struct Timings {
    double static_data = 0.0;
    double loading = 0.0;
    double compute = 0.0;
    double output = 0.0;
    double total = 0.0;
};

class ProgressReporter {
public:
    explicit ProgressReporter(std::uint64_t total)
        : total_(total), terminal_(::isatty(STDOUT_FILENO) != 0),
          step_(terminal_ ? 1 : 5) {}

    void update(std::uint64_t completed) {
        const unsigned percent = completed >= total_ ? 100U :
            static_cast<unsigned>(100.0L * static_cast<long double>(completed) /
                                  static_cast<long double>(total_));
        while (next_ <= percent) {
            if (terminal_) {
                std::cout << "\rProgress: " << std::setw(3) << next_ << "%" << std::flush;
            } else {
                std::cout << "Progress: " << next_ << "%\n";
            }
            next_ += step_;
        }
    }

    void finish() {
        update(total_);
        if (terminal_) std::cout << '\n';
    }

private:
    std::uint64_t total_;
    bool terminal_;
    unsigned step_;
    unsigned next_ = 0;
};

class AtomicOutputFile {
public:
    explicit AtomicOutputFile(const string &final_path)
        : final_path_(final_path),
          temporary_path_(final_path + ".tmp." +
                          std::to_string(static_cast<long long>(::getpid()))),
          stream_(temporary_path_) {
        if (!stream_) throw UkcError("Error opening " + temporary_path_ + " for output");
    }

    ~AtomicOutputFile() {
        if (!committed_) {
            stream_.close();
            std::remove(temporary_path_.c_str());
        }
    }

    std::ofstream &stream() { return stream_; }

    void commit() {
        stream_.close();
        if (!stream_) throw UkcError("Error writing " + temporary_path_);
        if (std::rename(temporary_path_.c_str(), final_path_.c_str()) != 0) {
            throw UkcError("Error replacing " + final_path_ + ": " + std::strerror(errno));
        }
        committed_ = true;
    }

private:
    string final_path_;
    string temporary_path_;
    std::ofstream stream_;
    bool committed_ = false;
};

struct TokenView {
    const char *begin = nullptr;
    const char *end = nullptr;
};

struct SelectedColumn {
    size_t source_index = 0;
    size_t destination_index = 0;
    string name;
};

struct InputRow {
    string cov_xy;
    string var_x;
    string meta;
    string x_missing;
    size_t line_number = 0;
};

struct RowResult {
    std::array<double, 5> summary{{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()}};
    int sample_size = 0;
    double quality_score = std::numeric_limits<double>::quiet_NaN();
    bool invalid = false;
    string error;
};

bool next_token(const char *&cursor, const char *end, TokenView &token) {
    while (cursor < end && std::isspace(static_cast<unsigned char>(*cursor))) ++cursor;
    if (cursor == end) return false;
    if (*cursor == '"') {
        ++cursor;
        token.begin = cursor;
        while (cursor < end && *cursor != '"') ++cursor;
        if (cursor == end) throw UkcError("Unterminated quoted field");
        token.end = cursor++;
        if (cursor < end && !std::isspace(static_cast<unsigned char>(*cursor))) {
            throw UkcError("Unexpected character after quoted field");
        }
        return true;
    }
    token.begin = cursor;
    while (cursor < end && !std::isspace(static_cast<unsigned char>(*cursor))) ++cursor;
    token.end = cursor;
    return true;
}

string token_string(const TokenView &token) { return string(token.begin, token.end); }

double parse_number(const TokenView &token, const string &filename,
                    size_t line_number, const string &field_name) {
    if (token.end - token.begin == 2 && token.begin[0] == 'N' && token.begin[1] == 'A') {
        return std::numeric_limits<double>::quiet_NaN();
    }
    errno = 0;
    char *parsed_end = nullptr;
    const double value = std::strtod(token.begin, &parsed_end);
    if (parsed_end != token.end || parsed_end == token.begin || errno == ERANGE) {
        throw UkcError(filename + ":" + std::to_string(line_number) +
                       ": invalid numeric value for " + field_name + " ('" +
                       token_string(token) + "')");
    }
    return value;
}

vector<string> parse_header(const string &line, const string &filename) {
    const char *cursor = line.data();
    const char *end = cursor + line.size();
    TokenView token;
    vector<string> names;
    std::unordered_set<string> seen;
    try {
        while (next_token(cursor, end, token)) {
            string name = token_string(token);
            if (!name.empty()) {
                if (!seen.insert(name).second) throw UkcError("Duplicate column '" + name + "'");
                names.push_back(std::move(name));
            }
        }
    } catch (const UkcError &error) {
        throw UkcError(filename + ":1: " + error.what());
    }
    if (names.empty()) throw UkcError(filename + ":1: empty header");
    return names;
}

vector<SelectedColumn> select_columns(const vector<string> &header,
                                      const vector<string> &requested,
                                      const string &filename) {
    std::unordered_map<string, size_t> indices;
    for (size_t i = 0; i < header.size(); ++i) indices.emplace(header[i], i);
    vector<SelectedColumn> selected;
    for (size_t i = 0; i < requested.size(); ++i) {
        const auto found = indices.find(requested[i]);
        if (found == indices.end()) {
            throw UkcError(filename + ": phenotype/covariate '" + requested[i] + "' not found");
        }
        selected.push_back({found->second, i, requested[i]});
    }
    std::sort(selected.begin(), selected.end(), [](const SelectedColumn &a,
                                                    const SelectedColumn &b) {
        return a.source_index < b.source_index;
    });
    return selected;
}

vector<double> parse_selected_values(const string &line,
                                     const vector<SelectedColumn> &selected,
                                     const string &filename,
                                     size_t line_number) {
    vector<double> values(selected.size(), std::numeric_limits<double>::quiet_NaN());
    const char *cursor = line.data();
    const char *end = cursor + line.size();
    TokenView token;
    try {
        if (!next_token(cursor, end, token)) throw UkcError("empty row");
        size_t source = 0;
        size_t wanted = 0;
        while (wanted < selected.size() && next_token(cursor, end, token)) {
            if (source == selected[wanted].source_index) {
                const auto &column = selected[wanted];
                values[column.destination_index] =
                    parse_number(token, filename, line_number, column.name);
                ++wanted;
            }
            ++source;
        }
        if (wanted != selected.size()) {
            throw UkcError("row ends before required column '" + selected[wanted].name + "'");
        }
    } catch (const UkcError &error) {
        const string prefix = filename + ":" + std::to_string(line_number) + ": ";
        if (string(error.what()).compare(0, prefix.size(), prefix) == 0) throw;
        throw UkcError(prefix + error.what());
    }
    return values;
}

double parse_second_value(const string &line, const string &filename,
                          size_t line_number, const string &field_name) {
    const char *cursor = line.data();
    const char *end = cursor + line.size();
    TokenView token;
    try {
        if (!next_token(cursor, end, token) || !next_token(cursor, end, token)) {
            throw UkcError("row does not contain " + field_name);
        }
        return parse_number(token, filename, line_number, field_name);
    } catch (const UkcError &error) {
        const string prefix = filename + ":" + std::to_string(line_number) + ": ";
        if (string(error.what()).compare(0, prefix.size(), prefix) == 0) throw;
        throw UkcError(prefix + error.what());
    }
}

string first_field(const string &line, const string &filename, size_t line_number) {
    const char *cursor = line.data();
    const char *end = cursor + line.size();
    TokenView token;
    try {
        if (!next_token(cursor, end, token)) throw UkcError("empty row");
        return token_string(token);
    } catch (const UkcError &error) {
        throw UkcError(filename + ":" + std::to_string(line_number) + ": " + error.what());
    }
}

std::ifstream open_input(const string &filename, bool binary = false) {
    std::ifstream stream(filename, binary ? std::ios::binary : std::ios::in);
    if (!stream) throw UkcError("Error opening " + filename + " for input");
    return stream;
}

vector<string> split_covariates(const string &input) {
    vector<string> values;
    size_t start = 0;
    while (start <= input.size()) {
        const size_t comma = input.find(',', start);
        const size_t stop = comma == string::npos ? input.size() : comma;
        size_t left = start;
        size_t right = stop;
        while (left < right && std::isspace(static_cast<unsigned char>(input[left]))) ++left;
        while (right > left && std::isspace(static_cast<unsigned char>(input[right - 1]))) --right;
        if (left == right) throw UkcError("--covar contains an empty name");
        values.emplace_back(input.substr(left, right - left));
        if (comma == string::npos) break;
        start = comma + 1;
    }
    return values;
}

void validate_unique_variables(const vector<string> &variables) {
    std::unordered_set<string> seen;
    for (const auto &name : variables) {
        if (name.empty()) throw UkcError("Phenotype/covariate names must not be empty");
        if (!seen.insert(name).second) {
            throw UkcError("Phenotype/covariate '" + name + "' was specified more than once");
        }
    }
}

MatrixXd read_legacy_covariance(const string &filename, const vector<string> &variables) {
    std::ifstream stream = open_input(filename);
    string line;
    if (!std::getline(stream, line)) throw UkcError(filename + ": missing header");
    const vector<SelectedColumn> selected = select_columns(parse_header(line, filename), variables, filename);
    MatrixXd result = MatrixXd::Zero(static_cast<Eigen::Index>(variables.size()),
                                     static_cast<Eigen::Index>(variables.size()));
    std::unordered_map<string, size_t> requested;
    for (size_t i = 0; i < variables.size(); ++i) requested.emplace(variables[i], i);
    vector<bool> found(variables.size(), false);
    size_t line_number = 1;
    while (std::getline(stream, line)) {
        ++line_number;
        const string row_name = first_field(line, filename, line_number);
        const auto match = requested.find(row_name);
        if (match == requested.end()) continue;
        if (found[match->second]) {
            throw UkcError(filename + ":" + std::to_string(line_number) +
                           ": duplicate row '" + row_name + "'");
        }
        const vector<double> row = parse_selected_values(line, selected, filename, line_number);
        for (size_t column = 0; column < row.size(); ++column) {
            result(static_cast<Eigen::Index>(match->second),
                   static_cast<Eigen::Index>(column)) = row[column];
        }
        found[match->second] = true;
    }
    for (size_t i = 0; i < found.size(); ++i) {
        if (!found[i]) throw UkcError(filename + ": missing covariance row '" + variables[i] + "'");
    }
    if (!result.allFinite()) throw UkcError(filename + ": selected covariance matrix contains NA or non-finite values");
    return result;
}

double read_legacy_y_quality(const string &filename, const vector<string> &variables) {
    std::ifstream stream = open_input(filename);
    std::unordered_map<string, size_t> requested;
    for (size_t i = 0; i < variables.size(); ++i) requested.emplace(variables[i], i);
    vector<bool> found(variables.size(), false);
    double quality = 1.0;
    string line;
    size_t line_number = 0;
    while (std::getline(stream, line)) {
        ++line_number;
        const string name = first_field(line, filename, line_number);
        const auto match = requested.find(name);
        if (match == requested.end()) continue;
        if (found[match->second]) throw UkcError(filename + ": duplicate missing-rate row '" + name + "'");
        const double missing = parse_second_value(line, filename, line_number, "missing_rate");
        if (!std::isfinite(missing) || missing < 0.0 || missing > 1.0) {
            throw UkcError(filename + ":" + std::to_string(line_number) +
                           ": missing_rate must be between 0 and 1");
        }
        quality *= 1.0 - missing;
        found[match->second] = true;
    }
    for (size_t i = 0; i < found.size(); ++i) {
        if (!found[i]) throw UkcError(filename + ": missing rate for '" + variables[i] + "'");
    }
    return quality;
}

class RegressionModel {
public:
    RegressionModel(const MatrixXd &phenotype_covariance, size_t covariate_count)
        : covariate_count_(covariate_count),
          theta_base_(MatrixXd::Zero(static_cast<Eigen::Index>(covariate_count + 2),
                                     static_cast<Eigen::Index>(covariate_count + 2))) {
        theta_base_(0, 0) = phenotype_covariance(0, 0);
        if (covariate_count_ == 0) return;
        theta_base_.block(0, 2, 1, static_cast<Eigen::Index>(covariate_count_)) =
            phenotype_covariance.block(0, 1, 1, static_cast<Eigen::Index>(covariate_count_));
        theta_base_.block(2, 0, static_cast<Eigen::Index>(covariate_count_), 1) =
            phenotype_covariance.block(1, 0, static_cast<Eigen::Index>(covariate_count_), 1);
        theta_base_.block(2, 2, static_cast<Eigen::Index>(covariate_count_),
                          static_cast<Eigen::Index>(covariate_count_)) =
            phenotype_covariance.block(1, 1, static_cast<Eigen::Index>(covariate_count_),
                                        static_cast<Eigen::Index>(covariate_count_));
        const MatrixXd covariance = theta_base_.block(
            2, 2, static_cast<Eigen::Index>(covariate_count_),
            static_cast<Eigen::Index>(covariate_count_));
        Eigen::FullPivLU<MatrixXd> decomposition(covariance);
        if (!decomposition.isInvertible()) {
            throw UkcError("Singularity detected in the covariance matrix; check whether covariates are linearly dependent");
        }
        covariate_inverse_ = decomposition.inverse();
    }

    int minimum_valid_sample_size() const { return static_cast<int>(covariate_count_ + 1); }

    std::array<double, 5> calculate(const vector<double> &cov_xy,
                                    double var_x, int sample_size) const {
        if (covariate_count_ == 0) {
            const double omega_inverse = 1.0 / var_x;
            const double beta = omega_inverse * cov_xy[0];
            const double rss = sample_size * (1.0 - beta * cov_xy[0]) /
                               (sample_size - 1.0);
            const double var_beta = rss * omega_inverse / sample_size;
            const double se = std::sqrt(var_beta);
            const double t = beta / se;
            const double log_p = cdf_t_log(-std::abs(t), sample_size - 2) + std::log(2.0);
            const double vif = var_beta / (rss / ((sample_size - 1.0) * var_x));
            return {{beta, se, t, -log_p * std::log10(std::exp(1.0)), vif}};
        }

        MatrixXd theta = theta_base_;
        theta(1, 0) = theta(0, 1) = cov_xy[0];
        theta(1, 1) = var_x;
        for (size_t i = 1; i < cov_xy.size(); ++i) {
            theta(1, static_cast<Eigen::Index>(i + 1)) = cov_xy[i];
            theta(static_cast<Eigen::Index>(i + 1), 1) = cov_xy[i];
        }
        const Eigen::Index dimension = theta.rows();
        const MatrixXd omega = theta.block(1, 1, dimension - 1, dimension - 1);
        MatrixXd omega_inverse = omega;
        const double q = 1.0 /
            (omega(0, 0) - (omega.block(0, 1, 1, dimension - 2) * covariate_inverse_ *
                             omega.block(1, 0, dimension - 2, 1)).value());
        omega_inverse(0, 0) = q;
        omega_inverse.block(0, 1, 1, dimension - 2) =
            -q * omega.block(0, 1, 1, dimension - 2) * covariate_inverse_;
        omega_inverse.block(1, 0, dimension - 2, 1) =
            -q * covariate_inverse_ * omega.block(1, 0, dimension - 2, 1);
        omega_inverse.block(1, 1, dimension - 2, dimension - 2) =
            covariate_inverse_ + q * covariate_inverse_ *
            omega.block(1, 0, dimension - 2, 1) *
            omega.block(0, 1, 1, dimension - 2) * covariate_inverse_;
        const VectorXd association = theta.block(0, 1, 1, dimension - 1).transpose();
        const VectorXd coefficients = omega_inverse * association;
        const double rss = sample_size *
            (1.0 - (coefficients.transpose() * association).value()) /
            (sample_size - (dimension - 2) - 1.0);
        const double var_beta = rss * omega_inverse(0, 0) / sample_size;
        const double se = std::sqrt(var_beta);
        const double t = coefficients[0] / se;
        const double log_p = cdf_t_log(-std::abs(t), sample_size - 2) + std::log(2.0);
        const double vif = var_beta / (rss / ((sample_size - 1.0) * theta(1, 1)));
        return {{coefficients[0], se, t,
                 -log_p * std::log10(std::exp(1.0)), vif}};
    }

private:
    size_t covariate_count_;
    MatrixXd theta_base_;
    MatrixXd covariate_inverse_;
};

bool any_non_finite(const vector<double> &values) {
    for (double value : values) if (!std::isfinite(value)) return true;
    return false;
}

void compute_numeric_batch(const vector<vector<double>> &cov_xy,
                           const vector<double> &var_x,
                           const vector<double> *x_missing,
                           vector<RowResult> &results,
                           const RegressionModel &model,
                           int total_sample_size,
                           double overall_non_missing_rate,
                           double y_quality_score,
                           bool use_missing_rate_estimate,
                           int threads) {
    const size_t count = var_x.size();
    results.assign(count, RowResult{});
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (long long raw = 0; raw < static_cast<long long>(count); ++raw) {
        const size_t i = static_cast<size_t>(raw);
        RowResult &result = results[i];
        vector<double> values(cov_xy.size());
        for (size_t column = 0; column < cov_xy.size(); ++column) values[column] = cov_xy[column][i];
        double quality = overall_non_missing_rate;
        if (use_missing_rate_estimate) {
            const double missing = (*x_missing)[i];
            if (!std::isfinite(missing) || missing < 0.0 || missing > 1.0) {
                result.error = "x_missing at variant " + std::to_string(i + 1) +
                               " must be between 0 and 1";
                continue;
            }
            quality = y_quality_score * (1.0 - missing);
        }
        result.quality_score = quality;
        result.sample_size = static_cast<int>(std::floor(total_sample_size * quality));
        if (any_non_finite(values) || !std::isfinite(var_x[i]) ||
            result.sample_size <= 2 ||
            result.sample_size <= model.minimum_valid_sample_size()) {
            result.invalid = true;
            continue;
        }
        result.summary = model.calculate(values, var_x[i], result.sample_size);
        for (double value : result.summary) {
            if (!std::isfinite(value)) result.invalid = true;
        }
    }
}

void write_result_rows(std::ofstream &output, const vector<string> &metadata,
                       const vector<RowResult> &results, size_t &invalid_rows) {
    for (size_t i = 0; i < results.size(); ++i) {
        if (!results[i].error.empty()) throw UkcError(results[i].error);
        if (results[i].invalid) ++invalid_rows;
        output << metadata[i] << ' ' << results[i].summary[0] << ' '
               << results[i].summary[1] << ' ' << results[i].summary[2] << ' '
               << results[i].summary[3] << ' ' << results[i].summary[4] << ' '
               << ' ' << results[i].sample_size << ' '
               << results[i].quality_score << '\n';
    }
}

MatrixXd select_v2_covariance(const ukc::Manifest &manifest,
                              const vector<size_t> &indices,
                              const vector<double> &full) {
    MatrixXd result(static_cast<Eigen::Index>(indices.size()),
                    static_cast<Eigen::Index>(indices.size()));
    const size_t width = static_cast<size_t>(manifest.phenotype_count);
    for (size_t row = 0; row < indices.size(); ++row) {
        for (size_t column = 0; column < indices.size(); ++column) {
            result(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                full[indices[row] * width + indices[column]];
        }
    }
    if (!result.allFinite()) throw UkcError("cov_yy.npy: selected covariance matrix contains NA or non-finite values");
    return result;
}

vector<size_t> select_v2_variables(const ukc::Manifest &manifest,
                                   const vector<string> &variables) {
    std::unordered_map<string, size_t> indices;
    for (const auto &phenotype : manifest.phenotypes) indices.emplace(phenotype.name, phenotype.index);
    vector<size_t> selected;
    for (const auto &name : variables) {
        const auto found = indices.find(name);
        if (found == indices.end()) throw UkcError("manifest.json: phenotype/covariate '" + name + "' not found");
        selected.push_back(found->second);
    }
    return selected;
}

double selected_y_quality(const ukc::Manifest &manifest,
                          const vector<size_t> &indices,
                          const vector<double> &missing) {
    double quality = 1.0;
    for (size_t index : indices) {
        const double value = missing[index];
        if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
            throw UkcError("y_missing.npy: missing rate for '" +
                           manifest.phenotypes[index].name + "' must be between 0 and 1");
        }
        quality *= 1.0 - value;
    }
    return quality;
}

void print_run_summary(const string &format, const string &io_mode,
                       int total_sample_size, const ukc::Manifest *manifest,
                       std::uint64_t variant_count, size_t selected_count,
                       const vector<string> &variables, int threads,
                       bool use_missing, std::uint64_t estimated_bytes) {
    std::cout << "Input format: " << format << '\n'
              << "I/O mode: " << io_mode << '\n'
              << "Sample count: " << total_sample_size << '\n'
              << "Variants: " << variant_count << '\n';
    if (manifest != nullptr) {
        std::cout << "Phenotypes in NSS: " << manifest->phenotype_count << '\n'
                  << "Numeric storage: cov_xy=<f4; cov_yy/var_x/missing=<f8\n";
    }
    std::cout << "Selected variables (" << selected_count << "): ";
    for (size_t i = 0; i < variables.size(); ++i) std::cout << (i ? "," : "") << variables[i];
    std::cout << '\n' << "Missing-rate sample estimation: " << (use_missing ? "yes" : "no")
              << '\n' << "Estimated working memory: " << std::fixed << std::setprecision(2)
              << static_cast<double>(estimated_bytes) / (1024.0 * 1024.0) << " MiB\n"
              << "Threads: " << threads << '\n' << std::defaultfloat;
}

void report_timings(const Timings &timings, const string &io_mode) {
    std::cout << std::fixed << std::setprecision(3)
              << "Timing summary (seconds):\n"
              << "  manifest/static data: " << timings.static_data << '\n'
              << (io_mode == "memory" ? "  data loading:         " : "  data reading:         ")
              << timings.loading << '\n'
              << "  computation:          " << timings.compute << '\n'
              << "  result writing:       " << timings.output << '\n'
              << "  total:                " << timings.total << '\n'
              << "Peak RSS: " << static_cast<double>(ukc::peak_rss_kib()) / 1024.0 << " MiB\n"
              << std::defaultfloat;
}

void run_v2(const string &directory, const string &result_filename,
            const vector<string> &variables, bool total_size_supplied,
            int supplied_total_size,
            double overall_non_missing_rate, bool use_missing,
            const string &io_mode, int threads) {
    const auto total_start = Clock::now();
    Timings timings;
    const auto static_start = Clock::now();
    const ukc::Manifest manifest = ukc::read_manifest(directory);
    if (manifest.sample_count >
        static_cast<std::uint64_t>((std::numeric_limits<int>::max)())) {
        throw UkcError("manifest sample_count exceeds supported integer range");
    }
    const int total_sample_size = static_cast<int>(manifest.sample_count);
    if (total_size_supplied && supplied_total_size != total_sample_size) {
        throw UkcError("--totalsize " + std::to_string(supplied_total_size) +
                       " does not match manifest sample_count " +
                       std::to_string(total_sample_size));
    }
    const int fixed_sample_size = static_cast<int>(std::floor(
        total_sample_size * overall_non_missing_rate));
    if (!use_missing && (fixed_sample_size <= 2 ||
        fixed_sample_size <= static_cast<int>(variables.size()))) {
        throw UkcError("Effective sample size is too small for the requested model");
    }
    ukc::validate_manifest(directory, manifest, true);
    const vector<size_t> indices = select_v2_variables(manifest, variables);
    ukc::NpyReader cov_reader(ukc::join_path(directory, manifest.cov_yy.path),
                             "<f8", manifest.cov_yy.shape);
    const MatrixXd covariance = select_v2_covariance(manifest, indices, cov_reader.read_all());
    double y_quality = 1.0;
    if (use_missing) {
        if (!manifest.has_x_missing || !manifest.has_y_missing) {
            throw UkcError("NSS v2 manifest does not provide both x_missing and y_missing arrays");
        }
        ukc::NpyReader y_reader(ukc::join_path(directory, manifest.y_missing.path),
                               "<f8", manifest.y_missing.shape);
        y_quality = selected_y_quality(manifest, indices, y_reader.read_all());
    }
    RegressionModel model(covariance, variables.size() - 1);
    timings.static_data = elapsed_seconds(static_start);

    const std::uint64_t per_variant = static_cast<std::uint64_t>(variables.size() + 1 +
        (use_missing ? 1 : 0)) * sizeof(double);
    const std::uint64_t estimated = per_variant *
        (io_mode == "memory" ? manifest.variant_count : kBatchSize);
    print_run_summary("NSS v2", io_mode, total_sample_size, &manifest,
                      manifest.variant_count, variables.size(), variables,
                      threads, use_missing, estimated);

    AtomicOutputFile output_file(result_filename);
    std::ofstream &output = output_file.stream();
    const string meta_path = ukc::join_path(directory, manifest.meta_path);
    std::ifstream meta = open_input(meta_path);
    string meta_header;
    if (!std::getline(meta, meta_header)) throw UkcError(meta_path + ": missing header");
    output << meta_header << " BETA SE T-STAT -log10_P VIF nobs Quality-Score\n";
    ProgressReporter progress(manifest.variant_count);
    size_t invalid_rows = 0;

    vector<std::unique_ptr<ukc::NpyReader>> cov_readers;
    for (size_t index : indices) {
        const auto &descriptor = manifest.phenotypes[index].cov_xy;
        cov_readers.emplace_back(new ukc::NpyReader(
            ukc::join_path(directory, descriptor.path), "<f4", descriptor.shape));
    }
    ukc::NpyReader var_reader(ukc::join_path(directory, manifest.var_x.path),
                              "<f8", manifest.var_x.shape);
    std::unique_ptr<ukc::NpyReader> x_reader;
    if (use_missing) {
        x_reader.reset(new ukc::NpyReader(
            ukc::join_path(directory, manifest.x_missing.path), "<f8",
            manifest.x_missing.shape));
    }

    if (io_mode == "memory") {
        const auto load_start = Clock::now();
        vector<vector<double>> cov_xy;
        for (auto &reader : cov_readers) cov_xy.push_back(reader->read_all());
        vector<double> var_x = var_reader.read_all();
        vector<double> x_missing;
        if (x_reader) x_missing = x_reader->read_all();
        vector<string> metadata;
        metadata.reserve(static_cast<size_t>(manifest.variant_count));
        string line;
        while (std::getline(meta, line)) metadata.push_back(std::move(line));
        if (metadata.size() != manifest.variant_count) {
            throw UkcError(meta_path + ": metadata row-count mismatch while loading");
        }
        timings.loading += elapsed_seconds(load_start);

        vector<RowResult> results;
        const auto compute_start = Clock::now();
        compute_numeric_batch(cov_xy, var_x, x_reader ? &x_missing : nullptr,
                              results, model, total_sample_size,
                              overall_non_missing_rate, y_quality, use_missing,
                              threads);
        timings.compute += elapsed_seconds(compute_start);
        const auto output_start = Clock::now();
        write_result_rows(output, metadata, results, invalid_rows);
        timings.output += elapsed_seconds(output_start);
        progress.finish();
    } else {
        std::uint64_t completed = 0;
        while (completed < manifest.variant_count) {
            const size_t count = static_cast<size_t>(std::min<std::uint64_t>(
                kBatchSize, manifest.variant_count - completed));
            const auto load_start = Clock::now();
            vector<vector<double>> cov_xy(cov_readers.size(), vector<double>(count));
            for (size_t column = 0; column < cov_readers.size(); ++column) {
                cov_readers[column]->read(cov_xy[column].data(), count);
            }
            vector<double> var_x(count);
            var_reader.read(var_x.data(), count);
            vector<double> x_missing;
            if (x_reader) {
                x_missing.resize(count);
                x_reader->read(x_missing.data(), count);
            }
            vector<string> metadata(count);
            for (size_t i = 0; i < count; ++i) {
                if (!std::getline(meta, metadata[i])) {
                    throw UkcError(meta_path + ": metadata row-count mismatch near row " +
                                   std::to_string(completed + i + 1));
                }
            }
            timings.loading += elapsed_seconds(load_start);
            vector<RowResult> results;
            const auto compute_start = Clock::now();
            compute_numeric_batch(cov_xy, var_x, x_reader ? &x_missing : nullptr,
                                  results, model, total_sample_size,
                                  overall_non_missing_rate, y_quality,
                                  use_missing, threads);
            timings.compute += elapsed_seconds(compute_start);
            const auto output_start = Clock::now();
            write_result_rows(output, metadata, results, invalid_rows);
            timings.output += elapsed_seconds(output_start);
            completed += count;
            progress.update(completed);
        }
        progress.finish();
    }

    const auto commit_start = Clock::now();
    output_file.commit();
    timings.output += elapsed_seconds(commit_start);
    timings.total = elapsed_seconds(total_start);
    if (invalid_rows != 0) {
        std::cerr << "Warning: " << invalid_rows
                  << " SNP rows produced non-finite results because of NA or invalid covariance values\n";
    }
    report_timings(timings, io_mode);
}

void process_legacy_batch(const vector<InputRow> &batch,
                          vector<RowResult> &results,
                          const vector<SelectedColumn> &selected,
                          const RegressionModel &model,
                          const string &cov_xy_file,
                          const string &var_x_file,
                          const string &x_missing_file,
                          int total_sample_size,
                          double overall_non_missing_rate,
                          double y_quality,
                          bool use_missing, int threads) {
    results.assign(batch.size(), RowResult{});
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (long long raw = 0; raw < static_cast<long long>(batch.size()); ++raw) {
        const size_t i = static_cast<size_t>(raw);
        try {
            const vector<double> cov_xy = parse_selected_values(
                batch[i].cov_xy, selected, cov_xy_file, batch[i].line_number);
            const double var_x = parse_second_value(
                batch[i].var_x, var_x_file, batch[i].line_number, "var_x");
            double quality = overall_non_missing_rate;
            if (use_missing) {
                const double missing = parse_second_value(
                    batch[i].x_missing, x_missing_file, batch[i].line_number,
                    "x_missing");
                if (!std::isfinite(missing) || missing < 0.0 || missing > 1.0) {
                    throw UkcError(x_missing_file + ":" +
                                   std::to_string(batch[i].line_number) +
                                   ": x_missing must be between 0 and 1");
                }
                quality = y_quality * (1.0 - missing);
            }
            RowResult &result = results[i];
            result.quality_score = quality;
            result.sample_size = static_cast<int>(std::floor(total_sample_size * quality));
            if (any_non_finite(cov_xy) || !std::isfinite(var_x) ||
                result.sample_size <= 2 ||
                result.sample_size <= model.minimum_valid_sample_size()) {
                result.invalid = true;
                continue;
            }
            result.summary = model.calculate(cov_xy, var_x, result.sample_size);
            for (double value : result.summary) if (!std::isfinite(value)) result.invalid = true;
        } catch (const std::exception &error) {
            results[i].error = error.what();
        }
    }
}

void run_v1(const string &prefix, const string &result_filename,
            const vector<string> &variables, int total_sample_size,
            double overall_non_missing_rate, bool use_missing,
            const string &io_mode, int threads) {
    const auto total_start = Clock::now();
    Timings timings;
    const string cov_xy_file = prefix + "_cov_xy.table";
    const string cov_yy_file = prefix + "_cov_yy.table";
    const string var_x_file = prefix + "_var_x.table";
    const string meta_file = prefix + "_meta.table";
    const string x_missing_file = prefix + "_x_missing.table";
    const string y_missing_file = prefix + "_y_missing.table";
    const auto static_start = Clock::now();
    const MatrixXd covariance = read_legacy_covariance(cov_yy_file, variables);
    const double y_quality = use_missing ? read_legacy_y_quality(y_missing_file, variables) : 1.0;
    const std::uint64_t variants = ukc::count_tsv_rows(meta_file);
    if (variants == 0) throw UkcError(meta_file + ": input contains no data rows");
    RegressionModel model(covariance, variables.size() - 1);
    timings.static_data = elapsed_seconds(static_start);
    print_run_summary("NSS v1 legacy text", io_mode, total_sample_size, nullptr,
                      variants, variables.size(), variables, threads, use_missing,
                      (io_mode == "memory" ? variants : kBatchSize) *
                      static_cast<std::uint64_t>(variables.size() + 2) * sizeof(double));

    std::ifstream cov_xy = open_input(cov_xy_file);
    std::ifstream var_x = open_input(var_x_file);
    std::ifstream meta = open_input(meta_file);
    std::ifstream x_missing;
    if (use_missing) x_missing = open_input(x_missing_file);
    string cov_header, var_header, meta_header, x_header;
    if (!std::getline(cov_xy, cov_header)) throw UkcError(cov_xy_file + ": missing header");
    if (!std::getline(var_x, var_header)) throw UkcError(var_x_file + ": missing header");
    if (!std::getline(meta, meta_header)) throw UkcError(meta_file + ": missing header");
    if (use_missing && !std::getline(x_missing, x_header)) throw UkcError(x_missing_file + ": missing header");
    const vector<SelectedColumn> selected = select_columns(parse_header(cov_header, cov_xy_file), variables, cov_xy_file);
    parse_header(var_header, var_x_file);
    parse_header(meta_header, meta_file);
    if (use_missing) parse_header(x_header, x_missing_file);

    AtomicOutputFile output_file(result_filename);
    std::ofstream &output = output_file.stream();
    output << meta_header << " BETA SE T-STAT -log10_P VIF nobs Quality-Score\n";
    ProgressReporter progress(variants);
    size_t invalid_rows = 0;
    std::uint64_t completed = 0;
    while (completed < variants) {
        const size_t requested = io_mode == "memory"
            ? static_cast<size_t>(variants - completed)
            : static_cast<size_t>(std::min<std::uint64_t>(kBatchSize, variants - completed));
        const auto load_start = Clock::now();
        vector<InputRow> batch;
        batch.reserve(requested);
        for (size_t i = 0; i < requested; ++i) {
            InputRow row;
            const bool a = static_cast<bool>(std::getline(cov_xy, row.cov_xy));
            const bool b = static_cast<bool>(std::getline(var_x, row.var_x));
            const bool c = static_cast<bool>(std::getline(meta, row.meta));
            const bool d = use_missing ? static_cast<bool>(std::getline(x_missing, row.x_missing)) : a;
            if (!(a && b && c && d)) {
                throw UkcError("Input row-count mismatch near data row " +
                               std::to_string(completed + i + 1));
            }
            row.line_number = static_cast<size_t>(completed + i + 2);
            batch.push_back(std::move(row));
        }
        timings.loading += elapsed_seconds(load_start);
        vector<RowResult> results;
        const auto compute_start = Clock::now();
        process_legacy_batch(batch, results, selected, model, cov_xy_file,
                             var_x_file, x_missing_file, total_sample_size,
                             overall_non_missing_rate, y_quality, use_missing,
                             threads);
        timings.compute += elapsed_seconds(compute_start);
        vector<string> metadata(batch.size());
        for (size_t i = 0; i < batch.size(); ++i) metadata[i] = std::move(batch[i].meta);
        const auto output_start = Clock::now();
        write_result_rows(output, metadata, results, invalid_rows);
        timings.output += elapsed_seconds(output_start);
        completed += batch.size();
        progress.update(completed);
    }
    string extra;
    if (std::getline(cov_xy, extra) || std::getline(var_x, extra) ||
        std::getline(meta, extra) || (use_missing && std::getline(x_missing, extra))) {
        throw UkcError("Input row-count mismatch after data row " + std::to_string(variants));
    }
    progress.finish();
    const auto commit_start = Clock::now();
    output_file.commit();
    timings.output += elapsed_seconds(commit_start);
    timings.total = elapsed_seconds(total_start);
    if (invalid_rows != 0) {
        std::cerr << "Warning: " << invalid_rows
                  << " SNP rows produced non-finite results because of NA or invalid covariance values\n";
    }
    report_timings(timings, io_mode);
}

void print_banner() {
    std::cout << "-----------------------------------\n"
              << "UKBioCoin (UKC) engine Version: 2.0\n"
              << "Written by: Jing-cheng He\n"
              << "GitHub: https://github.com/Ttttt47/UKBioCoin\n"
              << "Email: jc_he@zju.edu.cn\n"
              << "-----------------------------------\n";
}

}  // namespace

int main(int argc, const char *argv[]) {
    print_banner();
    po::options_description options("\nUsage: Linear regression based on summary data.\n\nOptions:\n");
    options.add_options()
        ("help,h", "help")
        ("use-missing-rate-file", "Legacy compatibility option.")
        ("use-missing-rate-estimate", "Use missing-rate arrays/files to estimate sample size.")
        ("file", po::value<string>(), "NSS v2 directory or NSS v1 file prefix.")
        ("phe", po::value<string>(), "Phenotype name.")
        ("covar", po::value<string>(), "Comma-separated covariate names.")
        ("out", po::value<string>(), "Output file prefix.")
        ("totalsize", po::value<int>(), "Base sample size (v2 defaults to manifest; v1 defaults to 292216).")
        ("overall-non-missing-rate", po::value<double>()->default_value(0.9, "0.9"), "Overall non-missing rate.")
        ("threads", po::value<int>()->default_value(1), "Worker threads (default: 1).")
        ("io-mode", po::value<string>()->default_value("stream"), "I/O mode: stream or memory (default: stream).");

    try {
        po::variables_map arguments;
        po::store(po::parse_command_line(argc, argv, options), arguments);
        po::notify(arguments);
        if (argc == 1) {
            std::cout << "No parameter is given, exiting...\n";
            return EXIT_SUCCESS;
        }
        if (arguments.count("help")) {
            std::cout << options << '\n';
            return EXIT_SUCCESS;
        }
        if (!arguments.count("file") || !arguments.count("phe") ||
            !arguments.count("out")) {
            throw UkcError("Required options: --file, --phe, and --out");
        }
        const string input = arguments["file"].as<string>();
        const string phenotype = arguments["phe"].as<string>();
        const string output = arguments["out"].as<string>();
        const string io_mode = arguments["io-mode"].as<string>();
        const int threads = arguments["threads"].as<int>();
        const double overall = arguments["overall-non-missing-rate"].as<double>();
        const bool use_missing = arguments.count("use-missing-rate-estimate") != 0;
        if (input.empty() || phenotype.empty() || output.empty()) {
            throw UkcError("--file, --phe, and --out must not be empty");
        }
        if (io_mode != "stream" && io_mode != "memory") {
            throw UkcError("--io-mode must be 'stream' or 'memory'");
        }
        if (threads <= 0) throw UkcError("--threads must be greater than zero");
        if (!std::isfinite(overall) || overall <= 0.0 || overall > 1.0) {
            throw UkcError("--overall-non-missing-rate must be in (0, 1]");
        }
#ifndef _OPENMP
        if (threads != 1) throw UkcError("This executable was built without OpenMP; --threads must be 1");
#endif
        vector<string> variables{phenotype};
        if (arguments.count("covar") && !arguments["covar"].as<string>().empty()) {
            const vector<string> covariates = split_covariates(arguments["covar"].as<string>());
            variables.insert(variables.end(), covariates.begin(), covariates.end());
        }
        validate_unique_variables(variables);

        const bool v2 = ukc::directory_exists(input) &&
                        ukc::file_exists(ukc::join_path(input, "manifest.json"));
        const bool total_size_supplied = arguments.count("totalsize") != 0;
        const int supplied_total_size = total_size_supplied
            ? arguments["totalsize"].as<int>() : 292216;
        if (total_size_supplied && supplied_total_size <= 0) {
            throw UkcError("--totalsize must be greater than zero");
        }
        Eigen::setNbThreads(1);
        const string result_filename = output + "_results.table";
        if (v2) {
            run_v2(input, result_filename, variables, total_size_supplied,
                   supplied_total_size, overall, use_missing, io_mode, threads);
        } else {
            const int fixed = static_cast<int>(std::floor(supplied_total_size * overall));
            if (!use_missing && (fixed <= 2 ||
                fixed <= static_cast<int>(variables.size()))) {
                throw UkcError("Effective sample size is too small for the requested model");
            }
            run_v1(input, result_filename, variables, supplied_total_size, overall,
                   use_missing, io_mode, threads);
        }
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

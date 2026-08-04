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
#include <iostream>
#include <limits>
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

#include "tools/cpp_cdfs-master/cdf_chisqt/cdf_base.h"

namespace po = boost::program_options;
using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::size_t;
using std::string;
using std::vector;

namespace {

constexpr size_t kBatchSize = 4096;

class UkcError : public std::runtime_error {
public:
    explicit UkcError(const string &message) : std::runtime_error(message) {}
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

class AtomicOutputFile {
public:
    explicit AtomicOutputFile(const string &final_path)
        : final_path_(final_path),
          temporary_path_(final_path + ".tmp." + std::to_string(static_cast<long long>(::getpid()))),
          stream_(temporary_path_, std::ios::out) {
        if (!stream_) {
            throw UkcError("Error opening " + temporary_path_ + " for output");
        }
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
        if (!stream_) {
            throw UkcError("Error writing " + temporary_path_);
        }
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

bool next_token(const char *&cursor, const char *end, TokenView &token) {
    while (cursor < end && std::isspace(static_cast<unsigned char>(*cursor))) {
        ++cursor;
    }
    if (cursor == end) {
        return false;
    }

    if (*cursor == '"') {
        ++cursor;
        token.begin = cursor;
        while (cursor < end && *cursor != '"') {
            ++cursor;
        }
        if (cursor == end) {
            throw UkcError("Unterminated quoted field");
        }
        token.end = cursor;
        ++cursor;
        if (cursor < end && !std::isspace(static_cast<unsigned char>(*cursor))) {
            throw UkcError("Unexpected character after quoted field");
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

string token_string(const TokenView &token) {
    return string(token.begin, token.end);
}

double parse_number(const TokenView &token,
                    const string &filename,
                    size_t line_number,
                    const string &field_name) {
    if (token.end - token.begin == 2 && token.begin[0] == 'N' && token.begin[1] == 'A') {
        return std::numeric_limits<double>::quiet_NaN();
    }

    errno = 0;
    char *parsed_end = nullptr;
    const double value = std::strtod(token.begin, &parsed_end);
    if (parsed_end != token.end || token.begin == parsed_end || errno == ERANGE) {
        throw UkcError(filename + ":" + std::to_string(line_number) +
                       ": invalid numeric value for " + field_name +
                       " ('" + token_string(token) + "')");
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
            if (name.empty()) {
                continue;
            }
            if (!seen.insert(name).second) {
                throw UkcError("Duplicate column '" + name + "'");
            }
            names.push_back(std::move(name));
        }
    } catch (const UkcError &error) {
        throw UkcError(filename + ":1: " + error.what());
    }
    if (names.empty()) {
        throw UkcError(filename + ":1: empty header");
    }
    return names;
}

vector<SelectedColumn> select_columns(const vector<string> &header,
                                      const vector<string> &requested,
                                      const string &filename) {
    std::unordered_map<string, size_t> index;
    for (size_t i = 0; i < header.size(); ++i) {
        index.emplace(header[i], i);
    }

    vector<SelectedColumn> selected;
    selected.reserve(requested.size());
    vector<string> missing;
    for (size_t i = 0; i < requested.size(); ++i) {
        const auto found = index.find(requested[i]);
        if (found == index.end()) {
            missing.push_back(requested[i]);
        } else {
            selected.push_back(SelectedColumn{found->second, i, requested[i]});
        }
    }
    if (!missing.empty()) {
        std::ostringstream message;
        message << filename << ": phenotypes/covariates not found:";
        for (const auto &name : missing) {
            message << ' ' << name;
        }
        throw UkcError(message.str());
    }
    std::sort(selected.begin(), selected.end(), [](const SelectedColumn &left,
                                                    const SelectedColumn &right) {
        return left.source_index < right.source_index;
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
        if (!next_token(cursor, end, token)) {
            throw UkcError("empty row");
        }

        size_t source_index = 0;
        size_t selected_index = 0;
        while (selected_index < selected.size() && next_token(cursor, end, token)) {
            if (source_index == selected[selected_index].source_index) {
                const auto &column = selected[selected_index];
                values[column.destination_index] =
                    parse_number(token, filename, line_number, column.name);
                ++selected_index;
            }
            ++source_index;
        }
        if (selected_index != selected.size()) {
            throw UkcError("row ends before required column '" +
                           selected[selected_index].name + "'");
        }
    } catch (const UkcError &error) {
        const string prefix = filename + ":" + std::to_string(line_number) + ": ";
        const string message = error.what();
        if (message.compare(0, prefix.size(), prefix) == 0) {
            throw;
        }
        throw UkcError(prefix + message);
    }
    return values;
}

double parse_second_value(const string &line,
                          const string &filename,
                          size_t line_number,
                          const string &field_name) {
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
        const string message = error.what();
        if (message.compare(0, prefix.size(), prefix) == 0) {
            throw;
        }
        throw UkcError(prefix + message);
    }
}

string first_field(const string &line, const string &filename, size_t line_number) {
    const char *cursor = line.data();
    const char *end = cursor + line.size();
    TokenView token;
    try {
        if (!next_token(cursor, end, token)) {
            throw UkcError("empty row");
        }
        return token_string(token);
    } catch (const UkcError &error) {
        throw UkcError(filename + ":" + std::to_string(line_number) + ": " + error.what());
    }
}

std::ifstream open_input(const string &filename) {
    std::ifstream stream(filename, std::ios::in);
    if (!stream) {
        throw UkcError("Error opening " + filename + " for input");
    }
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
        while (left < right && std::isspace(static_cast<unsigned char>(input[left]))) {
            ++left;
        }
        while (right > left && std::isspace(static_cast<unsigned char>(input[right - 1]))) {
            --right;
        }
        if (left == right) {
            throw UkcError("--covar contains an empty name");
        }
        values.emplace_back(input.substr(left, right - left));
        if (comma == string::npos) {
            break;
        }
        start = comma + 1;
    }
    return values;
}

void validate_unique_variables(const vector<string> &variables) {
    std::unordered_set<string> seen;
    for (const auto &name : variables) {
        if (name.empty()) {
            throw UkcError("Phenotype/covariate names must not be empty");
        }
        if (!seen.insert(name).second) {
            throw UkcError("Phenotype/covariate '" + name + "' was specified more than once");
        }
    }
}

MatrixXd read_covariance_matrix(const string &filename, const vector<string> &variables) {
    std::ifstream stream = open_input(filename);
    string header_line;
    if (!std::getline(stream, header_line)) {
        throw UkcError(filename + ": missing header");
    }
    const vector<string> header = parse_header(header_line, filename);
    const vector<SelectedColumn> selected = select_columns(header, variables, filename);

    MatrixXd result = MatrixXd::Zero(static_cast<Eigen::Index>(variables.size()),
                                     static_cast<Eigen::Index>(variables.size()));
    std::unordered_map<string, size_t> requested_rows;
    for (size_t i = 0; i < variables.size(); ++i) {
        requested_rows.emplace(variables[i], i);
    }
    vector<bool> found_rows(variables.size(), false);

    string line;
    size_t line_number = 1;
    while (std::getline(stream, line)) {
        ++line_number;
        const string row_name = first_field(line, filename, line_number);
        const auto requested = requested_rows.find(row_name);
        if (requested == requested_rows.end()) {
            continue;
        }
        if (found_rows[requested->second]) {
            throw UkcError(filename + ":" + std::to_string(line_number) +
                           ": duplicate row '" + row_name + "'");
        }
        const vector<double> row =
            parse_selected_values(line, selected, filename, line_number);
        for (size_t column = 0; column < row.size(); ++column) {
            result(static_cast<Eigen::Index>(requested->second),
                   static_cast<Eigen::Index>(column)) = row[column];
        }
        found_rows[requested->second] = true;
    }

    for (size_t i = 0; i < found_rows.size(); ++i) {
        if (!found_rows[i]) {
            throw UkcError(filename + ": missing covariance row '" + variables[i] + "'");
        }
    }
    if (!result.allFinite()) {
        throw UkcError(filename + ": selected covariance matrix contains NA or non-finite values");
    }
    return result;
}

double read_y_quality_score(const string &filename, const vector<string> &variables) {
    std::ifstream stream = open_input(filename);
    std::unordered_map<string, size_t> requested;
    for (size_t i = 0; i < variables.size(); ++i) {
        requested.emplace(variables[i], i);
    }
    vector<bool> found(variables.size(), false);
    double quality_score = 1.0;
    string line;
    size_t line_number = 0;
    while (std::getline(stream, line)) {
        ++line_number;
        const string name = first_field(line, filename, line_number);
        const auto match = requested.find(name);
        if (match == requested.end()) {
            continue;
        }
        if (found[match->second]) {
            throw UkcError(filename + ":" + std::to_string(line_number) +
                           ": duplicate row '" + name + "'");
        }
        const double missing_rate =
            parse_second_value(line, filename, line_number, "missing_rate");
        if (!std::isfinite(missing_rate) || missing_rate < 0.0 || missing_rate > 1.0) {
            throw UkcError(filename + ":" + std::to_string(line_number) +
                           ": missing_rate must be between 0 and 1");
        }
        quality_score *= (1.0 - missing_rate);
        found[match->second] = true;
    }
    for (size_t i = 0; i < found.size(); ++i) {
        if (!found[i]) {
            throw UkcError(filename + ": missing rate for '" + variables[i] + "'");
        }
    }
    return quality_score;
}

class RegressionModel {
public:
    RegressionModel(const MatrixXd &phenotype_covariance, size_t covariate_count)
        : covariate_count_(covariate_count),
          theta_base_(MatrixXd::Zero(static_cast<Eigen::Index>(covariate_count + 2),
                                     static_cast<Eigen::Index>(covariate_count + 2))) {
        theta_base_(0, 0) = phenotype_covariance(0, 0);
        if (covariate_count_ > 0) {
            theta_base_.block(0, 2, 1, static_cast<Eigen::Index>(covariate_count_)) =
                phenotype_covariance.block(0, 1, 1,
                                             static_cast<Eigen::Index>(covariate_count_));
            theta_base_.block(2, 0, static_cast<Eigen::Index>(covariate_count_), 1) =
                phenotype_covariance.block(1, 0,
                                             static_cast<Eigen::Index>(covariate_count_), 1);
            theta_base_.block(2, 2,
                              static_cast<Eigen::Index>(covariate_count_),
                              static_cast<Eigen::Index>(covariate_count_)) =
                phenotype_covariance.block(1, 1,
                                             static_cast<Eigen::Index>(covariate_count_),
                                             static_cast<Eigen::Index>(covariate_count_));

            const MatrixXd covariance = theta_base_.block(
                2, 2, static_cast<Eigen::Index>(covariate_count_),
                static_cast<Eigen::Index>(covariate_count_));
            Eigen::FullPivLU<MatrixXd> decomposition(covariance);
            if (!decomposition.isInvertible()) {
                throw UkcError("Singularity detected in the covariance matrix; "
                               "check whether covariates are linearly dependent");
            }
            covariate_inverse_ = covariance.inverse();
        }
    }

    std::array<double, 5> calculate(const vector<double> &cov_xy,
                                    double var_x,
                                    int sample_size) const {
        if (covariate_count_ == 0) {
            return calculate_without_covariates(cov_xy[0], var_x, sample_size);
        }

        MatrixXd theta = theta_base_;
        theta(1, 0) = cov_xy[0];
        theta(0, 1) = cov_xy[0];
        theta(1, 1) = var_x;
        for (size_t i = 1; i < cov_xy.size(); ++i) {
            theta(1, static_cast<Eigen::Index>(i + 1)) = cov_xy[i];
            theta(static_cast<Eigen::Index>(i + 1), 1) = cov_xy[i];
        }
        return calculate_with_covariates(theta, sample_size);
    }

    int minimum_valid_sample_size() const {
        return static_cast<int>(covariate_count_ + 1);
    }

private:
    std::array<double, 5> calculate_without_covariates(double cov_xy,
                                                       double var_x,
                                                       int sample_size) const {
        const double omega_inverse = 1.0 / var_x;
        const double beta = omega_inverse * cov_xy;
        const double rss = sample_size * (1.0 - beta * cov_xy) /
                           (sample_size - 1.0);
        const double var_beta = rss * omega_inverse / sample_size;
        const double standard_error = std::sqrt(var_beta);
        const double t_stat = beta / standard_error;
        const double log_p = cdf_t_log(-std::abs(t_stat), sample_size - 2) + std::log(2.0);
        const double log10_p = log_p * std::log10(std::exp(1.0));
        const double vif = var_beta /
                           (rss / ((sample_size - 1.0) * var_x));
        return {{beta, standard_error, t_stat, -log10_p, vif}};
    }

    std::array<double, 5> calculate_with_covariates(const MatrixXd &theta,
                                                    int sample_size) const {
        const Eigen::Index dimension = theta.rows();
        MatrixXd omega = theta.block(1, 1, dimension - 1, dimension - 1);
        MatrixXd omega_inverse = omega;

        const double q = 1.0 /
            (omega(0, 0) -
             (omega.block(0, 1, 1, dimension - 2) * covariate_inverse_ *
              omega.block(1, 0, dimension - 2, 1)).value());
        omega_inverse(0, 0) = q;
        omega_inverse.block(0, 1, 1, dimension - 2) =
            -q * omega.block(0, 1, 1, dimension - 2) * covariate_inverse_;
        omega_inverse.block(1, 0, dimension - 2, 1) =
            -q * covariate_inverse_ * omega.block(1, 0, dimension - 2, 1);
        omega_inverse.block(1, 1, dimension - 2, dimension - 2) =
            covariate_inverse_ +
            q * covariate_inverse_ * omega.block(1, 0, dimension - 2, 1) *
                omega.block(0, 1, 1, dimension - 2) * covariate_inverse_;

        const VectorXd association = theta.block(0, 1, 1, dimension - 1).transpose();
        const VectorXd coefficients = omega_inverse * association;
        const double rss = sample_size *
                           (1.0 - (coefficients.transpose() * association).value()) /
                           (sample_size - (dimension - 2) - 1.0);
        const double var_beta = rss * omega_inverse(0, 0) / sample_size;
        const double standard_error = std::sqrt(var_beta);
        const double t_stat = coefficients[0] / standard_error;
        const double log_p = cdf_t_log(-std::abs(t_stat), sample_size - 2) +
                             std::log(2.0);
        const double log10_p = log_p * std::log10(std::exp(1.0));
        const double vif = var_beta /
                           (rss / ((sample_size - 1.0) * theta(1, 1)));
        return {{coefficients[0], standard_error, t_stat, -log10_p, vif}};
    }

    size_t covariate_count_;
    MatrixXd theta_base_;
    MatrixXd covariate_inverse_;
};

bool contains_nan(const vector<double> &values) {
    for (double value : values) {
        if (std::isnan(value)) {
            return true;
        }
    }
    return false;
}

void process_batch(const vector<InputRow> &batch,
                   vector<RowResult> &results,
                   const vector<SelectedColumn> &selected_columns,
                   const RegressionModel &model,
                   const string &cov_xy_filename,
                   const string &var_x_filename,
                   const string &x_missing_filename,
                   int total_sample_size,
                   double overall_non_missing_rate,
                   double y_quality_score,
                   bool use_missing_rate_estimate,
                   int threads) {
    results.resize(batch.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (long long raw_index = 0; raw_index < static_cast<long long>(batch.size()); ++raw_index) {
        const size_t index = static_cast<size_t>(raw_index);
        const InputRow &row = batch[index];
        RowResult &result = results[index];
        try {
            const vector<double> cov_xy = parse_selected_values(
                row.cov_xy, selected_columns, cov_xy_filename, row.line_number);
            const double var_x = parse_second_value(
                row.var_x, var_x_filename, row.line_number, "var_x");

            double quality_score = overall_non_missing_rate;
            int sample_size = static_cast<int>(std::floor(
                total_sample_size * overall_non_missing_rate));
            if (use_missing_rate_estimate) {
                const double x_missing = parse_second_value(
                    row.x_missing, x_missing_filename, row.line_number,
                    "x_missing");
                if (!std::isfinite(x_missing) || x_missing < 0.0 || x_missing > 1.0) {
                    throw UkcError(x_missing_filename + ":" +
                                   std::to_string(row.line_number) +
                                   ": x_missing must be between 0 and 1");
                }
                quality_score = y_quality_score * (1.0 - x_missing);
                sample_size = static_cast<int>(std::floor(
                    total_sample_size * quality_score));
            }

            result.sample_size = sample_size;
            result.quality_score = quality_score;
            result.invalid = contains_nan(cov_xy) || std::isnan(var_x);
            if (result.invalid || sample_size <= 2 ||
                sample_size <= model.minimum_valid_sample_size()) {
                result.invalid = true;
                continue;
            }
            result.summary = model.calculate(cov_xy, var_x, sample_size);
            for (double value : result.summary) {
                if (!std::isfinite(value)) {
                    result.invalid = true;
                    break;
                }
            }
        } catch (const std::exception &error) {
            result.error = error.what();
        }
    }
}

void calculate_and_save(const string &cov_xy_filename,
                        const string &var_x_filename,
                        const string &meta_filename,
                        const string &result_filename,
                        const string &x_missing_filename,
                        const MatrixXd &phenotype_covariance,
                        const vector<string> &variables,
                        int total_sample_size,
                        double overall_non_missing_rate,
                        double y_quality_score,
                        bool use_missing_rate_estimate,
                        int threads) {
    std::ifstream cov_xy = open_input(cov_xy_filename);
    std::ifstream var_x = open_input(var_x_filename);
    std::ifstream meta = open_input(meta_filename);
    std::ifstream x_missing;
    if (use_missing_rate_estimate) {
        x_missing = open_input(x_missing_filename);
    }

    string cov_xy_header;
    string var_x_header;
    string meta_header;
    string x_missing_header;
    if (!std::getline(cov_xy, cov_xy_header)) {
        throw UkcError(cov_xy_filename + ": missing header");
    }
    if (!std::getline(var_x, var_x_header)) {
        throw UkcError(var_x_filename + ": missing header");
    }
    if (!std::getline(meta, meta_header)) {
        throw UkcError(meta_filename + ": missing header");
    }
    if (use_missing_rate_estimate && !std::getline(x_missing, x_missing_header)) {
        throw UkcError(x_missing_filename + ": missing header");
    }

    const vector<string> cov_xy_columns = parse_header(cov_xy_header, cov_xy_filename);
    parse_header(var_x_header, var_x_filename);
    parse_header(meta_header, meta_filename);
    if (use_missing_rate_estimate) {
        parse_header(x_missing_header, x_missing_filename);
    }
    const vector<SelectedColumn> selected_columns =
        select_columns(cov_xy_columns, variables, cov_xy_filename);
    RegressionModel model(phenotype_covariance, variables.size() - 1);

    AtomicOutputFile output_file(result_filename);
    std::ofstream &output = output_file.stream();
    output << meta_header
           << " BETA SE T-STAT -log10_P VIF nobs Quality-Score\n";

    vector<InputRow> batch;
    batch.reserve(kBatchSize);
    vector<RowResult> results;
    size_t data_line = 0;
    size_t invalid_rows = 0;
    size_t next_progress = 10000;

    while (true) {
        batch.clear();
        for (size_t i = 0; i < kBatchSize; ++i) {
            InputRow row;
            const bool have_cov_xy = static_cast<bool>(std::getline(cov_xy, row.cov_xy));
            const bool have_var_x = static_cast<bool>(std::getline(var_x, row.var_x));
            const bool have_meta = static_cast<bool>(std::getline(meta, row.meta));
            const bool have_x_missing = use_missing_rate_estimate
                ? static_cast<bool>(std::getline(x_missing, row.x_missing))
                : have_cov_xy;

            if (!have_cov_xy && !have_var_x && !have_meta &&
                (!use_missing_rate_estimate || !have_x_missing)) {
                break;
            }
            if (!(have_cov_xy && have_var_x && have_meta && have_x_missing)) {
                throw UkcError("Input row-count mismatch near data row " +
                               std::to_string(data_line + 1));
            }
            ++data_line;
            row.line_number = data_line + 1;
            batch.push_back(std::move(row));
        }

        if (batch.empty()) {
            break;
        }

        process_batch(batch, results, selected_columns, model,
                      cov_xy_filename, var_x_filename, x_missing_filename,
                      total_sample_size, overall_non_missing_rate,
                      y_quality_score, use_missing_rate_estimate, threads);

        for (size_t i = 0; i < batch.size(); ++i) {
            if (!results[i].error.empty()) {
                throw UkcError(results[i].error);
            }
            if (results[i].invalid) {
                ++invalid_rows;
            }
            output << batch[i].meta << ' '
                   << results[i].summary[0] << ' '
                   << results[i].summary[1] << ' '
                   << results[i].summary[2] << ' '
                   << results[i].summary[3] << ' '
                   << results[i].summary[4] << ' '
                   << ' ' << results[i].sample_size << ' '
                   << results[i].quality_score << '\n';
        }

        while (data_line >= next_progress) {
            std::cout << next_progress << '\n';
            next_progress += 10000;
        }
    }

    if (data_line == 0) {
        throw UkcError("Input files contain no data rows");
    }
    output_file.commit();
    if (invalid_rows > 0) {
        std::cerr << "Warning: " << invalid_rows
                  << " SNP rows produced non-finite results because of NA or invalid covariance values\n";
    }
}

void print_banner() {
    std::cout << "-----------------------------------\n"
              << "UKBioCoin (UKC) engine Version: V1.4\n"
              << "Written by: Jing-cheng He\n"
              << "GitHub: https://github.com/Ttttt47/UKBioCoin\n"
              << "Email: jc_he@zju.edu.cn\n"
              << "-----------------------------------\n";
}

}  // namespace

int main(int argc, const char *argv[]) {
    print_banner();

    po::options_description options("\n Usage: Linear regression based on summary data. \n\n Options: \n");
    options.add_options()
        ("help,h", "help")
        ("use-missing-rate-file", "Whether missing rate files are available (legacy compatibility option).")
        ("use-missing-rate-estimate", "Use missing-rate files to estimate sample size.")
        ("file", po::value<string>(), "Prefix of the covariance matrix files.")
        ("phe", po::value<string>(), "Phenotype name matching the covariance matrix header.")
        ("covar", po::value<string>(), "Comma-separated covariate names.")
        ("out", po::value<string>(), "Output file prefix.")
        ("totalsize", po::value<int>()->default_value(292216), "Total sample size.")
        ("overall-non-missing-rate", po::value<double>()->default_value(0.9, "0.9"), "Overall non-missing rate.")
        ("threads", po::value<int>()->default_value(1), "Number of worker threads (default: 1).");

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

        const string input_prefix = arguments["file"].as<string>();
        const string phenotype = arguments["phe"].as<string>();
        const string output_prefix = arguments["out"].as<string>();
        const int total_sample_size = arguments["totalsize"].as<int>();
        const double overall_non_missing_rate =
            arguments["overall-non-missing-rate"].as<double>();
        const int threads = arguments["threads"].as<int>();
        const bool use_missing_rate_estimate =
            arguments.count("use-missing-rate-estimate") > 0;

        if (input_prefix.empty() || phenotype.empty() || output_prefix.empty()) {
            throw UkcError("--file, --phe, and --out must not be empty");
        }
        if (total_sample_size <= 0) {
            throw UkcError("--totalsize must be greater than zero");
        }
        if (!std::isfinite(overall_non_missing_rate) ||
            overall_non_missing_rate <= 0.0 || overall_non_missing_rate > 1.0) {
            throw UkcError("--overall-non-missing-rate must be in (0, 1]");
        }
        if (threads <= 0) {
            throw UkcError("--threads must be greater than zero");
        }
#ifndef _OPENMP
        if (threads != 1) {
            throw UkcError("This executable was built without OpenMP; --threads must be 1");
        }
#endif

        vector<string> variables{phenotype};
        if (arguments.count("covar")) {
            const string covariate_argument = arguments["covar"].as<string>();
            if (!covariate_argument.empty()) {
                const vector<string> covariates = split_covariates(covariate_argument);
                variables.insert(variables.end(), covariates.begin(), covariates.end());
            }
        }
        validate_unique_variables(variables);

        const int fixed_sample_size = static_cast<int>(std::floor(
            total_sample_size * overall_non_missing_rate));
        if (!use_missing_rate_estimate &&
            (fixed_sample_size <= 2 ||
             fixed_sample_size <= static_cast<int>(variables.size()))) {
            throw UkcError("Effective sample size is too small for the requested model");
        }

        std::cout << "Using parameter:\n"
                  << "file = " << input_prefix << '\n'
                  << "phe = " << phenotype << '\n';
        if (variables.size() == 1) {
            std::cout << "No covar is given.\n";
        } else {
            std::cout << "covar = ";
            for (size_t i = 1; i < variables.size(); ++i) {
                std::cout << variables[i] << ' ';
            }
            std::cout << '\n';
        }
        if (use_missing_rate_estimate) {
            std::cout << "using estimated Sample size.\n";
        } else {
            std::cout << "Sample size set as " << fixed_sample_size << '\n';
        }
        std::cout << "threads = " << threads << '\n'
                  << "out = " << output_prefix << '\n'
                  << "Start computing...\n";

        const string cov_xy_filename = input_prefix + "_cov_xy.table";
        const string var_x_filename = input_prefix + "_var_x.table";
        const string cov_yy_filename = input_prefix + "_cov_yy.table";
        const string meta_filename = input_prefix + "_meta.table";
        const string result_filename = output_prefix + "_results.table";
        const string x_missing_filename = input_prefix + "_x_missing.table";
        const string y_missing_filename = input_prefix + "_y_missing.table";

        const auto start = std::chrono::steady_clock::now();
        const MatrixXd phenotype_covariance =
            read_covariance_matrix(cov_yy_filename, variables);
        double y_quality_score = 1.0;
        if (use_missing_rate_estimate) {
            y_quality_score = read_y_quality_score(y_missing_filename, variables);
        }

        Eigen::setNbThreads(1);
        calculate_and_save(cov_xy_filename, var_x_filename, meta_filename,
                           result_filename, x_missing_filename,
                           phenotype_covariance, variables, total_sample_size,
                           overall_non_missing_rate, y_quality_score,
                           use_missing_rate_estimate, threads);
        const auto end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = end - start;
        std::cout << "Completed in " << elapsed.count() << " seconds.\n";
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

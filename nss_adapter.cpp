#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <sys/resource.h>
#include <unistd.h>

#include "nss_v2.hpp"

namespace fs = boost::filesystem;
namespace po = boost::program_options;
using std::size_t;
using std::string;
using std::vector;

namespace {

constexpr size_t kBatchSize = 4096;

struct Token {
    const char *begin = nullptr;
    const char *end = nullptr;
};

bool next_token(const char *&cursor, const char *end, Token &token) {
    while (cursor < end && std::isspace(static_cast<unsigned char>(*cursor))) ++cursor;
    if (cursor == end) return false;
    if (*cursor == '"') {
        token.begin = ++cursor;
        while (cursor < end && *cursor != '"') ++cursor;
        if (cursor == end) throw ukc::NssError("unterminated quoted field");
        token.end = cursor++;
        if (cursor < end && !std::isspace(static_cast<unsigned char>(*cursor))) {
            throw ukc::NssError("unexpected character after quoted field");
        }
        return true;
    }
    token.begin = cursor;
    while (cursor < end && !std::isspace(static_cast<unsigned char>(*cursor))) ++cursor;
    token.end = cursor;
    return true;
}

string token_string(const Token &token) { return string(token.begin, token.end); }

double token_number(const Token &token, const string &path, size_t line) {
    if (token.end - token.begin == 2 && token.begin[0] == 'N' && token.begin[1] == 'A') {
        return std::numeric_limits<double>::quiet_NaN();
    }
    errno = 0;
    char *parsed_end = nullptr;
    const double value = std::strtod(token.begin, &parsed_end);
    if (parsed_end == token.begin || parsed_end != token.end || errno == ERANGE) {
        throw ukc::NssError(path + ":" + std::to_string(line) +
                            ": invalid numeric value '" + token_string(token) + "'");
    }
    return value;
}

vector<string> tokens(const string &line, const string &path, size_t line_number) {
    vector<string> result;
    const char *cursor = line.data();
    const char *end = cursor + line.size();
    Token token;
    try {
        while (next_token(cursor, end, token)) result.push_back(token_string(token));
    } catch (const std::exception &error) {
        throw ukc::NssError(path + ":" + std::to_string(line_number) + ": " + error.what());
    }
    return result;
}

vector<string> read_header(std::ifstream &input, const string &path) {
    string line;
    if (!std::getline(input, line)) throw ukc::NssError(path + ": missing header");
    vector<string> names = tokens(line, path, 1);
    if (names.empty()) throw ukc::NssError(path + ": empty header");
    std::unordered_set<string> unique;
    for (const auto &name : names) {
        if (name.empty() || !unique.insert(name).second) {
            throw ukc::NssError(path + ": phenotype names must be non-empty and unique");
        }
    }
    return names;
}

std::ifstream input_file(const string &path) {
    std::ifstream input(path);
    if (!input) throw ukc::NssError("Error opening " + path + " for input");
    return input;
}

vector<double> read_second_column(const string &path) {
    std::ifstream input = input_file(path);
    string line;
    if (!std::getline(input, line)) throw ukc::NssError(path + ": missing header");
    vector<double> values;
    size_t line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        const char *cursor = line.data();
        const char *end = cursor + line.size();
        Token token;
        try {
            if (!next_token(cursor, end, token) || !next_token(cursor, end, token)) {
                throw ukc::NssError("row has fewer than two fields");
            }
            values.push_back(token_number(token, path, line_number));
            if (next_token(cursor, end, token)) {
                throw ukc::NssError("row has more than two fields");
            }
        } catch (const std::exception &error) {
            const string prefix = path + ":" + std::to_string(line_number) + ": ";
            if (string(error.what()).compare(0, prefix.size(), prefix) == 0) throw;
            throw ukc::NssError(prefix + error.what());
        }
    }
    if (values.empty()) throw ukc::NssError(path + ": no data rows");
    return values;
}

vector<double> read_cov_yy(const string &path, const vector<string> &names) {
    std::ifstream input = input_file(path);
    const vector<string> header = read_header(input, path);
    if (header != names) throw ukc::NssError(path + ": header does not match cov_xy phenotype order");
    const size_t width = names.size();
    vector<double> matrix(width * width);
    vector<bool> found(width, false);
    std::unordered_map<string, size_t> indices;
    for (size_t i = 0; i < width; ++i) indices.emplace(names[i], i);
    string line;
    size_t line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        const char *cursor = line.data();
        const char *end = cursor + line.size();
        Token token;
        if (!next_token(cursor, end, token)) throw ukc::NssError(path + ": empty row");
        const string row_name = token_string(token);
        const auto match = indices.find(row_name);
        if (match == indices.end()) throw ukc::NssError(path + ": unknown row '" + row_name + "'");
        if (found[match->second]) throw ukc::NssError(path + ": duplicate row '" + row_name + "'");
        for (size_t column = 0; column < width; ++column) {
            if (!next_token(cursor, end, token)) throw ukc::NssError(path + ": short covariance row '" + row_name + "'");
            matrix[match->second * width + column] = token_number(token, path, line_number);
        }
        if (next_token(cursor, end, token)) throw ukc::NssError(path + ": extra covariance field");
        found[match->second] = true;
    }
    for (size_t i = 0; i < width; ++i) {
        if (!found[i]) throw ukc::NssError(path + ": missing covariance row '" + names[i] + "'");
    }
    return matrix;
}

vector<double> read_y_missing(const string &path, const vector<string> &names) {
    std::ifstream input = input_file(path);
    string line;
    if (!std::getline(input, line)) throw ukc::NssError(path + ": missing header");
    std::unordered_map<string, size_t> indices;
    for (size_t i = 0; i < names.size(); ++i) indices.emplace(names[i], i);
    vector<double> values(names.size(), std::numeric_limits<double>::quiet_NaN());
    vector<bool> found(names.size(), false);
    size_t line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        const char *cursor = line.data();
        const char *end = cursor + line.size();
        Token name_token, value_token;
        if (!next_token(cursor, end, name_token) || !next_token(cursor, end, value_token)) {
            throw ukc::NssError(path + ": malformed missing-rate row");
        }
        const auto match = indices.find(token_string(name_token));
        if (match == indices.end()) continue;
        if (found[match->second]) throw ukc::NssError(path + ": duplicate missing-rate row");
        values[match->second] = token_number(value_token, path, line_number);
        found[match->second] = true;
    }
    for (size_t i = 0; i < names.size(); ++i) {
        if (!found[i]) throw ukc::NssError(path + ": missing rate for '" + names[i] + "'");
    }
    return values;
}

void copy_metadata_as_tsv(const string &source, const string &destination,
                          std::uint64_t expected_rows) {
    std::ifstream input = input_file(source);
    std::ofstream output(destination);
    if (!output) throw ukc::NssError("Error opening " + destination + " for output");
    string line;
    std::uint64_t rows = 0;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const vector<string> fields = tokens(line, source, line_number);
        if (fields.empty()) throw ukc::NssError(source + ": empty metadata row");
        for (size_t i = 0; i < fields.size(); ++i) {
            output << (i == 0 ? "" : "\t") << fields[i];
        }
        output << '\n';
        if (line_number > 1) ++rows;
    }
    if (rows != expected_rows) {
        throw ukc::NssError(source + ": metadata row count " + std::to_string(rows) +
                            " does not match var_x length " + std::to_string(expected_rows));
    }
}

bool host_little() {
    const std::uint16_t value = 1;
    return *reinterpret_cast<const unsigned char *>(&value) == 1;
}

class OutputPool {
public:
    explicit OutputPool(size_t capacity) : capacity_(std::max<size_t>(1, capacity)) {}

    std::ofstream &get(const string &path) {
        auto found = streams_.find(path);
        if (found != streams_.end()) {
            order_.erase(found->second.order);
            order_.push_front(path);
            found->second.order = order_.begin();
            return *found->second.stream;
        }
        if (streams_.size() == capacity_) {
            const string evicted = order_.back();
            order_.pop_back();
            streams_.erase(evicted);
        }
        std::unique_ptr<std::ofstream> output(new std::ofstream(
            path, std::ios::binary | std::ios::app));
        if (!*output) throw ukc::NssError("Error opening " + path + " for append");
        order_.push_front(path);
        Entry entry;
        entry.stream = std::move(output);
        entry.order = order_.begin();
        auto inserted = streams_.emplace(path, std::move(entry));
        return *inserted.first->second.stream;
    }

    void write(const string &path, const vector<float> &values) {
        std::ofstream &output = get(path);
        if (host_little()) {
            output.write(reinterpret_cast<const char *>(values.data()),
                         static_cast<std::streamsize>(values.size() * sizeof(float)));
        } else {
            for (float value : values) {
                unsigned char *bytes = reinterpret_cast<unsigned char *>(&value);
                std::reverse(bytes, bytes + sizeof(float));
                output.write(reinterpret_cast<const char *>(&value), sizeof(float));
            }
        }
        if (!output) throw ukc::NssError("Error writing " + path);
    }

    void close() { streams_.clear(); order_.clear(); }

private:
    struct Entry {
        std::unique_ptr<std::ofstream> stream;
        std::list<string>::iterator order;
    };
    size_t capacity_;
    std::list<string> order_;
    std::unordered_map<string, Entry> streams_;
};

size_t file_pool_capacity() {
    struct rlimit limit {};
    if (getrlimit(RLIMIT_NOFILE, &limit) != 0) return 64;
    const std::uint64_t safe = limit.rlim_cur == RLIM_INFINITY
        ? 128 : std::max<std::uint64_t>(1, limit.rlim_cur / 4);
    return static_cast<size_t>(std::min<std::uint64_t>(128, safe));
}

void parse_cov_row(const string &line, const string &path, size_t line_number,
                   size_t width, vector<float> &values, bool &label_mismatch) {
    const char *cursor = line.data();
    const char *end = cursor + line.size();
    Token token;
    try {
        if (!next_token(cursor, end, token)) throw ukc::NssError("empty row");
        char *label_end = nullptr;
        const string label = token_string(token);
        const unsigned long long row_label = std::strtoull(label.c_str(), &label_end, 10);
        label_mismatch = label_end == label.c_str() || *label_end != '\0' ||
                         row_label != line_number - 1;
        values.resize(width);
        for (size_t column = 0; column < width; ++column) {
            if (!next_token(cursor, end, token)) throw ukc::NssError("row ends before all phenotype fields");
            values[column] = static_cast<float>(token_number(token, path, line_number));
        }
        if (next_token(cursor, end, token)) throw ukc::NssError("row has extra phenotype fields");
    } catch (const std::exception &error) {
        const string prefix = path + ":" + std::to_string(line_number) + ": ";
        if (string(error.what()).compare(0, prefix.size(), prefix) == 0) throw;
        throw ukc::NssError(prefix + error.what());
    }
}

string indexed_cov_path(size_t index) {
    std::ostringstream name;
    name << "cov_xy/" << std::setw(8) << std::setfill('0') << index << ".npy";
    return name.str();
}

std::uint64_t directory_size(const fs::path &directory) {
    std::uint64_t bytes = 0;
    for (fs::recursive_directory_iterator it(directory), end; it != end; ++it) {
        if (fs::is_regular_file(it->path())) bytes += fs::file_size(it->path());
    }
    return bytes;
}

void safe_remove_output(const fs::path &path) {
    const fs::path absolute = fs::absolute(path).lexically_normal();
    if (absolute == absolute.root_path() || absolute.filename().empty()) {
        throw ukc::NssError("Refusing to replace unsafe output path " + absolute.string());
    }
    fs::remove_all(absolute);
}

void convert(const string &prefix, const string &output_path,
             std::uint64_t sample_count, int threads, bool force) {
    const auto start = std::chrono::steady_clock::now();
    const string cov_xy_path = prefix + "_cov_xy.table";
    const string cov_yy_path = prefix + "_cov_yy.table";
    const string var_x_path = prefix + "_var_x.table";
    const string meta_path = prefix + "_meta.table";
    const string x_missing_path = prefix + "_x_missing.table";
    const string y_missing_path = prefix + "_y_missing.table";
    const fs::path output(output_path);
    const fs::path temporary(output_path + ".tmp." +
                             std::to_string(static_cast<long long>(::getpid())));
    if (fs::exists(output) && !force) {
        throw ukc::NssError(output_path + " already exists; use --force to replace it");
    }
    if (fs::exists(temporary)) safe_remove_output(temporary);

    try {
        std::ifstream cov_xy = input_file(cov_xy_path);
        const vector<string> names = read_header(cov_xy, cov_xy_path);
        const vector<double> var_x = read_second_column(var_x_path);
        const std::uint64_t variants = var_x.size();
        const vector<double> cov_yy = read_cov_yy(cov_yy_path, names);
        const bool have_x_missing = ukc::file_exists(x_missing_path);
        const bool have_y_missing = ukc::file_exists(y_missing_path);
        if (have_x_missing != have_y_missing) {
            throw ukc::NssError("Legacy input must provide both x/y missing-rate files or neither");
        }
        vector<double> x_missing;
        vector<double> y_missing;
        if (have_x_missing) {
            x_missing = read_second_column(x_missing_path);
            if (x_missing.size() != variants) throw ukc::NssError(x_missing_path + ": row count does not match var_x");
            y_missing = read_y_missing(y_missing_path, names);
        }

        fs::create_directories(temporary / "cov_xy");
        copy_metadata_as_tsv(meta_path, (temporary / "meta.tsv").string(), variants);
        ukc::write_npy_f64((temporary / "var_x.npy").string(), var_x, {variants});
        ukc::write_npy_f64((temporary / "cov_yy.npy").string(), cov_yy,
                           {static_cast<std::uint64_t>(names.size()),
                            static_cast<std::uint64_t>(names.size())});
        if (have_x_missing) {
            ukc::write_npy_f64((temporary / "x_missing.npy").string(), x_missing, {variants});
            ukc::write_npy_f64((temporary / "y_missing.npy").string(), y_missing,
                               {static_cast<std::uint64_t>(names.size())});
        }
        for (size_t column = 0; column < names.size(); ++column) {
            ukc::create_npy_f32((temporary / indexed_cov_path(column)).string(), variants);
        }

        OutputPool pool(file_pool_capacity());
        std::uint64_t completed = 0;
        std::uint64_t mismatched_labels = 0;
        const bool terminal = ::isatty(STDOUT_FILENO) != 0;
        unsigned next_logged_percent = 0;
        while (completed < variants) {
            const size_t count = static_cast<size_t>(std::min<std::uint64_t>(
                kBatchSize, variants - completed));
            vector<string> lines(count);
            for (size_t i = 0; i < count; ++i) {
                if (!std::getline(cov_xy, lines[i])) {
                    throw ukc::NssError(cov_xy_path + ": row count is shorter than var_x at row " +
                                        std::to_string(completed + i + 1));
                }
            }
            vector<vector<float>> rows(count);
            vector<unsigned char> mismatch(count, 0);
            vector<string> errors(count);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
            for (long long raw = 0; raw < static_cast<long long>(count); ++raw) {
                const size_t i = static_cast<size_t>(raw);
                try {
                    bool bad_label = false;
                    parse_cov_row(lines[i], cov_xy_path,
                                  static_cast<size_t>(completed + i + 2),
                                  names.size(), rows[i], bad_label);
                    mismatch[i] = bad_label ? 1 : 0;
                } catch (const std::exception &error) {
                    errors[i] = error.what();
                }
            }
            for (size_t i = 0; i < count; ++i) {
                if (!errors[i].empty()) throw ukc::NssError(errors[i]);
                mismatched_labels += mismatch[i];
            }
            vector<float> column_values(count);
            for (size_t column = 0; column < names.size(); ++column) {
                for (size_t row = 0; row < count; ++row) column_values[row] = rows[row][column];
                pool.write((temporary / indexed_cov_path(column)).string(), column_values);
            }
            completed += count;
            const unsigned percent = completed >= variants ? 100U :
                static_cast<unsigned>(100.0L * static_cast<long double>(completed) /
                                      static_cast<long double>(variants));
            if (terminal) {
                std::cout << "\rConverting COVxy: " << std::setw(3) << percent
                          << "%" << std::flush;
            } else {
                while (next_logged_percent <= percent) {
                    std::cout << "Converting COVxy: " << next_logged_percent << "%\n";
                    next_logged_percent += 5;
                }
            }
        }
        pool.close();
        if (terminal) std::cout << '\n';
        string extra;
        if (std::getline(cov_xy, extra)) {
            throw ukc::NssError(cov_xy_path + ": row count is longer than var_x");
        }
        if (mismatched_labels != 0) {
            std::cerr << "Warning: ignored " << mismatched_labels
                      << " non-sequential/duplicate legacy COVxy row labels; file order was preserved\n";
        }

        ukc::Manifest manifest;
        manifest.sample_count = sample_count;
        manifest.variant_count = variants;
        manifest.phenotype_count = names.size();
        manifest.meta_path = "meta.tsv";
        manifest.cov_yy = {"cov_yy.npy", "<f8",
            {static_cast<std::uint64_t>(names.size()),
             static_cast<std::uint64_t>(names.size())}, "C"};
        manifest.var_x = {"var_x.npy", "<f8", {variants}, "C"};
        manifest.has_x_missing = have_x_missing;
        manifest.has_y_missing = have_y_missing;
        if (have_x_missing) {
            manifest.x_missing = {"x_missing.npy", "<f8", {variants}, "C"};
            manifest.y_missing = {"y_missing.npy", "<f8",
                {static_cast<std::uint64_t>(names.size())}, "C"};
        }
        for (size_t i = 0; i < names.size(); ++i) {
            ukc::PhenotypeDescriptor phenotype;
            phenotype.index = i;
            phenotype.name = names[i];
            phenotype.cov_xy = {indexed_cov_path(i), "<f4", {variants}, "C"};
            manifest.phenotypes.push_back(std::move(phenotype));
        }
        ukc::write_manifest(temporary.string(), manifest);
        const ukc::Manifest checked = ukc::read_manifest(temporary.string());
        ukc::validate_manifest(temporary.string(), checked, true);
        const std::uint64_t output_bytes = directory_size(temporary);
        std::uint64_t cov_xy_bytes = 0;
        for (size_t i = 0; i < names.size(); ++i) {
            cov_xy_bytes += fs::file_size(temporary / indexed_cov_path(i));
        }
        const std::uint64_t legacy_bytes = fs::file_size(cov_xy_path);

        const fs::path backup(output_path + ".backup." +
                              std::to_string(static_cast<long long>(::getpid())));
        bool moved_existing = false;
        if (fs::exists(backup)) safe_remove_output(backup);
        if (fs::exists(output)) {
            fs::rename(output, backup);
            moved_existing = true;
        }
        try {
            fs::rename(temporary, output);
        } catch (...) {
            if (moved_existing && !fs::exists(output)) fs::rename(backup, output);
            throw;
        }
        if (moved_existing) safe_remove_output(backup);
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        std::cout << "Converted NSS v1 -> v2\n"
                  << "  sample_count: " << sample_count << '\n'
                  << "  variant_count: " << variants << '\n'
                  << "  phenotype_count: " << names.size() << '\n'
                  << "  file-pool capacity: " << file_pool_capacity() << '\n'
                  << "  elapsed: " << std::fixed << std::setprecision(3) << seconds << " s\n"
                  << "  peak RSS: " << static_cast<double>(ukc::peak_rss_kib()) / 1024.0 << " MiB\n"
                  << "  COVxy storage / legacy COVxy: "
                  << (100.0 * static_cast<double>(cov_xy_bytes) / legacy_bytes) << "%\n"
                  << "  v2 total size: " << static_cast<double>(output_bytes) /
                        (1024.0 * 1024.0) << " MiB\n";
    } catch (...) {
        if (fs::exists(temporary)) safe_remove_output(temporary);
        throw;
    }
}

}  // namespace

int main(int argc, const char *argv[]) {
    po::options_description options("UKBioCoinNSSAdapter: convert legacy NSS v1 to NSS v2\n\nOptions");
    options.add_options()
        ("help,h", "help")
        ("input-prefix", po::value<string>(), "Legacy NSS file prefix.")
        ("output", po::value<string>(), "Output .nss directory.")
        ("sample-count", po::value<std::uint64_t>(), "Original base sample count (required).")
        ("threads", po::value<int>()->default_value(1), "Parser threads (default: 1).")
        ("force", "Replace an existing output directory.");
    try {
        po::variables_map arguments;
        po::store(po::parse_command_line(argc, argv, options), arguments);
        po::notify(arguments);
        if (arguments.count("help") || argc == 1) {
            std::cout << options << '\n';
            return EXIT_SUCCESS;
        }
        if (!arguments.count("input-prefix") || !arguments.count("output") ||
            !arguments.count("sample-count")) {
            throw ukc::NssError("Required options: --input-prefix, --output, and --sample-count");
        }
        const int threads = arguments["threads"].as<int>();
        const std::uint64_t sample_count = arguments["sample-count"].as<std::uint64_t>();
        if (threads <= 0) throw ukc::NssError("--threads must be greater than zero");
        if (sample_count == 0) throw ukc::NssError("--sample-count must be a positive integer");
#ifndef _OPENMP
        if (threads != 1) throw ukc::NssError("This executable was built without OpenMP; --threads must be 1");
#endif
        convert(arguments["input-prefix"].as<string>(),
                arguments["output"].as<string>(), sample_count, threads,
                arguments.count("force") != 0);
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

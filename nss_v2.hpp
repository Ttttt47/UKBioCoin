#ifndef UKBIOCOIN_NSS_V2_HPP
#define UKBIOCOIN_NSS_V2_HPP

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ukc {

class NssError : public std::runtime_error {
public:
    explicit NssError(const std::string &message) : std::runtime_error(message) {}
};

struct ArrayDescriptor {
    std::string path;
    std::string dtype;
    std::vector<std::uint64_t> shape;
    std::string order = "C";
};

struct PhenotypeDescriptor {
    std::size_t index = 0;
    std::string name;
    ArrayDescriptor cov_xy;
};

struct Manifest {
    std::uint64_t sample_count = 0;
    std::uint64_t variant_count = 0;
    std::uint64_t phenotype_count = 0;
    std::string meta_path = "meta.tsv";
    ArrayDescriptor cov_yy;
    ArrayDescriptor var_x;
    bool has_x_missing = false;
    ArrayDescriptor x_missing;
    bool has_y_missing = false;
    ArrayDescriptor y_missing;
    std::vector<PhenotypeDescriptor> phenotypes;
};

struct NpyHeader {
    std::string dtype;
    std::vector<std::uint64_t> shape;
    bool fortran_order = false;
    std::uint64_t element_count = 0;
    std::size_t item_size = 0;
    std::streamoff data_offset = 0;
};

bool file_exists(const std::string &path);
bool directory_exists(const std::string &path);
std::string join_path(const std::string &base, const std::string &relative);

Manifest read_manifest(const std::string &directory);
void write_manifest(const std::string &directory, const Manifest &manifest);
void validate_manifest(const std::string &directory,
                       const Manifest &manifest,
                       bool validate_all_cov_xy = true,
                       bool validate_metadata = true);

class NpyReader {
public:
    NpyReader(const std::string &path,
              const std::string &expected_dtype,
              const std::vector<std::uint64_t> &expected_shape);

    const NpyHeader &header() const { return header_; }
    void read(double *destination, std::size_t count);
    std::vector<double> read_all();
    std::uint64_t remaining() const { return header_.element_count - elements_read_; }

private:
    std::string path_;
    std::ifstream stream_;
    NpyHeader header_;
    std::uint64_t elements_read_ = 0;
};

void write_npy_f32(const std::string &path, const std::vector<float> &values);
void write_npy_f64(const std::string &path,
                   const std::vector<double> &values,
                   const std::vector<std::uint64_t> &shape);
void write_npy_f64_header(std::ofstream &output,
                          const std::vector<std::uint64_t> &shape,
                          const std::string &path);
void write_npy_f64_payload(std::ofstream &output,
                           const double *values,
                           std::size_t count,
                           const std::string &path);
void create_npy_f32(const std::string &path, std::uint64_t count);
void append_f32_payload(const std::string &path, const std::vector<float> &values);

std::uint64_t count_tsv_rows(const std::string &path);
std::uint64_t peak_rss_kib();

}  // namespace ukc

#endif

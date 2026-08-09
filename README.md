# UKBioCoin

UKBioCoin (UKC) is a computational device that enables GWAS for UKB but does not rely on UKB original data.

Previously, UKB genotype files (ukb_geno), UKB phenotype files (ukb_phe), and UKB covariate files (ukb_cov) are needed to run
```
plink --bfile ukb_geno --phe ukb_phe --cov ukc_cov 2 --out out1
plink --bfile ukb_geno --phe ukb_phe --cov ukb_cov 1 --out out2
plink --bfile ukb_geno --phe ukb_phe --cov ukb_cov 1 2 --out out1_2
```
Now, you don't need these three files but still have the flexibility to build your own UKB model with UKBioCoin! 

Obviously, these outputs can be integrated into the routines for GWAS, polygenic genetic score, Mendelian randomization, and SNP-heritability estimation.

As UKC is developed on summary statistics and advanced programming technique, it is about 100 times faster than the GWAS on the original UKBioBank data.


A step-by-step example can be found in `./demo`.




## Using UKBioCoin

### Docker images

We provide prebuilt Docker images that package the UKBioCoin computational engine together with the Naive Summary Statistics (NSS) generated from UK Biobank (UKB) and Westlake Biobank for Chinese (WBBC). These images allow users to run covariate-adjustable GWAS without direct access to individual-level data.

> **Docker compatibility notice:** The Docker images listed below still contain the old UKC core. They do not support NSS v2, `--io-mode`, `--threads`, or `--output-format`, and they require at least one covariate via `--covar`. NSS v2 and the current runtime options require UKC 2.0 built from source.

A full description of all UKB phenotypes included in the NSS is available in the file `Description of phenotype.xlsx` in this repository. We also provide the first 20 UKC principal components (accessible as `PC{n}`, where `n` ranges from 1–20) and the first 20 UKB PCs (`ukbPC{n}`). When specifying phenotypes/covariates, a UKB Field ID 34 should be written as X34.0.0 when used in the command line. 

For WBBC, the included NSS cover height, sex, age, and the first five principal components. When specifying WBBC phenotypes or covariates in UKBioCoin, simply use:
`height`, `sex`, `age`, and `PC1–PC5`.

To use the images, first log in to the repository using the read-only account below:

```bash
docker login -u cn-east-3@CFB4B8RZY5ZZ5RYWIHRI -p \
3f70cdbb9813cc57ad3a4b08de08217bcd5575ba09daef2b41b6388f0f336bc1 \
swr.cn-east-3.myhuaweicloud.com
```

Then pull the full image containing the UKBioCoin algorithm and the NSS data using:

```bash
docker pull swr.cn-east-3.myhuaweicloud.com/ukbiocoin/ukbiocoin_full:v1.3_sc02_wbbc
```

The software is stored in `/UKB` and the NSS is stored in `/NSS`. Note that the NSS data are compressed. You should first decompress them using:
```bash
tar -I 'zstd -v' -xvf XXX.tar.zst
```

To pull the image containing only the UKBioCoin algorithm, use:

```bash
docker pull swr.cn-east-3.myhuaweicloud.com/ukbiocoin/ukbiocoin_none:v1.3
```

#### Security note

The provided Docker login credentials are for a dedicated, read-only account created specifically for distribution. This account has no permission to modify or delete any repository content.

To ensure the exact version of the environment and NSS data used in our manuscript is reproducible, we provide the immutable image digest below. Users can pull this specific version to guarantee content integrity, as the digest will change if any part of the image is modified:

```text
swr.cn-east-3.myhuaweicloud.com/ukbiocoin/ukbiocoin_full@sha256:9ac4148b6ec997d312fb4d56918bde36c3cef473146fbcc307f3fd3457aa3f47
```

### Running the updated UKBioCoin core from the command line

The commands in this section require UKBioCoin to be built from the current source code; they are not supported by the Docker images above.

Example:

```bash
UKBioCoin --file test_data/sam.nss \
          --phe X31.0.0 \
          --covar X1160.0.0,X1200.0.0,X1289.0.0,PC1,PC2,PC3,PC4,PC5 \
          --out test_data/test \
          --threads 16 \
          --io-mode stream \
          --use-missing-rate-estimate
```

- `--file`: Specifies either an NSS v2 directory containing `manifest.json`, or a legacy NSS v1 file prefix. UKC detects the format automatically.

- `--phe`: Specifies the phenotype to be analyzed. In this example, it is set to `X31.0.0`.

- `--covar`: Specifies comma-separated covariates. In this example, they are `X1160.0.0`, `X1200.0.0`, `X1289.0.0`, and `PC1`–`PC5`. Omit this option to run without covariates.

- `--out`: Specifies the output prefix. Depending on `--output-format`, UKC writes either `<out>_results.table` or `<out>_results.npy`.

- `--totalsize`: For NSS v2 this is normally omitted because the base sample size is read from `manifest.json`. If supplied, it must exactly match `dimensions.sample_count`. Legacy NSS v1 retains the historical default of `292216`.

- `--use-missing-rate-estimate`: Uses x/y missing-rate arrays to estimate per-SNP sample size. For NSS v2, both optional arrays must be declared in the manifest; for legacy NSS, UKC reads `xxx_x_missing.table` and `xxx_y_missing.table`.

- `--threads`: Number of worker threads used for parsing, regression, and table-result formatting. The default is `1`, which preserves the resource usage of earlier releases. Results are written in input order and are deterministic across thread counts.

- `--io-mode`: `stream` (default) reads the selected vectors in fixed-size blocks and has bounded memory use. `memory` loads only the requested phenotype/covariate vectors and shared arrays before computation, then reports loading, computation, and writing time separately.

- `--output-format`: Result format, either `table` (default) or `npy`. The `npy`
  mode writes only the numeric regression results to `<out>_results.npy`; it does
  not read, validate, or join SNP metadata.

To run a regression without covariates, omit `--covar`:

```bash
UKBioCoin --file test_data/sam.nss \
          --phe X31.0.0 \
          --out test_data/test.no_covar \
          --threads 16
```


When `--use-missing-rate-estimate` is not used, the following option controls the fixed sample-size estimate:

- `--overall-non-missing-rate`: The overall non-missing rate of the data. If `--use-missing-rate-estimate` isn't specified, the sample size will be set to the value of `--totalsize` times value of parameter `--overall-non-missing-rate` (default is 0.9).

### Output formats

#### Table output

The default result is a table containing the estimated coefficients and associated statistics for every SNP. Its row order matches the input NSS.

Result rows are formatted in parallel into bounded per-batch buffers and then written sequentially, preserving the historical table bytes, metadata, numeric precision, and input order. The timing summary reports formatting and file-writing time separately while retaining the combined `result writing` value.

Example:

```text
#CHROM ID POS Effect_Allele REF_Allele ALT_Allele REF_FREQ BETA SE T-STAT -log10_P VIF nobs Quality-Score
1 rs568927457 13453 C C T 0.00119 0.066221 0.0394034 1.68059 1.03225 0.99 262994 0.9
1 rs199856693 14933 A A G 0.01211 0.0230924 0.0124185 1.85952 1.20097 0.99 262994 0.9
1 rs533630043 15585 A A G 0.00192 -0.00582162 0.0310293 -0.187617 0.0699801 0.99 262994 0.9
1 rs568149713 15777 G G A 0.00331 -0.0474425 0.0236265 -2.00802 1.35025 0.99 262994 0.9
...
```

The leading columns are copied from the input metadata file; their number and meaning depend on that file. UKC appends the following result columns:

- `BETA`: This column contains the estimated regression coefficients for each predictor variable.

- `SE`: This column contains the standard error estimates for each coefficient.

- `T-STAT`: This column contains the t-statistic values for each coefficient, which indicate the strength of the evidence against the null hypothesis that the true coefficient is zero.

- `-log10_P`: This column contains the negative base-10 logarithm of the p-value for each coefficient. Larger values indicate stronger evidence against the null hypothesis.

- `VIF`: This column contains the variance inflation factor for the SNP coefficient.

- `nobs`: This column contains the estimated number of observations used in the regression.

- `Quality-Score`: This column contains the quality score of the regression, which is the product of non-missing rates of SNP, covariates and phenotype.

#### Numeric NPY output

For faster numeric-only output, use:

```bash
UKBioCoin --file test_data/sam.nss \
          --phe X31.0.0 \
          --out test_data/test.numeric \
          --threads 16 \
          --output-format npy
```

This creates `test_data/test.numeric_results.npy`, a standard C-order,
little-endian `float64` matrix with shape `(number_of_variants, 7)`. Its columns
are, in order, `BETA`, `SE`, `T-STAT`, `-log10_P`, `VIF`, `nobs`, and
`Quality-Score`. No column-name row or variant metadata is stored in this
array. Invalid numeric result rows retain `NaN` values. Because metadata is not
accessed in this mode, the input `meta.tsv` or legacy `_meta.table` file may be
absent.

### Performance recommendation

For a single large analysis whose downstream workflow only needs numeric
results, prefer `--io-mode stream --output-format npy`. Stream mode reads the
selected vectors sequentially in bounded batches, while NPY output avoids both
metadata parsing and floating-point-to-text conversion. This combination keeps
memory use nearly constant and was the fastest tested configuration.

Use `table` when downstream tools require SNP IDs, alleles, positions, or other
metadata in the same result file. Use `memory` only when a specialized workflow
benefits from resident vectors: in the current single-pass engine it usually
provides no speed benefit and can consume hundreds of MiB or more.

In a synthetic benchmark with 4,000,000 variants, one phenotype, seven
covariates, and eight threads, `stream + npy` took 1.238 seconds and 4.6 MiB peak
RSS. In comparison, `stream + table` took 3.048 seconds, while `memory + npy`
took 2.340 seconds and 661 MiB peak RSS. Absolute timings depend on the hardware
and filesystem. See the full [benchmark results](tests/benchmark_results.md#numeric-npy-result-output-benchmark)
for the test design and single-thread measurements.




## Building UKBioCoin with your own dataset

### Building the UKC executable

To use the latest UKC core, build the `UKBioCoin` executable from the current source code.

#### Reproducible development environment

The supplied environment includes R, the R packages used by `demo/script.R`, NumPy (for interoperability tests), PLINK2, Boost, CMake, Ninja, and a C++ compiler. Eigen is intentionally not duplicated and is supplied to CMake as an include path.

```bash
mamba env create -f environment.yml
mamba activate ukbiocoin-dev
```

#### Building and testing UKBioCoin

UKBioCoin requires C++14. If Eigen is installed in a standard system location,
CMake will find it automatically. Otherwise, provide the directory containing
`Eigen/Dense` through `EIGEN3_INCLUDE_DIR`:

```bash
cmake -S . -B build -G Ninja \
  -DEIGEN3_INCLUDE_DIR=/path/to/eigen3
cmake --build build
ctest --test-dir build --output-on-failure
```

The executables are written to `build/UKBioCoin` and `build/UKBioCoinNSSAdapter`.

### Generating Naive Summary Statistics (NSS) for UKBioCoin

UKBioCoin 2.0 uses NSS v2 by default. It stores standard little-endian, C-order NumPy arrays and splits COVxy into one float32 vector per phenotype. COVyy, var_x, and missing-rate arrays remain float64. This avoids parsing unrelated phenotype values and substantially reduces COVxy storage. Legacy text NSS v1 remains a supported input format.

We provide a script at `demo/script.R` for you to generate the NSS files from the original PLINK's bfile/pfile and the phenotype(covariates) data. To use this script, you should make sure that you have PLINK/PLINK2 installed in your system and the `plink`/`plink2` command is in your `PATH`.

You may first test the script with the provided demo data at `demo/` folder with the following command:

```bash
Rscript script.R \
  --bfile euro_10Ksnp \
  --pheno euro_10Ksnp.all.phe \
  --threads 16 \
  --memory 8000 \
  --prefix euro_10Ksnp \
  --ukc /path/to/UKBioCoin \
  --novisualize
```

The default output is `2.nss/euro_10Ksnp.nss`. Use `--nss-format legacy` to generate the old files under `2.matrix`, and use `--force` only when an existing NSS output should be replaced. The generator obtains the base sample count from the input `.fam` or `.psam`, validates the SNP ID/order across PLINK outputs, and writes each COVxy vector immediately instead of building a complete SNP-by-phenotype matrix in memory.

For additional options, run:

```bash
Rscript script.R --help
```

### NSS v2 layout and manifest

An NSS v2 directory has the following layout:

```text
prefix.nss/
├── manifest.json
├── meta.tsv
├── cov_yy.npy
├── var_x.npy
├── x_missing.npy          # optional
├── y_missing.npy          # optional
└── cov_xy/
    ├── 00000000.npy
    ├── 00000001.npy
    └── ...
```

`manifest.json` uses schema `org.ukbiocoin.nss`, schema version `2`, and records `sample_count`, `variant_count`, `phenotype_count`, every dtype/shape/path, and the stable zero-based index and original name of every phenotype. All three dimensions are required positive integers. UKC validates the manifest, all NPY headers and payload sizes, every selected variable, and optional missing-rate arrays. The metadata file and its row count are additionally validated for the default table output; numeric-only NPY output does not access it.

NPY files are uncompressed standard NumPy files. Each COVxy vector has dtype `<f4` and shape `[variant_count]`; COVyy, var_x, x-missing, and y-missing use `<f8`. Missing values are IEEE NaN. The regression itself always computes in double precision.

### Converting legacy NSS v1

The adapter converts an existing text NSS without loading its full COVxy matrix. `--sample-count` is required because this value cannot be inferred from legacy NSS files.

```bash
UKBioCoinNSSAdapter \
  --input-prefix /path/to/legacy_prefix \
  --output /path/to/prefix.nss \
  --sample-count 292216 \
  --threads 16
```

The adapter writes to a temporary directory, validates the completed manifest and all arrays, and atomically renames it. It refuses to overwrite an output unless `--force` is given. Legacy row labels are treated as advisory so files created by older versions with duplicated labels remain convertible; physical row order is preserved and a warning is reported.

### Legacy NSS v1 input formats

You may also generate the NSS files by yourself. The following is the detailed description of the input formats.

#### xxx_meta.table

For the `xxx_meta.table` format, the first row lists all the metadata fields' names. These metadata are information about the SNPs to be added to the collunms of the result file, So the number of rows of it should be the same with `xxx_var_x.table`(which is equal to the number of SNPs). You can tailor them with your interests.

For example:

```
#CHROM ID POS REF_Allele ALT_Allele REF_FREQ
1 rs568927457 13453 C T 0.001190
1 rs554760071 13483 C G 0.001235
1 rs199856693 14933 A G 0.012111
1 rs533630043 15585 A G 0.001922
...
```

#### xxx_cov_xy.table

For the `xxx_cov_xy.table` format, the first row lists all the phenotypes, including principal components, as column names. Each subsequent row represents a single SNP and contains covariances between that SNP and all the phenotypes. The first column in each row is the row number, enclosed in quotes, and the remaining columns contain the covariance values.

For example:

```
"phenotype1" "phenotype2" "phenotype3" "PC1" "PC2"
"1" 0.245 0.168 0.322 0.054 0.003
"2" 0.082 0.316 0.154 0.092 0.007
"3" -0.015 -0.062 -0.042 -0.012 0.001
...
```

#### xxx_cov_yy.table

For the `xxx_cov_yy.table` format, the first row lists all the phenotypes, including principal components, as column headers. Each subsequent row represents a single phenotype and contains covariance values between that phenotype and itself. The first column in each row lists the name of the phenotype, enclosed in quotes, and the remaining columns contain the covariance values.

For example:

```
"phenotype1" "phenotype2" "phenotype3" "PC1" "PC2"
"phenotype1" 0.374 0.129 0.521 0.045 0.003
"phenotype2" 0.129 0.238 0.187 0.009 0.002
"phenotype3" 0.521 0.187 0.606 0.027 0.004
"PC1" 0.045 0.009 0.027 0.121 0.002
"PC2" 0.003 0.002 0.004 0.002 0.041
```

#### xxx_var_x.table

For the `xxx_var_x.table` format, the first row is an arbitrary column name, such as "`var_x`". The second row lists the variance of each SNP, with each row representing a single SNP. The first column in each row is the row number, enclosed in quotes, and the second column contains the variance value.

For example:
```
"var_x"
"1" 0.153
"2" 0.214
"3" 0.082
...
```

## Access statement
The NSS provided with UKC consist solely of summary statistics (SNP–phenotype and SNP–PC covariances, phenotype covariance matrices, variances, and missingness rates) and contain no individual-level UK Biobank data. They comply with UKB’s policy on derived summary statistics and are therefore made openly accessible for all users.

## License
This project is covered under the GNU General Public License v2.0.

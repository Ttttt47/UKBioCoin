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


Examples can be found in `./demo`.




# How to use UKBioCoin within docker images

We have packed the UKBioCoin algorithm with the NSS calculated from UKBioBank data into a seamless docker image for convenient use. 



Description of all the phenotype included in our data is given in `./description of phenotype0422.csv` in the respotory. We also provide the first 20 PC as covarites, use them as `PC{n}` with `{n}` replaced with integer ranging from 1 to 20 in the commandline. If you are using phenotype with FieldID `34`, replace it with `X34.0.0` in the commandline calling UKBioCoin algorithm.

To use the image, you should first log in to the repository using:

```
docker login -u cn-east-3@CFB4B8RZY5ZZ5RYWIHRI -p \
3f70cdbb9813cc57ad3a4b08de08217bcd5575ba09daef2b41b6388f0f336bc1 \
swr.cn-east-3.myhuaweicloud.com
```

Then pull the full image containing the UKBioCoin algorithm and the UKB NSS data using:

```
docker pull swr.cn-east-3.myhuaweicloud.com/ukbiocoin/ukbiocoin_full:v1.3
```

The data and the software are stored in /UKB. Note that the NSS data are compressed. You should first decompress them using:
```
tar -I 'zstd -v' -xvf 10M_0107_15.tar.zst
```

To pull the image containing only the UKBioCoin algorithm, use:

```bash
docker pull swr.cn-east-3.myhuaweicloud.com/ukbiocoin/ukbiocoin_none:v1.3
```


## Running UKBioCoin in command line
Here is an example using UKBioCoin in command line.
```{bash}
UKBioCoin --file test_data/sam \ 
            --phe X31.0.0 \
            --covar X1160.0.0,X1200.0.0,X1289.0.0,PC1,PC2,PC3,PC4,PC5 \
            --out test_data/test \
            --totalsize 292216 \
            --use-missing-rate-estimate 
```
- `--file`: Specifies the input file prefix. In this example, it is set to `test_data/sam`, the software will then try to find `test_data/sam_cov_xy.table`, `test_data/sam_cov_yy.table`, `test_data/sam_var_x.table` and `test_data/sam_meta.table`.

- `--phe`: Specifies the phenotype to be analyzed. In this example, it is set to `X31.0.0`.

- `--covar`: Specifies the covariates to be included in the analysis. Multiple covariates can be included by separating them with commas and no space in between. In this example, the covariates included are `X1160.0.0`, `X1200.0.0`, `X1289.0.0`, `PC1`, `PC2`, `PC3`, `PC4`, and `PC5`. If you don't want any covariates, skip these parameter.

- `--out`: Specifies the output file prefix. In this example, it is set to `test_data/test`. The output files will have this prefix with `_results.table` ended, namely, `test_data/test_results.table`.

- `--totalsize`: Total sample size of the regression, default is 292216 (sample size of our working UKB data.)

- `--use-missing-rate-estimate`: Whether use missing_rate files to estimate sample size. If so, the software will try to find missing rates file `xxx_x_missing.table` and `xxx_y_missing.table` and use them to estimate the sample size.


**PS**: Additionally, you can also set the following parameter to estimate the sample size instead of using `--use-missing-rate-estimate`.

- `--overall-non-missing-rate`: The overall non-missing rate of the data. If `--use-missing-rate-estimate` isn't specified, the sample size will be set to the value of `--totalsize` times value of parameter `--overall-non-missing-rate` (default is 0.9).

## Output format
The regression results file is typically a tabular format that presents the estimated coefficients and their associated statistical information for every SNPs, the ordering of the results are the same with the input NSS files.

for example:
```
#CHROM ID POS Effect_Allele REF_Allele ALT_Allele REF_FREQ BETA SE T-STAT -log10_P nobs Quality-Score
1 rs568927457 13453 C C T 0.00119 0.066221 0.0394034 1.68059 1.03225 262994 0.9
1 rs199856693 14933 A A G 0.01211 0.0230924 0.0124185 1.85952 1.20097 262994 0.9
1 rs533630043 15585 A A G 0.00192 -0.00582162 0.0310293 -0.187617 0.0699801 262994 0.9
1 rs568149713 15777 G G A 0.00331 -0.0474425 0.0236265 -2.00802 1.35025 262994 0.9
...
```

The first 7 columns is the same with the input meta file (`xxx_meta.table`) describing info of the SNPs, and the rest columns are the regression results. The meaning of the columns are as follows:

- `BETA`: This column contains the estimated regression coefficients for each predictor variable.

- `SE`: This column contains the standard error estimates for each coefficient.

- `T-STAT`: This column contains the t-statistic values for each coefficient, which indicate the strength of the evidence against the null hypothesis that the true coefficient is zero.

- `-log10_P`: This column contains the negative base-10 logarithm of the p-value for each coefficient, which measures the statistical significance of each predictor variable. The more postive the value, the more significant the variable is.

- `nobs`: This column contains the estimated number of observations used in the regression.

- `Quality-Score`: This column contains the quality score of the regression, which is the product of non-missing rates of SNP, covariates and phenotype.




# How to build a `UKBioCoin' device with your own dataset

## Building UKC excutable `UKBioCoin`

If you are intend to build a UKBioCoin device with your own dataset, you should first build UKC excutable `UKBioCoin`. Beside using the docker image mentioned above, you can also build the UKBioCoin excutable from the source code `main.cpp`. 

The UKBioCoin excutable depends on the following libraries:
- Eigen
- Boost
- cpp_cdfs (already included in `tools` folder)

### Eigen

Eigen is a C++ template library for linear algebra. It provides a wide range of matrix and vector operations, making it useful for various numerical computations.

To install Eigen, you can follow these steps:

1. Download the latest stable release of [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download).
2. Extract the downloaded archive.
3. Copy the `Eigen` directory to a location where your compiler can find it. For example, you can copy it to `/usr/local/include`.

### Boost

Boost is a collection of peer-reviewed C++ libraries that provide support for various tasks and functionalities, such as string manipulation, file system operations, and more.

To install Boost, you can use the following commands:

```bash

sudo apt-get install libboost-all-dev
    
```

### Building UKBioCoin
To build UKBioCoin, make sure you have installed Eigen, Boost, and a C++ compiler that supports C++11. Also, make sure that `tools` is in the same directory with `main.cpp`. You can use the following commands:

```bash
g++ -std=c++11 main.cpp -I "/path-to-eigen-X.X.X/" -lboost_program_options -o UKBioCoin 
```

## Generating Naive Summary Statistics (NSS) for UKBioCoin

UKBioCoin uses Naive Summary Statistics (NSS) as input, which includes files with the same prefix and three different suffixes: `_meta.table`, `_cov_xy.table`, `_cov_yy.table`, and `_var_x.table`. Additionaly, if you want to estimate the sample size, you should also provide the missing rate files `xxx_x_missing.table` and `xxx_y_missing.table`.

We provide a script at `demo/script.R` for you to generate the NSS files from the original PLINK's bfile/pfile and the phenotype(covariates) data. To use this script, you should make sure that you have PLINK/PLINK2 installed in your system and the `plink`/`plink2` command is in your `PATH`.

You may first test the script with the provided demo data at `demo/` folder with the following command:

```bash
Rscript script.R --bfile  euro_10Ksnp --pheno euro_10Ksnp.all.phe --threads 2 --memory 8000 --prefix euro_10Ksnp --ukc path-to-UKC-executable
```

For more usage of the script, you can use the following command to get the help information:

```bash
Rscript script.R --help
```

## Input Formats as Naive Summary Statistics for UKBioCoin

You may also generate the NSS files by yourself. The following is the detailed description of the input formats.

## xxx_meta.table

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

## xxx_cov_xy.table

For the `xxx_cov_xy.table` format, the first row lists all the phenotypes, including principal components, as column names. Each subsequent row represents a single SNP and contains covariances between that SNP and all the phenotypes. The first column in each row is the row number, enclosed in quotes, and the remaining columns contain the covariance values.

For example:

```
"phenotype1" "phenotype2" "phenotype3" "PC1" "PC2"
"1" 0.245 0.168 0.322 0.054 0.003
"2" 0.082 0.316 0.154 0.092 0.007
"3" -0.015 -0.062 -0.042 -0.012 0.001
...
```

## xxx_cov_yy.table

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

## xxx_var_x.table

For the `xxx_var_x.table` format, the first row is an arbitrary column name, such as "`var_x`". The second row lists the variance of each SNP, with each row representing a single SNP. The first column in each row is the row number, enclosed in quotes, and the second column contains the variance value.

For example:
```
"var_x"
"1" 0.153
"2" 0.214
"3" 0.082
...
```

# Access statement
The NSS provided with UKC consist solely of summary statistics (SNP–phenotype and SNP–PC covariances, phenotype covariance matrices, variances, and missingness rates) and contain no individual-level UK Biobank data. They comply with UKB’s policy on derived summary statistics and are therefore made openly accessible for all users.

# License
This project is covered under the GNU General Public License v2.0.



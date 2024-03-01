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


Examples can be found in `./test_data`

# How to use UKBioCoin within docker images
We have packed the UKBioCoin algorithm with the NSS calculated from UKBioBank data into a seamless docker image for convenient use. 
Description of all the phenotype included in our data is given in `./description of phenotype0327.csv` in the respotory. We also provide the first 20 PC as covarites, use them as `PC{n}` with `{n}` replaced with integer ranging from 1 to 20 in the commandline. If you are using phenotype with FieldID `34`, replace it with `X34.0.0` in the commandline calling UKBioCoin algorithm.

To use the image, you should first log in to the repository using:

```
docker login -u cn-east-3@CFB4B8RZY5ZZ5RYWIHRI -p \
3f70cdbb9813cc57ad3a4b08de08217bcd5575ba09daef2b41b6388f0f336bc1 \
swr.cn-east-3.myhuaweicloud.com
```

Then pull the image using:

```
docker pull swr.cn-east-3.myhuaweicloud.com/ukbiocoin/ukbiocoin:v6
```

The data and the software are stored in /UKB. Note that the NSS data are compressed. You should first decompress them using:
```
tar -I 'zstd -v' -xvf 14M_15.tar.zst
```


## Running UKBioCoin in command line
Here is an example using UKBioCoin in command line.
```{bash}
UKBioCoin --file test_data/sam \ 
            --phe X31.0.0 \
            --covar X1160.0.0,X1200.0.0,X1289.0.0,PC1,PC2,PC3,PC4,PC5 \
            --out test_data/test \
            --size 270000
```
`--file`: Specifies the input file prefix. In this example, it is set to `test_data/sam`, the software will then try to find `test_data/sam_cov_xy.table`, `test_data/sam_cov_yy.table`, `test_data/sam_var_x.table` and `test_data/sam_meta.table`.

`--phe`: Specifies the phenotype to be analyzed. In this example, it is set to `X31.0.0`.

`--covar`: Specifies the covariates to be included in the analysis. Multiple covariates can be included by separating them with commas and no space in between. In this example, the covariates included are `X1160.0.0`, `X1200.0.0`, `X1289.0.0`, `PC1`, `PC2`, `PC3`, `PC4`, and `PC5`. If you don't want any covariates, skip these parameter.

`--out`: Specifies the output file prefix. In this example, it is set to `test_data/test`. The output files will have this prefix with `_results.table` ended, namely, `test_data/test_results.table`.

`--size`: Sample size of the regression, default is 270000 (0.9 times 300000, an approximation of the sample size of the used UKB data.)

## Output format
The regression results file is typically a tabular format that presents the estimated coefficients and their associated statistical information for every SNPs, the ordering of the results are the same with the input NSS files.

for example:
```
#CHROM ID POS REF_Allele ALT_Allele REF_FREQ BETA SE T-STAT -log10_P
1 rs568927457 13453 C T 0.00119 0.00880088 0.0389788 0.225786 0.0854624
1 rs554760071 13483 C G 0.001235 -0.10489 0.0382549 -2.74188 2.21401
1 rs199856693 14933 A G 0.012111 -0.0148937 0.0122842 -1.21243 0.647142
1 rs533630043 15585 A G 0.001922 -0.0610412 0.0306952 -1.98862 1.33028
...
```

`BETA`: This column contains the estimated regression coefficients for each predictor variable.

`SE`: This column contains the standard error estimates for each coefficient.

`T-STAT`: This column contains the t-statistic values for each coefficient, which indicate the strength of the evidence against the null hypothesis that the true coefficient is zero.

`-log10_P`: This column contains the negative base-10 logarithm of the p-value for each coefficient, which measures the statistical significance of each predictor variable. The more negative the value, the more significant the variable is.


# How to build a `UKBioCoin' device with your own dataset
UKBioCoin uses Naive Summary Statistics as input, which includes files with the same prefix and three different suffixes: `_cov_xy.table`, `_cov_yy.table`, and `_var_x`.table. For example, `sam_cov_xy.table`, `sam_cov_yy.table`, and `sam_var_x.table`. These three files represent the covariance matrix between SNPs and phenotypes, the covariance matrix of phenotypes, and the covariance matrix of var_x.

## Input Formats as Naive Summary Statistics for UKBioCoin

UKBioCoin uses Naive Summary Statistics as input, which includes files with the same prefix and three different suffixes: `_cov_xy.table`, `_cov_yy.table`, and `_var_x`.table. For example, `sam_cov_xy.table`, `sam_cov_yy.table`, and `sam_var_x.table`. These three files represent the covariance matrix between SNPs and phenotypes, the covariance matrix of phenotypes, and the variance of all the SNPs.

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

# UKBioCoin

This project is currently developling...

Examples can be found in `./test_data`

# How to use UKBioCoin within docker images
We have packed the UKBioCoin algorithm with the NSS calculated from UKBioBank data into a seamless docker image for convenient use. To use the image, you should first log in to the repository using:

```
docker login -u cn-east-3@CFB4B8RZY5ZZ5RYWIHRI -p \
3f70cdbb9813cc57ad3a4b08de08217bcd5575ba09daef2b41b6388f0f336bc1 \
swr.cn-east-3.myhuaweicloud.com
```

Then pull the image using:

```
docker pull swr.cn-east-3.myhuaweicloud.com/ukbiocoin/ukbiocoin:v3
```

The data and the software are stored in /UKB. Note that the NSS data are compressed. You should first decompress them using:
```
tar -I 'zstd -v' -xvf 14M_15.tar.zst
```

## Running UKBiocoin in command line
Here is an example using UKBioCoin in command line.
```{bash}
.\main.exe --file test_data/sam \ 
            --phe X31.0.0 \
            --covar X1160.0.0,X1200.0.0,X1289.0.0,PC1,PC2,PC3,PC4,PC5 \
            --out test_data/test
```
`--file`: Specifies the input file prefix. In this example, it is set to `test_data/sam`.

`--phe`: Specifies the phenotype to be analyzed. In this example, it is set to `X31.0.0`.

`--covar`: Specifies the covariates to be included in the analysis. Multiple covariates can be included by separating them with commas and no space in between. In this example, the covariates included are `X1160.0.0`, `X1200.0.0`, `X1289.0.0`, `PC1`, `PC2`, `PC3`, `PC4`, and `PC5`.

`--out`: Specifies the output file prefix. In this example, it is set to `test_data/test`. The output files will have this prefix.

## Output format
The regression results file is typically a tabular format that presents the estimated coefficients and their associated statistical information for every SNPs, the ordering of the results are the same with the input NSS files.

for example:
```
BETA SE T-STAT log10_P
2.58495 0.317914 8.13097 -15.3689
-0.0746259 0.0386125 -1.93269 -1.27347
0.130353 0.0612135 2.12948 -1.47866
0.388451 0.390303 0.995254 -0.495375
...
```

`BETA`: This column contains the estimated regression coefficients for each predictor variable.

`SE`: This column contains the standard error estimates for each coefficient.

`T-STAT`: This column contains the t-statistic values for each coefficient, which indicate the strength of the evidence against the null hypothesis that the true coefficient is zero.

`log10_P`: This column contains the negative base-10 logarithm of the p-value for each coefficient, which measures the statistical significance of each predictor variable. The more negative the value, the more significant the variable is.


# How to use UKBioCoin with your own dataset
UKBioCoin uses Naive Summary Statistics as input, which includes files with the same prefix and three different suffixes: `_cov_xy.table`, `_cov_yy.table`, and `_var_x`.table. For example, `sam_cov_xy.table`, `sam_cov_yy.table`, and `sam_var_x.table`. These three files represent the covariance matrix between SNPs and phenotypes, the covariance matrix of phenotypes, and the covariance matrix of var_x.

## Input Formats as Naive Summary Statistics for UKBioCoin

UKBioCoin uses Naive Summary Statistics as input, which includes files with the same prefix and three different suffixes: `_cov_xy.table`, `_cov_yy.table`, and `_var_x`.table. For example, `sam_cov_xy.table`, `sam_cov_yy.table`, and `sam_var_x.table`. These three files represent the covariance matrix between SNPs and phenotypes, the covariance matrix of phenotypes, and the variance of all the SNPs.

## xxx_cov_xy.table

For the `xxx_cov_xy.table` format, the first row lists all the phenotypes, including principal components, as column names. Each subsequent row represents a single SNP and contains covariances between that SNP and all the phenotypes. The first column in each row is the row number, enclosed in quotes, and the remaining columns contain the covariance values.

For example:

```
"phenotype1" "phenotype2" "phenotype3" "PC1" "PC2"
"1" 0.245 0.168 0.322 0.054 0.003
"2" 0.082 0.316 0.154 0.092 0.007
"3" -0.015 -0.062 -0.042 -0.012 0.001
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
```

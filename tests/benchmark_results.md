# UKC 2.0 acceptance results

Tested on 2026-08-07 with the `ukbiocoin-dev` environment, Eigen from
`/home/jingcheng/fsl/include/eigen3`, and 16 threads where noted. Wall time and
peak RSS were recorded with GNU `time`; timing values are medians of three runs.

## Demo (`demo/euro_10Ksnp`)

The analysis selected `Trait1_H2_0.25` with `PC1`-`PC5` and missing-rate sample
size estimation over 100,000 variants.

| Input/mode | Threads | Median wall | Median peak RSS |
|---|---:|---:|---:|
| NSS v1 stream | 1 | 0.35 s | 6.1 MiB |
| NSS v2 stream | 1 | 0.22 s | 6.1 MiB |
| NSS v2 stream | 16 | 0.13 s | 6.1 MiB |
| NSS v2 memory | 1 | 0.23 s | 27.1 MiB |
| NSS v2 memory | 16 | 0.14 s | 27.1 MiB |

- R-direct v2, adapter v2, v2 stream/memory, and 1/16-thread outputs had
  identical 100,000-row SNP ordering. R-direct and adapter result files were
  byte-identical.
- Against legacy v1, all displayed statistical fields had maximum absolute
  difference `1e-5` or less and correlation above `0.999999`; `nobs` and
  `Quality-Score` were exact.
- No-covariate BETA matched PLINK2 exactly after allele alignment and exclusion
  of 267 frequency-0.5 ambiguous variants. SE correlation was `0.999969`,
  T-statistic correlation was `0.999994`, and P-value correlation was
  `0.999992` under UKC's preserved statistical definitions.
- Legacy COVxy was 20,992,911 bytes; v2 COVxy was 4,800,960 bytes (`22.87%`).
- The R v2 generator peaked at 304.1 MiB versus 315.5 MiB for legacy generation,
  consistent with eliminating the 100,000-by-12 in-memory COVxy matrix.
- Adapter conversion took 0.14 s and peaked at 7.6 MiB.

## Synthetic wide-phenotype benchmark

The synthetic input contained 100,000 variants and 128 phenotypes; `W127` was
selected so the legacy parser had to traverse the full text row.

| Input/mode | Threads | Median wall | Median peak RSS |
|---|---:|---:|---:|
| NSS v1 stream | 1 | 0.86 s | 14.2 MiB |
| NSS v2 stream | 1 | 0.12 s | 4.5 MiB |
| NSS v2 stream | 16 | 0.10 s | 4.5 MiB |
| NSS v2 memory | 1 | 0.13 s | 18.1 MiB |
| NSS v2 memory | 16 | 0.11 s | 18.1 MiB |

All v2 stream/memory and 1/16-thread result files were byte-identical. Wide v2
COVxy occupied `25.72%` of the legacy text COVxy size.

## Automated checks

Regular CTest and the AddressSanitizer/UndefinedBehaviorSanitizer build passed.
The regression suite covers NumPy interoperability, manifest dimensions and
sample count, stream/memory and thread determinism, no-covariate computation,
NA preservation, duplicate legacy row labels, missing files, invalid dtype and
shape, truncated NPY payloads, singular covariance, row-count mismatch, and
invalid CLI values.

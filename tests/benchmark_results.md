# UKC 2.0 acceptance results

Tested on 2026-08-07 with the `ukbiocoin-dev` environment, locally installed
Eigen headers, and 16 threads where noted. Wall time and
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

## Buffered text-output benchmark

Tested on 2026-08-08 after replacing per-value iostream insertion with bounded
batch buffers, Boost.Charconv numeric formatting, and ordered bulk writes. The
output remains the historical `_results.table` format. A synthetic NSS with
1,000,000 variants and 64 phenotypes selected one phenotype and seven
covariates. Values below are medians of three runs; every optimized output was
byte-identical to the pre-optimization executable.

| Implementation | Threads | Total | Result writing | Output speedup | Total speedup |
|---|---:|---:|---:|---:|---:|
| Previous iostream | 1 | 2.099 s | 1.045 s | 1.00x | 1.00x |
| Buffered/charconv | 1 | 1.302 s | 0.254 s | 4.11x | 1.61x |
| Previous iostream | 8 | 1.341 s | 1.022 s | 1.00x | 1.00x |
| Buffered/charconv | 8 | 0.670 s | 0.236 s | 4.33x | 2.00x |

The result table was 86.6 MiB. A second test with 4,000,000 variants produced a
352.6 MiB table: 8-thread result writing fell from 4.037 s to 1.012 s (3.99x),
and total time fell from 5.523 s to 2.891 s (1.91x). The optimized stream output
used 6.0 MiB peak RSS. The optimized memory and stream modes produced
byte-identical results; memory mode used 981.7 MiB peak RSS.

## Numeric NPY result-output benchmark

Tested on 2026-08-09 with a synthetic NSS v2 containing 4,000,000 variants and
16 phenotypes. One phenotype and seven covariates were selected. Table and NPY
runs were alternated to reduce cache and transient-load bias; values are medians
of three runs in stream mode.

| Format | Threads | Reading | Computation | Result writing | Total | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| table | 1 | 0.230 s | 3.987 s | 1.406 s | 5.781 s | 6.1 MiB |
| npy | 1 | 0.092 s | 3.853 s | 0.321 s | 4.274 s | 4.6 MiB |
| table | 8 | 0.400 s | 0.862 s | 1.614 s | 3.048 s | 4.6 MiB |
| npy | 8 | 0.187 s | 0.672 s | 0.366 s | 1.238 s | 4.6 MiB |

NPY reduced result-writing time by 4.38x at one thread and 4.41x at eight
threads. Total time improved by 1.35x and 2.46x, respectively. The numeric NPY
was 213.6 MiB versus 267.1 MiB for the table because metadata and decimal text
were omitted. Single- and eight-thread NPY files were byte-identical. Sampled
rows agreed with the table values within its six-significant-digit text
precision.

The eight-thread NPY memory mode took a median 2.340 s and 661 MiB peak RSS,
compared with 1.238 s and 4.6 MiB for stream mode. Loading all selected vectors
therefore provided no benefit for this single-pass analysis.

## Automated checks

Regular CTest and the AddressSanitizer/UndefinedBehaviorSanitizer build passed.
The regression suite covers NumPy interoperability, manifest dimensions and
sample count, stream/memory and thread determinism, no-covariate computation,
numeric result NPY shape/dtype and metadata independence, NA preservation,
duplicate legacy row labels, missing files, invalid dtype and shape, truncated
NPY payloads, singular covariance, row-count mismatch, and invalid CLI values.

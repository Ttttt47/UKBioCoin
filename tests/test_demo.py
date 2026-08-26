#!/usr/bin/env python3
"""Manual end-to-end validation on the repository's 100,000-variant demo.

Usage:
    python3 tests/test_demo.py UKBioCoin UKBioCoinNSSAdapter Rscript

The test intentionally is not registered as a routine CTest: it runs PLINK2
twice to generate fresh v2 and legacy NSS data, creates a validation plot, and
executes a matrix of full 100,000-row regressions.  ``plink2`` is expected next
to the supplied Rscript or otherwise on PATH.
"""

import argparse
import hashlib
import json
import math
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

import numpy as np


VARIANT_COUNT = 100_000
PHENOTYPE = "Trait1_H2_0.25"
CANONICAL_METADATA = ["#CHROM", "POS", "ID", "REF", "ALT1", "A1", "REF_FREQ"]
RESULT_COLUMNS = [
    "BETA", "SE", "T_STAT", "NEG_LOG10_P", "VIF", "OBS_CT", "QUALITY_SCORE"
]

# The runtime table uses six significant digits.  Comparing it with its binary
# NPY counterpart therefore needs a six-digit text-rounding tolerance.
TABLE_RTOL = 6e-6

# COVxy is intentionally float32 in NSS v2.  Legacy NSS stores decimal text,
# and its cov_yy/var_x writers also round-trip through R text formatting before
# adapter conversion.  These bounds are tight enough to detect row/direction or
# formula errors while allowing those documented representation differences.
NSS_FORMAT_RTOL = 2e-4
NSS_FORMAT_ATOL = 2e-6

# The demo scales phenotypes with R's sample-SD convention, while UKC rebuilds
# OLS quantities from normalized covariances and its explicit residual DF.
# Hence PLINK BETA is very close, but SE/T/P can differ by a few per mille.
PLINK_BETA_RTOL = 3e-4
PLINK_BETA_ATOL = 2e-6
PLINK_SE_T_RTOL = 8e-2
PLINK_SE_T_ATOL = 2e-4
PLINK_P_RTOL = 1.5e-1
PLINK_P_ATOL = 2e-2


def resolve_executable(value):
    candidate = pathlib.Path(value).expanduser()
    if candidate.is_file():
        return str(candidate.resolve())
    resolved = shutil.which(value)
    if resolved is None:
        raise SystemExit(f"executable not found: {value}")
    return str(pathlib.Path(resolved).resolve())


def display_command(command):
    return " ".join(str(item) for item in command)


def run(command, env, label, timeout=600):
    command = [str(item) for item in command]
    print(f"[{label}] {display_command(command)}", flush=True)
    result = subprocess.run(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env=env, timeout=timeout,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"{label} failed ({result.returncode})\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def file_digest(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as stream:
        while True:
            block = stream.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def metadata_digest_update(digest, fields):
    digest.update("\t".join(fields).encode("utf-8"))
    digest.update(b"\n")


def parse_number(text):
    return math.nan if text in ("NA", "NaN", "nan") else float(text)


def read_metadata(path, expected_rows=VARIANT_COUNT):
    path = pathlib.Path(path)
    digest = hashlib.sha256()
    ids = []
    a1 = []
    with path.open(encoding="utf-8") as stream:
        header_line = stream.readline()
        if not header_line:
            raise AssertionError(f"metadata is empty: {path}")
        header = header_line.rstrip("\r\n").split()
        if header != CANONICAL_METADATA:
            raise AssertionError(f"noncanonical metadata header in {path}: {header}")
        for line_number, line in enumerate(stream, start=2):
            fields = line.rstrip("\r\n").split()
            if len(fields) != len(CANONICAL_METADATA):
                raise AssertionError(
                    f"{path}:{line_number}: expected 7 metadata fields, got {len(fields)}"
                )
            if fields[5] not in (fields[3], fields[4]):
                raise AssertionError(f"{path}:{line_number}: A1 is not REF or ALT1")
            metadata_digest_update(digest, fields)
            ids.append(fields[2])
            a1.append(fields[5])
    if len(ids) != expected_rows:
        raise AssertionError(f"{path}: expected {expected_rows} rows, got {len(ids)}")
    if len(set(ids)) != len(ids):
        raise AssertionError(f"{path}: duplicate variant IDs")
    return digest.hexdigest(), ids, a1


def pvar_ids(path):
    ids = []
    with pathlib.Path(path).open(encoding="utf-8") as stream:
        header = stream.readline().rstrip("\r\n").split("\t")
        try:
            id_index = header.index("ID")
        except ValueError as error:
            raise AssertionError(f"pvar has no ID column: {path}") from error
        for line in stream:
            fields = line.rstrip("\r\n").split("\t")
            ids.append(fields[id_index])
    return ids


def verify_glm(path, expected_ids, expected_a1, collect_statistics=False):
    path = pathlib.Path(path)
    statistics = np.empty((len(expected_ids), 4), dtype=np.float64) \
        if collect_statistics else None
    with path.open(encoding="utf-8") as stream:
        header = stream.readline().rstrip("\r\n").split()
        required = ["ID", "A1"]
        if collect_statistics:
            required += ["BETA", "SE", "T_STAT", "P"]
        missing = [name for name in required if name not in header]
        if missing:
            raise AssertionError(f"{path}: missing GLM columns {missing}")
        indices = {name: header.index(name) for name in required}
        test_index = header.index("TEST") if "TEST" in header else None
        row = 0
        for line_number, line in enumerate(stream, start=2):
            fields = line.rstrip("\r\n").split()
            if test_index is not None and fields[test_index] != "ADD":
                continue
            if row >= len(expected_ids):
                raise AssertionError(f"{path}:{line_number}: extra additive row")
            if fields[indices["ID"]] != expected_ids[row]:
                raise AssertionError(f"{path}:{line_number}: variant order changed")
            if fields[indices["A1"]] != expected_a1[row]:
                raise AssertionError(f"{path}:{line_number}: A1 direction changed")
            if collect_statistics:
                statistics[row] = [
                    parse_number(fields[indices[name]])
                    for name in ("BETA", "SE", "T_STAT", "P")
                ]
            row += 1
    if row != len(expected_ids):
        raise AssertionError(f"{path}: expected {len(expected_ids)} ADD rows, got {row}")
    return statistics


def read_result_table(path, expected_rows=VARIANT_COUNT):
    path = pathlib.Path(path)
    statistics = np.empty((expected_rows, 7), dtype=np.float64)
    digest = hashlib.sha256()
    with path.open(encoding="utf-8") as stream:
        header = stream.readline().rstrip("\r\n").split("\t")
        expected_header = CANONICAL_METADATA + RESULT_COLUMNS
        if header != expected_header:
            raise AssertionError(f"unexpected result header in {path}: {header}")
        row = 0
        for line_number, line in enumerate(stream, start=2):
            if row >= expected_rows:
                raise AssertionError(f"{path}:{line_number}: extra result row")
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(expected_header):
                raise AssertionError(f"{path}:{line_number}: result row-width mismatch")
            metadata_digest_update(digest, fields[:7])
            statistics[row] = [parse_number(value) for value in fields[7:]]
            row += 1
    if row != expected_rows:
        raise AssertionError(f"{path}: expected {expected_rows} rows, got {row}")
    return digest.hexdigest(), statistics


def check_npy_result(path, expected_rows=VARIANT_COUNT):
    values = np.load(path, allow_pickle=False)
    if values.dtype.str != "<f8" or values.shape != (expected_rows, 7):
        raise AssertionError(
            f"{path}: expected <f8 {(expected_rows, 7)}, got "
            f"{values.dtype.str} {values.shape}"
        )
    return values


def assert_nan_pattern(actual, expected, label):
    if not np.array_equal(np.isnan(actual), np.isnan(expected)):
        raise AssertionError(f"{label}: NaN pattern differs")


def assert_close(actual, expected, label, rtol, atol):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    if actual.shape != expected.shape:
        raise AssertionError(f"{label}: shape {actual.shape} != {expected.shape}")
    assert_nan_pattern(actual, expected, label)
    np.testing.assert_allclose(
        actual, expected, rtol=rtol, atol=atol, equal_nan=True,
        err_msg=label,
    )
    finite = np.isfinite(actual) & np.isfinite(expected)
    max_absolute = float(np.max(np.abs(actual[finite] - expected[finite]))) \
        if finite.any() else 0.0
    print(f"  {label}: max |delta|={max_absolute:.3g}")


def verify_v2_layout(path, phenotype_count, expected_metadata_digest):
    path = pathlib.Path(path)
    manifest = json.loads((path / "manifest.json").read_text(encoding="utf-8"))
    expected_dimensions = {
        "sample_count": 200,
        "variant_count": VARIANT_COUNT,
        "phenotype_count": phenotype_count,
    }
    if manifest.get("dimensions") != expected_dimensions:
        raise AssertionError(f"{path}: wrong manifest dimensions")
    if manifest.get("schema") != "org.ukbiocoin.nss" or \
            manifest.get("schema_version") != 2:
        raise AssertionError(f"{path}: wrong manifest schema/version")
    metadata_hash, _, _ = read_metadata(path / "meta.tsv")
    if metadata_hash != expected_metadata_digest:
        raise AssertionError(f"{path}: metadata differs from fresh canonical metadata")

    expected_arrays = (
        ("cov_yy.npy", "<f8", (phenotype_count, phenotype_count)),
        ("var_x.npy", "<f8", (VARIANT_COUNT,)),
        ("x_missing.npy", "<f8", (VARIANT_COUNT,)),
        ("y_missing.npy", "<f8", (phenotype_count,)),
    )
    for relative, dtype, shape in expected_arrays:
        values = np.load(path / relative, mmap_mode="r", allow_pickle=False)
        if values.dtype.str != dtype or values.shape != shape:
            raise AssertionError(f"{path / relative}: bad dtype/shape")
    phenotypes = manifest.get("phenotypes", [])
    if len(phenotypes) != phenotype_count:
        raise AssertionError(f"{path}: wrong phenotype descriptor count")
    for index, descriptor in enumerate(phenotypes):
        expected_relative = f"cov_xy/{index:08d}.npy"
        if descriptor.get("index") != index or \
                descriptor.get("cov_xy", {}).get("path") != expected_relative:
            raise AssertionError(f"{path}: malformed phenotype descriptor {index}")
        values = np.load(path / expected_relative, mmap_mode="r", allow_pickle=False)
        if values.dtype.str != "<f4" or values.shape != (VARIANT_COUNT,):
            raise AssertionError(f"{path / expected_relative}: bad dtype/shape")
    return manifest


def compare_v2_payloads(fresh, adapted, phenotype_count):
    fresh = pathlib.Path(fresh)
    adapted = pathlib.Path(adapted)
    for relative in ("cov_yy.npy", "var_x.npy", "x_missing.npy", "y_missing.npy"):
        assert_close(
            np.load(adapted / relative, allow_pickle=False),
            np.load(fresh / relative, allow_pickle=False),
            f"adapter payload {relative}", NSS_FORMAT_RTOL, NSS_FORMAT_ATOL,
        )
    for index in range(phenotype_count):
        relative = f"cov_xy/{index:08d}.npy"
        assert_close(
            np.load(adapted / relative, allow_pickle=False),
            np.load(fresh / relative, allow_pickle=False),
            f"adapter payload {relative}", NSS_FORMAT_RTOL, NSS_FORMAT_ATOL,
        )


def regression_command(binary, nss, prefix, model, io_mode, threads,
                       output_format, legacy=False):
    command = [
        binary, "--file", nss, "--phe", PHENOTYPE,
        "--use-missing-rate-estimate", "--threads", str(threads),
        "--io-mode", io_mode, "--output-format", output_format,
        "--out", prefix,
    ]
    if model[1] is not None:
        command += ["--covar", model[1]]
    if legacy:
        command += ["--totalsize", "200"]
    return command


def run_regression(binary, nss, prefix, model, io_mode, threads,
                   output_format, legacy, env, label):
    run(
        regression_command(binary, nss, prefix, model, io_mode, threads,
                           output_format, legacy),
        env, label,
    )
    suffix = "_results.npy" if output_format == "npy" else "_results.table"
    output = pathlib.Path(str(prefix) + suffix)
    if not output.is_file():
        raise AssertionError(f"{label}: missing output {output}")
    return output


def compare_table_with_npy(table_values, npy_values, label):
    # Six significant digits imply a relative rounding error below 5e-6 for
    # every nonzero finite value; use a small margin for parser round trips.
    assert_close(table_values, npy_values, label, TABLE_RTOL, 1e-300)


def compare_format_result(actual, reference, label):
    # Observation counts must remain exact even when covariance storage changes.
    np.testing.assert_array_equal(actual[:, 5], reference[:, 5])
    assert_close(
        actual[:, :5], reference[:, :5], label + " statistics",
        NSS_FORMAT_RTOL, NSS_FORMAT_ATOL,
    )
    assert_close(
        actual[:, 6], reference[:, 6], label + " quality",
        2e-12, 2e-12,
    )


def plink_numeric_check(ukc_values, plink_values, label,
                        covariate_approximation=False):
    beta, se, t_statistic, p_value = plink_values.T
    neg_log10_p = np.full(p_value.shape, math.nan)
    valid_p = np.isfinite(p_value) & (p_value > 0.0)
    neg_log10_p[valid_p] = -np.log10(p_value[valid_p])

    def require_correlation(actual, expected, name, minimum=0.999):
        finite = np.isfinite(actual) & np.isfinite(expected)
        correlation = float(np.corrcoef(actual[finite], expected[finite])[0, 1])
        print(f"  {name}: correlation={correlation:.9f}")
        if correlation < minimum:
            raise AssertionError(
                f"{name}: correlation {correlation} is below {minimum}"
            )

    if covariate_approximation:
        # The NSS joint model is reconstructed from separately estimated
        # phenotype/SNP and covariate/SNP vectors plus pairwise cov_yy entries.
        # With missing genotypes, that is not algebraically identical to
        # PLINK's per-SNP complete-case joint regression.  Test direction and a
        # bounded absolute-error profile instead of pretending this is an
        # exact file-format comparison.
        comparisons = (
            (ukc_values[:, 0], beta, "BETA", 0.075),
            (ukc_values[:, 1], se, "SE", 0.075),
            (ukc_values[:, 2], t_statistic, "T_STAT", 0.75),
            (ukc_values[:, 3], neg_log10_p, "NEG_LOG10_P", 1.0),
        )
        for actual, expected, statistic, maximum in comparisons:
            name = f"{label} {statistic}"
            require_correlation(actual, expected, name, minimum=0.99)
            finite = np.isfinite(actual) & np.isfinite(expected)
            differences = np.abs(actual[finite] - expected[finite])
            q50, q99, observed_maximum = np.quantile(
                differences, [0.5, 0.99, 1.0]
            )
            print(
                f"  {name}: |delta| median={q50:.3g}, "
                f"p99={q99:.3g}, max={observed_maximum:.3g}"
            )
            if observed_maximum > maximum:
                raise AssertionError(
                    f"{name}: maximum absolute difference {observed_maximum} "
                    f"exceeds {maximum}"
                )
        return

    assert_close(
        ukc_values[:, 0], beta, label + " BETA",
        PLINK_BETA_RTOL, PLINK_BETA_ATOL,
    )
    assert_close(
        ukc_values[:, 1], se, label + " SE",
        PLINK_SE_T_RTOL, PLINK_SE_T_ATOL,
    )
    require_correlation(ukc_values[:, 1], se, label + " SE")
    assert_close(
        ukc_values[:, 2], t_statistic, label + " T_STAT",
        PLINK_SE_T_RTOL, PLINK_SE_T_ATOL,
    )
    require_correlation(ukc_values[:, 2], t_statistic, label + " T_STAT")
    assert_close(
        ukc_values[:, 3], neg_log10_p, label + " NEG_LOG10_P",
        PLINK_P_RTOL, PLINK_P_ATOL,
    )
    require_correlation(
        ukc_values[:, 3], neg_log10_p, label + " NEG_LOG10_P"
    )


def exercise_model(binary, model, fresh_v2, legacy_prefix, adapted_v2,
                   fresh_metadata_digest, generated_visual_result,
                   plink_statistics, root, env):
    name = model[0]
    model_root = root / name
    model_root.mkdir(parents=True)
    print(f"\n[model] {name}: covariates={model[1] or '(none)'}", flush=True)

    # v2: stream/memory, 1/16 threads, and table/NPY must all agree.
    v2_table_paths = {}
    if generated_visual_result is not None:
        v2_table_paths["stream16"] = pathlib.Path(generated_visual_result)
    else:
        v2_table_paths["stream16"] = run_regression(
            binary, fresh_v2, model_root / "v2_stream_table_16", model,
            "stream", 16, "table", False, env, f"{name} v2 stream table t16",
        )
    v2_table_paths["stream1"] = run_regression(
        binary, fresh_v2, model_root / "v2_stream_table_1", model,
        "stream", 1, "table", False, env, f"{name} v2 stream table t1",
    )
    v2_table_paths["memory16"] = run_regression(
        binary, fresh_v2, model_root / "v2_memory_table_16", model,
        "memory", 16, "table", False, env, f"{name} v2 memory table t16",
    )
    reference_table_digest = file_digest(v2_table_paths["stream1"])
    for mode, path in v2_table_paths.items():
        if file_digest(path) != reference_table_digest:
            raise AssertionError(f"{name}: v2 table differs for {mode}")
    table_metadata_hash, v2_table_values = read_result_table(
        v2_table_paths["stream1"]
    )
    if table_metadata_hash != fresh_metadata_digest:
        raise AssertionError(f"{name}: v2 result metadata changed")

    v2_npy_paths = {}
    for mode, io_mode, threads in (
        ("stream1", "stream", 1),
        ("stream16", "stream", 16),
        ("memory16", "memory", 16),
    ):
        v2_npy_paths[mode] = run_regression(
            binary, fresh_v2, model_root / f"v2_{mode}_npy", model,
            io_mode, threads, "npy", False, env, f"{name} v2 {mode} npy",
        )
    v2_npy_reference = check_npy_result(v2_npy_paths["stream1"])
    for mode, path in v2_npy_paths.items():
        values = check_npy_result(path)
        if not np.array_equal(values, v2_npy_reference, equal_nan=True):
            raise AssertionError(f"{name}: v2 NPY differs for {mode}")
    compare_table_with_npy(v2_table_values, v2_npy_reference,
                           f"{name} v2 table/NPY")

    # Direct legacy and adapter-converted v2 paths must preserve metadata and
    # agree numerically with the freshly generated v2 NSS.
    legacy_table = run_regression(
        binary, legacy_prefix, model_root / "legacy_table", model,
        "stream", 16, "table", True, env, f"{name} legacy table",
    )
    legacy_npy = run_regression(
        binary, legacy_prefix, model_root / "legacy_npy", model,
        "stream", 16, "npy", True, env, f"{name} legacy npy",
    )
    adapted_table = run_regression(
        binary, adapted_v2, model_root / "adapted_table", model,
        "stream", 16, "table", False, env, f"{name} adapted v2 table",
    )
    adapted_npy = run_regression(
        binary, adapted_v2, model_root / "adapted_npy", model,
        "stream", 16, "npy", False, env, f"{name} adapted v2 npy",
    )

    legacy_metadata_hash, legacy_table_values = read_result_table(legacy_table)
    adapted_metadata_hash, adapted_table_values = read_result_table(adapted_table)
    if legacy_metadata_hash != fresh_metadata_digest or \
            adapted_metadata_hash != fresh_metadata_digest:
        raise AssertionError(f"{name}: legacy/adapter result metadata changed")
    legacy_npy_values = check_npy_result(legacy_npy)
    adapted_npy_values = check_npy_result(adapted_npy)
    compare_table_with_npy(legacy_table_values, legacy_npy_values,
                           f"{name} legacy table/NPY")
    compare_table_with_npy(adapted_table_values, adapted_npy_values,
                           f"{name} adapted table/NPY")
    compare_format_result(legacy_npy_values, v2_npy_reference,
                          f"{name} legacy vs fresh v2")
    compare_format_result(adapted_npy_values, v2_npy_reference,
                          f"{name} adapter vs fresh v2")
    plink_numeric_check(
        v2_npy_reference, plink_statistics, f"{name} UKC vs PLINK",
        covariate_approximation=model[1] is not None,
    )
    return v2_table_paths["stream1"], v2_npy_paths["stream1"]


def adjustment_table(path, metadata_columns, expected_rows=VARIANT_COUNT):
    path = pathlib.Path(path)
    values = np.empty((expected_rows, 2), dtype=np.float64)
    digest = hashlib.sha256()
    with path.open(encoding="utf-8") as stream:
        header = stream.readline().rstrip("\r\n").split("\t")
        expected_header = metadata_columns + ["NEG_LOG10_UNADJ", "NEG_LOG10_GC"]
        if header != expected_header:
            raise AssertionError(f"{path}: wrong adjustment header {header}")
        row = 0
        for line_number, line in enumerate(stream, start=2):
            if row >= expected_rows:
                raise AssertionError(f"{path}:{line_number}: extra adjustment row")
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(expected_header):
                raise AssertionError(f"{path}:{line_number}: adjustment row-width mismatch")
            if metadata_columns == CANONICAL_METADATA:
                metadata_digest_update(digest, fields[:7])
            else:
                if fields[0] != str(row):
                    raise AssertionError(f"{path}:{line_number}: wrong ROW_INDEX")
            values[row] = [parse_number(value) for value in fields[-2:]]
            row += 1
    if row != expected_rows:
        raise AssertionError(f"{path}: expected {expected_rows} rows, got {row}")
    return digest.hexdigest(), values


def logged_lambda(result, key):
    prefix = key + ":"
    values = [
        float(line[len(prefix):].strip()) for line in result.stdout.splitlines()
        if line.startswith(prefix)
    ]
    if len(values) != 1:
        raise AssertionError(f"missing/duplicate {key} in adjustment log")
    return values[0]


def exercise_gc(binary, main_table, main_npy, metadata_digest, root, env):
    print("\n[GC] comparing table/NPY adjustment on the full demo result", flush=True)
    root.mkdir(parents=True, exist_ok=True)
    cases = {}
    for input_name, input_path in (("table", main_table), ("npy", main_npy)):
        for output_format in ("table", "npy"):
            prefix = root / f"gc_{input_name}_to_{output_format}"
            command = [binary, "--adjust-file", input_path, "--out", prefix,
                       "--output-format", output_format]
            result = run(command, env, f"GC {input_name}->{output_format}")
            suffix = ".adjusted.npy" if output_format == "npy" else ".adjusted"
            cases[(input_name, output_format)] = (result, pathlib.Path(str(prefix) + suffix))

    table_main_metadata, table_main_values = read_result_table(main_table)
    if table_main_metadata != metadata_digest:
        raise AssertionError("main table metadata changed before GC")
    npy_main_values = check_npy_result(main_npy)

    table_hash, table_to_table = adjustment_table(
        cases[("table", "table")][1], CANONICAL_METADATA
    )
    if table_hash != metadata_digest:
        raise AssertionError("table GC did not preserve canonical metadata")
    _, npy_to_table = adjustment_table(
        cases[("npy", "table")][1], ["#ROW_INDEX"]
    )
    table_to_npy = np.load(cases[("table", "npy")][1], allow_pickle=False)
    npy_to_npy = np.load(cases[("npy", "npy")][1], allow_pickle=False)
    for label, values in (("table->NPY", table_to_npy), ("NPY->NPY", npy_to_npy)):
        if values.dtype.str != "<f8" or values.shape != (VARIANT_COUNT, 2):
            raise AssertionError(f"GC {label} violates the <f8 N x 2 contract")

    if not np.array_equal(table_to_table, table_to_npy, equal_nan=True):
        raise AssertionError("table GC output differs between table and NPY destinations")
    if not np.array_equal(npy_to_table, npy_to_npy, equal_nan=True):
        raise AssertionError("NPY GC output differs between table and NPY destinations")
    if not np.array_equal(table_to_table[:, 0], table_main_values[:, 3],
                          equal_nan=True):
        raise AssertionError("table GC changed NEG_LOG10_UNADJ")
    if not np.array_equal(npy_to_npy[:, 0], npy_main_values[:, 3],
                          equal_nan=True):
        raise AssertionError("NPY GC changed NEG_LOG10_UNADJ")

    # Table input has six-digit T/P values while NPY retains full float64, so
    # their separately estimated lambdas/results are close rather than bitwise equal.
    table_raw = logged_lambda(cases[("table", "table")][0], "lambda_raw")
    npy_raw = logged_lambda(cases[("npy", "table")][0], "lambda_raw")
    if not math.isclose(table_raw, npy_raw, rel_tol=2e-5, abs_tol=2e-6):
        raise AssertionError(f"table/NPY GC lambda differs: {table_raw} vs {npy_raw}")
    assert_close(
        table_to_table, npy_to_npy, "demo table/NPY GC",
        4e-5, 2e-5,
    )


def generate_demo(rscript, demo_script, pfile, pheno, ukc, output, prefix,
                  nss_format, visualize, env):
    command = [
        rscript, demo_script, "--pfile", pfile, "--pheno", pheno,
        "--threads", "16", "--memory", "4000", "--prefix", prefix,
        "--ukc", ukc, "--out", output, "--nss-format", nss_format,
    ]
    if not visualize:
        command.append("--novisualize")
    return run(command, env, f"generate {nss_format} demo", timeout=900)


def make_covariate_plink_reference(plink2, pfile, scaled_pheno, frequency,
                                   output_prefix, env):
    pathlib.Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)
    command = [
        plink2, "--pfile", pfile, "--silent", "--glm", "hide-covar",
        "skip-invalid-pheno", "--no-input-missing-phenotype",
        "--pheno", scaled_pheno, "--pheno-name", PHENOTYPE,
        "--covar", scaled_pheno, "--covar-name", "PC1,PC2,PC3",
        "--read-freq", frequency, "--threads", "16", "--memory", "4000",
        "--out", output_prefix,
    ]
    run(command, env, "PLINK PC1-PC3 reference")
    return pathlib.Path(str(output_prefix) + f".{PHENOTYPE}.glm.linear")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ukc", help="UKBioCoin executable")
    parser.add_argument("adapter", help="UKBioCoinNSSAdapter executable")
    parser.add_argument("rscript", help="Rscript from the UKC demo environment")
    args = parser.parse_args()

    ukc = resolve_executable(args.ukc)
    adapter = resolve_executable(args.adapter)
    rscript = resolve_executable(args.rscript)
    repository = pathlib.Path(__file__).resolve().parents[1]
    demo = repository / "demo"
    demo_script = demo / "script.R"
    pfile = demo / "euro_10Ksnp"
    pheno = demo / "euro_10Ksnp.all.phe"
    for required in (demo_script, pheno, pathlib.Path(str(pfile) + ".pgen"),
                     pathlib.Path(str(pfile) + ".pvar"),
                     pathlib.Path(str(pfile) + ".psam")):
        if not required.is_file():
            raise SystemExit(f"required demo input is missing: {required}")

    env = os.environ.copy()
    env["PATH"] = str(pathlib.Path(rscript).parent) + os.pathsep + env.get("PATH", "")
    env["OMP_THREAD_LIMIT"] = "16"
    plink2 = shutil.which("plink2", path=env["PATH"])
    if plink2 is None:
        raise SystemExit("plink2 was not found next to Rscript or on PATH")

    with tempfile.TemporaryDirectory(prefix="ukc-real-demo-test-") as temporary:
        root = pathlib.Path(temporary)
        v2_root = root / "generated_v2"
        legacy_root = root / "generated_legacy"
        generate_demo(
            rscript, demo_script, pfile, pheno, ukc, v2_root, "fresh",
            "v2", True, env,
        )
        generate_demo(
            rscript, demo_script, pfile, pheno, ukc, legacy_root, "legacy",
            "legacy", False, env,
        )

        fresh_v2 = v2_root / "2.nss" / "fresh.nss"
        legacy_prefix = legacy_root / "2.matrix" / "legacy"
        legacy_metadata = pathlib.Path(str(legacy_prefix) + "_meta.table")
        adapted_v2 = root / "adapted.nss"
        run(
            [adapter, "--input-prefix", legacy_prefix, "--output", adapted_v2,
             "--sample-count", "200", "--threads", "16"],
            env, "adapt legacy demo",
        )

        fresh_hash, variant_ids, a1 = read_metadata(fresh_v2 / "meta.tsv")
        legacy_hash, legacy_ids, legacy_a1 = read_metadata(legacy_metadata)
        adapted_hash, adapted_ids, adapted_a1 = read_metadata(adapted_v2 / "meta.tsv")
        if not (fresh_hash == legacy_hash == adapted_hash):
            raise AssertionError("fresh, legacy, and adapted metadata differ")
        if variant_ids != legacy_ids or variant_ids != adapted_ids or \
                a1 != legacy_a1 or a1 != adapted_a1:
            raise AssertionError("metadata ID/A1 rows differ between NSS formats")
        source_ids = pvar_ids(pathlib.Path(str(pfile) + ".pvar"))
        if len(source_ids) != VARIANT_COUNT or variant_ids != source_ids:
            raise AssertionError("metadata row order differs from the source pvar")

        phenotype_names = pheno.read_text(encoding="utf-8").splitlines()[0].split("\t")[2:]
        if len(phenotype_names) < 4 or not all(
                name in phenotype_names for name in (PHENOTYPE, "PC1", "PC2", "PC3")):
            raise AssertionError("demo phenotype header is unexpected")
        fresh_manifest = verify_v2_layout(fresh_v2, len(phenotype_names), fresh_hash)
        verify_v2_layout(adapted_v2, len(phenotype_names), fresh_hash)
        compare_v2_payloads(fresh_v2, adapted_v2, len(phenotype_names))

        glm_directory = v2_root / "1.plink_temp"
        no_cov_plink = None
        for phenotype_name in phenotype_names:
            glm_path = glm_directory / f"single_reg.{phenotype_name}.glm.linear"
            collected = verify_glm(
                glm_path, variant_ids, a1,
                collect_statistics=phenotype_name == PHENOTYPE,
            )
            if phenotype_name == PHENOTYPE:
                no_cov_plink = collected
        if no_cov_plink is None:
            raise AssertionError("did not collect the primary PLINK reference")

        plot = v2_root / "3.analysis" / f"Validation.{PHENOTYPE}.png"
        visual_result = v2_root / "3.analysis" / f"test.{PHENOTYPE}_results.table"
        if not plot.is_file() or plot.stat().st_size < 10_000 or \
                plot.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
            raise AssertionError("demo visualization PNG was not generated correctly")
        if not visual_result.is_file():
            raise AssertionError("visual demo did not write its UKC result table")
        if (legacy_root / "3.analysis").exists():
            raise AssertionError("--novisualize legacy generation unexpectedly made plots")

        covariate_glm = make_covariate_plink_reference(
            plink2, pfile, v2_root / "scaled_pheno.table",
            v2_root / "0.basic" / "0.freq.afreq",
            root / "plink_reference" / "pc1_pc3", env,
        )
        pc_plink = verify_glm(
            covariate_glm, variant_ids, a1, collect_statistics=True
        )

        runs_root = root / "regressions"
        no_cov_table, no_cov_npy = exercise_model(
            ukc, ("no_covariates", None), fresh_v2, legacy_prefix, adapted_v2,
            fresh_hash, visual_result, no_cov_plink, runs_root, env,
        )
        exercise_model(
            ukc, ("pc1_pc3", "PC1,PC2,PC3"), fresh_v2, legacy_prefix,
            adapted_v2, fresh_hash, None, pc_plink, runs_root, env,
        )
        exercise_gc(
            ukc, no_cov_table, no_cov_npy, fresh_hash, root / "adjustment", env
        )

        # The manifest's declared names/order should also match the phenotype table.
        manifest_names = [item["name"] for item in fresh_manifest["phenotypes"]]
        if manifest_names != phenotype_names:
            raise AssertionError("manifest phenotype order differs from phenotype input")

    print(
        "\nAll 100,000-variant demo tests passed. "
        "Cross-format tolerance allows float32 COVxy/text round trips; "
        "PLINK SE/T/P tolerance allows the demo's sample-SD covariance normalization.",
        flush=True,
    )


if __name__ == "__main__":
    main()

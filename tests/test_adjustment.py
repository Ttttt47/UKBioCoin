#!/usr/bin/env python3
"""End-to-end tests for UKBioCoin's genomic-control adjustment mode.

Usage:
    python3 tests/test_adjustment.py /path/to/UKBioCoin

The numeric reference implementation below uses the identity
P(ChiSq(1) >= z^2) = erfc(z / sqrt(2)).  For extreme z, it evaluates the
normal-tail asymptotic expansion in log space so the reference itself does not
underflow.
"""

import math
import pathlib
import subprocess
import sys
import tempfile

import numpy as np


CANONICAL_METADATA = [
    "#CHROM", "POS", "ID", "REF", "ALT1", "A1", "REF_FREQ", "EXTRA"
]
CANONICAL_STATISTICS = [
    "BETA", "SE", "T_STAT", "NEG_LOG10_P", "VIF", "OBS_CT",
    "QUALITY_SCORE",
]
OLD_TO_CANONICAL = {
    "CHROM": "#CHROM",
    "REF_Allele": "REF",
    "ALT_Allele": "ALT1",
    "Effect_Allele": "A1",
    "T-STAT": "T_STAT",
    "-log10_P": "NEG_LOG10_P",
    "nobs": "OBS_CT",
    "Quality-Score": "QUALITY_SCORE",
}


def write_text(path, text):
    pathlib.Path(path).write_text(text, encoding="utf-8")


def combined_output(result):
    return result.stdout + result.stderr


def expect_failure(command, contains=None):
    result = subprocess.run(
        [str(item) for item in command], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode == 0:
        raise AssertionError(
            "command unexpectedly succeeded: " + " ".join(map(str, command))
        )
    if contains is not None and contains.lower() not in combined_output(result).lower():
        raise AssertionError(
            f"failure did not contain {contains!r}:\n{combined_output(result)}"
        )
    return result


def invoke_adjust(binary, input_path, output_prefix, output_format=None,
                  lambda_value=None, extra=(), expect_ok=True,
                  error_contains=None):
    command = [
        binary, "--adjust-file", str(input_path), "--out", str(output_prefix)
    ]
    if output_format is not None:
        command += ["--output-format", output_format]
    if lambda_value is not None:
        command += ["--lambda", str(lambda_value)]
    command += list(extra)
    if not expect_ok:
        return expect_failure(command, error_contains)
    result = subprocess.run(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        raise AssertionError(combined_output(result))
    return result


def output_path(prefix, output_format="table"):
    suffix = ".adjusted.npy" if output_format == "npy" else ".adjusted"
    return pathlib.Path(str(prefix) + suffix)


def assert_no_temporary_output(prefix):
    prefix = pathlib.Path(prefix)
    for final_suffix in (".adjusted", ".adjusted.npy"):
        final_name = prefix.name + final_suffix
        matches = list(prefix.parent.glob(final_name + ".tmp.*"))
        if matches:
            raise AssertionError(f"temporary adjustment output was left behind: {matches}")


def parse_log_value(result, key):
    prefix = key + ":"
    values = [
        line[len(prefix):].strip()
        for line in result.stdout.splitlines()
        if line.startswith(prefix)
    ]
    if len(values) != 1:
        raise AssertionError(f"expected one {key!r} log line, got {values}")
    return float(values[0])


def assert_logged_lambda(result, raw, used):
    # The command-line log intentionally uses a compact human-readable format.
    if not math.isclose(parse_log_value(result, "lambda_raw"), raw,
                        rel_tol=1e-5, abs_tol=1e-7):
        raise AssertionError("logged lambda_raw is incorrect")
    if not math.isclose(parse_log_value(result, "lambda_used"), used,
                        rel_tol=1e-5, abs_tol=1e-7):
        raise AssertionError("logged lambda_used is incorrect")


def chi1_neg_log10_from_t(t_statistic, lambda_used=1.0):
    z = abs(float(t_statistic)) / math.sqrt(lambda_used)
    if z < 26.0:
        return -math.log10(math.erfc(z / math.sqrt(2.0)))

    # 2 Phi(-z) = 2 phi(z) / z *
    #   (1 - z^-2 + 3 z^-4 - 15 z^-6 + ...).
    inverse_square = 1.0 / (z * z)
    term = 1.0
    series = 1.0
    previous = abs(term)
    for order in range(1, 10000):
        term *= -(2.0 * order - 1.0) * inverse_square
        magnitude = abs(term)
        if magnitude >= previous:
            break
        series += term
        if magnitude <= sys.float_info.epsilon * abs(series):
            break
        previous = magnitude
    log_survival = (
        math.log(2.0) - 0.5 * z * z - 0.5 * math.log(2.0 * math.pi)
        - math.log(z) + math.log(series)
    )
    return -log_survival / math.log(10.0)


def value_text(value):
    if isinstance(value, str):
        return value
    return format(float(value), ".17g")


def make_record(index, t_statistic, neg_log10_p):
    ref = "A" if index % 2 == 0 else "C"
    alt = "G" if index % 2 == 0 else "T"
    record = {
        "#CHROM": str(1 + index // 4),
        "POS": str(1000 + index),
        "ID": f"rs{index + 1}",
        "REF": ref,
        "ALT1": alt,
        "A1": alt,
        "REF_FREQ": format(0.1 + index * 0.01, ".17g"),
        "EXTRA": f"tag {index}",
        "BETA": format((index + 1) * 0.01, ".17g"),
        "SE": format(0.1 + index * 0.001, ".17g"),
        "T_STAT": value_text(t_statistic),
        "NEG_LOG10_P": value_text(neg_log10_p),
        "VIF": "0.99",
        "OBS_CT": str(1000 - index),
        "QUALITY_SCORE": "1",
    }
    return record


def write_result_table(path, records, legacy=False):
    if legacy:
        names = [
            "EXTRA", "Effect_Allele", "ID", "CHROM", "REF_FREQ",
            "ALT_Allele", "POS", "REF_Allele", "Quality-Score", "nobs",
            "VIF", "-log10_P", "T-STAT", "SE", "BETA",
        ]
    else:
        names = CANONICAL_METADATA + CANONICAL_STATISTICS

    def canonical(name):
        return OLD_TO_CANONICAL.get(name, name)

    lines = ["\t".join(names)]
    for record in records:
        lines.append("\t".join(record[canonical(name)] for name in names))
    write_text(path, "\n".join(lines) + "\n")


def table_values(path):
    lines = pathlib.Path(path).read_text(encoding="utf-8").splitlines()
    if not lines:
        raise AssertionError("adjustment table is empty")
    header = lines[0].split("\t")
    if len(header) < 2 or header[-2:] != ["NEG_LOG10_UNADJ", "NEG_LOG10_GC"]:
        raise AssertionError(f"unexpected adjustment header: {header}")
    values = np.empty((len(lines) - 1, 2), dtype=np.dtype("<f8"))
    rows = []
    for row_index, line in enumerate(lines[1:]):
        fields = line.split("\t")
        if len(fields) != len(header):
            raise AssertionError(f"row width differs from header: {fields}")
        rows.append(fields)
        for column in range(2):
            text = fields[len(fields) - 2 + column]
            values[row_index, column] = math.nan if text == "NA" else float(text)
    return header, rows, values


def input_number(value):
    if not isinstance(value, str):
        return float(value)
    lowered = value.lower()
    if lowered in ("na", "nan"):
        return math.nan
    if lowered in ("inf", "+inf"):
        return math.inf
    if lowered == "-inf":
        return -math.inf
    return float(value)


def valid_pair(t_statistic, neg_log10_p):
    return (
        math.isfinite(t_statistic) and math.isfinite(neg_log10_p)
        and neg_log10_p >= 0.0
    )


def expected_matrix(pairs, lambda_used):
    result = np.full((len(pairs), 2), math.nan, dtype=np.dtype("<f8"))
    for row, (t_value, p_value) in enumerate(pairs):
        t_statistic = input_number(t_value)
        neg_log10_p = input_number(p_value)
        if valid_pair(t_statistic, neg_log10_p):
            result[row, 0] = neg_log10_p
            result[row, 1] = chi1_neg_log10_from_t(t_statistic, lambda_used)
    return result


def make_result_npy(path, pairs):
    array = np.zeros((len(pairs), 7), dtype=np.dtype("<f8"))
    for row, (t_value, p_value) in enumerate(pairs):
        array[row] = [0.01 * (row + 1), 0.1, input_number(t_value),
                      input_number(p_value), 0.99, 1000 - row, 1.0]
    np.save(path, array, allow_pickle=False)
    return array


def assert_adjusted_values(actual, expected):
    if actual.shape != expected.shape:
        raise AssertionError(f"shape mismatch: {actual.shape} != {expected.shape}")
    np.testing.assert_array_equal(np.isnan(actual), np.isnan(expected))
    np.testing.assert_allclose(actual, expected, rtol=3e-13, atol=3e-13,
                               equal_nan=True)


def test_header_compatibility_and_four_formats(binary, root):
    pairs = [
        (0.3, chi1_neg_log10_from_t(0.3)),
        (-1.0, chi1_neg_log10_from_t(1.0)),
        (2.0, chi1_neg_log10_from_t(2.0)),
        (-3.0, chi1_neg_log10_from_t(3.0)),
        ("NA", 1.0),
        ("nan", 1.0),
        ("inf", 1.0),
        (2.0, "NA"),
        (2.0, -0.1),
    ]
    records = [make_record(index, *pair) for index, pair in enumerate(pairs)]
    canonical_input = root / "canonical.table"
    legacy_input = root / "legacy.table"
    npy_input = root / "input.npy"
    write_result_table(canonical_input, records)
    write_result_table(legacy_input, records, legacy=True)
    make_result_npy(npy_input, pairs)

    valid_chi_squares = [
        input_number(t) ** 2 for t, p in pairs
        if valid_pair(input_number(t), input_number(p))
    ]
    lambda_raw = float(np.median(np.asarray(valid_chi_squares))) / 0.456
    lambda_used = max(1.0, lambda_raw)
    expected = expected_matrix(pairs, lambda_used)

    canonical_prefix = root / "canonical_to_table"
    legacy_prefix = root / "legacy_to_table"
    table_npy_prefix = root / "table_to_npy"
    npy_table_prefix = root / "npy_to_table"
    npy_npy_prefix = root / "npy_to_npy"

    canonical_run = invoke_adjust(binary, canonical_input, canonical_prefix)
    legacy_run = invoke_adjust(binary, legacy_input, legacy_prefix)
    table_npy_run = invoke_adjust(
        binary, canonical_input, table_npy_prefix, output_format="npy"
    )
    npy_table_run = invoke_adjust(binary, npy_input, npy_table_prefix)
    npy_npy_run = invoke_adjust(
        binary, npy_input, npy_npy_prefix, output_format="npy"
    )
    for result in (canonical_run, legacy_run, table_npy_run, npy_table_run,
                   npy_npy_run):
        assert_logged_lambda(result, lambda_raw, lambda_used)

    canonical_bytes = output_path(canonical_prefix).read_bytes()
    legacy_bytes = output_path(legacy_prefix).read_bytes()
    if canonical_bytes != legacy_bytes:
        raise AssertionError("new and legacy table aliases produced different output")

    header, rows, canonical_values = table_values(output_path(canonical_prefix))
    expected_header = CANONICAL_METADATA + ["NEG_LOG10_UNADJ", "NEG_LOG10_GC"]
    if header != expected_header:
        raise AssertionError(f"metadata was not canonicalized: {header}")
    for index, row in enumerate(rows):
        if row[:8] != [records[index][name] for name in CANONICAL_METADATA]:
            raise AssertionError("metadata values or row order changed")
    # All seven regression columns must be removed before adjusted columns are added.
    for name in CANONICAL_STATISTICS:
        if name in header:
            raise AssertionError(f"input regression column was retained: {name}")

    npy_header, npy_rows, npy_table_values = table_values(
        output_path(npy_table_prefix)
    )
    if npy_header != ["#ROW_INDEX", "NEG_LOG10_UNADJ", "NEG_LOG10_GC"]:
        raise AssertionError(f"unexpected NPY-to-table header: {npy_header}")
    if [row[0] for row in npy_rows] != [str(i) for i in range(len(pairs))]:
        raise AssertionError("NPY-to-table row indices are incorrect")

    table_npy_values = np.load(
        output_path(table_npy_prefix, "npy"), allow_pickle=False
    )
    npy_npy_values = np.load(
        output_path(npy_npy_prefix, "npy"), allow_pickle=False
    )
    for array in (table_npy_values, npy_npy_values):
        if array.dtype.str != "<f8" or array.shape != (len(pairs), 2):
            raise AssertionError(
                f"adjusted NPY violates the <f8 N x 2 contract: "
                f"{array.dtype.str}, {array.shape}"
            )

    for actual in (canonical_values, npy_table_values, table_npy_values,
                   npy_npy_values):
        assert_adjusted_values(actual, expected)
    if not np.array_equal(canonical_values, npy_table_values, equal_nan=True):
        raise AssertionError("table and NPY input paths produced different numbers")
    if not np.array_equal(table_npy_values, npy_npy_values, equal_nan=True):
        raise AssertionError("table and NPY output paths produced different numbers")

    # The five deliberately invalid rows remain in place and are excluded from
    # the lambda median.  Table uses NA; NPY uses NaN.
    invalid_indices = [4, 5, 6, 7, 8]
    for index in invalid_indices:
        if rows[index][-2:] != ["NA", "NA"]:
            raise AssertionError("invalid table row was not retained as NA")
    if not np.isnan(table_npy_values[invalid_indices]).all():
        raise AssertionError("invalid NPY rows were not retained as NaN")


def write_simple_table(path, t_values, p_values=None):
    if p_values is None:
        p_values = [chi1_neg_log10_from_t(value) for value in t_values]
    lines = ["ID\tT_STAT\tNEG_LOG10_P"]
    for index, (t_value, p_value) in enumerate(zip(t_values, p_values)):
        lines.append(f"rs{index + 1}\t{value_text(t_value)}\t{value_text(p_value)}")
    write_text(path, "\n".join(lines) + "\n")


def test_lambda_medians_clamp_and_manual(binary, root):
    odd_input = root / "odd.table"
    odd_t = [3.0, -1.0, 2.0]
    write_simple_table(odd_input, odd_t)
    odd_run = invoke_adjust(binary, odd_input, root / "odd")
    odd_raw = 4.0 / 0.456
    assert_logged_lambda(odd_run, odd_raw, odd_raw)
    _, _, odd_values = table_values(output_path(root / "odd"))
    # Rounded reference from PLINK v2.0.0-a.6.9LM (29 Jan 2025), using an
    # input P consistent with the one-df chi-square tail of T^2. PLINK sorts
    # its report by P; these values are restored to the input rs1/rs2/rs3 order.
    plink_reference = np.asarray([
        [2.56867, 0.5071],
        [0.498516, 0.133337],
        [1.34199, 0.301465],
    ])
    np.testing.assert_allclose(odd_values, plink_reference,
                               rtol=2e-5, atol=2e-6)
    if np.any(odd_values[:, 1] > odd_values[:, 0] + 1e-13):
        raise AssertionError(
            "automatic lambda inflated a chi-square-consistent P value")

    even_input = root / "even.table"
    even_t = [4.0, -1.0, 3.0, 2.0]
    write_simple_table(even_input, even_t)
    even_run = invoke_adjust(binary, even_input, root / "even")
    even_raw = ((4.0 + 9.0) / 2.0) / 0.456
    assert_logged_lambda(even_run, even_raw, even_raw)

    clamp_input = root / "clamp.table"
    clamp_t = [0.1, 0.3, 0.2]
    write_simple_table(clamp_input, clamp_t)
    clamp_prefix = root / "clamp"
    clamp_run = invoke_adjust(binary, clamp_input, clamp_prefix)
    clamp_raw = 0.04 / 0.456
    assert_logged_lambda(clamp_run, clamp_raw, 1.0)
    _, _, clamp_values = table_values(output_path(clamp_prefix))
    expected_clamp = expected_matrix(
        list(zip(clamp_t, [chi1_neg_log10_from_t(t) for t in clamp_t])), 1.0
    )
    assert_adjusted_values(clamp_values, expected_clamp)
    np.testing.assert_allclose(clamp_values[:, 1], clamp_values[:, 0],
                               rtol=3e-13, atol=3e-13)

    # A manually supplied value is used directly, even when it is below one.
    manual_lambda = 0.25
    manual_prefix = root / "manual"
    manual_run = invoke_adjust(
        binary, clamp_input, manual_prefix, lambda_value=manual_lambda
    )
    assert_logged_lambda(manual_run, clamp_raw, manual_lambda)
    _, _, manual_values = table_values(output_path(manual_prefix))
    expected_manual = expected_matrix(
        list(zip(clamp_t, [chi1_neg_log10_from_t(t) for t in clamp_t])),
        manual_lambda,
    )
    assert_adjusted_values(manual_values, expected_manual)
    if not np.all(manual_values[:, 1] > clamp_values[:, 1]):
        raise AssertionError("manual lambda below one was apparently clamped")


def test_extreme_log_space_values(binary, root):
    t_values = [10.0, 20.0, 40.0, 50.0, 100.0]
    reference = np.asarray([
        22.81702340982209,
        88.25906534741159,
        349.13597646368186,
        544.6653058663327,
        2173.5705128733703,
    ])
    # Cross-check the fixed references against the independent implementation.
    np.testing.assert_allclose(
        reference, [chi1_neg_log10_from_t(value) for value in t_values],
        rtol=2e-15, atol=2e-13,
    )
    extreme_input = root / "extreme.table"
    write_simple_table(extreme_input, t_values, reference)
    extreme_prefix = root / "extreme"
    result = invoke_adjust(
        binary, extreme_input, extreme_prefix, lambda_value=1.0
    )
    expected_raw = float(np.median(np.square(t_values))) / 0.456
    assert_logged_lambda(result, expected_raw, 1.0)
    _, _, values = table_values(output_path(extreme_prefix))
    if not np.isfinite(values).all():
        raise AssertionError("extreme log-space GC values underflowed or overflowed")
    np.testing.assert_allclose(values[:, 0], reference, rtol=2e-15, atol=2e-13)
    np.testing.assert_allclose(values[:, 1], reference, rtol=3e-13, atol=3e-13)


def save_npy(path, array):
    np.save(path, array, allow_pickle=False)
    if not pathlib.Path(path).exists():
        raise AssertionError(f"numpy unexpectedly changed output path {path}")


def test_strict_npy_validation(binary, root):
    base = np.zeros((2, 7), dtype=np.dtype("<f8"))
    base[:, 2] = [1.0, 2.0]
    base[:, 3] = [chi1_neg_log10_from_t(1.0), chi1_neg_log10_from_t(2.0)]

    cases = []
    wrong_dtype = root / "wrong_dtype.npy"
    save_npy(wrong_dtype, base.astype(np.float32))
    cases.append((wrong_dtype, "dtype"))

    big_endian = root / "big_endian.npy"
    save_npy(big_endian, base.astype(np.dtype(">f8")))
    cases.append((big_endian, "dtype"))

    one_dimensional = root / "one_dimensional.npy"
    save_npy(one_dimensional, np.arange(7, dtype=np.dtype("<f8")))
    cases.append((one_dimensional, "two-dimensional"))

    three_dimensional = root / "three_dimensional.npy"
    save_npy(three_dimensional, np.zeros((1, 1, 7), dtype=np.dtype("<f8")))
    cases.append((three_dimensional, "two-dimensional"))

    six_columns = root / "six_columns.npy"
    save_npy(six_columns, np.zeros((2, 6), dtype=np.dtype("<f8")))
    cases.append((six_columns, "7 columns"))

    eight_columns = root / "eight_columns.npy"
    save_npy(eight_columns, np.zeros((2, 8), dtype=np.dtype("<f8")))
    cases.append((eight_columns, "7 columns"))

    fortran_order = root / "fortran_order.npy"
    save_npy(fortran_order, np.asfortranarray(base))
    cases.append((fortran_order, "fortran"))

    good_for_damage = root / "good_for_damage.npy"
    save_npy(good_for_damage, base)
    good_bytes = good_for_damage.read_bytes()
    truncated = root / "truncated.npy"
    truncated.write_bytes(good_bytes[:-8])
    cases.append((truncated, "file size"))

    trailing = root / "trailing.npy"
    trailing.write_bytes(good_bytes + b"trailing")
    cases.append((trailing, "file size"))

    corrupt_magic = root / "corrupt_magic.npy"
    corrupt = bytearray(good_bytes)
    corrupt[0] = 0
    corrupt_magic.write_bytes(corrupt)
    cases.append((corrupt_magic, "npy"))

    for index, (path, error_text) in enumerate(cases):
        prefix = root / f"strict_failure_{index}"
        invoke_adjust(
            binary, path, prefix, expect_ok=False, error_contains=error_text
        )
        if output_path(prefix).exists() or output_path(prefix, "npy").exists():
            raise AssertionError(f"invalid NPY produced an output: {path}")
        assert_no_temporary_output(prefix)

    all_invalid = root / "all_invalid.npy"
    invalid = base.copy()
    invalid[:, 2] = math.nan
    invalid[:, 3] = math.nan
    save_npy(all_invalid, invalid)
    prefix = root / "no_valid_npy"
    invoke_adjust(
        binary, all_invalid, prefix, expect_ok=False,
        error_contains="no valid tests",
    )
    if output_path(prefix).exists():
        raise AssertionError("all-invalid NPY unexpectedly produced an output")
    assert_no_temporary_output(prefix)


def test_table_errors_cli_and_atomicity(binary, root):
    valid_input = root / "cli_valid.table"
    write_simple_table(valid_input, [1.0, 2.0, 3.0])

    expect_failure(
        [binary, "--adjust-file", valid_input], "requires --out"
    )
    for option, value in (("--file", "nss"), ("--phe", "Y"),
                          ("--covar", "C1")):
        prefix = root / ("conflict_" + option[2:])
        invoke_adjust(
            binary, valid_input, prefix, extra=(option, value),
            expect_ok=False, error_contains="mutually exclusive",
        )
        assert_no_temporary_output(prefix)

    regression_option_prefix = root / "regression_option"
    invoke_adjust(
        binary, valid_input, regression_option_prefix,
        extra=("--threads", "2"), expect_ok=False,
        error_contains="regression-only",
    )
    assert_no_temporary_output(regression_option_prefix)

    invalid_format_prefix = root / "invalid_format"
    invoke_adjust(
        binary, valid_input, invalid_format_prefix, output_format="csv",
        expect_ok=False, error_contains="output-format",
    )
    assert_no_temporary_output(invalid_format_prefix)

    for index, value in enumerate(("0", "-1", "nan", "inf", "1e309", "bad")):
        prefix = root / f"invalid_lambda_{index}"
        invoke_adjust(
            binary, valid_input, prefix, lambda_value=value,
            expect_ok=False,
        )
        if output_path(prefix).exists() or output_path(prefix, "npy").exists():
            raise AssertionError(f"invalid lambda {value!r} produced output")
        assert_no_temporary_output(prefix)

    expect_failure(
        [binary, "--lambda", "1", "--file", "unused", "--phe", "Y",
         "--out", root / "lambda_without_adjust"],
        "only valid with --adjust-file",
    )

    missing_t = root / "missing_t.table"
    write_text(missing_t, "ID\tNEG_LOG10_P\nrs1\t1\n")
    missing_p = root / "missing_p.table"
    write_text(missing_p, "ID\tT_STAT\nrs1\t1\n")
    duplicate_t = root / "duplicate_t.table"
    write_text(
        duplicate_t,
        "ID\tT_STAT\tT-STAT\tNEG_LOG10_P\nrs1\t1\t1\t1\n",
    )
    bad_width = root / "bad_width.table"
    write_text(
        bad_width,
        "ID\tT_STAT\tNEG_LOG10_P\nrs1\t1\n",
    )
    for index, (path, error_text) in enumerate((
        (missing_t, "requires"),
        (missing_p, "requires"),
        (duplicate_t, "both map"),
        (bad_width, "field count"),
    )):
        prefix = root / f"bad_table_{index}"
        invoke_adjust(
            binary, path, prefix, expect_ok=False, error_contains=error_text
        )
        assert_no_temporary_output(prefix)

    # A failed scan must leave a pre-existing regular final byte-for-byte intact.
    malformed = root / "malformed_for_atomic.table"
    write_text(
        malformed,
        "ID\tT_STAT\tNEG_LOG10_P\nrs1\tbad\t1\n",
    )
    atomic_prefix = root / "atomic_regular"
    final = output_path(atomic_prefix)
    sentinel = b"pre-existing adjusted result\n"
    final.write_bytes(sentinel)
    invoke_adjust(
        binary, malformed, atomic_prefix, expect_ok=False,
        error_contains="invalid numeric",
    )
    if final.read_bytes() != sentinel:
        raise AssertionError("failed adjustment damaged the old final output")
    assert_no_temporary_output(atomic_prefix)

    # Force rename(2) to fail after the temporary file has been fully written.
    # The existing final directory and its marker stand in for an old artifact.
    commit_prefix = root / "atomic_commit"
    commit_final = output_path(commit_prefix)
    commit_final.mkdir()
    marker = commit_final / "sentinel"
    marker.write_bytes(sentinel)
    invoke_adjust(
        binary, valid_input, commit_prefix, expect_ok=False,
        error_contains="replacing",
    )
    if not commit_final.is_dir() or marker.read_bytes() != sentinel:
        raise AssertionError("failed atomic commit damaged the old final artifact")
    assert_no_temporary_output(commit_prefix)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_adjustment.py /path/to/UKBioCoin")
    binary = str(pathlib.Path(sys.argv[1]).resolve())
    if not pathlib.Path(binary).is_file():
        raise SystemExit(f"UKBioCoin binary does not exist: {binary}")

    with tempfile.TemporaryDirectory(prefix="ukc-adjustment-tests-") as tmp:
        root = pathlib.Path(tmp)
        test_header_compatibility_and_four_formats(binary, root)
        test_lambda_medians_clamp_and_manual(binary, root)
        test_extreme_log_space_values(binary, root)
        test_strict_npy_validation(binary, root)
        test_table_errors_cli_and_atomicity(binary, root)
    print("All UKBioCoin adjustment tests passed.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import math
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile

import numpy as np


def write(path, text):
    pathlib.Path(path).write_text(text, encoding="utf-8")


def make_fixture(root):
    prefix = root / "fixture"
    write(str(prefix) + "_cov_xy.table", '''"Y" "C1" "junk"
"1" 0.2 0.1 9
"2" -0.1 0.05 8
"3" 0.03 -0.02 7
''')
    write(str(prefix) + "_cov_yy.table", '''"Y" "C1" "junk"
"Y" 1 0.2 0
"C1" 0.2 1 0
"junk" 0 0 1
''')
    write(str(prefix) + "_var_x.table", '''"var_x"
"1" 0.5
"2" 0.4
"3" 0.2
''')
    write(str(prefix) + "_meta.table", '''ID A
rs1 A
rs2 C
rs3 G
''')
    write(str(prefix) + "_y_missing.table", '''"" "missing_rate"
"Y" 0.1
"C1" 0.2
"junk" 0
''')
    write(str(prefix) + "_x_missing.table", '''"" "x_missing"
"1" 0.01
"2" 0.02
"3" 0.03
''')
    return prefix


def invoke(binary, prefix, output, extra, expect_ok=True, totalsize=True):
    command = [binary, "--file", str(prefix), "--phe", "Y",
               "--out", str(output)]
    if totalsize:
        command += ["--totalsize", "100"]
    command += extra
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    if expect_ok and result.returncode != 0:
        raise AssertionError(result.stdout + result.stderr)
    if not expect_ok and result.returncode == 0:
        raise AssertionError("command unexpectedly succeeded: " + " ".join(command))
    return result


def result_rows(output):
    lines = pathlib.Path(str(output) + "_results.table").read_text().splitlines()
    return [line.split() for line in lines]


def close(actual, expected, tolerance=5e-6):
    if not math.isclose(float(actual), expected, rel_tol=tolerance, abs_tol=tolerance):
        raise AssertionError(f"{actual} != {expected}")


def main():
    binary = sys.argv[1]
    adapter = sys.argv[2]
    no_arguments = subprocess.run([binary], text=True, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    if no_arguments.returncode != 0 or "No parameter is given" not in no_arguments.stdout:
        raise AssertionError("empty command line is not handled cleanly")

    with tempfile.TemporaryDirectory(prefix="ukc-tests-") as tmp:
        root = pathlib.Path(tmp)
        prefix = make_fixture(root)

        no_cov_1 = root / "no_cov_1"
        no_cov_16 = root / "no_cov_16"
        invoke(binary, prefix, no_cov_1,
               ["--overall-non-missing-rate", "1", "--threads", "1"])
        invoke(binary, prefix, no_cov_16,
               ["--overall-non-missing-rate", "1", "--threads", "16"])
        if pathlib.Path(str(no_cov_1) + "_results.table").read_bytes() != \
                pathlib.Path(str(no_cov_16) + "_results.table").read_bytes():
            raise AssertionError("1-thread and 16-thread no-covariate outputs differ")
        expected_no_covariate = (
            "ID A BETA SE T-STAT -log10_P VIF nobs Quality-Score\n"
            "rs1 A 0.4 0.13633 2.93406 2.3802 0.99  100 1\n"
            "rs2 C -0.25 0.156911 -1.59326 0.941865 0.99  100 1\n"
            "rs3 G 0.15 0.224227 0.668965 0.296631 0.99  100 1\n")
        if pathlib.Path(str(no_cov_1) + "_results.table").read_text() != \
                expected_no_covariate:
            raise AssertionError("optimized output formatting changed the legacy table bytes")
        rows = result_rows(no_cov_1)
        if len(rows) != 4:
            raise AssertionError("no-covariate output lost rows")
        beta = 0.2 / 0.5
        rss = 100.0 * (1.0 - beta * 0.2) / 99.0
        var_beta = rss * (1.0 / 0.5) / 100.0
        close(rows[1][2], beta)
        close(rows[1][3], math.sqrt(var_beta))
        close(rows[1][4], beta / math.sqrt(var_beta))
        close(rows[1][6], 0.99)

        with_cov_1 = root / "with_cov_1"
        with_cov_16 = root / "with_cov_16"
        invoke(binary, prefix, with_cov_1,
               ["--covar", "C1", "--use-missing-rate-estimate", "--threads", "1"])
        invoke(binary, prefix, with_cov_16,
               ["--covar", "C1", "--use-missing-rate-estimate", "--threads", "16"])
        if pathlib.Path(str(with_cov_1) + "_results.table").read_bytes() != \
                pathlib.Path(str(with_cov_16) + "_results.table").read_bytes():
            raise AssertionError("1-thread and 16-thread covariate outputs differ")

        wide_names = " ".join(f'\"W{i}\"' for i in range(20))
        wide_values = " ".join(str(i) for i in range(20))
        write(str(prefix) + "_cov_xy.table",
              f'\"junk\" \"Y\" {wide_names} \"C1\"\n'
              f'\"1\" 9 0.2 {wide_values} 0.1\n'
              f'\"2\" 8 -0.1 {wide_values} 0.05\n'
              f'\"3\" 7 0.03 {wide_values} -0.02\n')
        reordered = root / "reordered_wide"
        invoke(binary, prefix, reordered,
               ["--covar", "C1", "--use-missing-rate-estimate", "--threads", "16"])
        if pathlib.Path(str(with_cov_1) + "_results.table").read_bytes() != \
                pathlib.Path(str(reordered) + "_results.table").read_bytes():
            raise AssertionError("reordered/wide COVXY columns changed results")

        write(str(prefix) + "_cov_xy.table", '''"Y" "C1" "junk"
"1" NA 0.1 9
"2" -0.1 0.05 8
"3" 0.03 -0.02 7
''')
        na_output = root / "na"
        invoke(binary, prefix, na_output,
               ["--overall-non-missing-rate", "1", "--threads", "1"])
        if "nan" not in pathlib.Path(str(na_output) + "_results.table").read_text().lower():
            raise AssertionError("NA row was not preserved as a non-finite result")

        prefix = make_fixture(root)
        missing = invoke(binary, prefix, root / "missing_col",
                         ["--covar", "ABSENT"], expect_ok=False)
        if "not found" not in missing.stderr:
            raise AssertionError("missing-column error lacks context")

        write(str(prefix) + "_cov_xy.table", '''"Y" "C1" "junk"
"1" bad 0.1 9
"2" -0.1 0.05 8
"3" 0.03 -0.02 7
''')
        malformed = invoke(binary, prefix, root / "malformed", [], expect_ok=False)
        if "invalid numeric" not in malformed.stderr:
            raise AssertionError("malformed-number error lacks context")

        prefix = make_fixture(root)
        write(str(prefix) + "_meta.table", "ID A\nrs1 A\nrs2 C\n")
        mismatch = invoke(binary, prefix, root / "mismatch", [], expect_ok=False)
        if "row-count mismatch" not in mismatch.stderr:
            raise AssertionError("row mismatch was not detected")

        prefix = make_fixture(root)
        write(str(prefix) + "_cov_xy.table", '''"Y" "C1" "C2"
"1" 0.2 0.1 0.1
"2" -0.1 0.05 0.05
"3" 0.03 -0.02 -0.02
''')
        write(str(prefix) + "_cov_yy.table", '''"Y" "C1" "C2"
"Y" 1 0.2 0.2
"C1" 0.2 1 1
"C2" 0.2 1 1
''')
        singular = invoke(binary, prefix, root / "singular",
                          ["--covar", "C1,C2"], expect_ok=False)
        if "Singularity" not in singular.stderr:
            raise AssertionError("singular covariance matrix was not detected")

        prefix = make_fixture(root)
        invalid_threads = invoke(binary, prefix, root / "threads",
                                 ["--threads", "0"], expect_ok=False)
        if "--threads" not in invalid_threads.stderr:
            raise AssertionError("invalid thread count was not rejected")

        # Convert the fixture with the real C++ adapter, then exercise both v2 I/O modes.
        prefix = make_fixture(root)
        v2 = root / "fixture.nss"
        converted = subprocess.run(
            [adapter, "--input-prefix", str(prefix), "--output", str(v2),
             "--sample-count", "100", "--threads", "16"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if converted.returncode != 0:
            raise AssertionError(converted.stdout + converted.stderr)
        manifest = json.loads((v2 / "manifest.json").read_text())
        if manifest["dimensions"] != {
                "sample_count": 100, "variant_count": 3,
                "phenotype_count": 3}:
            raise AssertionError("adapter wrote incorrect manifest dimensions")
        if np.load(v2 / "cov_xy" / "00000000.npy").dtype != np.dtype("<f4"):
            raise AssertionError("COVxy is not interoperable float32 NPY")
        if np.load(v2 / "var_x.npy").dtype != np.dtype("<f8"):
            raise AssertionError("var_x is not interoperable float64 NPY")
        if np.load(v2 / "cov_yy.npy").shape != (3, 3):
            raise AssertionError("cov_yy NPY has the wrong shape")

        v2_stream_1 = root / "v2_stream_1"
        v2_stream_16 = root / "v2_stream_16"
        v2_memory_16 = root / "v2_memory_16"
        invoke(binary, v2, v2_stream_1,
               ["--covar", "C1", "--use-missing-rate-estimate",
                "--threads", "1"], totalsize=False)
        invoke(binary, v2, v2_stream_16,
               ["--covar", "C1", "--use-missing-rate-estimate",
                "--threads", "16"], totalsize=False)
        memory_run = invoke(binary, v2, v2_memory_16,
                            ["--covar", "C1", "--use-missing-rate-estimate",
                             "--threads", "16", "--io-mode", "memory"],
                            totalsize=False)
        reference_bytes = pathlib.Path(str(v2_stream_1) + "_results.table").read_bytes()
        for output in (v2_stream_16, v2_memory_16):
            if pathlib.Path(str(output) + "_results.table").read_bytes() != reference_bytes:
                raise AssertionError("v2 stream/memory or 1/16-thread outputs differ")
        if "data loading" not in memory_run.stdout or "computation" not in memory_run.stdout:
            raise AssertionError("memory mode did not report separate load/compute timing")

        v1_reference = root / "v1_reference"
        invoke(binary, prefix, v1_reference,
               ["--covar", "C1", "--use-missing-rate-estimate", "--threads", "1"])
        v1_rows = result_rows(v1_reference)
        v2_rows = result_rows(v2_stream_1)
        if [row[:2] for row in v1_rows] != [row[:2] for row in v2_rows]:
            raise AssertionError("v1/v2 metadata or SNP order differs")
        for old, new in zip(v1_rows[1:], v2_rows[1:]):
            for column in range(2, len(old)):
                close(new[column], float(old[column]), tolerance=2e-5)

        explicit = root / "v2_explicit"
        invoke(binary, v2, explicit, ["--threads", "1"], totalsize=True)
        v2_no_cov_16 = root / "v2_no_cov_16"
        v2_no_cov_memory = root / "v2_no_cov_memory"
        invoke(binary, v2, v2_no_cov_16,
               ["--overall-non-missing-rate", "1", "--threads", "16"],
               totalsize=False)
        invoke(binary, v2, v2_no_cov_memory,
               ["--overall-non-missing-rate", "1", "--threads", "16",
                "--io-mode", "memory"], totalsize=False)
        if pathlib.Path(str(v2_no_cov_16) + "_results.table").read_bytes() != \
                pathlib.Path(str(v2_no_cov_memory) + "_results.table").read_bytes():
            raise AssertionError("v2 no-covariate stream/memory outputs differ")
        rows = result_rows(v2_no_cov_16)
        quantized_cov = float(np.float32(0.2))
        beta = quantized_cov / 0.5
        rss = 100.0 * (1.0 - beta * quantized_cov) / 99.0
        expected_se = math.sqrt(rss * 2.0 / 100.0)
        close(rows[1][2], beta)
        close(rows[1][3], expected_se)
        mismatch = subprocess.run(
            [binary, "--file", str(v2), "--phe", "Y", "--out",
             str(root / "bad_n"), "--totalsize", "101"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if mismatch.returncode == 0 or "manifest sample_count 100" not in mismatch.stderr:
            raise AssertionError("v2 totalsize mismatch was not rejected clearly")

        missing_sample = subprocess.run(
            [adapter, "--input-prefix", str(prefix), "--output",
             str(root / "missing_sample.nss")], text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if missing_sample.returncode == 0 or "--sample-count" not in missing_sample.stderr:
            raise AssertionError("adapter did not require --sample-count")
        overwrite = subprocess.run(
            [adapter, "--input-prefix", str(prefix), "--output", str(v2),
             "--sample-count", "100"], text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        if overwrite.returncode == 0 or "--force" not in overwrite.stderr:
            raise AssertionError("adapter overwrote output without --force")

        def corrupt(name, callback):
            target = root / name
            shutil.copytree(v2, target)
            callback(target)
            result = invoke(binary, target, root / (name + "_out"), [],
                            expect_ok=False, totalsize=False)
            return result

        bad_sample = corrupt("bad_sample", lambda path: (
            (lambda doc: (doc["dimensions"].__setitem__("sample_count", 0),
                          (path / "manifest.json").write_text(json.dumps(doc))))(
                json.loads((path / "manifest.json").read_text()))))
        if "positive integer" not in bad_sample.stderr:
            raise AssertionError("invalid manifest sample_count was not rejected")

        def remove_sample_count(path):
            doc = json.loads((path / "manifest.json").read_text())
            del doc["dimensions"]["sample_count"]
            (path / "manifest.json").write_text(json.dumps(doc))
        missing_count = corrupt("missing_count", remove_sample_count)
        if "sample_count" not in missing_count.stderr:
            raise AssertionError("missing manifest sample_count was not rejected")

        def change_dtype(path):
            doc = json.loads((path / "manifest.json").read_text())
            doc["arrays"]["var_x"]["dtype"] = ">f8"
            (path / "manifest.json").write_text(json.dumps(doc))
        bad_dtype = corrupt("bad_dtype", change_dtype)
        if "dtype" not in bad_dtype.stderr:
            raise AssertionError("non-little-endian manifest dtype was not rejected")

        def change_shape(path):
            doc = json.loads((path / "manifest.json").read_text())
            doc["arrays"]["var_x"]["shape"] = [2]
            (path / "manifest.json").write_text(json.dumps(doc))
        bad_shape = corrupt("bad_shape", change_shape)
        if "shape" not in bad_shape.stderr:
            raise AssertionError("manifest shape mismatch was not rejected")

        def truncate(path):
            array = path / "var_x.npy"
            array.write_bytes(array.read_bytes()[:-1])
        truncated = corrupt("truncated", truncate)
        if "file size" not in truncated.stderr:
            raise AssertionError("truncated NPY was not rejected")

        def remove_vector(path):
            (path / "cov_xy" / "00000002.npy").unlink()
        absent_vector = corrupt("absent_vector", remove_vector)
        if "Error opening" not in absent_vector.stderr:
            raise AssertionError("missing phenotype vector was not rejected")

        nan_v2 = root / "nan_v2"
        shutil.copytree(v2, nan_v2)
        nan_values = np.load(nan_v2 / "cov_xy" / "00000000.npy")
        nan_values[0] = np.nan
        np.save(nan_v2 / "cov_xy" / "00000000.npy", nan_values,
                allow_pickle=False)
        nan_output = root / "nan_v2_output"
        invoke(binary, nan_v2, nan_output,
               ["--overall-non-missing-rate", "1"], totalsize=False)
        if "nan" not in pathlib.Path(str(nan_output) + "_results.table").read_text().lower():
            raise AssertionError("v2 NaN row was not preserved")

        duplicate_root = root / "duplicate"
        duplicate_root.mkdir()
        duplicate_prefix = make_fixture(duplicate_root)
        duplicate_cov = pathlib.Path(str(duplicate_prefix) + "_cov_xy.table")
        duplicate_cov.write_text(duplicate_cov.read_text().replace('"2" -0.1',
                                                                    '"1" -0.1'))
        duplicate_v2 = root / "duplicate.nss"
        duplicate_result = subprocess.run(
            [adapter, "--input-prefix", str(duplicate_prefix), "--output",
             str(duplicate_v2), "--sample-count", "100", "--threads", "16"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if duplicate_result.returncode != 0 or "row labels" not in duplicate_result.stderr:
            raise AssertionError("adapter did not tolerate/warn about duplicate row labels")

        invalid_mode = invoke(binary, v2, root / "mode", ["--io-mode", "cache"],
                              expect_ok=False, totalsize=False)
        if "--io-mode" not in invalid_mode.stderr:
            raise AssertionError("invalid io mode was not rejected")

    print("UKC v1/v2 and adapter tests passed")


if __name__ == "__main__":
    main()

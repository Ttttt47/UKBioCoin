#!/usr/bin/env python3
import math
import pathlib
import subprocess
import sys
import tempfile


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


def invoke(binary, prefix, output, extra, expect_ok=True):
    command = [binary, "--file", str(prefix), "--phe", "Y",
               "--totalsize", "100", "--out", str(output)] + extra
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

    print("candidate tests passed")


if __name__ == "__main__":
    main()

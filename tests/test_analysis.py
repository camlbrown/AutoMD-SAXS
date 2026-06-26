"""Unit tests for ``automd_saxs.analysis`` and the ``analyze`` CLI command.

No ATSAS/numpy/pytest required; reads captured fixtures. Runs standalone.
"""

import csv
import io
import json
import os
import shutil
import sys
import tempfile
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs import analysis, cli  # noqa: E402

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")


def _almost(a, b, tol=1e-6):
    assert abs(a - b) <= tol, "expected {0}, got {1}".format(b, a)


def _make_analysis_dir(tmp):
    """Populate a dir with the SAXS fixtures under the names analysis expects."""
    shutil.copy(os.path.join(FIX, "crysol_summary.txt"),
                os.path.join(tmp, "crysol_summary.txt"))
    shutil.copy(os.path.join(FIX, "profiles.fit"), os.path.join(tmp, "profiles_001_1.fit"))
    shutil.copy(os.path.join(FIX, "Rg_distr.txt"), os.path.join(tmp, "Rg_distr_001_1.txt"))
    return tmp


def test_collect_saxs_results_finds_all():
    with tempfile.TemporaryDirectory() as tmp:
        _make_analysis_dir(tmp)
        res = analysis.collect_saxs_results(tmp)
        assert len(res["records"]) == 3
        _almost(res["summary"]["bestChi2"], 0.98)
        assert res["summary"]["bestFrame"] == 1
        assert res["fit"] is not None and len(res["fit"]) == 3
        assert res["rgDistribution"] is not None and len(res["rgDistribution"]) == 3
        assert len(res["filesFound"]) == 3


def test_collect_handles_missing_files():
    with tempfile.TemporaryDirectory() as tmp:
        res = analysis.collect_saxs_results(tmp)
        assert res["records"] == []
        assert res["summary"]["bestChi2"] is None
        assert res["fit"] is None


def test_write_summary_csv_sorted_by_chi2():
    with tempfile.TemporaryDirectory() as tmp:
        _make_analysis_dir(tmp)
        res = analysis.collect_saxs_results(tmp)
        out = os.path.join(tmp, "s.csv")
        analysis.write_summary_csv(res["records"], out)
        rows = list(csv.reader(open(out)))
        assert rows[0] == ["rank", "frame", "model", "rg", "chi2"]
        # best chi2 (frame 1) ranked first
        assert rows[1][1] == "1"
        assert rows[1][4] == "0.98"


def test_build_analysis_manifest_metrics_and_outputs():
    with tempfile.TemporaryDirectory() as tmp:
        _make_analysis_dir(tmp)
        res = analysis.collect_saxs_results(tmp)
        m = analysis.build_analysis_manifest(res, summary_csv="s.csv", saxs_dat="exp.dat")
        data = m.to_dict()
        _almost(data["metrics"]["bestChi2"], 0.98)
        assert data["metrics"]["bestFrame"] == 1
        assert "s.csv" in data["outputs"]["summaryTables"]
        assert any(f.endswith(".fit") for f in data["outputs"]["saxsFits"])
        # ATSAS commands recorded for traceability
        assert [s["name"] for s in data["steps"]] == ["shanum", "crysol", "gajoe"]


def test_analyze_cli_end_to_end():
    with tempfile.TemporaryDirectory() as tmp:
        adir = os.path.join(tmp, "ana")
        os.makedirs(adir, exist_ok=True)
        _make_analysis_dir(adir)
        out = os.path.join(tmp, "out")
        buf = io.StringIO()
        with redirect_stdout(buf):
            code = cli.main(["analyze", "--analysis-dir", adir, "--out-dir", out,
                             "--saxs", "exp.dat"])
        assert code == 0
        assert os.path.isfile(os.path.join(out, "saxs_summary.csv"))
        data = json.load(open(os.path.join(out, "analysis_manifest.json")))
        _almost(data["metrics"]["bestChi2"], 0.98)
        assert "best chi^2    : 0.98" in buf.getvalue()


def _run_standalone():
    funcs = sorted((n, o) for n, o in globals().items()
                   if n.startswith("test_") and callable(o))
    failures = []
    for name, fn in funcs:
        try:
            fn()
            print("PASS", name)
        except Exception as exc:  # noqa: BLE001
            failures.append((name, exc))
            print("FAIL", name, "->", exc)
    print("\n{0}/{1} passed".format(len(funcs) - len(failures), len(funcs)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_run_standalone())

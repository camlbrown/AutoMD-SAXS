"""Unit tests for ``automd_saxs.structural`` and the CLoNe lazy-dep guard.

Pure-stdlib helpers and parsers are tested directly. The numerical CLoNe path is
tested only for its clean failure without numpy/scipy/sklearn (and, if those are
importable, a tiny end-to-end clustering sanity check). Runs standalone.
"""

import csv
import io
import os
import sys
import tempfile
from contextlib import redirect_stderr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs import structural, cli  # noqa: E402
from automd_saxs.clone import CLoNe  # noqa: E402
from automd_saxs.command_runner import MissingDependencyError  # noqa: E402


def _almost(a, b, tol=1e-9):
    assert abs(a - b) <= tol, "expected {0}, got {1}".format(b, a)


# --- percentile / summary ---------------------------------------------------

def test_percentile_matches_numpy_default():
    _almost(structural._percentile([0, 1], 25), 0.25)
    _almost(structural._percentile([0, 1], 50), 0.5)
    _almost(structural._percentile([1, 2, 3, 4], 25), 1.75)
    _almost(structural._percentile([5], 50), 5.0)


def test_cluster_summary_two_clusters():
    data = [[0, 0], [1, 0], [10, 0], [11, 0]]
    labels = [0, 0, 1, 1]
    centers = [0, 2]  # data indices of cluster centres
    summ = structural.cluster_summary(data, labels, centers, headers=["x", "y"])
    assert len(summ) == 2
    c0 = summ[0]
    assert c0["n_elements"] == 2
    _almost(c0["dimensions"][0]["center"], 0.0)     # data[0][0]
    _almost(c0["dimensions"][0]["median"], 0.5)     # median of [0,1]
    _almost(c0["dimensions"][0]["iqr"], 0.5)
    _almost(summ[1]["dimensions"][0]["center"], 10.0)
    _almost(summ[1]["dimensions"][0]["median"], 10.5)


def test_cluster_summary_uses_labels_all_for_count():
    data = [[0], [1], [2]]
    labels = [0, 0, -1]          # one point became an outlier
    labels_all = [0, 0, 0]       # all three originally in cluster 0
    summ = structural.cluster_summary(data, labels, [0], labels_all=labels_all)
    assert summ[0]["n_elements"] == 3
    assert summ[0]["n_after_outliers"] == 2


def test_write_cluster_summary_csv():
    data = [[0, 0], [1, 0], [10, 0], [11, 0]]
    summ = structural.cluster_summary(data, [0, 0, 1, 1], [0, 2], headers=["x", "y"])
    with tempfile.TemporaryDirectory() as tmp:
        path = structural.write_cluster_summary_csv(summ, os.path.join(tmp, "s.csv"))
        rows = list(csv.reader(open(path)))
        assert rows[0] == ["cluster", "center_index", "n_elements", "n_after_outliers",
                           "dimension", "center", "median", "iqr"]
        # 2 clusters x 2 dimensions = 4 data rows
        assert len(rows) == 1 + 4
        assert rows[1][0] == "1"  # cluster index is 1-based in the CSV


# --- parsers ----------------------------------------------------------------

def test_parse_pca_coords():
    text = "PC1(0.43) PC2(0.21)\n0.1 0.2\n0.3 0.4\n"
    headers, rows = structural.parse_pca_coords(text)
    assert headers == ["PC1(0.43)", "PC2(0.21)"]
    assert len(rows) == 2
    _almost(rows[1][0], 0.3)


def test_parse_cluster_results():
    text = "0 2\n0 5 1.20\n0 4 1.10\n1 6 0.90\n"
    res = structural.parse_cluster_results(text)
    assert res["centers"] == [0, 2]
    assert len(res["points"]) == 3
    assert res["points"][0] == (0, 5, 1.20)


# --- CLoNe lazy-dep behaviour ----------------------------------------------

def test_clone_fit_requires_numerics_or_clusters():
    """Without numpy/scipy/sklearn, fit raises a clear error; with them, it clusters."""
    try:
        import numpy as np  # noqa: F401
        import scipy  # noqa: F401
        import sklearn  # noqa: F401
    except ImportError:
        # deps absent (the dev machine): must fail clearly, not with an obscure error
        try:
            CLoNe().fit([[0.0, 0.0], [1.0, 1.0]])
        except MissingDependencyError:
            return
        raise AssertionError("expected MissingDependencyError without numerics")
    # deps present: a trivial two-blob dataset should produce a result with labels
    import numpy as np
    rng = [[0.0, 0.0], [0.1, 0.1], [0.0, 0.1], [10.0, 10.0], [10.1, 10.0], [10.0, 10.1]]
    clone = CLoNe(pdc=4, n_resize=1, filt=0)
    clone.fit(np.array(rng))
    assert hasattr(clone, "labels_") and len(clone.labels_) == len(rng)


def test_cluster_cli_errors_cleanly_without_mdtraj():
    try:
        import mdtraj  # noqa: F401
        import sklearn  # noqa: F401
        return  # deps present: skip the failure-path assertion
    except ImportError:
        pass
    with tempfile.TemporaryDirectory() as tmp:
        err = io.StringIO()
        with redirect_stderr(err):
            code = cli.main(["cluster", "--traj", "x.xtc", "--topo", "x.pdb",
                             "--out-dir", tmp])
        assert code == 1
        assert "requires mdtraj" in err.getvalue()


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

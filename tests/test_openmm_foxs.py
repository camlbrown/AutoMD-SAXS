"""Unit tests for ``automd_saxs.openmm.foxs``. No FoXS needed; runs standalone."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs.openmm import foxs  # noqa: E402

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")


def _almost(a, b, tol=1e-6):
    assert abs(a - b) <= tol, "expected {0}, got {1}".format(b, a)


def test_command_builders():
    assert foxs.foxs_profile_command("m.pdb") == ["foxs", "-p", "m.pdb"]
    assert foxs.foxs_fit_command("exp.dat", "m.pdb") == ["foxs", "exp.dat", "m.pdb"]
    mf = foxs.multifoxs_command("exp.dat", ["a.dat", "b.dat"], output="ens.dat")
    assert mf == ["multi_foxs", "exp.dat", "-o", "ens.dat", "a.dat", "b.dat"]


def test_parse_foxs_fit():
    with open(os.path.join(FIX, "foxs_fit.dat")) as fh:
        parsed = foxs.parse_foxs_fit(fh.read())
    _almost(parsed["chi2"], 1.23)
    _almost(parsed["c1"], 1.01)
    _almost(parsed["c2"], 2.50)
    assert len(parsed["data"]) == 3
    _almost(parsed["data"][0][0], 0.0)        # q
    _almost(parsed["data"][0][1], 1000.0)     # I_exp
    assert len(parsed["data"][0]) == 4


def test_parse_foxs_log():
    log = "foo\nc1 = 0.99 c2 = -1.5\nChi^2 = 3.14\n"
    out = foxs.parse_foxs_log(log)
    _almost(out["c1"], 0.99)
    _almost(out["c2"], -1.5)
    _almost(out["chi2"], 3.14)


def test_frame_index():
    assert foxs.frame_index("md_000000012.pdb") == 12
    assert foxs.frame_index("structure_7.pdb") == 7
    assert foxs.frame_index("model.pdb") is None


def test_best_fit_and_summary():
    records = [
        {"model": "structure_0.pdb", "frame": 0, "chi2": 1.45, "rg": 25.3},
        {"model": "structure_1.pdb", "frame": 1, "chi2": 0.98, "rg": 24.8},
        {"model": "structure_2.pdb", "frame": 2, "chi2": None, "rg": 26.1},  # no fit
    ]
    best = foxs.best_fit(records)
    assert best["frame"] == 1
    summ = foxs.summarize_fits(records)
    _almost(summ["bestChi2"], 0.98)
    assert summ["bestFrame"] == 1
    _almost(summ["rgMean"], (25.3 + 24.8 + 26.1) / 3)


def test_summary_empty():
    assert foxs.summarize_fits([]) == {"bestChi2": None, "bestFrame": None, "rgMean": None}


def _run_standalone():
    funcs = sorted((n, o) for n, o in globals().items()
                   if n.startswith("test_") and callable(o))
    failures = []
    for name, fn in funcs:
        try:
            fn(); print("PASS", name)
        except Exception as exc:  # noqa: BLE001
            failures.append((name, exc)); print("FAIL", name, "->", exc)
    print("\n{0}/{1} passed".format(len(funcs) - len(failures), len(funcs)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_run_standalone())

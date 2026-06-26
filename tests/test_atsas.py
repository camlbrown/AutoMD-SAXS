"""Unit tests for ``automd_saxs.atsas`` parsers and command builders.

No ATSAS/numpy/pytest required; parses captured fixtures. Runs standalone.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs import atsas  # noqa: E402

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")


def _read(name):
    with open(os.path.join(FIX, name)) as fh:
        return fh.read()


def _almost(a, b, tol=1e-6):
    assert abs(a - b) <= tol, "expected {0}, got {1}".format(b, a)


# --- command builders -------------------------------------------------------

def test_crysol_command_flags():
    argv = atsas.crysol_command("exp.dat", smax=0.3)
    assert argv[0] == "crysol"
    assert "-lm" in argv and argv[argv.index("-lm") + 1] == "30"
    assert "-cst" in argv
    assert argv[argv.index("-sm") + 1] == "0.3"


def test_shanum_and_gajoe_commands():
    assert atsas.shanum_command("exp.dat") == ["shanum", "exp.dat"]
    assert atsas.gajoe_command("profiles.txt") == ["gajoe", "-p", "profiles.txt"]


# --- crysol_summary ---------------------------------------------------------

def test_parse_crysol_summary():
    records = atsas.parse_crysol_summary(_read("crysol_summary.txt"))
    assert len(records) == 3
    assert records[0]["model"] == "structure_0.pdb"
    assert records[0]["frame"] == 0
    _almost(records[1]["rg"], 24.80)
    _almost(records[1]["chi2"], 0.980)


def test_best_fit_and_summary():
    records = atsas.parse_crysol_summary(_read("crysol_summary.txt"))
    best = atsas.best_fit(records)
    assert best["frame"] == 1            # lowest chi2 (0.98)
    summary = atsas.summarize_fits(records)
    _almost(summary["bestChi2"], 0.980)
    assert summary["bestFrame"] == 1
    _almost(summary["rgMean"], (25.30 + 24.80 + 26.10) / 3)


def test_summary_empty():
    s = atsas.summarize_fits([])
    assert s == {"bestChi2": None, "bestFrame": None, "rgMean": None}


# --- fit / Rg distribution / xvg --------------------------------------------

def test_parse_fit_skips_header():
    rows = atsas.parse_fit(_read("profiles.fit"))
    assert len(rows) == 3
    _almost(rows[0][0], 0.01)      # q
    _almost(rows[0][1], 1000.0)    # I(q)
    assert len(rows[0]) == 4       # q, I, err, fit


def test_parse_rg_distribution_skips_five():
    rows = atsas.parse_rg_distribution(_read("Rg_distr.txt"))
    assert len(rows) == 3
    _almost(rows[1][0], 24.0)
    _almost(rows[1][1], 0.40)
    _almost(rows[1][2], 0.55)


def test_parse_xvg_metadata_and_rows():
    meta, rows = atsas.parse_xvg(_read("gyrate.xvg"))
    assert meta.get("title") == "Radius of gyration"
    assert len(rows) == 2
    _almost(rows[0][1], 2.345)


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

"""Unit tests for ``automd_saxs.postprocess`` and the ``postprocess`` CLI.

No GROMACS needed: command plans (argv + stdin) are checked directly, and the
CLI is exercised in dry-run. Runs standalone.
"""

import io
import json
import os
import sys
import tempfile
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs import postprocess as pp  # noqa: E402
from automd_saxs import cli  # noqa: E402


# --- per-repeat plan --------------------------------------------------------

def test_repeat_plan_order_and_groups():
    steps = pp.plan_repeat_processing("md.gro", "md.tpr", "md.xtc", "GMX.gro", skip=2)
    labels = [s[0] for s in steps]
    assert labels == [
        "editconf_center", "make_ndx", "editconf_nosolv", "convert_tpr",
        "trjconv_nosolv", "trjconv_nojump", "trjconv_fit", "editconf_final_pdb",
        "extract_frames",
    ]
    by = {s[0]: s for s in steps}
    # make_ndx selects Protein then quits
    assert by["make_ndx"][2] == "1\nq\n"
    # nojump uses System (0); fit uses System twice
    assert by["trjconv_nojump"][2] == "0\n"
    assert by["trjconv_fit"][2] == "0\n0\n"
    # frame extraction is per-frame separated with the configured stride
    assert "-sep" in by["extract_frames"][1]
    assert by["extract_frames"][1][by["extract_frames"][1].index("-skip") + 1] == "2"


def test_repeat_plan_threads_source_files():
    steps = pp.plan_repeat_processing("/r/md.gro", "/r/md.tpr", "/r/md.xtc", "/p/GMX.gro")
    by = {s[0]: s for s in steps}
    assert "/r/md.tpr" in by["convert_tpr"][1]
    assert "/r/md.xtc" in by["trjconv_nosolv"][1]
    assert "/p/GMX.gro" in by["editconf_final_pdb"][1]


# --- combine plan -----------------------------------------------------------

def test_combine_plan_trjcat_and_align():
    steps = pp.plan_combine(["a/final.xtc", "b/final.xtc", "c/final.xtc"], "ref.pdb", skip=2)
    labels = [s[0] for s in steps]
    assert labels == ["trjcat", "align_combined"]
    trjcat = steps[0][1]
    # all three inputs concatenated with -cat
    assert trjcat.count("final.xtc") == 0  # they're separate argv tokens
    assert "a/final.xtc" in trjcat and "c/final.xtc" in trjcat
    assert "-cat" in trjcat
    align = steps[1]
    assert "combined_aligned.xtc" in align[1]
    assert align[2] == "3\n1\n"   # fit on C-alpha (3), output Protein (1)


# --- CLI dry-run ------------------------------------------------------------

def test_postprocess_cli_dry_run_plans_all_steps():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = os.path.join(tmp, "job.json")
        with open(cfg, "w") as fh:
            json.dump({"protein_file": "tiny.pdb", "system": "Protein",
                       "force_field": "amber14sb", "n_repeats": 3,
                       "ionic_concentration_M": 0.15, "simulation_time_ns": 5}, fh)
        buf = io.StringIO()
        with redirect_stdout(buf):
            code = cli.main(["postprocess", "--config", cfg, "--work-dir", tmp, "--dry-run"])
        assert code == 0
        out = buf.getvalue()
        # 3 repeats x 9 steps + 2 combine steps = 29 planned commands
        assert "29 planned command(s)" in out
        assert "combined_aligned.xtc" in out
        # dry-run writes nothing
        assert not os.path.isdir(os.path.join(tmp, "tiny_simulation", "production"))


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

"""Unit tests for ``automd_saxs.cli``. No external deps; runs standalone.

Exercises the CLI in-process (no subprocess), capturing stdout, so it works
without GROMACS/ATSAS/Slurm/pytest.
"""

import io
import json
import os
import sys
import tempfile
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs import cli  # noqa: E402
from automd_saxs import mdp  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG = os.path.join(REPO, "tests", "fixtures", "configurations.txt")
PDB = os.path.join(REPO, "tests", "fixtures", "tiny.pdb")
MD_TEMPLATE = os.path.join(REPO, "mdp_files", "md.mdp")


def _run(argv):
    buf = io.StringIO()
    with redirect_stdout(buf):
        code = cli.main(argv)
    return code, buf.getvalue()


def test_plan_exits_zero_and_reports_genion_fix():
    code, text = _run([
        "plan", "--config", CONFIG, "--pdb", PDB,
        "--work-dir", "/tmp/clitest",
        "--slurm-dir", "slurms", "--mdp-dir", "mdp_files",
    ])
    assert code == 0
    assert "genion -conc       : 0.2 M" in text       # not the legacy 0.15
    assert "box padding (-d)  : 3.500 nm" in text
    assert "afterok:${JOBID_setup}" in text           # DAG rendered


def test_plan_writes_manifest():
    with tempfile.TemporaryDirectory() as tmp:
        manifest_path = os.path.join(tmp, "m.json")
        code, _ = _run([
            "plan", "--config", CONFIG, "--pdb", PDB,
            "--work-dir", tmp, "--manifest", manifest_path,
        ])
        assert code == 0
        data = json.load(open(manifest_path))
        assert data["parameters"]["genionConcentrationM"] == 0.2
        assert data["parameters"]["boxPaddingNm"] == 3.5


def test_dmax_subcommand_output():
    code, text = _run(["dmax", "--pdb", PDB, "--experimental", "6.5"])
    assert code == 0
    assert "MODEL_DMAX=5.000" in text
    assert "BOX_PADDING=3.500" in text


def test_render_mdp_subcommand_leaves_template_untouched():
    if not os.path.isfile(MD_TEMPLATE):
        print("skip: md.mdp not found")
        return
    before = open(MD_TEMPLATE).read()
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "md.mdp")
        code, _ = _run(["render-mdp", "--template", MD_TEMPLATE, "--out", out,
                        "--config", CONFIG])
        assert code == 0
        assert mdp.parse_mdp(open(out).read())["nsteps"] == str(50 * 500000)
    assert open(MD_TEMPLATE).read() == before  # template byte-for-byte unchanged


def test_no_command_returns_2():
    code, _ = _run([])
    assert code == 2


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

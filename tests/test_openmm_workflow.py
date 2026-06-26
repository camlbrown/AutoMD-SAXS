"""Unit tests for ``automd_saxs.openmm.workflow`` and the OpenMM CLI.

No OpenMM/FoXS needed: stage planning, manifest shaping, and the CLI dry-run are
pure. The executing Workflow.run is checked for its clean not-yet-implemented /
missing-dependency behaviour. Runs standalone.
"""

import io
import json
import os
import sys
import tempfile
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs.openmm import cli, workflow  # noqa: E402
from automd_saxs.openmm.schema import OpenMMConfig  # noqa: E402


def test_plan_stages_protein_no_saxs():
    cfg = OpenMMConfig(pdb="m.pdb", n_repeats=2)
    names = [n for n, _ in workflow.plan_stages(cfg)]
    assert names[:5] == ["validate_inputs", "prepare_structure", "solvate",
                         "minimize", "equilibrate"]
    assert "production_rep1" in names and "production_rep2" in names
    assert "production_rep3" not in names
    # no SAXS -> no foxs/multifoxs, but clustering still runs
    assert "foxs" not in names and "multifoxs" not in names
    assert names[-1] == "cluster"


def test_plan_stages_with_saxs_adds_foxs():
    cfg = OpenMMConfig(pdb="m.pdb", saxs="d.dat", n_repeats=1)
    names = [n for n, _ in workflow.plan_stages(cfg)]
    assert "foxs" in names and "multifoxs" in names


def test_new_manifest_pipeline_and_inputs():
    cfg = OpenMMConfig(pdb="m.pdb", saxs="d.dat", simulation_time_ns=20)
    m = workflow.new_manifest(cfg).to_dict()
    assert m["pipeline"] == "automd-saxs-openmm"
    assert m["inputs"] == {"pdb": "m.pdb", "saxs": "d.dat"}
    assert m["parameters"]["productionSteps"] == 20 * 500000


def test_workflow_run_dry_run_returns_planned():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        cfg = OpenMMConfig(pdb="m.pdb", job_name="wfjob")
        result = workflow.Workflow(cfg, tmp, dry_run=True).run()
        assert result["status"] == "planned"
        assert result["job_dir"].endswith("wfjob")


def test_cli_plan_end_to_end():
    with tempfile.TemporaryDirectory() as tmp:
        cfg_path = os.path.join(tmp, "job.json")
        with open(cfg_path, "w") as fh:
            json.dump({"pdb": "m.pdb", "saxs": "d.dat", "simulation_time_ns": 10,
                       "force_field": "amber14", "water_model": "tip3p", "n_repeats": 2}, fh)
        manifest = os.path.join(tmp, "m.json")
        buf = io.StringIO()
        with redirect_stdout(buf):
            code = cli.main(["plan", "--config", cfg_path, "--manifest", manifest])
        assert code == 0
        out = buf.getvalue()
        assert "Pipeline stages" in out and "production_rep2" in out
        data = json.load(open(manifest))
        assert data["pipeline"] == "automd-saxs-openmm"
        assert any(s["name"] == "foxs" for s in data["steps"])


def test_cli_validate_rejects_bad_config():
    with tempfile.TemporaryDirectory() as tmp:
        cfg_path = os.path.join(tmp, "bad.json")
        with open(cfg_path, "w") as fh:
            json.dump({"pdb": "m.gro"}, fh)  # not a .pdb
        code = cli.main(["validate", "--config", cfg_path])
        assert code == 1


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

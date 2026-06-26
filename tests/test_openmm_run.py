"""Tests for OpenMM paths, frame extraction guard, and Workflow.run wiring.

No OpenMM/mdtraj/FoXS needed: the dry-run path is exercised fully; the real path
is checked for its clean missing-dependency failure (which is what happens until
the BilboMD image provides the binaries). Runs standalone.
"""

import io
import json
import os
import sys
import tempfile
from contextlib import redirect_stderr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs.command_runner import CommandRunner, MissingDependencyError  # noqa: E402
from automd_saxs.openmm import cli, foxs, frames, workflow  # noqa: E402
from automd_saxs.openmm.paths import OpenMMPaths  # noqa: E402
from automd_saxs.openmm.schema import OpenMMConfig  # noqa: E402


def _cfg(tmp, **kw):
    base = dict(pdb="m.pdb", n_repeats=2)
    base.update(kw)
    return OpenMMConfig(**base)


# --- paths ------------------------------------------------------------------

def test_paths_layout_and_create():
    with tempfile.TemporaryDirectory() as tmp:
        p = OpenMMPaths(tmp, "job1", n_repeats=2)
        assert p.job_dir.endswith("job1")
        assert p.repeat_trajectory(1).endswith("production/rep1/production.dcd")
        made = p.create()
        for d in made:
            assert os.path.isdir(d)


def test_paths_explicit_job_dir_no_nesting():
    # The --out contract: results land directly in the given dir (no job_name subdir)
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "results")
        p = OpenMMPaths(job_dir=out, n_repeats=1)
        assert p.job_dir == os.path.abspath(out)
        assert p.manifest_path == os.path.join(os.path.abspath(out), "manifest.json")


def test_run_out_dir_writes_manifest_directly():
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "jobout")
        cfg = OpenMMConfig(pdb="m.pdb", job_name="ignored_for_out", n_repeats=1,
                           simulation_time_ns=5)
        result = workflow.Workflow(cfg, out_dir=out, dry_run=True).run()
        assert result["status"] == "planned"
        # manifest is directly in --out, not nested under job_name
        assert os.path.isfile(os.path.join(out, "manifest.json"))
        assert result["job_dir"] == os.path.abspath(out)


# --- frames lazy guard ------------------------------------------------------

def test_extract_frames_requires_mdtraj():
    try:
        import mdtraj  # noqa: F401
        return  # present: skip
    except ImportError:
        pass
    try:
        frames.extract_frames("t.dcd", "top.pdb", "/tmp/out", 2)
    except MissingDependencyError:
        return
    raise AssertionError("expected MissingDependencyError")


# --- FoXS orchestration in dry-run ------------------------------------------

def test_run_foxs_fits_dry_run_records_commands():
    runner = CommandRunner(dry_run=True)
    recs = foxs.run_foxs_fits(["frames/structure_0.pdb", "frames/structure_1.pdb"],
                              "exp.dat", runner, cwd="saxs")
    assert recs == []  # nothing parsed in dry-run
    assert [p.label for p in runner.planned] == ["foxs:structure_0.pdb", "foxs:structure_1.pdb"]
    assert runner.planned[0].argv == ["foxs", "exp.dat", "frames/structure_0.pdb"]


# --- Workflow.run dry-run ---------------------------------------------------

def test_workflow_run_dry_run_prepares_and_plans():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = OpenMMConfig(pdb="m.pdb", saxs="d.dat", job_name="jobx", n_repeats=2,
                           simulation_time_ns=5)
        result = workflow.Workflow(cfg, tmp, dry_run=True).run()
        assert result["status"] == "planned"
        job_dir = os.path.join(tmp, "jobx")
        assert os.path.isfile(os.path.join(job_dir, "job.json"))
        assert os.path.isdir(os.path.join(job_dir, "production", "rep2"))
        data = json.load(open(os.path.join(job_dir, "manifest.json")))
        assert data["pipeline"] == "automd-saxs-openmm"
        names = [s["name"] for s in data["steps"]]
        assert "production_rep2" in names and "foxs" in names and "multifoxs" in names
        # the recorded multifoxs step carries a command
        mf = [s for s in data["steps"] if s["name"] == "multifoxs" and "command" in s]
        assert mf and mf[0]["command"][0] == "multi_foxs"


# --- Workflow.run real path fails cleanly without OpenMM --------------------

def test_workflow_run_real_fails_cleanly_without_openmm():
    try:
        import openmm  # noqa: F401
        return  # present: skip (would actually run)
    except ImportError:
        pass
    with tempfile.TemporaryDirectory() as tmp:
        cfg = OpenMMConfig(pdb="m.pdb", job_name="jobr", n_repeats=1)
        wf = workflow.Workflow(cfg, tmp, dry_run=False)
        try:
            wf.run()
        except MissingDependencyError:
            data = json.load(open(os.path.join(tmp, "jobr", "manifest.json")))
            assert data["status"] == "failed"
            assert data["notes"]
            return
        raise AssertionError("expected MissingDependencyError")


def test_cli_run_dry_run_and_real():
    with tempfile.TemporaryDirectory() as tmp:
        cfg_path = os.path.join(tmp, "job.json")
        with open(cfg_path, "w") as fh:
            json.dump({"pdb": "m.pdb", "job_name": "clijob", "n_repeats": 1,
                       "simulation_time_ns": 5}, fh)
        assert cli.main(["run", "--config", cfg_path, "--work-dir", tmp, "--dry-run"]) == 0
        # real run without openmm -> exit 1 (skip if openmm present)
        try:
            import openmm  # noqa: F401
            return
        except ImportError:
            pass
        err = io.StringIO()
        with redirect_stderr(err):
            code = cli.main(["run", "--config", cfg_path, "--work-dir", tmp])
        assert code == 1
        assert "OpenMM" in err.getvalue()


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

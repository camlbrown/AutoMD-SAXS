"""Unit tests for the ``run`` subcommand. No external binaries; runs standalone.

The dry-run path is exercised fully. The real path is exercised only up to the
first missing binary, which must fail cleanly with a failed manifest -- exactly
the situation on a machine without GROMACS.
"""

import io
import json
import os
import sys
import tempfile
from contextlib import redirect_stdout, redirect_stderr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs import cli  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PDB = os.path.join(REPO, "tests", "fixtures", "tiny.pdb")
MDP_DIR = os.path.join(REPO, "mdp_files")
SLURM_DIR = os.path.join(REPO, "slurms")


def _write_cfg(path, system="Protein"):
    cfg = {
        "protein_file": "tiny.pdb", "system": system, "force_field": "amber14sb",
        "box_shape": "dodecahedron", "ionic_concentration_M": 0.2,
        "simulation_time_ns": 10, "n_repeats": 2,
    }
    with open(path, "w") as fh:
        json.dump(cfg, fh)
    return path


def _run(argv):
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        code = cli.main(argv)
    return code, out.getvalue(), err.getvalue()


def test_run_dry_run_prepares_job():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _write_cfg(os.path.join(tmp, "job.json"))
        code, _, _ = _run(["run", "--config", cfg, "--pdb", PDB, "--work-dir", tmp,
                           "--mdp-dir", MDP_DIR, "--slurm-dir", SLURM_DIR, "--dry-run"])
        assert code == 0
        sim = os.path.join(tmp, "tiny_simulation")
        # directories created
        assert os.path.isdir(os.path.join(sim, "pdb2gmx"))
        # job.json written
        assert os.path.isfile(os.path.join(sim, "job.json"))
        # production mdp rendered into per-job mdp dir with overridden nsteps
        rendered = os.path.join(sim, "mdp_files", "md.mdp")
        assert os.path.isfile(rendered)
        from automd_saxs import mdp
        assert mdp.parse_mdp(open(rendered).read())["nsteps"] == str(10 * 500000)
        # manifest written, status planned, genion step present
        data = json.load(open(os.path.join(sim, "manifest.json")))
        assert data["status"] == "planned"
        assert any(s["name"] == "genion" for s in data["steps"])
        assert any(s["name"].startswith("sbatch:") for s in data["steps"])


def test_run_real_fails_cleanly_without_binaries():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _write_cfg(os.path.join(tmp, "job.json"))
        code, _, err = _run(["run", "--config", cfg, "--pdb", PDB, "--work-dir", tmp,
                             "--mdp-dir", MDP_DIR, "--slurm-dir", SLURM_DIR])
        assert code == 1
        assert "not found on PATH" in err
        # a failed manifest is still written for traceability
        data = json.load(open(os.path.join(tmp, "tiny_simulation", "manifest.json")))
        assert data["status"] == "failed"
        assert data["notes"]


def test_run_requires_pdb_when_executing():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _write_cfg(os.path.join(tmp, "job.json"))
        code, _, err = _run(["run", "--config", cfg, "--work-dir", tmp,
                             "--mdp-dir", MDP_DIR, "--slurm-dir", SLURM_DIR])
        assert code == 2
        assert "--pdb is required" in err


def test_run_protein_ligand_execution_blocked():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _write_cfg(os.path.join(tmp, "job.json"), system="Protein-ligand")
        code, _, err = _run(["run", "--config", cfg, "--pdb", PDB, "--work-dir", tmp,
                             "--mdp-dir", MDP_DIR, "--slurm-dir", SLURM_DIR])
        assert code == 2
        assert "protein-ligand execution is not yet ported" in err


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

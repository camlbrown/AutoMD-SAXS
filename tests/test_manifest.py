"""Unit tests for ``automd_saxs.manifest``. No external deps; runs standalone."""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs import manifest as mf  # noqa: E402
from automd_saxs.config import ForceField, JobConfig, SystemType  # noqa: E402


def _cfg(**kw):
    base = dict(protein_file="1AKI.pdb", saxs_file="1AKI.dat",
                ionic_concentration_M=0.2, simulation_time_ns=50)
    base.update(kw)
    return JobConfig(**base)


def test_manifest_core_schema_keys():
    m = mf.Manifest(_cfg())
    data = m.to_dict()
    for key in ("pipeline", "status", "inputs", "outputs", "metrics"):
        assert key in data
    assert data["pipeline"] == "automd-saxs-gromacs"
    assert data["status"] == mf.STATUS_PLANNED
    assert data["inputs"] == {"pdb": "1AKI.pdb", "saxs": "1AKI.dat"}
    # all output categories present and empty
    assert set(data["outputs"]) == {
        "trajectories", "structures", "saxsFits", "summaryTables", "plots", "logs"
    }
    assert all(v == [] for v in data["outputs"].values())
    assert data["metrics"] == {"bestChi2": None, "bestFrame": None, "rgMean": None}


def test_parameters_derived_from_config():
    m = mf.Manifest(_cfg(simulation_time_ns=50))
    params = m.to_dict()["parameters"]
    assert params["productionNsteps"] == 50 * 500000
    assert params["waterModel"] == "tip3p"
    assert params["ionicConcentrationM"] == 0.2


def test_add_output_validates_category():
    m = mf.Manifest(_cfg())
    m.add_output("structures", "best.pdb")
    assert m.to_dict()["outputs"]["structures"] == ["best.pdb"]
    try:
        m.add_output("nonsense", "x")
    except ValueError:
        return
    raise AssertionError("expected ValueError for bad category")


def test_set_metric_validates_key():
    m = mf.Manifest(_cfg())
    m.set_metric("bestChi2", 1.23)
    assert m.to_dict()["metrics"]["bestChi2"] == 1.23
    try:
        m.set_metric("notAMetric", 1)
    except ValueError:
        return
    raise AssertionError("expected ValueError for bad metric")


def test_json_is_parseable_and_includes_config():
    m = mf.Manifest(_cfg())
    data = json.loads(m.to_json())
    assert data["config"]["protein_file"] == "1AKI.pdb"


def test_build_plan_manifest_captures_box_padding_and_steps():
    cfg = _cfg()
    steps = [("pdb2gmx", ["gmx_mpi", "pdb2gmx"]), ("genion", ["gmx_mpi", "genion"])]
    m = mf.build_plan_manifest(cfg, box_padding_nm=3.5, steps=steps)
    data = m.to_dict()
    assert data["parameters"]["boxPaddingNm"] == 3.5
    assert data["parameters"]["genionConcentrationM"] == 0.2
    assert [s["name"] for s in data["steps"]] == ["pdb2gmx", "genion"]
    assert data["steps"][0]["command"] == ["gmx_mpi", "pdb2gmx"]


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

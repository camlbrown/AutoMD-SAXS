"""Unit tests for ``automd_saxs.openmm.schema``. No OpenMM needed; runs standalone."""

import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs.config import SystemType  # noqa: E402
from automd_saxs.openmm.schema import (  # noqa: E402
    OpenMMConfig,
    OpenMMForceField,
    WaterModel,
)


def _expect_error(fn, *a, **k):
    try:
        fn(*a, **k)
    except (ValueError, TypeError):
        return
    raise AssertionError("expected error")


def test_minimal_valid():
    c = OpenMMConfig(pdb="m.pdb")
    assert c.system is SystemType.PROTEIN
    assert c.force_field is OpenMMForceField.AMBER14
    assert c.uses_saxs is False


def test_forcefield_files():
    c = OpenMMConfig(pdb="m.pdb", force_field=OpenMMForceField.AMBER14,
                     water_model=WaterModel.TIP3P)
    files = c.forcefield_files()
    assert files[0] == "amber14-all.xml"
    assert files[1] == "amber14/tip3p.xml"


def test_production_and_equilibration_steps():
    c = OpenMMConfig(pdb="m.pdb", simulation_time_ns=50, timestep_fs=2.0, equilibration_ns=0.2)
    assert c.production_steps() == 50 * 500000
    assert c.equilibration_steps() == 100000


def test_fractional_ns_allowed():
    # sub-ns / fractional ns is valid as long as it yields whole steps
    c = OpenMMConfig(pdb="m.pdb", simulation_time_ns=0.02, timestep_fs=2.0)
    assert c.production_steps() == 10000
    _expect_error(OpenMMConfig, pdb="m.pdb", simulation_time_ns=0)


def test_rejects_bad_inputs():
    _expect_error(OpenMMConfig, pdb="m.gro")
    _expect_error(OpenMMConfig, pdb="m.pdb", saxs="d.txt")
    _expect_error(OpenMMConfig, pdb="m.pdb", ph=20)
    _expect_error(OpenMMConfig, pdb="m.pdb", temperature_K=0)
    _expect_error(OpenMMConfig, pdb="m.pdb", n_repeats=0)


def test_unsupported_ff_water_combo():
    # charmm36 + tip3pfb is not in the supported map
    _expect_error(OpenMMConfig, pdb="m.pdb", force_field=OpenMMForceField.CHARMM36,
                  water_model=WaterModel.TIP3PFB)


def test_json_round_trip():
    c = OpenMMConfig(pdb="m.pdb", saxs="d.dat", system=SystemType.PROTEIN_LIGAND,
                     simulation_time_ns=20, seed=42, ionic_concentration_M=0.2)
    r = OpenMMConfig.from_json(c.to_json())
    assert r.to_dict() == c.to_dict()
    assert r.system is SystemType.PROTEIN_LIGAND
    assert r.seed == 42


def test_from_dict_unknown_to_extra():
    c = OpenMMConfig.from_dict({"pdb": "m.pdb", "future_knob": 7})
    assert c.extra["future_knob"] == 7


def test_from_file_json():
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump({"pdb": "m.pdb", "simulation_time_ns": 10, "force_field": "amber14",
                   "water_model": "tip3p"}, fh)
        path = fh.name
    try:
        c = OpenMMConfig.from_file(path)
        assert c.simulation_time_ns == 10
    finally:
        os.unlink(path)


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

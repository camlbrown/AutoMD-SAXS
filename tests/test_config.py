"""Unit tests for ``automd_saxs.config`` and ``automd_saxs.legacy``.

No GROMACS/ATSAS/Slurm/numpy/pytest required. Runs under pytest and standalone
(``python tests/test_config.py``).
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs import legacy  # noqa: E402
from automd_saxs.config import (  # noqa: E402
    BoxShape,
    ForceField,
    JobConfig,
    SystemType,
)

FIXTURE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "fixtures", "configurations.txt"
)


def _almost(a, b, tol=1e-9):
    assert abs(a - b) <= tol, "expected {0}, got {1}".format(b, a)


def _expect_value_error(fn, *args, **kwargs):
    try:
        fn(*args, **kwargs)
    except ValueError:
        return
    raise AssertionError("expected ValueError")


# --- enums ------------------------------------------------------------------

def test_system_from_str_aliases():
    assert SystemType.from_str("Protein") is SystemType.PROTEIN
    assert SystemType.from_str("protein-ligand") is SystemType.PROTEIN_LIGAND
    assert SystemType.from_str("Protein_ligand") is SystemType.PROTEIN_LIGAND
    _expect_value_error(SystemType.from_str, "membrane")


def test_force_field_from_path_basename():
    assert ForceField.from_str("amber14sb") is ForceField.AMBER14SB
    # legacy FORCE_FIELD was an absolute path ending in the name
    assert ForceField.from_str("/x/ff_files/charmm36m") is ForceField.CHARMM36M
    _expect_value_error(ForceField.from_str, "opls")


def test_box_shape_rectangular_alias_is_triclinic():
    assert BoxShape.from_str("rectangular") is BoxShape.TRICLINIC
    assert BoxShape.from_str("Dodecahedron") is BoxShape.DODECAHEDRON
    assert BoxShape.from_str("auto") is BoxShape.AUTO


# --- JobConfig validation ---------------------------------------------------

def test_minimal_valid_config():
    cfg = JobConfig(protein_file="prot.pdb")
    assert cfg.system is SystemType.PROTEIN
    assert cfg.uses_saxs is False
    assert cfg.water_model() == "tip3p"


def test_rejects_non_pdb_protein():
    _expect_value_error(JobConfig, protein_file="prot.gro")


def test_rejects_non_dat_saxs():
    _expect_value_error(JobConfig, protein_file="prot.pdb", saxs_file="curve.txt")


def test_rejects_nonpositive_ionic_concentration():
    _expect_value_error(JobConfig, protein_file="prot.pdb", ionic_concentration_M=0)


def test_protein_ligand_requires_amber():
    _expect_value_error(
        JobConfig,
        protein_file="prot.pdb",
        system=SystemType.PROTEIN_LIGAND,
        force_field=ForceField.CHARMM36M,
    )
    # amber is fine, and switches the water model to spce
    cfg = JobConfig(
        protein_file="prot.pdb",
        system=SystemType.PROTEIN_LIGAND,
        force_field=ForceField.AMBER14SB,
    )
    assert cfg.water_model() == "spce"


# --- derived quantities -----------------------------------------------------

def test_number_of_steps_matches_legacy_for_2fs():
    # legacy: SIMULATION_TIME * 500000 at 2 fs
    cfg = JobConfig(protein_file="p.pdb", simulation_time_ns=50, timestep_fs=2.0)
    assert cfg.number_of_steps() == 50 * 500000


def test_number_of_steps_other_timestep():
    cfg = JobConfig(protein_file="p.pdb", simulation_time_ns=1, timestep_fs=4.0)
    assert cfg.number_of_steps() == 250000


def test_genion_concentration_uses_chosen_value():
    # The legacy run_MD.sh hardcoded 0.15; the model must honour the choice.
    cfg = JobConfig(protein_file="p.pdb", ionic_concentration_M=0.2)
    _almost(cfg.genion_concentration(), 0.2)


def test_partition_flag():
    assert JobConfig(protein_file="p.pdb", partition="batch").partition_flag() == "--partition=batch"
    assert JobConfig(protein_file="p.pdb").partition_flag() == ""


# --- serialisation round-trip ----------------------------------------------

def test_json_round_trip():
    cfg = JobConfig(
        protein_file="p.pdb",
        system=SystemType.PROTEIN_LIGAND,
        force_field=ForceField.AMBER14SB,
        saxs_file="d.dat",
        box_shape=BoxShape.TRICLINIC,
        ionic_concentration_M=0.15,
        simulation_time_ns=20,
        dmax_nm=7.0,
        partition="gpu",
    )
    restored = JobConfig.from_json(cfg.to_json())
    assert restored.to_dict() == cfg.to_dict()
    assert restored.system is SystemType.PROTEIN_LIGAND
    assert restored.box_shape is BoxShape.TRICLINIC


def test_from_dict_keeps_unknown_keys_in_extra():
    cfg = JobConfig.from_dict({"protein_file": "p.pdb", "future_field": 123})
    assert cfg.extra["future_field"] == 123


# --- legacy configurations.txt parsing -------------------------------------

def test_parse_configurations_skips_functions_and_comments():
    values = legacy.parse_configurations(FIXTURE)
    # scientific keys captured
    assert values["SYSTEM"] == "Protein"
    assert values["GMX_MODULE"] == "gromacs/2021.2/intel"
    assert values["PARTITION"] == "--partition=batch"
    # function definitions must NOT appear as keys
    assert "load_gmx" not in values
    assert "unload_gmx" not in values


def test_config_from_legacy_file_maps_all_fields():
    cfg = legacy.config_from_legacy_file(FIXTURE)
    assert cfg.protein_file == "1AKI.pdb"
    assert cfg.system is SystemType.PROTEIN
    assert cfg.force_field is ForceField.AMBER14SB
    assert cfg.saxs_file == "1AKI.dat"
    assert cfg.box_shape is BoxShape.DODECAHEDRON
    _almost(cfg.ionic_concentration_M, 0.2)
    assert cfg.simulation_time_ns == 50
    assert cfg.disulfide is True
    _almost(cfg.dmax_nm, 6.5)
    assert cfg.gmx_module == "gromacs/2021.2/intel"
    assert cfg.email == "scientist@example.org"
    assert cfg.partition == "batch"
    # the genion bug is fixed: planning would use 0.2, not 0.15
    _almost(cfg.genion_concentration(), 0.2)


def test_legacy_none_tokens_become_none():
    text = (
        "PROTEIN_FILE=x.pdb\n"
        "SYSTEM=Protein\n"
        "FORCE_FIELD=/a/b/amber14sb\n"
        "BOX_SHAPE=dodecahedron\n"
        "IONIC_CONCENTRATION=0.15\n"
        "SIMULATION_TIME=10\n"
        "DISULFIDE=no\n"
        "SAXS_FILE=None\n"
        "DMAX=Model\n"
        "PARTITION=\n"
        "EMAIL_ADDR=\n"
        "GMX_MODULE=\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write(text)
        path = fh.name
    try:
        cfg = legacy.config_from_legacy_file(path)
        assert cfg.saxs_file is None
        assert cfg.uses_saxs is False
        assert cfg.dmax_nm is None           # "Model" sentinel -> None
        assert cfg.partition is None
        assert cfg.email is None
        assert cfg.gmx_module is None
        assert cfg.disulfide is False
    finally:
        os.unlink(path)


def _run_standalone():
    funcs = sorted(
        (name, obj)
        for name, obj in globals().items()
        if name.startswith("test_") and callable(obj)
    )
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

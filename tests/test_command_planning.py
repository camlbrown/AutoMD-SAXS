"""Unit tests for command planning: gromacs_commands, slurm, command_runner.

No GROMACS/ATSAS/Slurm/numpy/pytest needed. Runs standalone.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs import gromacs_commands as gmx  # noqa: E402
from automd_saxs import slurm as slurm_mod  # noqa: E402
from automd_saxs.command_runner import (  # noqa: E402
    CommandRunner,
    MissingDependencyError,
)
from automd_saxs.config import BoxShape, ForceField, JobConfig, SystemType  # noqa: E402


def _protein_cfg(**kw):
    base = dict(protein_file="p.pdb", ionic_concentration_M=0.2, disulfide=True)
    base.update(kw)
    return JobConfig(**base)


# --- gromacs_commands -------------------------------------------------------

def test_pdb2gmx_includes_ff_water_and_ss():
    cfg = _protein_cfg(force_field=ForceField.CHARMM36M)
    argv = gmx.pdb2gmx_command(cfg)
    assert argv[:2] == ["gmx_mpi", "pdb2gmx"]
    assert "-ff" in argv and argv[argv.index("-ff") + 1] == "charmm36m"
    assert "-water" in argv and argv[argv.index("-water") + 1] == "tip3p"
    assert "-ss" in argv  # disulfide=True


def test_pdb2gmx_omits_ss_when_no_disulfide():
    argv = gmx.pdb2gmx_command(_protein_cfg(disulfide=False))
    assert "-ss" not in argv


def test_editconf_box_maps_shape_and_padding():
    cfg = _protein_cfg(box_shape=BoxShape.TRICLINIC)
    argv = gmx.editconf_box_command(cfg, 2.5)
    assert argv[argv.index("-bt") + 1] == "triclinic"
    assert argv[argv.index("-d") + 1] == "2.5"


def test_editconf_box_rejects_auto():
    cfg = _protein_cfg(box_shape=BoxShape.AUTO)
    try:
        gmx.editconf_box_command(cfg, 2.0)
    except ValueError:
        return
    raise AssertionError("expected ValueError for AUTO box shape")


def test_genion_uses_config_concentration_not_hardcoded():
    cfg = _protein_cfg(ionic_concentration_M=0.2)
    argv = gmx.genion_command(cfg)
    assert argv[argv.index("-conc") + 1] == "0.2"   # legacy hardcoded 0.15
    assert "-neutral" in argv


def test_setup_pipeline_order_and_no_wildcards():
    cfg = _protein_cfg()
    steps = gmx.setup_pipeline(cfg, 2.0, "/mdp/ions.mdp")
    labels = [label for label, _ in steps]
    assert labels == ["pdb2gmx", "editconf_box", "editconf_center",
                      "solvate", "ions_grompp", "genion"]
    # no argument should be a bare wildcard
    for _, argv in steps:
        for token in argv:
            assert "*" not in token


# --- slurm ------------------------------------------------------------------

def test_pipeline_dag_protein():
    cfg = _protein_cfg(n_repeats=3)
    jobs = {j.name: j for j in slurm_mod.build_pipeline(cfg)}
    assert jobs["minim1"].depends_on == []
    assert jobs["minim2"].depends_on == ["minim1"]
    assert jobs["nvt"].depends_on == ["minim2"]
    assert jobs["npt"].depends_on == ["nvt"]
    assert jobs["setup"].depends_on == ["npt"]
    assert jobs["rep1"].depends_on == ["setup"]
    assert jobs["rep3"].depends_on == ["setup"]
    # protein path uses the plain script names
    assert jobs["nvt"].script == "nvt/nvt.slurm"
    assert jobs["rep2"].script == "production/rep2/rep2.slurm"


def test_pipeline_dag_protein_ligand_uses_lig_variants():
    cfg = JobConfig(protein_file="p.pdb", system=SystemType.PROTEIN_LIGAND,
                    force_field=ForceField.AMBER14SB, n_repeats=2)
    jobs = {j.name: j for j in slurm_mod.build_pipeline(cfg)}
    assert jobs["nvt"].script == "nvt/lig_nvt.slurm"
    assert jobs["npt"].script == "npt/lig_npt.slurm"
    assert jobs["setup"].script == "production/md_setup/lig_md_setup.slurm"
    assert jobs["rep1"].script == "production/rep1/lig_rep1.slurm"
    assert "rep3" not in jobs  # n_repeats=2


def test_postprocess_node_depends_on_all_repeats():
    cfg = _protein_cfg(n_repeats=3)
    jobs = {j.name: j for j in slurm_mod.build_pipeline(cfg)}
    assert "postprocess" in jobs
    assert jobs["postprocess"].script == "analysis/postprocess.slurm"
    assert jobs["postprocess"].depends_on == ["rep1", "rep2", "rep3"]


def test_analyze_and_cluster_depend_on_postprocess():
    cfg = _protein_cfg(saxs_file="d.dat", n_repeats=3)
    jobs = {j.name: j for j in slurm_mod.build_pipeline(cfg)}
    assert jobs["analyze"].depends_on == ["postprocess"]
    assert jobs["cluster"].depends_on == ["postprocess"]


def test_no_analyze_node_without_saxs():
    cfg = _protein_cfg()  # no saxs_file
    jobs = {j.name: j for j in slurm_mod.build_pipeline(cfg)}
    assert "analyze" not in jobs
    # cluster (and therefore postprocess) still present
    assert "cluster" in jobs and "postprocess" in jobs


def test_cluster_node_present_even_without_saxs():
    cfg = _protein_cfg(n_repeats=3)  # no saxs_file
    jobs = {j.name: j for j in slurm_mod.build_pipeline(cfg)}
    assert "cluster" in jobs
    assert jobs["cluster"].script == "analysis/cluster.slurm"
    assert jobs["cluster"].depends_on == ["postprocess"]


def test_no_postprocess_when_all_analysis_disabled():
    cfg = _protein_cfg()
    jobs = {j.name: j for j in
            slurm_mod.build_pipeline(cfg, include_cluster=False, include_analyze=False)}
    assert "cluster" not in jobs and "postprocess" not in jobs


def test_postprocess_sbatch_has_multi_repeat_dependency():
    cfg = _protein_cfg(saxs_file="d.dat", n_repeats=2)
    subs = dict((j.name, argv) for j, argv in
                slurm_mod.plan_submissions(cfg, "/slurms", "/sim"))
    assert any(a == "--dependency=afterok:${JOBID_rep1}:${JOBID_rep2}"
               for a in subs["postprocess"])
    assert any(a == "--dependency=afterok:${JOBID_postprocess}" for a in subs["analyze"])


def test_sbatch_command_builds_dependency_flag():
    cfg = _protein_cfg()
    subs = dict((j.name, argv) for j, argv in
                slurm_mod.plan_submissions(cfg, "/slurms", "/sim", "--partition=batch"))
    min1 = subs["minim1"]
    assert "sbatch" == min1[0] and "--parsable" in min1
    assert "--partition=batch" in min1
    assert "/sim" == min1[-1]
    # minim2 depends on minim1 via symbolic token
    assert any(a == "--dependency=afterok:${JOBID_minim1}" for a in subs["minim2"])
    # rep1 depends on setup
    assert any(a == "--dependency=afterok:${JOBID_setup}" for a in subs["rep1"])


# --- command_runner ---------------------------------------------------------

def test_runner_dry_run_records_without_executing():
    r = CommandRunner(dry_run=True)
    planned = r.run(["gmx_mpi", "pdb2gmx"], label="pdb2gmx", stdin="SOL")
    assert planned.label == "pdb2gmx"
    assert r.planned_labels() == ["pdb2gmx"]
    assert "SOL" in planned.as_shell()


def test_runner_raises_for_missing_binary_when_not_dry():
    r = CommandRunner(dry_run=False)
    try:
        r.run(["definitely-not-a-real-binary-xyz", "--help"])
    except MissingDependencyError:
        return
    raise AssertionError("expected MissingDependencyError")


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

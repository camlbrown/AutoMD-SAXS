"""Slurm job-pipeline planning (no submission).

The legacy ``run_MD.sh`` submitted a fixed dependency chain with ``sbatch
--dependency=afterok:<id>``::

    minim1 -> minim2 -> nvt -> npt -> md_setup -> {rep1, rep2, ... repN}

with separate script variants for the protein and protein-ligand paths. This
module reproduces that DAG as inspectable data and builds the ``sbatch`` argv for
each job, without ever calling ``sbatch``. Real job IDs are only known at submit
time, so for planning each dependency is rendered as a symbolic token
(``${JOBID_minim2}``) that mirrors the shell variables the legacy script used.
"""

from typing import Dict, List, Optional

from .config import JobConfig, SystemType


class SlurmJob:
    """One node in the submission DAG."""

    def __init__(self, name, script, depends_on=None):
        self.name = name
        self.script = script                      # path relative to slurms/
        self.depends_on = list(depends_on or [])  # job names this waits for

    def __repr__(self):
        return "SlurmJob(name={0!r}, script={1!r}, depends_on={2})".format(
            self.name, self.script, self.depends_on
        )


def _variant(config: JobConfig, base: str) -> str:
    """Prefix a script basename with ``lig_`` on the protein-ligand path."""
    if config.system is SystemType.PROTEIN_LIGAND:
        return "lig_" + base
    return base


def build_pipeline(config: JobConfig, include_analyze=None, include_cluster=True) -> List[SlurmJob]:
    """Return the ordered list of Slurm jobs with their dependency edges.

    Matches the legacy chain and the per-system script selection. The number of
    production repeats follows :attr:`JobConfig.n_repeats` (legacy fixed it at 3).

    When any post-MD analysis is requested, a ``postprocess`` job is appended that
    depends on *all* production repeats (trajectory cleanup + trjcat/align into the
    combined trajectory). The analysis jobs then depend on ``postprocess``:

    * ``cluster`` (when ``include_cluster``) -- CLoNe/PCA structural clustering of
      the combined trajectory; run for every job, as the legacy post-processing
      did, independent of SAXS.
    * ``analyze`` (default: whenever the job has SAXS data, matching the legacy
      ``SAXS_FILE != None`` guard) -- SAXS chi^2/Rg fit summary.
    """
    if include_analyze is None:
        include_analyze = config.uses_saxs

    jobs = [
        SlurmJob("minim1", "minim1/minim1.slurm"),
        SlurmJob("minim2", "minim2/minim2.slurm", depends_on=["minim1"]),
        SlurmJob("nvt", "nvt/" + _variant(config, "nvt.slurm"), depends_on=["minim2"]),
        SlurmJob("npt", "npt/" + _variant(config, "npt.slurm"), depends_on=["nvt"]),
        SlurmJob(
            "setup",
            "production/md_setup/" + _variant(config, "md_setup.slurm"),
            depends_on=["npt"],
        ),
    ]
    repeat_names = []
    for i in range(1, config.n_repeats + 1):
        name = "rep{0}".format(i)
        repeat_names.append(name)
        jobs.append(
            SlurmJob(
                name,
                "production/rep{0}/{1}".format(i, _variant(config, "rep{0}.slurm".format(i))),
                depends_on=["setup"],
            )
        )

    if include_cluster or include_analyze:
        jobs.append(SlurmJob("postprocess", "analysis/postprocess.slurm",
                             depends_on=list(repeat_names)))
    if include_cluster:
        jobs.append(SlurmJob("cluster", "analysis/cluster.slurm", depends_on=["postprocess"]))
    if include_analyze:
        jobs.append(SlurmJob("analyze", "analysis/analyze.slurm", depends_on=["postprocess"]))
    return jobs


def sbatch_command(
    job: SlurmJob,
    slurm_dir: str,
    simulation_dir: str,
    dependency_ids: Optional[Dict[str, str]] = None,
    partition_flag: str = "",
    export_all: bool = True,
) -> List[str]:
    """Build the ``sbatch`` argv for one job.

    ``dependency_ids`` maps a dependency job name to its (symbolic or real) job
    ID. When planning, pass :func:`symbolic_jobids` to render
    ``--dependency=afterok:${JOBID_<name>}`` like the legacy script.
    """
    import os

    dependency_ids = dependency_ids or {}
    argv = ["sbatch", "--parsable", "-J", job.name]
    if partition_flag:
        argv.append(partition_flag)
    if export_all:
        argv.append("--export=ALL")
    if job.depends_on:
        ids = ":".join(dependency_ids.get(dep, dep) for dep in job.depends_on)
        argv.append("--dependency=afterok:{0}".format(ids))
    argv.append(os.path.join(slurm_dir, job.script))
    argv.append(simulation_dir)
    return argv


def symbolic_jobids(jobs: List[SlurmJob]) -> Dict[str, str]:
    """Map each job name to the shell-style token used during planning."""
    return {job.name: "${{JOBID_{0}}}".format(job.name) for job in jobs}


def plan_submissions(
    config: JobConfig,
    slurm_dir: str,
    simulation_dir: str,
    partition_flag: str = "",
):
    """Return ``[(job, sbatch_argv), ...]`` for the whole pipeline (planning only)."""
    jobs = build_pipeline(config)
    ids = symbolic_jobids(jobs)
    return [
        (job, sbatch_command(job, slurm_dir, simulation_dir, ids, partition_flag))
        for job in jobs
    ]

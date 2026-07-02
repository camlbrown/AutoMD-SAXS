"""OpenMM job directory layout (derived data; I/O only on ``create``).

A flat, worker-friendly layout (no Slurm stage dirs): one job directory holding
the prepared/solvated/minimised structures, per-repeat trajectories, extracted
frames, SAXS outputs, and clustering results. Designed to map cleanly onto a
BilboMD worker job directory in Phase 3.
"""

import os
from typing import List


class OpenMMPaths:
    def __init__(self, work_dir=None, job_name=None, n_repeats=3, job_dir=None):
        if n_repeats < 1:
            raise ValueError("n_repeats must be >= 1")
        self.n_repeats = n_repeats
        if job_dir is not None:
            # Explicit output directory (the worker contract: --out OUTPUT_DIR).
            # Results + manifest land directly here, with no job_name nesting.
            self.job_dir = os.path.abspath(job_dir)
            self.work_dir = os.path.dirname(self.job_dir)
            self.job_name = os.path.basename(self.job_dir)
        else:
            if work_dir is None or job_name is None:
                raise ValueError("provide either job_dir, or both work_dir and job_name")
            self.work_dir = os.path.abspath(work_dir)
            self.job_name = job_name
            self.job_dir = os.path.join(self.work_dir, job_name)

    # structures
    @property
    def prepared_pdb(self):
        return os.path.join(self.job_dir, "prepared.pdb")

    @property
    def solvated_pdb(self):
        return os.path.join(self.job_dir, "solvated.pdb")

    @property
    def minimized_pdb(self):
        return os.path.join(self.job_dir, "minimized.pdb")

    @property
    def equilibrated_state(self):
        return os.path.join(self.job_dir, "equilibrated.xml")

    # per-repeat
    def repeat_dir(self, i):
        if not 1 <= i <= self.n_repeats:
            raise ValueError("repeat index {0} out of range".format(i))
        return os.path.join(self.job_dir, "production", "rep{0}".format(i))

    def repeat_trajectory(self, i):
        return os.path.join(self.repeat_dir(i), "production.dcd")

    def repeat_log(self, i):
        return os.path.join(self.repeat_dir(i), "production.log")

    # downstream
    @property
    def frames_dir(self):
        return os.path.join(self.job_dir, "frames")

    @property
    def saxs_dir(self):
        return os.path.join(self.job_dir, "saxs")

    @property
    def ensemble_dir(self):
        return os.path.join(self.job_dir, "ensemble")

    @property
    def clustering_dir(self):
        return os.path.join(self.job_dir, "clustering")

    @property
    def combined_trajectory(self):
        return os.path.join(self.job_dir, "production", "combined.dcd")

    @property
    def combined_topology(self):
        # Solute-only reference written next to combined.dcd by combine_trajectories.
        return os.path.join(self.job_dir, "production", "combined.pdb")

    @property
    def solute_topology(self):
        # Solute-only single-frame topology for reading the protein-only
        # production DCDs (written before the production loop).
        return os.path.join(self.job_dir, "production", "solute_top.pdb")

    @property
    def config_path(self):
        return os.path.join(self.job_dir, "job.json")

    @property
    def manifest_path(self):
        return os.path.join(self.job_dir, "manifest.json")

    @property
    def progress_path(self):
        # Live stage/progress file polled by the BilboMD backend while running.
        return os.path.join(self.job_dir, "progress.json")

    def all_dirs(self) -> List[str]:
        dirs = [self.job_dir, os.path.join(self.job_dir, "production"),
                self.frames_dir, self.saxs_dir, self.ensemble_dir, self.clustering_dir]
        for i in range(1, self.n_repeats + 1):
            dirs.append(self.repeat_dir(i))
        return dirs

    def create(self) -> List[str]:
        made = self.all_dirs()
        for d in made:
            os.makedirs(d, exist_ok=True)
        return made

"""Job directory layout.

The legacy ``simulation_setup.sh`` built a fixed directory tree under
``<protein>_simulation/`` and then wrote every path back into
``configurations.txt`` (``MINIM1_DIR``, ``REPEAT_DIR1`` ...). That coupled the
layout to shell state and made it untestable.

:class:`JobPaths` reproduces the same layout as *derived data* from a base
directory and the protein name, so it can be computed and asserted without
touching the filesystem. ``mkdir`` only happens if you call :meth:`create`.
"""

import os
from typing import List

# Stage subdirectories created directly under the simulation directory, matching
# the legacy ``STAGE_DIRS`` array (repeat dirs are generated per ``n_repeats``).
_STAGE_SUBDIRS = (
    "pdb2gmx",
    "solvate",
    "genion",
    "minim1",
    "minim2",
    "nvt",
    "npt",
    "production",
    "SAXS",
    "ligand_setup",
)


class JobPaths:
    """Resolved filesystem layout for one job. No I/O until :meth:`create`."""

    def __init__(self, base_dir, protein_name, n_repeats=3):
        if n_repeats < 1:
            raise ValueError("n_repeats must be >= 1, got {0}".format(n_repeats))
        self.base_dir = os.path.abspath(base_dir)
        self.protein_name = protein_name
        self.n_repeats = n_repeats
        self.simulation_dir = os.path.join(
            self.base_dir, "{0}_simulation".format(protein_name)
        )

    # --- single-stage directories ---------------------------------------
    def stage_dir(self, stage):
        return os.path.join(self.simulation_dir, stage)

    @property
    def pdb2gmx_dir(self):
        return self.stage_dir("pdb2gmx")

    @property
    def solvate_dir(self):
        return self.stage_dir("solvate")

    @property
    def genion_dir(self):
        return self.stage_dir("genion")

    @property
    def minim1_dir(self):
        return self.stage_dir("minim1")

    @property
    def minim2_dir(self):
        return self.stage_dir("minim2")

    @property
    def nvt_dir(self):
        return self.stage_dir("nvt")

    @property
    def npt_dir(self):
        return self.stage_dir("npt")

    @property
    def production_dir(self):
        return self.stage_dir("production")

    @property
    def saxs_dir(self):
        return self.stage_dir("SAXS")

    @property
    def ligand_setup_dir(self):
        return self.stage_dir("ligand_setup")

    # --- per-repeat directories -----------------------------------------
    def repeat_dir(self, index):
        """Production directory for repeat ``index`` (1-based), e.g. ``production/rep1``."""
        if not 1 <= index <= self.n_repeats:
            raise ValueError(
                "repeat index {0} out of range 1..{1}".format(index, self.n_repeats)
            )
        return os.path.join(self.production_dir, "rep{0}".format(index))

    def processed_dir(self, index):
        return os.path.join(self.repeat_dir(index), "processed")

    def extract_frames_dir(self, index):
        return os.path.join(self.processed_dir(index), "extract_frames")

    # --- aggregate -------------------------------------------------------
    def all_dirs(self) -> List[str]:
        """Every directory the job needs, in creation order (parents first)."""
        dirs = [self.simulation_dir]
        for sub in _STAGE_SUBDIRS:
            dirs.append(self.stage_dir(sub))
        for i in range(1, self.n_repeats + 1):
            dirs.append(self.repeat_dir(i))
            dirs.append(self.processed_dir(i))
            dirs.append(self.extract_frames_dir(i))
        return dirs

    def create(self) -> List[str]:
        """Create the directory tree (idempotent). Returns the created paths."""
        made = self.all_dirs()
        for d in made:
            os.makedirs(d, exist_ok=True)
        return made

    @property
    def config_path(self) -> str:
        """Where the machine-readable job config is written (replaces configurations.txt)."""
        return os.path.join(self.simulation_dir, "job.json")

    @property
    def mdp_dir(self) -> str:
        """Per-job copy of the .mdp templates (so rendering never mutates originals)."""
        return os.path.join(self.simulation_dir, "mdp_files")

    @property
    def manifest_path(self) -> str:
        """Where the run manifest is written."""
        return os.path.join(self.simulation_dir, "manifest.json")

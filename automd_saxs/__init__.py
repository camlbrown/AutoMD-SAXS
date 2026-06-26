"""AutoMD-SAXS workflow package (Phase 1 refactor).

This package is being grown incrementally out of the legacy GROMACS/ATSAS/Slurm
shell pipeline (``simulation_setup.sh`` + ``run_MD.sh`` + ``slurms/``). The goal
of Phase 1 is to lift the *scientific logic* out of bash heredocs and ad-hoc
shell into typed, unit-testable Python while preserving behaviour.

Nothing in this package runs GROMACS, ATSAS, or Slurm. External scientific
binaries are treated as unavailable by default; modules here only *plan* and
*compute*.
"""

__version__ = "0.0.1"

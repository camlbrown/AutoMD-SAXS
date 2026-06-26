"""OpenMM/FoXS pipeline for AutoMD-SAXS (Phase 2).

The higher-accuracy, explicit-solvent all-atom MD + SAXS refinement path, ported
away from GROMACS/Slurm/ATSAS. It reuses the tested Phase 1 foundation
(``manifest``, ``command_runner``, ``dmax``, ``clone``, ``structural``) and adds
an OpenMM MD stack and FoXS/MultiFoXS SAXS analysis.

Design mirrors Phase 1's conventions:

* heavy scientific deps (``openmm``, ``pdbfixer``, ``foxs``/``multi_foxs``) are
  imported lazily and raise a clear ``MissingDependencyError`` if invoked without
  them, so the package imports and its pure logic stays testable on a bare
  interpreter;
* a JSON/YAML job config (:mod:`automd_saxs.openmm.schema`);
* a dry-run ``plan`` that lists the stages and commands without executing.

OpenMM coding patterns (PDBFixer + Modeller.addHydrogens(pH), addSolvent with
padding + ionic strength, LangevinMiddleIntegrator, reporters, energy gating,
platform handling) and the FoXS/MultiFoXS invocation follow the existing BilboMD
worker scripts so the eventual BilboMD integration is natural.
"""

__all__ = ["schema", "foxs", "workflow"]

"""AutoMD-SAXS workflow package (OpenMM branch).

The higher-accuracy, explicit-solvent all-atom MD + SAXS refinement pipeline,
built on OpenMM and FoXS/MultiFoXS. The OpenMM workflow lives in
:mod:`automd_saxs.openmm`; the top-level modules here are the shared, tested
foundation it reuses (``config``, ``manifest``, ``command_runner``, ``dmax``,
``clone``, ``structural``).

The GROMACS/ATSAS/Slurm origin of this project is preserved on the
``AutoMD-SAXs-GROMACS`` branch. This branch ports the workflow away from those
tools. Heavy scientific dependencies (openmm, pdbfixer, mdtraj, FoXS) are
imported lazily and run inside the BilboMD image; the pure logic here is testable
without them.

Entry point: ``python -m automd_saxs.openmm`` (``plan`` / ``run`` / ``validate``).
"""

__version__ = "0.1.0"

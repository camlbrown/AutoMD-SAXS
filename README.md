# AutoMD-SAXS (OpenMM branch)

AutoMD-SAXS is an automated, all-atom molecular dynamics + SAXS refinement
workflow for protein (and, later, protein–ligand) systems.

This is the **OpenMM branch**: the higher-accuracy, explicit-solvent refinement
path, ported away from the original GROMACS/Slurm/ATSAS pipeline onto **OpenMM**
for MD and **FoXS/MultiFoXS** for SAXS fitting. It is designed for clean,
worker-style execution (JSON/YAML config + a machine-readable manifest) so it can
run as a BilboMD worker job.

> The original GROMACS/Slurm/ATSAS workflow is preserved on the
> [`AutoMD-SAXs-GROMACS`](https://github.com/camlbrown/AutoMD-SAXS/tree/AutoMD-SAXs-GROMACS)
> branch.

## Pipeline

```
validate → prepare (PDBFixer + addHydrogens@pH) → solvate (explicit water + ions)
  → minimize → equilibrate (NVT/NPT) → production repeats → frame extraction
  → FoXS (per-frame χ²) → MultiFoXS (ensemble) → CLoNe/PCA clustering → manifest
```

## Layout

```text
automd_saxs/
  config.py, manifest.py, command_runner.py, dmax.py,   # shared foundation
  clone.py, structural.py                               # CLoNe/PCA clustering
  openmm/
    schema.py        # typed OpenMMConfig (JSON/YAML)
    workflow.py      # plan_stages + Workflow.run orchestration
    paths.py         # job directory layout
    prepare.py       # PDBFixer / addHydrogens / addSolvent
    md.py            # minimise / equilibrate / production
    frames.py        # DCD frame extraction + combine
    foxs.py          # FoXS / MultiFoXS command builders + parsers
    cli.py           # `plan` / `run` / `validate`
tests/               # stdlib-only; no OpenMM/FoXS needed
```

## Usage

```bash
# Validate and dry-run plan (no OpenMM/FoXS required)
python -m automd_saxs.openmm plan --config job.json
python -m automd_saxs.openmm validate --config job.json

# Run (executes inside an environment providing OpenMM + FoXS, e.g. the BilboMD image)
python -m automd_saxs.openmm run --config job.json --work-dir ./jobs
# add --dry-run to record the planned stages/commands without executing
```

Example `job.json`:

```json
{
  "job_name": "lyz_refine",
  "pdb": "lysozyme.pdb",
  "saxs": "lysozyme.dat",
  "system": "Protein",
  "force_field": "amber14",
  "water_model": "tip3p",
  "simulation_time_ns": 50,
  "n_repeats": 3,
  "ionic_concentration_M": 0.15,
  "ph": 7.0,
  "temperature_K": 300
}
```

## Dependencies

The MD/SAXS binaries (OpenMM, PDBFixer, mdtraj, FoXS/MultiFoXS) are provided by
the BilboMD Podman image and are imported lazily — the package imports and its
pure logic (config, planning, parsers, summaries) is testable without them.

```bash
pytest          # or: for t in tests/test_*.py; do python "$t"; done
```

## Status

Phase 2 of the project (see `OPENMM_NOTES.md`). Planning and orchestration are
implemented and unit-tested; real end-to-end execution and final validation are
performed in the BilboMD image during worker integration (Phase 3).

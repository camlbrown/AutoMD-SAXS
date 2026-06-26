# Phase 2 notes (OpenMM/FoXS branch)

The OpenMM-compatible all-atom MD + SAXS refinement path, on branch
`AutoMD-SAXS-OpenMM` (branched from the Phase 1 GROMACS branch to reuse its
tested foundation). It ports the workflow away from GROMACS/Slurm/ATSAS to
OpenMM + FoXS/MultiFoXS. ATSAS is no longer required.

## Why a subpackage

The OpenMM stack lives under `automd_saxs/openmm/` so it can reuse the tested
Phase 1 modules (`manifest`, `command_runner`, `dmax`, `clone`, `structural`)
without disturbing the GROMACS modules (kept as legacy reference on this branch).
The CLI is `python -m automd_saxs.openmm`.

## Scientific approach (vs BilboMD)

BilboMD's OpenMM path uses implicit solvent + Rg-restrained conformational
sampling. AutoMD-SAXS is the **higher-accuracy explicit-solvent equilibrium**
refinement path: PME, LangevinMiddleIntegrator, MonteCarloBarostat (NPT),
minimise → equilibrate → production repeats → frame extraction → FoXS fit →
MultiFoXS ensemble → clustering. It reuses BilboMD's *coding patterns* (PDBFixer +
`Modeller.addHydrogens(pH)`, `addSolvent` with padding + ionic strength, energy
gating after minimisation, reporters, platform selection) and the FoXS/MultiFoXS
invocation, so eventual BilboMD integration is natural.

## Modules (`automd_saxs/openmm/`)

| Module | Responsibility | Status |
| --- | --- | --- |
| `schema.py` | Typed `OpenMMConfig` (JSON/YAML), enums, validation, derived steps | done + tested |
| `foxs.py` | FoXS/MultiFoXS argv builders + `.fit` parser + chi^2/Rg summary | done + tested |
| `workflow.py` | `plan_stages` (pure) + manifest + `Workflow.run` orchestration | done; dry-run tested, real unrun |
| `paths.py` | OpenMM job directory layout (worker-friendly) | done + tested |
| `prepare.py` | PDBFixer + addHydrogens + addSolvent (lazy OpenMM) | implemented, unrun |
| `md.py` | minimise / equilibrate / production (lazy OpenMM) | implemented, unrun |
| `frames.py` | DCD frame extraction + combine (lazy mdtraj) | implemented, guard tested |
| `cli.py` | `plan` / `run` / `validate` subcommands | done + tested |

`Workflow.run` is the end-to-end orchestrator: create dirs + job.json + manifest,
then prepare → solvate → minimize → equilibrate → production repeats → frame
extraction → (FoXS per frame + MultiFoXS if SAXS) → clustering, folding outputs
and chi^2/Rg metrics into the manifest. `--dry-run` records stages + the
FoXS/MultiFoXS command plan without importing OpenMM; a real run executes the lazy
OpenMM/mdtraj/FoXS stages and (here) fails cleanly with a `failed` manifest since
those binaries live in the BilboMD image.

## Pipeline (dry-run `plan`)

```
validate_inputs → prepare_structure → solvate → minimize → equilibrate
  → production_rep1..N → extract_frames → [foxs → multifoxs (if SAXS)] → cluster
```

## Testing

Pure-stdlib, no OpenMM/FoXS/numpy/pytest required. Heavy deps (`openmm`,
`pdbfixer`, `foxs`/`multi_foxs`) are imported lazily and raise
`MissingDependencyError` if invoked without them. Phase 2 adds 32 tests
(schema 8, foxs 6, workflow/CLI 6, run/paths/frames 6, plus updates); full repo
suite is 138.

```bash
python -m automd_saxs.openmm plan --config job.json   # dry-run plan
for t in tests/test_openmm_*.py; do python "$t"; done
```

## Pending (Phase 2) — needs the BilboMD image / installed binaries to validate

* Real end-to-end execution of `Workflow.run` inside the BilboMD Podman image
  (OpenMM + mdtraj + FoXS/MultiFoXS). Only dry-run plans + pure logic verified here.
* FoXS fit-file naming: `run_foxs_fits` parses `<frame>*.fit`/`*.dat` defensively;
  confirm the exact FoXS output name in the image and tighten parsing.
* Protein-ligand parameterisation in OpenMM (protein-only is the robust path;
  mark protein-ligand experimental until tested).
* A short CPU smoke test once OpenMM is available.
* Reproducible-seed plumbing is in the schema (`seed`) and integrators; verify
  determinism in the image.

## Phase 3 (later)

Integrate as a BilboMD worker job type. The flat OpenMM job layout
(`paths.OpenMMPaths`), JSON config, manifest, and `python -m automd_saxs.openmm
run` CLI are designed for worker execution. Read `/home/kri42825/bilbomd/CLAUDE.md`
before any BilboMD changes.

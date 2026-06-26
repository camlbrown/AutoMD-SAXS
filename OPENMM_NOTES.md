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
| `workflow.py` | `plan_stages` (pure) + manifest shaping + `Workflow` orchestration | plan done + tested; run skeleton |
| `prepare.py` | PDBFixer + addHydrogens + addSolvent (lazy OpenMM) | implemented, unrun |
| `md.py` | minimise / equilibrate / production (lazy OpenMM) | implemented, unrun |
| `cli.py` | `plan` / `validate` subcommands | done + tested |

## Pipeline (dry-run `plan`)

```
validate_inputs → prepare_structure → solvate → minimize → equilibrate
  → production_rep1..N → extract_frames → [foxs → multifoxs (if SAXS)] → cluster
```

## Testing

Pure-stdlib, no OpenMM/FoXS/numpy/pytest required. Heavy deps (`openmm`,
`pdbfixer`, `foxs`/`multi_foxs`) are imported lazily and raise
`MissingDependencyError` if invoked without them. Phase 2 adds 20 tests
(schema 8, foxs 6, workflow/CLI 6); full repo suite is 132.

```bash
python -m automd_saxs.openmm plan --config job.json   # dry-run plan
for t in tests/test_openmm_*.py; do python "$t"; done
```

## Pending (Phase 2)

* End-to-end `Workflow.run` execution (wire prepare → md → frames → foxs →
  multifoxs → cluster through a runner; needs OpenMM/FoXS to validate).
* Per-frame FoXS orchestration + MultiFoXS ensemble collection into the manifest
  (reuse `analysis.py`/`foxs.summarize_fits`).
* Frame extraction from DCD (mdtraj) + combined-trajectory clustering reuse.
* Protein-ligand parameterisation in OpenMM (start protein-only, mark
  protein-ligand experimental until tested).
* A short CPU smoke test once OpenMM is available.
* Reproducible-seed plumbing is in the schema (`seed`) and integrators.

# CLAUDE.md

This file provides guidance to Claude Code when working in the AutoMD-SAXS repository and when preparing it for later BilboMD integration.

## Project purpose

AutoMD-SAXS is an all-atom molecular dynamics and SAXS-analysis workflow originally written around GROMACS, Slurm, ATSAS tools, shell scripts, Python utilities, and local HPC assumptions.

The current development goal is **not** merely to tidy the old scripts. The goal is to turn the scientific workflow into a clean, testable, BilboMD-compatible all-atom MD refinement pipeline that can eventually run as a new BilboMD worker pathway and optionally accept structures produced by the BilboMD Carbonara pipeline.

AutoMD-SAXS should become the higher-accuracy all-atom refinement path in BilboMD, distinct from existing coarse-grained or structure-generation workflows.

## Key local repositories and paths

Expected local paths on the development machine:

```bash
/home/kri42825/AutoMD-SAXS              # AutoMD-SAXS source repository
/home/kri42825/bilbomd                 # BilboMD repository
/home/kri42825/bilbomd/CLAUDE.md       # BilboMD-specific Claude Code guidance
/home/kri42825/carbonara-pseudoWaxsis  # Carbonara source repository used for BilboMD Carbonara integration
```

Before any BilboMD integration work, read:

```bash
/home/kri42825/bilbomd/CLAUDE.md
```

The BilboMD repository uses local Podman-based development, prebuilt images for normal local operation, and source-run development for newly added worker/backend/UI functionality when images do not yet contain local changes. Follow the BilboMD `CLAUDE.md` rather than inventing new BilboMD infrastructure.

## Current AutoMD-SAXS repository structure

The current repository is a legacy GROMACS/Slurm/ATSAS-oriented workflow. Important files and directories include:

```text
simulation_setup.sh     # interactive/front-end setup script for user choices
run_MD.sh               # main execution script that consumes configurations.txt and runs MD setup/submission
ff_convert/             # PDB atom/residue conversion scripts for AMBER/CHARMM compatibility
ff_files/               # bundled force-field files used by the legacy GROMACS workflow
mdp_files/              # GROMACS .mdp parameter templates
slurms/                 # Slurm job scripts and analysis scripts
examples/               # example workflows and test inputs
AutoMD_SAXS_Manual.pdf  # legacy user-facing manual
workflow.svg            # legacy workflow schematic
automdsaxs.yml          # legacy conda environment definition
```

The legacy public README describes AutoMD-SAXS as a framework for automated molecular dynamics simulation setup and SAXS-based analysis of protein and protein-ligand systems, built for HPC with Slurm. It currently lists GROMACS, Slurm, Conda, and ATSAS as requirements. It is important to not this code is not robust and some code functions could be implemented incorrectly - the original code was writted before agentic coding agents were avialble. The legacy documented user flow is:

```bash
sh simulation_setup.sh -p *Protein*.pdb -s *SAXS*.dat
sh run_MD.sh *Protein*_simulation/
```

## Development phases

There are three major phases. Do not collapse them into one uncontrolled rewrite.

### Phase 1 — Improve the existing GROMACS/ATSAS branch

The user has created a GitHub branch called:

```text
AutoMD-SAXs-GROMACS
```

Branch names are case-sensitive. Inspect available branches with:

```bash
git branch -a
```

Use the exact existing branch name. Do not invent a new branch name without checking.

My git branch is available to access, as we have done for the bilboMD developement code. The repo is found at https://github.com/camlbrown/AutoMD-SAXS.
The first development branch which I cloned locally is https://github.com/camlbrown/AutoMD-SAXS/tree/AutoMD-SAXs-GROMACS.

Phase 1 is allowed to refactor the existing AutoMD-SAXS code substantially, but it should preserve the scientific intent and make the old workflow easier to understand, test, and port.

Important constraints for Phase 1:

- Do not install GROMACS locally.
- Do not install ATSAS locally.
- Do not require Slurm locally.
- Do not run real GROMACS, ATSAS, or Slurm jobs unless the user explicitly says the environment is ready.
- Add dry-run, validation, command-planning, and unit-testable layers wherever useful.
- Treat external scientific commands as unavailable by default and mock/wrap them in tests.
- Prefer Python modules and explicit config files over long shell scripts.
- Keep compatibility wrappers for `simulation_setup.sh` and `run_MD.sh` during early refactoring unless there is a clear reason to replace them.
- Do not remove legacy functionality without documenting what replaced it.

The main Phase 1 task is to transform the code from “interactive shell pipeline” into “well-structured workflow code with clear inputs, outputs, command planning, and testable units”.

### Phase 2 — Create an OpenMM/FoXS/MultiFoXS branch

After the legacy GROMACS branch is cleaned up, create a new branch for the OpenMM-compatible version. Suggested names are:


```text
AutoMD-SAXS-OpenMM
```

Do not create the branch until explicitly instructed.

Phase 2 should port the workflow away from GROMACS, Slurm, and ATSAS.

The target stack is:

- OpenMM for all-atom MD.
- FoXS for SAXS curve fitting where possible.
- MultiFoXS for ensemble selection where possible.
- Open-source Python analysis packages where useful, such as NumPy, SciPy, pandas, MDAnalysis, MDTraj, scikit-learn, Biopython, matplotlib, or existing BilboMD utilities.
- Existing BilboMD OpenMM patterns where applicable.

Avoid dependencies that require separate academic-only downloads, manual licence-gated installation, or local institutional modules. ATSAS must not be required in the OpenMM branch.

Phase 2 does not need to preserve the old two-script shell interface. The OpenMM version should be designed for clean worker execution. A Python package + CLI with JSON/YAML config is preferred.

### Phase 3 — Integrate the OpenMM AutoMD-SAXS pipeline into BilboMD

Only after the OpenMM/FoXS version has a clean execution contract should it be integrated into BilboMD.

The BilboMD integration should add a new worker pipeline/job type for all-atom MD refinement. It should not modify unrelated BilboMD workers, pipelines, backend routes, schemas, UI pages, or infrastructure except where required to register, submit, execute, monitor, and display the new AutoMD-SAXS job.

The future BilboMD pathway should support two starting modes:

1. User uploads an input PDB and SAXS data directly into an AutoMD-SAXS job.
2. User starts from one or more all-atom Carbonara outputs, after Carbonara has produced cg2all-reconstructed PDBs.

The second mode should become an optional Carbonara -> AutoMD-SAXS refinement handoff, not a hard dependency between the two pipelines.

## Git and publishing rules

Unless the user explicitly instructs otherwise:

- Do not push to GitHub.
- Do not open pull requests.
- Do not merge branches.
- Do not commit unless explicitly instructed.
- Do not rewrite shared history.
- Do not force-push.
- Work only in the repository and branch the user has selected.
- For broad refactors, report the intended files and plan before editing.

For AutoMD-SAXS, broad refactoring is allowed once the user has approved the phase, but changes must remain reviewable and scientifically traceable.

For BilboMD, changes must be conservative and scoped. BilboMD is owned by another maintainer. Treat it as a downstream integration target with stricter boundaries.

Carbonara, also owned by another maintainer, was merged into bilboMD previously using claude code. AutoMD-SAXS is the next logical step. 

## Scientific behaviour to preserve from the legacy workflow

The legacy workflow contains scientific decisions that should be preserved or explicitly replaced:

- Protein and protein-ligand setup modes.
- Force-field-aware input preparation for AMBER14SB and CHARMM36m-style workflows.
- Handling of disulfide-bond choices.
- Solvated all-atom simulation with explicit water and ions.
- User-defined ionic concentration, commonly including physiological salt around 0.15 M.
- Dmax/model-size-aware box setup and padding logic.
- Minimisation, equilibration, and production MD stages.
- Multiple production repeats, commonly three repeats.
- Trajectory frame extraction for SAXS analysis.
- Rg/compactness analysis.
- Theoretical SAXS fitting and chi-squared ranking.
- Ensemble selection/reweighting where available.
- Output summaries suitable for publication-quality interpretation.

Do not blindly preserve brittle implementation details. Preserve scientific intent, parameters, and expected outputs.

## Known weak points in the legacy implementation

The following areas should be treated as refactor targets:

- Long interactive shell scripts mixing UI, filesystem layout, scientific logic, command execution, and job submission.
- `configurations.txt` as shell-sourceable state.
- Direct `sed` edits of `.mdp` files.
- Hard-coded Slurm submission assumptions.
- Hard-coded GROMACS command names and module-loading logic.
- ATSAS-specific analysis assumptions.
- `ff_convert` scripts that try to recast atom/residue names in a brittle way. For instance BioExcel provide more robust conversion approaches.
- Force-field conversion logic that is difficult to test.
- Protein-ligand setup via shell loops and ad hoc file movement.
- Lack of dry-run mode.
- Lack of unit tests for parsing, config validation, command generation, Dmax calculation, and output collection.
- Mixed handling of local paths, generated files, and static template files.

The OpenMM branch may replace the entire workflow style if that produces a cleaner and more BilboMD-compatible pipeline.

## Recommended Phase 1 refactor architecture

For the GROMACS/ATSAS branch, a suggested move toward this structure without breaking everything at once is below, but if you have a better approach then do that:

```text
automd_saxs/
  __init__.py
  cli.py                    # command-line entry points
  config.py                 # typed config models and validation
  paths.py                  # job directory and file layout
  inputs.py                 # PDB/SAXS input validation and copying
  dmax.py                   # model Dmax and box-padding calculations
  forcefield.py             # force-field choice and input-prep planning
  command_runner.py         # subprocess wrapper and dry-run support
  gromacs_commands.py       # command construction only, not raw shell blobs
  slurm.py                  # Slurm command planning/submission wrapper
  atsas.py                  # ATSAS command planning/parsing wrapper
  analysis.py               # trajectory/SAXS analysis orchestration
  manifest.py               # output manifest creation
  legacy.py                 # compatibility helpers for old scripts

tests/
  test_config.py
  test_dmax.py
  test_paths.py
  test_command_planning.py
  test_manifest.py
```

Use this as a direction, not a rigid requirement. Inspect the repository first and make a concrete plan before restructuring.

### Dry-run contract

Add or preserve a way to run the pipeline without external scientific binaries:

```bash
python -m automd_saxs.cli plan --config job.json --dry-run
```

A dry run should validate inputs and emit the planned scientific steps/commands without running GROMACS, ATSAS, Slurm, or FoXS. This allows Claude Code to test logic on the local machine even when scientific executables are unavailable.

## Recommended Phase 2 OpenMM architecture

The OpenMM/FoXS branch should be designed as a worker-friendly Python package rather than a set of Slurm scripts.

Suggested structure:

```text
automd_saxs/
  __init__.py
  cli.py
  schema.py                 # JSON/YAML job schema
  workflow.py               # top-level orchestration
  prepare/
    structure.py            # PDB validation, cleanup, protonation assumptions
    forcefield.py           # OpenMM force-field selection
    solvate.py              # water box/ions/barostat setup
    ligands.py              # future ligand support, initially limited if needed
  md/
    minimize.py
    equilibrate.py
    production.py
    reporters.py
  saxs/
    frames.py               # frame extraction
    foxs.py                 # FoXS execution/parsing wrapper
    multifoxs.py            # MultiFoXS execution/parsing wrapper
  analysis/
    rg.py
    rmsd.py
    pca.py
    clustering.py
    summaries.py
  outputs/
    manifest.py
    package_results.py
  utils/
    subprocess.py
    logging.py
    units.py

tests/
```

The first OpenMM version should prioritise a robust protein-only path before attempting full protein-ligand parity. If ligand support is retained or added, it must be modular and honestly marked as supported only when tested.

## OpenMM development guidance

When porting to OpenMM:

- Inspect existing BilboMD OpenMM scripts and utilities before writing new code.
- Reuse BilboMD patterns where reasonable, especially for minimisation, heating/equilibration, production, reporters, and logging.
- Prefer explicit JSON/YAML config over interactive prompts.
- Prefer reproducible seeds where possible.
- Will want a nice UI for the MD setup choices, based on existing bilboMD pipeline
- Keep each MD stage independently testable at the orchestration level.
- Do not require GPUs for basic logic tests.
- Keep CPU-compatible short smoke-test options for development.
- Expose production settings cleanly but provide safe defaults.
- Record all effective parameters in the result manifest.

Potential MD stages:

```text
input validation -> structure preparation -> solvation/ions -> minimisation -> equilibration -> production repeats -> frame extraction -> SAXS analysis -> ensemble/summary outputs
```

## SAXS analysis replacement strategy

The OpenMM branch should remove ATSAS as a required dependency.

Replace ATSAS/CRYSOL-style dependencies with:

- FoXS for theoretical SAXS fits and chi-squared scoring where possible.
- MultiFoXS for ensemble selection where possible.
- Existing BilboMD FoXS parsing/rendering concepts where possible.
- Open-source Python analysis for trajectory and structural summaries e.g pca and clustering.

Expected analysis outputs include, where implemented:

- per-frame/time/per-model chi-squared and rg values;
- ranked best-fit structures;
- ensemble selection output;
- clustering and pca of combined trajectories
- fit files suitable for BilboMD plotting;
- summary tables as CSV/JSON;
- analysis manifest describing what was run and what failed/skipped.

If FoXS or MultiFoXS (available through bilboMD) is not locally available during early development, wrap it behind a command interface and test parsing with fixtures.

## Force-field/input preparation rewrite

The old `ff_convert` directory should be treated as a legacy compatibility layer, not as the desired long-term design.

For the OpenMM branch:

- Design a clean structure-preparation abstraction.
- Validate PDB input early and report actionable errors.
- Avoid silent atom/residue renaming.
- Keep a clear audit trail of any changes made to the input structure.
- Consider existing open-source structure-preparation utilities where appropriate, including tools already present in BilboMD/OpenMM workflows.
- Do not introduce licence-gated or manually downloaded tools as mandatory dependencies.

The first robust implementation can support a narrower set of inputs if it reports limitations clearly.

## BilboMD integration contract

When integrating into BilboMD, the AutoMD-SAXS pipeline should behave like a normal worker-executed job.

Expected contract:

```text
BilboMD UI
  -> backend validates job request and uploads
  -> backend stores job in MongoDB and queues it
  -> worker creates isolated job directory
  -> worker writes AutoMD-SAXS config JSON
  -> worker runs AutoMD-SAXS Python CLI or module
  -> worker streams/records progress
  -> worker collects outputs into a result manifest
  -> backend/UI display status and downloadable results
```

The AutoMD-SAXS CLI should be usable by the worker with a command similar to:

```bash
python -m automd_saxs run --config /path/to/job.json --work-dir /path/to/jobdir --out-dir /path/to/results
```

It should exit non-zero on failure and write clear logs. Avoid shell interpolation of user inputs. Use subprocess argument arrays, not string-built shell commands, unless there is no safe alternative.

### Output manifest

The pipeline should produce a machine-readable manifest, for example:

```json
{
  "pipeline": "automd-saxs-openmm",
  "status": "completed",
  "inputs": {
    "pdb": "input.pdb",
    "saxs": "input.dat"
  },
  "outputs": {
    "trajectories": [],
    "structures": [],
    "saxsFits": [],
    "summaryTables": [],
    "plots": [],
    "logs": []
  },
  "metrics": {
    "bestChi2": null,
    "bestFrame": null,
    "rgMean": null
  }
}
```

The exact schema can evolve, but BilboMD should not need to guess filenames from a loose directory tree.

## Carbonara -> AutoMD-SAXS refinement path

Carbonara and AutoMD-SAXS should remain separate pipelines, but BilboMD should eventually allow a user to refine Carbonara output structures using AutoMD-SAXS.

Expected staged design:

1. AutoMD-SAXS accepts any valid uploaded PDB as input.
2. Carbonara produces all-atom PDBs using cg2all.
3. A user can manually download/upload a Carbonara PDB into AutoMD-SAXS.
4. BilboMD later adds a convenience action such as “Refine with AutoMD-SAXS” from a Carbonara results page.
5. The backend records job lineage from the Carbonara job to the AutoMD-SAXS job.

Do not tightly couple the AutoMD-SAXS worker to Carbonara internals. The handoff should be through explicit files and metadata.

## BilboMD modification boundaries

When working in `/home/kri42825/bilbomd`:

- Read `/home/kri42825/bilbomd/CLAUDE.md` before editing.
- Do not alter unrelated workers or pipelines.
- Do not refactor shared infrastructure unless required for the AutoMD-SAXS job.
- Do not change the Carbonara worker unless implementing the explicit Carbonara -> AutoMD-SAXS handoff.
- Do not change deployment/Podman infrastructure unless necessary and approved.
- Keep UI/backend/schema changes minimal and focused on registering, submitting, monitoring, and displaying AutoMD-SAXS jobs.

Likely BilboMD areas, to be confirmed by inspection:

```text
packages/bilbomd-types
packages/mongodb-schema
apps/backend/src/controllers/jobs
apps/backend/src/validation
apps/backend/src/routes
apps/worker/src/services/functions
apps/worker/src/services/pipelines
apps/ui/src/features
apps/ui/src/schemas
```

Do not assume file names. Inspect existing patterns first.

## User-facing BilboMD workflow target

The eventual BilboMD AutoMD-SAXS UI should be simpler than the legacy terminal prompts.

Expected first UI fields:

- Job name.
- Input PDB.
- SAXS `.dat` file.
- Starting source: user upload or Carbonara output, once handoff exists.
- System type: protein initially; protein-ligand later
- Force field choice, if more than one is truly supported - use bilboMD MD pipelines for inspiration.
- Simulation length.
- Number of repeats.
- Ionic concentration.
- pH
- Disulfide handling.
- Box/padding mode, with sensible automatic defaults.
- Advanced settings collapsed by default.

Results should include:

- final structures;
- trajectories or reduced trajectory files;
- extracted SAXS-analysis frames;
- FoXS fit files;
- MultiFoXS ensemble files where applicable;
- chi-squared/Rg summary tables;
- structural analysis 
- logs;
- plots or data suitable for BilboMD plotting components.

## Testing guidance

### In AutoMD-SAXS

Use tests that do not depend on GROMACS, ATSAS, Slurm, or GPUs by default.

Recommended tests:

```bash
pytest
```

If shell scripts remain, consider:

```bash
bash -n simulation_setup.sh
bash -n run_MD.sh
```

Add tests for:

- config parsing and validation;
- input file validation;
- simulation directory creation;
- Dmax calculation from PDB coordinates;
- box-padding calculation;
- simulation-step calculation from ns and timestep;
- command planning;
- dry-run output;
- output manifest generation;
- FoXS/MultiFoXS parser fixtures;
- error messages for missing external dependencies.

Do not mark a real MD run as tested unless it actually ran.

### In BilboMD

Follow the BilboMD `CLAUDE.md` and use package-specific checks where possible, for example:

```bash
pnpm -F @bilbomd/worker test
pnpm -F @bilbomd/backend test
pnpm -F @bilbomd/ui test
pnpm -F @bilbomd/worker build
pnpm -F @bilbomd/backend build
pnpm -F @bilbomd/ui build
```

Do not claim full BilboMD integration is validated until a browser-driven job has run through the actual local stack.

## Reporting expectations for Claude Code

After each meaningful step, report:

- files inspected;
- files changed;
- why they were changed;
- what was preserved from the legacy workflow;
- what was intentionally changed or removed;
- what was tested;
- what could not be tested because external scientific tools were unavailable;
- next recommended step.

Be explicit about whether a result is:

```text
unit-tested / dry-run tested / parsed from fixture / syntax checked / not tested against real MD
```

## First task for Claude Code in AutoMD-SAXS

Start with inspection and planning only.

Suggested first prompt:

```text
Read CLAUDE.md in this AutoMD-SAXS repository and inspect the current repository structure.

We are in /home/kri42825/AutoMD-SAXS. The current goal is Phase 1: improve the existing GROMACS/ATSAS branch before later creating an OpenMM/FoXS/MultiFoXS branch.

Do not edit files yet.
Do not push, commit, merge, or open a PR.
Do not install GROMACS, ATSAS, or Slurm.
Do not try to run real MD.

Inspect simulation_setup.sh, run_MD.sh, ff_convert, ff_files, mdp_files, slurms, examples, and automdsaxs.yml.

Produce:
1. a repository map;
2. the current workflow as executed by the legacy scripts;
3. fragile areas and redundant code;
4. which parts should be preserved scientifically;
5. a Phase 1 refactor plan that can be tested without GROMACS/ATSAS/Slurm;
6. a proposed file/module structure;
7. the first small implementation task.
```

## Model-use guidance

Use Opus for architecture, repository mapping, and major design decisions.
Use Sonnet subagents for routine implementation, UI work, tests, small refactors, and build/lint fixing once the plan is agreed.

Do not use model capability as a reason to make broad uncontrolled changes. Keep changes staged and reviewable.

## Non-goals for now

Do not do the following unless explicitly instructed:

- Do not install GROMACS.
- Do not install ATSAS.
- Do not install Slurm.
- Do not build a full BilboMD worker before the OpenMM/FoXS CLI contract is stable.
- Do not add a BilboMD UI before the worker contract exists.
- Do not make AutoMD-SAXS depend directly on Carbonara internals.
- Do not remove the ability to understand the legacy workflow.
- Do not claim industry compatibility without checking dependencies.
- Do not silently drop protein-ligand support; mark it as unsupported/experimental if it is not carried through the OpenMM port.

## Success criteria

### End of Phase 1

- Legacy workflow is mapped and documented.
- Major shell-script responsibilities are separated into testable modules or clearly planned for separation.
- Dry-run/command-planning mode exists or is clearly implemented.
- Core calculations and config handling have tests.
- No GROMACS/ATSAS/Slurm execution is required for tests.
- Scientific behaviour preserved or explicitly documented where changed.

### End of Phase 2

- OpenMM-based all-atom MD pipeline can run from a clean CLI/config contract.
- ATSAS is no longer required.
- FoXS/MultiFoXS analysis is wrapped or integrated where possible.
- Outputs are described by a machine-readable manifest.
- A short CPU-compatible smoke test exists if practical.
- The code is ready to be called by a BilboMD worker.

### End of Phase 3

- BilboMD has a new AutoMD-SAXS worker pathway.
- The job can be submitted from the browser.
- The worker runs the OpenMM AutoMD-SAXS pipeline.
- Results are collected and displayed/downloadable.
- Carbonara outputs can be used as AutoMD-SAXS inputs, at least through an explicit handoff or upload path.
- Unrelated BilboMD workers and pipelines remain unchanged.

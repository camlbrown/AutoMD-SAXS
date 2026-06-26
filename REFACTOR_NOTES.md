# Phase 1 refactor notes (GROMACS/ATSAS branch)

This documents the incremental Phase 1 refactor: the legacy shell pipeline is
being lifted into a testable Python package (`automd_saxs/`) while the original
entry points (`simulation_setup.sh`, `run_MD.sh`) remain as thin wrappers that
delegate to it. Nothing here runs GROMACS, ATSAS, or Slurm; external scientific
binaries are treated as unavailable by default.

## Package layout

| Module | Responsibility |
| --- | --- |
| `automd_saxs/config.py` | Typed `JobConfig` + enums; validation; JSON round-trip |
| `automd_saxs/legacy.py` | Parse old `configurations.txt` (as data, never sourced) → `JobConfig` |
| `automd_saxs/dmax.py` | Model Dmax + box-padding calculators |
| `automd_saxs/paths.py` | `JobPaths` — the job directory tree as derived data |
| `automd_saxs/gromacs_commands.py` | argv builders for the setup pipeline |
| `automd_saxs/ligand.py` | Protein-ligand: structure split, topology splicing, tool plan |
| `automd_saxs/slurm.py` | Submission dependency DAG + `sbatch` argv |
| `automd_saxs/mdp.py` | `.mdp` templating (render, never `sed -i`) |
| `automd_saxs/atsas.py` | ATSAS `shanum`/`crysol`/`gajoe` command builders + output parsers |
| `automd_saxs/analysis.py` | SAXS result collection, ranked CSV summary, metrics → manifest |
| `automd_saxs/postprocess.py` | Trajectory cleanup + trjcat/align command plans (per-repeat + combine) |
| `automd_saxs/clone.py` | CLoNe clustering algorithm (faithful port; lazy numpy/scipy/sklearn) |
| `automd_saxs/structural.py` | Clustering/PCA orchestration + pure-stdlib summary stats + parsers |
| `automd_saxs/manifest.py` | Machine-readable run manifest |
| `automd_saxs/command_runner.py` | Subprocess wrapper + dry-run + missing-binary handling |
| `automd_saxs/cli.py` | `plan`, `run`, `postprocess`, `analyze`, `cluster`, `dmax`, `render-mdp` subcommands |

## Legacy → new mapping

| Legacy behaviour | Replacement |
| --- | --- |
| `configurations.txt` (sourced shell + functions) | `JobConfig` / `job.json`; `legacy.py` reads old files without executing them |
| Inline numpy/scipy Dmax heredoc + `bc` padding in `run_MD.sh` | `automd_saxs dmax` (pure-Python, numpy-optional) |
| `sed -i` on tracked `mdp_files/md.mdp` | `automd_saxs render-mdp` into a per-job mdp copy |
| Submission chain in `run_MD.sh` | `automd_saxs slurm` DAG (`plan`) |
| Hand-built directory tree | `automd_saxs.paths.JobPaths` |

## Verifiable bug fixes

1. **Ionic concentration ignored.** `run_MD.sh` hardcoded `genion ... -conc 0.15`,
   discarding the user's `IONIC_CONCENTRATION`. Now `-conc "${IONIC_CONCENTRATION:-0.15}"`,
   and `JobConfig.genion_concentration()` is the source of truth in planning.
2. **In-place template mutation.** `sed -i` edited the git-tracked `md.mdp` every
   run. `simulation_setup.sh` now copies `mdp_files/*.mdp` into
   `<job>/mdp_files/` and points `MDP_DIR` there; rendering targets that copy.
   Slurm scripts read `$MDP_DIR`, so they transparently use the per-job copy.
3. **`configurations.txt` as executable shell.** Replaced by data-only parsing.

## Wrapper delegation

* `run_MD.sh <dir> --dry-run` → `automd_saxs plan` (validate + print plan, no execution).
* `run_MD.sh` step count / box padding → `render-mdp` / `dmax` subcommands.
* `simulation_setup.sh` runs a non-fatal `plan` validation preview after writing the config.

## Testing

Pure-stdlib, no GROMACS/ATSAS/Slurm/numpy/pytest required. Each test file runs
under pytest and standalone (`python tests/<file>.py`). Current suite: 112 tests.

Modules with heavy scientific dependencies (numpy/scipy/scikit-learn/mdtraj for
`clone.py`/`structural.run_clustering`) import those lazily and raise a clear
`MissingDependencyError` if invoked without them, so the package still imports and
its pure helpers/parsers stay testable on a bare interpreter.

```bash
for t in tests/test_*.py; do python tests/$(basename "$t"); done
# or, in an environment with pytest:
pytest
```

## Run vs. plan

* `automd_saxs plan` — validate + print the plan (always dry).
* `automd_saxs run` — prepare the job dir (create tree, write `job.json`, render
  the production mdp, write the manifest) and then execute the GROMACS setup +
  `sbatch` submission through `CommandRunner`. `--dry-run` records commands
  instead of executing; a real run on a machine without the binaries fails
  cleanly with a `failed` manifest. Protein-ligand *execution* is still gated
  (planning works); protein execution and Slurm submission are wired but
  untested here (no binaries).

## Not yet ported / pending

* Protein-ligand **execution** (the `run` command plans it but refuses to execute
  it; the per-ligand file-staging between stage dirs is not yet reproduced).
* Inter-stage file staging for the protein `run` path (works as a command plan;
  real multi-binary execution unverified without GROMACS).
* SAXS analysis is wired: `automd_saxs analyze --analysis-dir <dir>` collects
  `crysol_summary.txt` / `*.fit` / `Rg_distr*.txt`, writes a chi^2-ranked
  `saxs_summary.csv`, and folds `bestChi2`/`bestFrame`/`rgMean` into a manifest.
  Phase 2 will swap ATSAS for FoXS/MultiFoXS behind the same `analysis.py` shape.
* Post-MD chain is wired into the Slurm DAG:
  `reps → postprocess → {cluster, analyze}`.
  - `postprocess` (`slurms/analysis/postprocess.slurm`) `afterok` all repeats:
    per-repeat trajectory cleanup (strip solvent, pbc nojump, fit, extract
    frames) then `trjcat` + align into `combined_aligned.xtc`. Ported in
    `automd_saxs/postprocess.py`; run via `automd_saxs postprocess`.
  - `cluster` (`slurms/analysis/cluster.slurm`, always; CLoNe pdc sweep 1,3,5,7)
    and `analyze` (`slurms/analysis/analyze.slurm`, SAXS jobs only) each depend
    on `postprocess`.
* CLoNe/PCA structural analysis is ported: `automd_saxs/clone.py` (faithful
  algorithm) + `automd_saxs/structural.py` (orchestration, pure-stdlib cluster
  summary stats, and parsers for `PCA_coords.txt`/`results.txt`), exposed as
  `automd_saxs cluster`. The numerical path needs numpy/scipy/sklearn/mdtraj;
  the summary/parsing helpers do not.
* Still pending: actually *running* the GROMACS post-processing, the ATSAS tools,
  and the numerical clustering (verified here only via dry-run plans, the
  lazy-dep error path, and pure helpers, as the dev machine lacks the stack).

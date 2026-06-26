# Phase 3 scope — integrate AutoMD-SAXS as a BilboMD worker

Goal: add a new **all-atom MD refinement** job type to BilboMD that runs the
OpenMM-branch pipeline (`python -m automd_saxs.openmm run`) as a worker job,
following existing BilboMD patterns conservatively. No changes to the 7 existing
job types (`pdb`, `crd`, `auto`, `alphafold`, `sans`, `scoper`, `multi`) or the
in-progress `carbonara` worker except where strictly required to register the new
job.

This is a **scope/plan only**. No BilboMD code is changed yet. Per
`/home/kri42825/bilbomd/CLAUDE.md`: work locally, do not commit/push/PR/merge or
create branches without explicit instruction; do not work on `main` or directly
on `feature/carbonara-worker`.

## Integration contract

BilboMD treats AutoMD-SAXS as a normal worker job. The worker's only coupling to
this package is one command + one output file:

```
worker creates job dir
  -> writes job.json (from the submitted form params)
  -> runs:  python -m automd_saxs.openmm run --config job.json --work-dir <jobdir>
  -> reads <jobdir>/<job_name>/manifest.json   (status, outputs, metrics)
  -> streams status; surfaces outputs/metrics to the UI
```

The OpenMM branch already provides exactly this contract: a JSON config
(`OpenMMConfig`), a CLI `run`, a machine-readable `manifest.json`
(`status`, `inputs`, `outputs.{trajectories,structures,saxsFits,summaryTables,…}`,
`metrics.{bestChi2,bestFrame,rgMean}`), and a flat job layout.

## Two starting modes (per project brief)

1. **Upload** — user uploads a PDB + SAXS `.dat` → AutoMD-SAXS job. (Phase 3a)
2. **Carbonara handoff** — start from a cg2all-reconstructed all-atom PDB produced
   by the Carbonara worker; user picks "Refine with AutoMD-SAXS" from a Carbonara
   result. The handoff is via explicit files + job lineage metadata, not a code
   dependency on Carbonara internals. (Phase 3c, after upload mode works.)

## Prerequisite: package availability in the worker image

OpenMM/PDBFixer/mdtraj/FoXS/MultiFoXS already live in the BilboMD worker image.
The `automd_saxs` package must also be importable there. Options (decision needed):
- vendor/copy the package into the worker image build, or
- `pip install` it from the repo/a wheel during image build.
Either way the worker calls it as `python -m automd_saxs.openmm`.

## Files to create/modify (mirrors the existing OpenMM `auto` job)

Job-type name proposal: string `automd-saxs`, discriminator `BilboMdAutoMDSAXS`,
UI route `jobs/automd-saxs`. (Names are a decision point.)

**packages/bilbomd-types** (1)
- `src/jobs/jobs.ts` — add `'automd-saxs'` to `JobType`; add `BilboMDAutoMDSAXSDTO`;
  add it to the `BilboMDMongoJobDTO` union. Reuse/adapt `OpenMMParametersDTO`.

**packages/mongodb-schema** (4)
- `src/interfaces/jobInterface.ts` — `IBilboMDAutoMDSAXSJob` (`__t: 'BilboMdAutoMDSAXS'`).
- `src/models/Job.ts` — schema + `Job.discriminator('BilboMdAutoMDSAXS', …)`.
- `src/models/index.ts` — export the model.
- `src/constants/jobConstants.ts` — display name.

**apps/backend/src** (5–6)
- `validation/automdSaxsJobSchema.ts` (+ export in `validation/index.ts`) — Yup
  schema: PDB required, SAXS `.dat` optional, MD params (sim length, repeats,
  ionic strength, pH, temperature, force field, disulfide).
- `controllers/jobs/handleBilboMDAutoMDSAXS.ts` — copy `handleBilboMDAutoJob.ts`:
  accept files, validate, build params, create job doc, `writeJobParams`,
  `queueJob({type:'automd-saxs', …})`.
- `controllers/jobs/createJob.ts` — add a case in `dispatchBilboMDJob()`.
- `routes/jobs.ts` — optional dedicated route (generic POST also works).
- multer fields: ensure `pdb_file` + `dat_file` accepted (already are for auto).

**apps/worker/src** (2–3)
- `services/pipelines/bilbomd-automd-saxs.ts` — pipeline: init job dir → write
  `job.json` from the stored params → run the AutoMD-SAXS CLI (reuse the spawn
  pattern in `services/functions/openmm-functions.ts`) → parse `manifest.json` →
  collect outputs/metrics → status + email.
- `workerHandlers/bilboMdHandler.ts` — add entry in `getPipelineExecutor()`.
- (optional) a small `services/functions/automd-saxs-functions.ts` for the spawn +
  manifest-parse helpers.

**apps/ui/src** (5)
- `features/automd-saxsjob/NewAutoMDSAXSJobForm.tsx` — form fields from the brief:
  job name, PDB, SAXS `.dat`, system type (protein first), force field, sim
  length, repeats, ionic concentration, pH, disulfide, box/padding mode
  (advanced collapsed).
- `schemas/BilboMDAutoMDSAXSJobSchema.ts` — Formik/Yup.
- `layout/MainLayout/index.tsx` — sidebar entry (consider a feature flag like the
  existing `enableBilboMd*`).
- `routes/MainRoutes.tsx` — lazy import + route.
- results: reuse existing results components; map manifest `metrics`/`outputs`
  (chi²/Rg, best structure, FoXS fits, plots) onto them.

Total ≈ 15 files, mostly new, following the `auto`/OpenMM template.

## Staged plan

- **3a — worker happy path:** package in image; minimal types+schema+backend
  handler+worker pipeline; submit a protein-only upload job, run it in the image,
  collect the manifest. Start with a tiny CPU smoke config.
- **3b — UI:** add the form, route, sidebar entry (feature-flagged), results
  mapping; browser-driven job.
- **3c — Carbonara handoff:** "Refine with AutoMD-SAXS" action from a Carbonara
  result; record job lineage; reuse the same worker.

## Job lifecycle / results

Submitted → queued (BullMQ `bilbomd` queue) → running (worker) →
completed/failed. The worker maps the AutoMD-SAXS `manifest.status` to BilboMD's
`JobStatus`; `manifest.outputs`/`metrics` populate the results view. UI polls via
the existing job API slices.

## Testing (BilboMD side)

Narrowest first, per BilboMD CLAUDE.md (Vitest, `__tests__/`):
`pnpm -F @bilbomd/worker {lint,build,test}`, then backend, then UI. Real
end-to-end MD validation happens in the Podman image (where OpenMM/FoXS exist) —
this is where the AutoMD-SAXS pipeline finally runs for real.

## Open decisions (need user input before 3a)

1. Job-type naming (`automd-saxs` / `BilboMdAutoMDSAXS` / `jobs/automd-saxs`)?
2. How to ship `automd_saxs` into the worker image (vendor vs pip install)?
3. New BilboMD branch name (e.g. `feature/automd-saxs-worker`) — requires explicit
   go-ahead since CLAUDE.md forbids creating branches without it; must NOT reuse
   `feature/carbonara-worker`.
4. Protein-only for v1 (recommended), protein-ligand marked experimental?
5. CPU smoke defaults vs production defaults for the first image run.

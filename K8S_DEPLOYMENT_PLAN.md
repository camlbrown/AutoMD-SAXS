# Kubernetes deployment plan — BilboMD + Carbonara + AutoMD-SAXS at Diamond

A planning document to brief Diamond scientific computing on deploying the
BilboMD platform (with our Carbonara and AutoMD-SAXS worker pipelines) onto a
Kubernetes (k8s) cluster. It records the current architecture, what already
exists, the genuinely hard problems, what we should change in our own code
first, and the questions to put to sci comp.

This spans two repos:
- `bilbomd` (worktrees: `feature/carbonara-worker`, `feature/automd-saxs-worker`)
- `AutoMD-SAXS` (the external `automd-saxs` CLI we pip-install into the worker)

---

## 0. TL;DR for the sci-comp conversation

- BilboMD is a multi-service app: **UI + backend (Express) + worker (Node/BullMQ)
  + MongoDB + Redis**, plus GPU MD tooling (OpenMM/CHARMM/FoXS) in the worker.
- **It already runs on k8s at NERSC** (a Rancher/“Spin” cluster) and a **Helm
  chart already exists** in `bilbomd/infra/helm/`. So the web tier is a solved
  problem; we adapt that chart for Diamond.
- **The one big architectural decision:** at NERSC the k8s worker is a *thin
  client* — it does **not** run MD in the pod; it submits Slurm jobs to the
  Perlmutter HPC cluster via an API. The GPU MD runs on HPC, not in k8s. For
  Diamond we must choose: **(A) run MD inside GPU k8s pods**, or **(B) submit MD
  to a separate Diamond HPC/Slurm cluster** like NERSC does. This choice drives
  everything else.
- **The hardest technical problem for our additions:** both Carbonara and
  AutoMD-SAXS currently run by having the worker **launch another container at
  runtime** (`podman run …`). That "container-launching-a-container" pattern is
  an anti-pattern on k8s. The clean fix is to **bake the tools into the worker
  image and call them in-process** — which AutoMD-SAXS already supports and
  Carbonara can be moved to.

---

## 1. Current architecture (what we run today, locally)

Six services on a Docker/Podman network (`bilbomd/infra/docker-compose.local*.yml`):

| Service | Image | Role | GPU |
| --- | --- | --- | --- |
| mongodb | mongo:8.3.x | job/user state | no |
| redis | redis:8.x | BullMQ job queue | no |
| backend | ghcr.io/bl1231/bilbomd-backend | REST API, enqueues jobs | no |
| worker | ghcr.io/bl1231/bilbomd-worker | runs the pipelines | **yes** |
| ui | ghcr.io/bl1231/bilbomd-ui | React SPA | no |
| colabfold / of3 | (prod) | structure-prediction sidecars | yes |

Data flow: **UI → backend → (Mongo doc + Redis/BullMQ enqueue) → worker → results
on a shared volume → backend serves results → UI**.

Shared state that must persist:
- **MongoDB** (jobs, users), **Redis** (the queue),
- **a shared job/results volume** mounted at `DATA_VOL=/bilbomd/uploads`
  (`infra/uploads-dev` locally). Every job is a directory `…/uploads/<uuid>/`.

The worker base image (`bilbomd-worker-base.dockerfile`) is built **FROM
`nvidia/cuda:12.9.2-devel`** and bakes in CHARMM, OpenMM (in `/opt/envs/openmm`),
and IMP/FoXS/MultiFoXS (`/usr/bin/foxs`, `/usr/bin/multi_foxs`). GPU is requested
locally via `--device nvidia.com/gpu=all` (podman CDI) or compose
`deploy.resources.reservations.devices`.

---

## 2. What already exists toward k8s (reuse, don't reinvent)

`bilbomd/infra/helm/` is a working Helm chart used at NERSC:
- `Chart.yaml`, `values.yaml` + `values-dev.yaml` / `values-prod.yaml`
- `templates/`: deployments for backend, worker, ui, redis, mongo; an ingress;
  Services; ConfigMaps/Secrets.
- `infra/HELM_NOTES.md` documents the NERSC/Rancher (Spin) deploy: kubeconfig
  setup, the required Secrets (`bilbomd-secrets`, `mongo-secrets`, `ui-tls`,
  `ghcr`, `registry-nersc`, `sfapi-priv-key`), and `helm install/upgrade`.

**For Diamond, the bulk of the work is adapting this chart** (registry, storage
class, ingress/TLS, GPU) — the manifests for the web tier mostly exist.

---

## 3. The key architectural fork: where does the MD actually run?

This is the most important thing to settle with sci comp.

**NERSC model (what upstream does):** the k8s ("Spin") cluster runs Mongo, Redis,
backend, UI, and a **thin worker** with `USE_NERSC=true`. That worker does **not**
run OpenMM/CHARMM; it submits Slurm jobs to **Perlmutter** (a separate GPU HPC
cluster) through the Superfacility API, and a separate `bilbomd-perlmutter-worker`
image runs the real MD on HPC compute nodes. k8s here is just the web/orchestration
tier; **the GPU compute is on HPC.**

**Two options for Diamond:**

- **Option A — MD inside k8s GPU pods.** The worker pod itself requests
  `nvidia.com/gpu: 1` and runs OpenMM/FoXS in-image. Simpler conceptually (one
  cluster), but requires the Diamond k8s cluster to have GPU nodes + the NVIDIA
  device plugin/operator, and we must solve the nested-container problem (§5).
- **Option B — k8s web tier + submit MD to Diamond HPC (Wilson/Slurm).** Mirrors
  NERSC: k8s runs the web stack, the worker submits MD to a Diamond HPC/Slurm
  cluster. Reuses the upstream NERSC submission pattern, keeps heavy GPU compute
  off k8s, but needs a Diamond equivalent of the Superfacility API / job-submit
  bridge and is more moving parts.

**Recommendation to discuss:** start with **Option A** for Carbonara + AutoMD-SAXS
(both are self-contained container workloads already), because Diamond may not
have a NERSC-style API bridge; revisit Option B if GPU-in-k8s capacity is limited.
Ask sci comp which Diamond clusters have GPUs and whether k8s GPU scheduling is
available.

---

## 4. The hard problems (and recommended solutions)

### 4.1 Nested containers — the #1 problem
Today the worker runs each Carbonara/AutoMD-SAXS task by spawning **another
container** from inside the worker:
- Carbonara: `podman run --rm -v <jobdir>:/job <carbonara-image> python runner …`
- AutoMD-SAXS: the `infra/automd-saxs/automd-saxs-run` shim does
  `podman run --rm -v $DATA_VOL:$DATA_VOL <automd-saxs-image> automd-saxs run …`

On k8s this "Docker-out-of-Docker" needs the container runtime socket or a
privileged pod — both discouraged/blocked on shared clusters, and the bind-mount
host-path semantics break (the spawned container sees the *node* FS, not the
pod's PVC).

**Solution — eliminate the nesting:**
- **AutoMD-SAXS: already solved by design.** The shim is a *dev-only* convenience.
  In a real deployment we **pip-install `automd-saxs` into the worker image** (the
  worker base already has OpenMM/FoXS) and set
  `AUTOMD_SAXS_BIN=/opt/envs/openmm/bin/automd-saxs`. The worker then calls the
  CLI **as a subprocess in its own pod** — no nested container. (Our
  `infra/automd-saxs/Dockerfile` is exactly this image layer.)
- **Carbonara: bake the runtime into the worker image too** (or a sibling image)
  and call its Python entry points directly, instead of `podman run`. Carbonara's
  runtime is larger (~5 GB: cg2all, torch, pyfoxs) so this inflates the worker
  image; alternatively run each Carbonara task as a **k8s Job** (one pod per task)
  that the worker creates via the k8s API. Baking-in is simpler operationally;
  k8s-Jobs is cleaner isolation. Decide per sci-comp guidance.

### 4.2 Shared storage
Every job is a directory under `DATA_VOL=/bilbomd/uploads`, written by the worker
and read by the backend (to serve results) — so **backend and worker need the
same volume**. On k8s that is a **ReadWriteMany (RWX) PersistentVolumeClaim**
(e.g. NFS/CephFS). NERSC uses a hostPath onto its CFS shared filesystem; Diamond
will need an RWX storage class. **Ask sci comp:** what RWX storage class is
available, and expected capacity (MD trajectories are large — see §6.2).

### 4.3 GPU scheduling
- Request `nvidia.com/gpu: 1` in the worker (and AutoMD-SAXS/Carbonara Job) pods;
  the NVIDIA device plugin sets `CUDA_VISIBLE_DEVICES`.
- The worker already round-robins GPUs from `CUDA_VISIBLE_DEVICES` into `OMM_GPU_ID`
  for OpenMM (`apps/worker/src/services/functions/openmm-functions.ts`). AutoMD-SAXS
  picks the fastest platform (CUDA) automatically and accepts an explicit platform.
- The current Helm worker template has **no GPU resource block** (NERSC offloads
  to HPC); for Option A we must add `resources.limits."nvidia.com/gpu": 1`.
- **Ask sci comp:** GPU models, count per node, device plugin / GPU operator
  present, and whether MPS (multi-process service, for >1 job per GPU) is allowed.

### 4.4 Images & registry
Web images are on `ghcr.io/bl1231/*`. Carbonara + AutoMD-SAXS runtime images are
built locally (`localhost/…`). For Diamond we must **push all images to a registry
the cluster can pull from** (Diamond-internal registry or a pull-through cache).
Pin **versioned tags** (AutoMD-SAXS is already tagged `v0.1.0`) so worker code and
runtime image versions stay in lockstep.

### 4.5 Secrets / config
Split the env (already enumerated in `infra/.env.example`) into:
- **ConfigMap** (non-secret): `MONGO_HOSTNAME/PORT/DB`, `REDIS_HOST/PORT`,
  `DATA_VOL`, `EXAMPLE_DATA`, feature flags, `AUTOMD_SAXS_BIN`, etc.
- **Secret**: `*_TOKEN_SECRET`, `SESSION_SECRET`, `HASH_IP_SALT`,
  `MONGO_USERNAME/PASSWORD`, ORCID secrets.

---

## 5. What we should change in OUR code BEFORE k8s

These are the actionable prep items on our side (mostly already aligned):

1. **AutoMD-SAXS: make the in-image CLI the default execution path.**
   - Production: worker image installs `automd-saxs` (our `infra/automd-saxs/Dockerfile`
     layer, or fold into `bilbomd-worker-base`), and `AUTOMD_SAXS_BIN` points at
     `/opt/envs/openmm/bin/automd-saxs`. The worker already shells out to
     `$AUTOMD_SAXS_BIN run --config … --out …` — so **no worker code change**, just
     image + env. The podman shim stays dev-only.
   - Add `scikit-learn` to the worker base (we currently add it in our layer) so
     clustering works without a separate layer.
2. **Paths must be PVC-relative, not host-absolute.** The worker already uses
   `DATA_VOL` for the job dir and passes absolute paths into the CLI; as long as the
   CLI runs **in the same pod/PVC mount** (not a nested container), the absolute
   paths resolve. Removing the nested container (item 1) fixes this automatically.
3. **Disk discipline (already improved).** AutoMD-SAXS now writes **protein-only**
   frames/trajectories (not the full water box) — ~75× smaller — and clustering is
   non-fatal. This matters because MD output is large; see §6.2.
4. **Carbonara: plan the de-nesting** (bake-in vs k8s-Job). No change needed for a
   first AutoMD-SAXS-only k8s trial, but required before Carbonara runs on k8s.
5. **Health/liveness:** the worker already exposes a config/health server (port
   3000, `/config`); good for k8s liveness/readiness probes.

---

## 6. Notes & gotchas

### 6.1 One queue, one consumer semantics
All workers drain a single Redis `bilbomd` queue. On k8s, scaling the worker to N
replicas means N consumers — fine (BullMQ handles it), but each replica needs a
GPU and shared PVC. (Locally we hit a related issue: two host workers on one queue
caused "Unknown job type" — on k8s, keep image/queue versions consistent.)

### 6.2 Storage sizing
Explicit-solvent MD produces large trajectories. Even protein-only frames + DCDs
add up across repeats and jobs. The local dev volume is only an 8 GB home quota —
**too small**; this caused a real job failure. Diamond PVC must be sized for many
multi-GB jobs, with a retention/cleanup policy (the worker has a cleanup step).

### 6.3 Job longevity
MD jobs run for minutes to hours. k8s worker pods must not be evicted mid-job;
set appropriate resource requests, `terminationGracePeriod`, and avoid aggressive
autoscaling on the worker. If using k8s-Jobs per task, set `backoffLimit`/`ttl`.

---

## 7. Questions to put to Diamond sci comp

1. Does the Diamond k8s cluster have **GPU nodes** + the **NVIDIA device plugin /
   GPU operator**? Which GPU models, how many per node? (Decides Option A vs B.)
2. Is there a **NERSC-style HPC submission bridge** (Slurm + API) if we prefer
   Option B (k8s web tier + HPC compute)?
3. What **RWX storage class** is available for the shared job/results volume, and
   what capacity/quota can we get?
4. What **container registry** should images be pushed to / pulled from? Any image
   size limits (Carbonara runtime is ~5 GB)?
5. **Ingress/TLS**: ingress controller, cert mechanism, DNS for the UI.
6. **Privileges**: are privileged pods / runtime sockets ever allowed (they
   shouldn't be needed once we de-nest), and what PodSecurity level applies?
7. Can we run **MPS** (multiple MD jobs per GPU) or is it one-job-per-GPU?

---

## 8. Suggested phased plan

1. **Phase 0 — facts:** answer §7 with sci comp; pick Option A or B.
2. **Phase 1 — web tier on Diamond k8s:** adapt the existing Helm chart (registry,
   RWX PVC, ConfigMap/Secret, ingress/TLS); deploy Mongo, Redis, backend, UI, and
   a non-GPU worker that runs the simple (non-MD) job types. Validate the platform.
3. **Phase 2 — AutoMD-SAXS on k8s (de-nested):** worker image includes the
   `automd-saxs` CLI + scikit-learn; `AUTOMD_SAXS_BIN` set to the in-image CLI; add
   `nvidia.com/gpu: 1` to the worker (Option A) or wire HPC submission (Option B).
   Run a small protein job end-to-end. **AutoMD-SAXS is the natural first GPU
   pipeline to land** because it's already a clean CLI + manifest contract.
4. **Phase 3 — Carbonara on k8s:** de-nest (bake-in or k8s-Job) and deploy.
5. **Phase 4 — hardening:** storage retention/cleanup, autoscaling policy,
   monitoring, MPS if permitted.

---

## Appendix — key files

- Helm chart: `bilbomd/infra/helm/` (+ `infra/HELM_NOTES.md`)
- Worker base image: `bilbomd/apps/worker/bilbomd-worker-base.dockerfile`
- Compose topology: `bilbomd/infra/docker-compose.local*.yml`
- Env keys: `bilbomd/infra/.env.example`, `infra/worker.env`
- AutoMD-SAXS runtime image + shim: `bilbomd/infra/automd-saxs/{Dockerfile,automd-saxs-run}`
- AutoMD-SAXS CLI contract: `AutoMD-SAXS` repo, `automd_saxs/openmm/cli.py`, tag `v0.1.0`
- Carbonara runtime image: `bilbomd/infra/carbonara/Dockerfile.carbonara-allatom-runtime`
- Worker GPU selection: `bilbomd/apps/worker/src/services/functions/openmm-functions.ts`

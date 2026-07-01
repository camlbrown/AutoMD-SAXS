# Route to k8s — from current BilboMD + our pipelines to running on the Diamond cluster

A concrete, sequenced action plan for the Friday meeting. Companion to the fuller
`K8S_DEPLOYMENT_PLAN.md`; this one is "what do we actually do, in order".

## The one-paragraph answer to "can't we just containerise + add to the images + edit the helm and deploy?"

**Yes — that is exactly the route.** BilboMD already ships Docker images + a Helm
chart and was deployed on k8s at NERSC. To add our two pipelines we (1) get their
tools *into the worker image* so the worker runs them as ordinary subprocesses in
its own pod (no "container launching a container"), (2) make three small edits to
the existing Helm chart for Diamond (registry, storage, GPU), (3) push the images
to a registry the cluster can pull from, and (4) `helm install`. **We are not
rewriting BilboMD or the pipelines** — Carbonara and AutoMD-SAXS keep the exact
same code and behaviour; we only change *how the worker invokes them* and *where
the images/volumes live.*

## The mental model (say this in the meeting)

- **k8s runs images, not source.** "Deploying new code" = build a new image →
  push it to a registry → tell Helm to use that image tag (`helm upgrade`). The
  cluster never compiles our code.
- BilboMD = 5 images (mongo, redis, backend, worker, ui). Only the **worker**
  needs GPU + the science tools. Our two pipelines are extra work the **worker**
  does — so they attach to the worker image.

## The only real technical obstacle, and how we remove it

Today, for local dev, the worker runs each Carbonara/AutoMD-SAXS task by launching
*another* container (`podman run …`). That "Docker-out-of-Docker" is the one thing
that doesn't belong on k8s. **Fix = bake the tools into the worker image and call
them in-process** (a normal subprocess), so there is no nested container:

- **AutoMD-SAXS** — trivial: `pip install` the `automd-saxs` CLI into the worker
  image's OpenMM env and set `AUTOMD_SAXS_BIN=/opt/envs/openmm/bin/automd-saxs`.
  The worker already calls `$AUTOMD_SAXS_BIN run --config … --out …` as a
  subprocess — **no worker code change**, just the image + one env var. (Our
  `infra/automd-saxs/Dockerfile` already builds this layer; fold it into the
  worker image / `bilbomd-worker-base`.)
- **Carbonara** — its runtime (cg2all + torch + FoXS, ~5 GB) already exists as a
  standalone image. Bake that runtime into the worker image and change the 6
  `podman run <carbonara-image> python <script> …` calls in
  `carbonara-functions.ts` to plain `spawn(python, [script, …])` (same scripts,
  same args, minus the `podman run -v … <image>` wrapper; the bind-mount is moot
  because it's the same pod filesystem). Behaviour is identical.

> Alternative (cleaner long-term, more moving parts): instead of baking Carbonara
> into the worker, have the worker create a **k8s Job per task** from the existing
> Carbonara image (needs a ServiceAccount + RBAC). Recommend **bake-in first** to
> get running, evolve to k8s-Jobs later if image size / isolation warrants it.

## Three edits to the existing Helm chart (that's all Diamond needs)

The chart is `infra/helm/` (worker template: `templates/worker-deployment.yaml`).
It was written for NERSC; for Diamond change only:

1. **Registry** — point `values` image `registry/owner/…:tag` at the Diamond
   registry (or a pull-through cache of ghcr) where we push our images.
2. **Storage** — NERSC used `hostPath` mounts onto its shared filesystem for
   `/bilbomd/uploads` (+ logs). Diamond: replace with a **ReadWriteMany PVC** that
   both backend and worker mount (both must see the same job dir). One value swap
   in the deployment templates.
3. **GPU** — the worker template currently has `resources: {}`. Add
   `resources.limits."nvidia.com/gpu": 1` so the worker pod gets a GPU (assuming
   Diamond's k8s has GPU nodes + the NVIDIA device plugin — a key meeting question).

Everything else (Deployments, Services, ConfigMaps/Secrets, ingress) is reusable
largely as-is.

## The route, in order

**Before/at the meeting — get these facts (they decide a few branches):**
- Does the Diamond k8s cluster have **GPU nodes + the NVIDIA device plugin**?
  (If not: MD must run on Diamond HPC instead — see Option B in the full plan.)
- Which **container registry** do we push to, and any image size limit? (Carbonara
  image is ~5 GB.)
- What **ReadWriteMany storage class** is available, and what capacity/quota?
  (MD output is large; the local 8 GB dev quota is far too small.)
- **Ingress/TLS** mechanism + DNS for the UI; and the **PodSecurity** level
  (confirm no privileged pods needed — once we de-nest, they aren't).
- Do they offer **GitOps (ArgoCD/Flux)** or do we `helm upgrade` manually?

**Implementation steps (mostly doable now, independent of the cluster):**
1. **De-nest AutoMD-SAXS** into the worker image (pip install + env). *(smallest)*
2. **De-nest Carbonara** into the worker image (bake runtime + `spawn` instead of
   `podman run`). *(the main engineering task; Carbonara is otherwise finished)*
3. **Build the worker image** with both baked in; build backend/ui unchanged.
4. **Push images** to the Diamond registry with **versioned tags**.
5. **Fork the Helm values** into `values-diamond.yaml`: registry, RWX PVC, GPU on
   the worker, ingress/TLS, ConfigMap/Secret (split per the full plan).
6. **`helm install`** the web tier first (mongo, redis, backend, ui, non-GPU
   worker) → validate the platform with a simple job type.
7. **Enable the GPU worker** → run a **Carbonara** job end-to-end (it's user-ready
   → first pipeline to land), then an **AutoMD-SAXS** job.
8. **Hardening:** storage retention/cleanup, Mongo + PVC backups, monitoring.

**Minimum viable first deploy to aim for:** web tier + GPU worker on Diamond k8s,
running one real Carbonara job from the browser. Everything else builds on that.

## Ongoing development after it's on k8s

You keep developing **locally exactly as now** (the host-source dev stack + login
link). To ship a change: commit → CI (or a manual step) builds & pushes the image →
bump the tag in `values-diamond.yaml` → `helm upgrade` (or GitOps auto-syncs).
Local stays the fast loop; k8s is for shared/production. Versioned tags make
rollbacks one command (`helm rollback`). Secrets live as k8s Secrets, never in git;
Mongo + the PVC persist across rollouts.

## What we already have (leverage points to state plainly)

- A working **Helm chart** proven on k8s at NERSC — we adapt, not author.
- **Carbonara** functionally complete + already containerised (image exists).
- **AutoMD-SAXS** with a clean CLI + JSON-config + manifest contract, already
  designed to run as an in-image subprocess (the dev-only podman shim is discarded
  in production).
- A **local dev workflow** (host-source stack) that stays unchanged for iteration.

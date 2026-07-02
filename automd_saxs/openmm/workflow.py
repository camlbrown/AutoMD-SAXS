"""Top-level OpenMM workflow orchestration.

:func:`plan_stages` is pure and testable -- it describes the pipeline for a given
config without importing OpenMM. :class:`Workflow` runs the stages; it imports the
MD/prep layers lazily so a real run requires OpenMM but planning does not.
"""

import glob
import os
from typing import List, Tuple

from ..command_runner import CommandRunner, MissingDependencyError
from ..manifest import Manifest, STATUS_PLANNED, STATUS_COMPLETED, STATUS_RUNNING, STATUS_FAILED
from .paths import OpenMMPaths
from .schema import OpenMMConfig

PIPELINE_NAME = "automd-saxs-openmm"

Stage = Tuple[str, str]


def plan_stages(config: OpenMMConfig) -> List[Stage]:
    """Ordered ``(name, detail)`` stages for this job (no execution)."""
    stages = [
        ("validate_inputs", "validate PDB/SAXS inputs and config"),
        ("prepare_structure",
         "PDBFixer: missing atoms/residues + addHydrogens(pH={0:g}){1}".format(
             config.ph, ", keep disulfides" if config.disulfide else "")),
        ("solvate",
         "explicit {0} water + ions @ {1:g} M, box padding {2}".format(
             config.water_model.value, config.ionic_concentration_M,
             "{0:g} nm".format(config.box_padding_nm) if config.box_padding_nm
             else "auto (from model Dmax)")),
        ("minimize",
         "energy minimisation (maxIterations={0})".format(config.minimize_max_iterations)),
        ("equilibrate",
         "NVT+NPT equilibration {0:g} ns ({1} steps)".format(
             config.equilibration_ns, config.equilibration_steps())),
    ]
    for i in range(1, config.n_repeats + 1):
        stages.append((
            "production_rep{0}".format(i),
            "{0} ns at {1:g} fs, {2} K ({3} steps)".format(
                config.simulation_time_ns, config.timestep_fs,
                config.temperature_K, config.production_steps())))
    stages.append(("extract_frames", "frame stride {0}".format(config.frame_stride)))
    if config.uses_saxs:
        stages.append(("foxs", "per-frame FoXS chi^2 vs experimental SAXS"))
        stages.append(("multifoxs", "MultiFoXS ensemble selection"))
    stages.append(("cluster", "CLoNe/PCA structural clustering of combined trajectory"))
    return stages


def new_manifest(config: OpenMMConfig, status: str = STATUS_PLANNED) -> Manifest:
    """A manifest pre-populated for the OpenMM pipeline."""
    manifest = Manifest(config=None, status=status)
    manifest.pipeline = PIPELINE_NAME
    manifest.inputs = {"pdb": config.pdb, "saxs": config.saxs}
    manifest.parameters = {
        "system": config.system.value,
        "forceField": config.force_field.value,
        "waterModel": config.water_model.value,
        "simulationTimeNs": config.simulation_time_ns,
        "timestepFs": config.timestep_fs,
        "productionSteps": config.production_steps(),
        "nRepeats": config.n_repeats,
        "temperatureK": config.temperature_K,
        "ionicConcentrationM": config.ionic_concentration_M,
        "pH": config.ph,
        "disulfide": config.disulfide,
        "usesSaxs": config.uses_saxs,
        "seed": config.seed,
    }
    manifest.config = None
    manifest.extra_config = config.to_dict()
    return manifest


class Workflow:
    """Runs the OpenMM pipeline. MD layers are imported lazily on ``run``."""

    def __init__(self, config: OpenMMConfig, work_dir=None, dry_run=True, out_dir=None):
        self.config = config
        self.work_dir = work_dir
        self.out_dir = out_dir
        self.dry_run = dry_run

    def _paths(self) -> OpenMMPaths:
        """Resolve the job layout: ``out_dir`` (results land there directly) wins;
        otherwise ``work_dir``/<job_name>."""
        if self.out_dir:
            return OpenMMPaths(job_dir=self.out_dir, n_repeats=self.config.n_repeats)
        return OpenMMPaths(self.work_dir or ".", self.config.job_name, self.config.n_repeats)

    def plan(self) -> List[Stage]:
        return plan_stages(self.config)

    def run(self):
        """Prepare the job and run (or, in dry-run, plan) the OpenMM pipeline.

        Always performs local orchestration (create dirs, write job.json + a
        manifest). In dry-run it records the planned stages and the FoXS/MultiFoXS
        commands without importing OpenMM. In a real run it executes the OpenMM MD
        stages, frame extraction, SAXS fitting, and clustering -- these import
        OpenMM/mdtraj/FoXS lazily and run inside the BilboMD image in Phase 3.
        Returns a result dict; on failure writes a ``failed`` manifest and raises.
        """
        cfg = self.config
        paths = self._paths()
        paths.create()
        with open(paths.config_path, "w") as handle:
            handle.write(cfg.to_json())

        manifest = new_manifest(cfg, STATUS_PLANNED if self.dry_run else STATUS_RUNNING)
        runner = CommandRunner(dry_run=self.dry_run)

        if self.dry_run:
            from . import foxs
            for name, _ in plan_stages(cfg):
                if name == "foxs":
                    # attach a representative command (frames not yet present)
                    manifest.add_step("foxs", command=foxs.foxs_fit_command(
                        cfg.saxs, os.path.join(paths.frames_dir, "structure_<i>.pdb")))
                elif name == "multifoxs":
                    manifest.add_step("multifoxs", command=foxs.multifoxs_command(
                        cfg.saxs, [os.path.join(paths.frames_dir, "structure_<i>.pdb.dat")],
                        output=os.path.join(paths.ensemble_dir, "ensemble.dat")))
                else:
                    manifest.add_step(name)
            manifest.write(paths.manifest_path)
            return {"status": "planned", "job_dir": paths.job_dir,
                    "manifest": paths.manifest_path, "stages": [n for n, _ in plan_stages(cfg)]}

        try:
            return self._execute(cfg, paths, manifest, runner)
        except Exception as exc:  # noqa: BLE001 - record then re-raise
            manifest.set_status(STATUS_FAILED).add_note(str(exc))
            manifest.write(paths.manifest_path)
            raise

    def _write_progress(self, paths, cfg, stage, current_repeat=0):
        """Write a small live-progress file the BilboMD backend polls while running.

        Stage-level signal; the per-repeat ns counter is derived by the backend
        from the StateDataReporter step column in each production log.
        """
        import json
        import time

        try:
            data = {
                "stage": stage,
                "nRepeats": cfg.n_repeats,
                "currentRepeat": current_repeat,
                "productionSteps": cfg.production_steps(),
                "equilibrationSteps": cfg.equilibration_steps(),
                "timestepFs": cfg.timestep_fs,
                "simulationTimeNs": cfg.simulation_time_ns,
                "reportIntervalSteps": cfg.report_interval_steps,
                "updatedAt": time.time(),
            }
            with open(paths.progress_path, "w") as handle:
                json.dump(data, handle)
        except Exception:  # noqa: BLE001 - progress reporting must never fail a run
            pass

    def _execute(self, cfg, paths, manifest, runner):
        from . import prepare, md, frames, foxs
        from .. import structural

        # 1. structure preparation + solvation
        self._write_progress(paths, cfg, "prepare_structure")
        audit = prepare.prepare_structure(cfg, cfg.pdb, paths.prepared_pdb)
        manifest.add_step("prepare_structure", status=STATUS_COMPLETED)
        manifest.set_parameter("prepAudit", audit)
        padding = prepare.resolve_box_padding(cfg, cfg.pdb)
        manifest.set_parameter("boxPaddingNm", padding)
        self._write_progress(paths, cfg, "solvate")
        prepare.solvate(cfg, paths.prepared_pdb, paths.solvated_pdb, padding)
        manifest.add_step("solvate", status=STATUS_COMPLETED)

        # 2. minimise + equilibrate
        self._write_progress(paths, cfg, "minimize")
        min_result = md.minimize(cfg, paths.solvated_pdb, paths.minimized_pdb)
        manifest.add_step("minimize", status=STATUS_COMPLETED)
        manifest.set_parameter("minimizedEnergyKJ", min_result["potentialEnergyKJ"])
        self._write_progress(paths, cfg, "equilibrate")
        md.equilibrate(cfg, paths.minimized_pdb, paths.equilibrated_state)
        manifest.add_step("equilibrate", status=STATUS_COMPLETED)

        # 3. production repeats. Fan out across GPUs when more than one is
        # visible: each repeat is pinned to its own CUDA device and they run
        # concurrently (concurrency = min(n_repeats, n_gpus)). On a single GPU
        # this collapses to the original sequential loop (concurrency 1).
        #
        # Production DCDs are written PROTEIN-ONLY (atom subset) so they are
        # ~75x smaller than the full solvated box -- this avoids disk-quota
        # blow-ups on long/large runs and cuts I/O. Downstream frame extraction
        # therefore reads them with the solute-only topology, not minimized_pdb.
        try:
            solute_top = frames.write_solute_topology(
                paths.minimized_pdb, paths.solute_topology)
            solute_idx = frames.solute_indices(paths.minimized_pdb)
        except Exception as exc:  # noqa: BLE001 - fall back to full-system DCDs
            manifest.add_note("protein-only DCD disabled (non-fatal): {0}".format(exc))
            solute_top, solute_idx = None, None
        # Topology used to read the production DCDs back: solute-only when we
        # wrote protein-only DCDs, otherwise the full solvated minimized PDB.
        prod_topology = solute_top if solute_idx else paths.minimized_pdb

        trajectories = [paths.repeat_trajectory(i) for i in range(1, cfg.n_repeats + 1)]
        n_gpus = md.gpu_count()
        concurrency = max(1, min(cfg.n_repeats, n_gpus))

        def _run_repeat(i):
            md.production(cfg, paths.equilibrated_state, paths.minimized_pdb,
                          paths.repeat_trajectory(i), paths.repeat_log(i),
                          device_index=((i - 1) % n_gpus), atom_subset=solute_idx)
            return i

        if concurrency <= 1:
            for i in range(1, cfg.n_repeats + 1):
                self._write_progress(paths, cfg, "production_rep{0}".format(i),
                                     current_repeat=i)
                _run_repeat(i)
                manifest.add_step("production_rep{0}".format(i), status=STATUS_COMPLETED)
                manifest.add_output("trajectories", paths.repeat_trajectory(i))
        else:
            from concurrent.futures import ThreadPoolExecutor, as_completed
            self._write_progress(paths, cfg, "production_rep1",
                                 current_repeat=cfg.n_repeats)
            manifest.add_note(
                "production: {0} repeats across {1} GPU(s)".format(
                    cfg.n_repeats, n_gpus))
            with ThreadPoolExecutor(max_workers=concurrency) as pool:
                futures = {pool.submit(_run_repeat, i): i
                           for i in range(1, cfg.n_repeats + 1)}
                for fut in as_completed(futures):
                    fut.result()  # propagate any repeat failure
            for i in range(1, cfg.n_repeats + 1):
                manifest.add_step("production_rep{0}".format(i), status=STATUS_COMPLETED)
                manifest.add_output("trajectories", paths.repeat_trajectory(i))

        # 4. frame extraction + combine. Track which repeat each frame came from
        # and its simulation time (ns) so the dashboard can toggle x-axis
        # (frame/ns) and per-repeat series.
        self._write_progress(paths, cfg, "extract_frames", current_repeat=cfg.n_repeats)
        # ns per extracted frame: DCD saves every report_interval_steps, then we
        # sub-sample by frame_stride.
        ns_per_frame = (cfg.frame_stride * cfg.report_interval_steps
                        * cfg.timestep_fs / 1.0e6)
        frame_pdbs = []
        frame_meta = []  # aligned with frame_pdbs: {repeat, frameInRep, timeNs}
        for i, traj in enumerate(trajectories, 1):
            pdbs = frames.extract_frames(
                traj, prod_topology,
                os.path.join(paths.frames_dir, "rep{0}".format(i)), cfg.frame_stride)
            for k, pdb in enumerate(pdbs):
                frame_pdbs.append(pdb)
                frame_meta.append({"repeat": i, "frameInRep": k,
                                   "timeNs": round(k * ns_per_frame, 4)})
        combined = frames.combine_trajectories(
            trajectories, prod_topology, paths.combined_trajectory, cfg.frame_stride)
        manifest.add_step("extract_frames", status=STATUS_COMPLETED)

        # Per-repeat structural time-series (Rg / Cα-RMSD / SASA over time) for the
        # Structural Analysis tab. Computed on the solute; advisory (never fatal).
        time_series = []
        for i in range(1, cfg.n_repeats + 1):
            ts = frames.structural_timeseries(
                paths.repeat_trajectory(i), prod_topology, stride=cfg.frame_stride)
            # Total Energy comes from the StateDataReporter log (one row per DCD
            # frame); the time-series is strided, so frame k maps to log row
            # k*frame_stride. Attach it as the solvent-free energy trace.
            energies = frames.read_energy_series(paths.repeat_log(i))
            for row in ts:
                log_idx = row["frame"] * cfg.frame_stride
                row["energy"] = (energies[log_idx]
                                 if 0 <= log_idx < len(energies) else None)
                row["repeat"] = i
                row["timeNs"] = round(row["frame"] * ns_per_frame, 4)
                time_series.append(row)
        if time_series:
            manifest.set_parameter("timeSeries", time_series)

        # 5. SAXS analysis (FoXS per frame + MultiFoXS ensemble)
        if cfg.uses_saxs:
            self._write_progress(paths, cfg, "foxs", current_repeat=cfg.n_repeats)
            records = foxs.run_foxs_fits(frame_pdbs, cfg.saxs, runner, cwd=paths.saxs_dir)
            # Attach a per-frame Rg (computed from the protein frame) so the
            # dashboard can chart chi^2 and Rg, and rgMean is populated. Records
            # are aligned with frame_pdbs order; use a global frame index so
            # frames from different repeats don't collide.
            per_frame = []
            for gi, (pdb, rec) in enumerate(zip(frame_pdbs, records)):
                rec["frame"] = gi
                rec["rg"] = frames.rg_of_pdb(pdb)
                meta = frame_meta[gi] if gi < len(frame_meta) else {}
                per_frame.append({
                    "frame": gi,
                    "repeat": meta.get("repeat"),
                    "frameInRep": meta.get("frameInRep"),
                    "timeNs": meta.get("timeNs"),
                    "chi2": rec.get("chi2"),
                    "rg": rec["rg"],
                })
            manifest.set_parameter("perFrame", per_frame)
            summary = foxs.summarize_fits(records)
            for key in ("bestChi2", "bestFrame", "rgMean"):
                if summary.get(key) is not None:
                    manifest.set_metric(key, summary[key])
            # MultiFoXS ensemble fit over the MD frames (the FoXS-based analogue
            # of the legacy GAJOE ensemble step). Pass the frame PDBs so multi_foxs
            # computes partial profiles and fits c1/c2, matching the Carbonara
            # invocation. Lenient: never fail the job on an ensemble-step error.
            if frame_pdbs:
                # Map a member filename multi_foxs echoes back to the GLOBAL frame
                # index (frame_pdbs order == perFrame order). Prefer a full-path
                # match (basenames collide across repeats: structure_5 in each).
                basename_to_gi = {}
                for gi, p in enumerate(frame_pdbs):
                    basename_to_gi.setdefault(os.path.basename(p), gi)
                path_to_gi = {p: gi for gi, p in enumerate(frame_pdbs)}

                def _member_frame(name):
                    if name in path_to_gi:
                        return path_to_gi[name]
                    return basename_to_gi.get(os.path.basename(name))

                try:
                    mf = foxs.run_multifoxs(
                        frame_pdbs, cfg.saxs, runner, paths.ensemble_dir,
                        num_states=5, frame_of=_member_frame)
                    if mf:
                        manifest.set_parameter("multifoxs", mf)
                    for f in sorted(glob.glob(
                            os.path.join(paths.ensemble_dir, "ensembles_size_*.txt"))
                            + glob.glob(os.path.join(
                                paths.ensemble_dir, "multi_state_model_*"))):
                        manifest.add_output("saxsFits", f)
                except Exception as exc:  # noqa: BLE001 - ensemble step is advisory
                    manifest.add_note("MultiFoXS ensemble failed (non-fatal): {0}".format(exc))
            manifest.add_step("foxs", status=STATUS_COMPLETED)
            manifest.add_step("multifoxs", status=STATUS_COMPLETED)

        # 6. structural clustering of the combined trajectory. Lenient: clustering
        # is a secondary analysis (needs scikit-learn), so a failure here records
        # a note and is skipped rather than failing an otherwise-successful job.
        self._write_progress(paths, cfg, "cluster", current_repeat=cfg.n_repeats)
        try:
            # The combined trajectory is solute-only, so cluster against the
            # solute-only topology written beside it (not the solvated structure).
            # Sweep a range of CLoNe pdc values (granularity), like the legacy
            # run_CLoNe.sh; PCA is shared across pdc so it is computed once.
            pdc_values = [1, 2, 3, 4, 5, 6, 7]
            default_pdc = 4  # matches the previous single-run default
            cluster = structural.run_clustering_sweep(
                combined, paths.combined_topology, paths.clustering_dir,
                at_sel="name CA", pca=2, pdc_values=pdc_values)
            if cluster.get("skipped"):
                manifest.add_step("cluster", status=STATUS_COMPLETED)
                manifest.add_note("clustering skipped: {0}".format(cluster["skipped"]))
            else:
                manifest.add_step("cluster", status=STATUS_COMPLETED)
                sweep = cluster.get("sweep", [])
                for entry in sweep:
                    for f in entry.get("outputFiles", []):
                        manifest.add_output("summaryTables", f)
                # Per-pdc clusterings for the UI toggle: {pdc, nClusters, labels}.
                clusterings = [{"pdc": e["pdc"], "nClusters": e["nClusters"],
                                "labels": e["labels"]} for e in sweep]
                manifest.set_parameter("clusterings", clusterings)
                default_entry = next(
                    (e for e in sweep if e["pdc"] == default_pdc),
                    (sweep[0] if sweep else None))
                if default_entry:
                    manifest.set_parameter("nClusters", default_entry["nClusters"])
                    manifest.set_parameter("defaultPdc", default_pdc)
                # PCA scatter (chi^2-coloured) for the Structural Analysis tab.
                # PCA row i aligns with combined frame i == perFrame[i]. Cluster
                # label uses the default pdc; the UI can switch via `clusterings`.
                pca_rows = cluster.get("pcaCoords", [])
                if pca_rows:
                    pf = manifest.parameters.get("perFrame", []) if cfg.uses_saxs else []
                    labels = default_entry["labels"] if default_entry else []
                    pca = []
                    for i, row in enumerate(pca_rows):
                        if len(row) < 2:
                            continue
                        pca.append({
                            "frame": i,
                            "pc1": row[0],
                            "pc2": row[1],
                            "chi2": (pf[i].get("chi2") if i < len(pf) else None),
                            "cluster": (int(labels[i]) if i < len(labels) else None),
                        })
                    manifest.set_parameter("pca", pca)
        except MissingDependencyError as exc:
            manifest.add_step("cluster", status=STATUS_FAILED)
            manifest.add_note("clustering skipped: {0}".format(exc))
        except Exception as exc:  # noqa: BLE001 - clustering must not fail the job
            manifest.add_step("cluster", status=STATUS_FAILED)
            manifest.add_note("clustering failed (non-fatal): {0}".format(exc))

        manifest.set_status(STATUS_COMPLETED)
        manifest.write(paths.manifest_path)
        return {"status": "completed", "job_dir": paths.job_dir,
                "manifest": paths.manifest_path}

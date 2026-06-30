"""Top-level OpenMM workflow orchestration.

:func:`plan_stages` is pure and testable -- it describes the pipeline for a given
config without importing OpenMM. :class:`Workflow` runs the stages; it imports the
MD/prep layers lazily so a real run requires OpenMM but planning does not.
"""

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

        # 3. production repeats
        trajectories = []
        for i in range(1, cfg.n_repeats + 1):
            self._write_progress(paths, cfg, "production_rep{0}".format(i), current_repeat=i)
            md.production(cfg, paths.equilibrated_state, paths.minimized_pdb,
                          paths.repeat_trajectory(i), paths.repeat_log(i))
            trajectories.append(paths.repeat_trajectory(i))
            manifest.add_step("production_rep{0}".format(i), status=STATUS_COMPLETED)
            manifest.add_output("trajectories", paths.repeat_trajectory(i))

        # 4. frame extraction + combine
        self._write_progress(paths, cfg, "extract_frames", current_repeat=cfg.n_repeats)
        frame_pdbs = []
        for i, traj in enumerate(trajectories, 1):
            frame_pdbs += frames.extract_frames(
                traj, paths.minimized_pdb,
                os.path.join(paths.frames_dir, "rep{0}".format(i)), cfg.frame_stride)
        combined = frames.combine_trajectories(
            trajectories, paths.minimized_pdb, paths.combined_trajectory, cfg.frame_stride)
        manifest.add_step("extract_frames", status=STATUS_COMPLETED)

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
                per_frame.append({"frame": gi, "chi2": rec.get("chi2"), "rg": rec["rg"]})
            manifest.set_parameter("perFrame", per_frame)
            summary = foxs.summarize_fits(records)
            for key in ("bestChi2", "bestFrame", "rgMean"):
                if summary.get(key) is not None:
                    manifest.set_metric(key, summary[key])
            profiles = [r["fitFile"] for r in records if r.get("fitFile")]
            if profiles:
                ensemble = os.path.join(paths.ensemble_dir, "ensemble.dat")
                foxs.run_multifoxs(profiles, cfg.saxs, runner, ensemble, cwd=paths.ensemble_dir)
                manifest.add_output("saxsFits", ensemble)
            manifest.add_step("foxs", status=STATUS_COMPLETED)
            manifest.add_step("multifoxs", status=STATUS_COMPLETED)

        # 6. structural clustering of the combined trajectory. Lenient: clustering
        # is a secondary analysis (needs scikit-learn), so a failure here records
        # a note and is skipped rather than failing an otherwise-successful job.
        self._write_progress(paths, cfg, "cluster", current_repeat=cfg.n_repeats)
        try:
            # The combined trajectory is solute-only, so cluster against the
            # solute-only topology written beside it (not the solvated structure).
            cluster = structural.run_clustering(
                combined, paths.combined_topology, paths.clustering_dir,
                at_sel="name CA", pca=2)
            if cluster.get("skipped"):
                manifest.add_step("cluster", status=STATUS_COMPLETED)
                manifest.add_note("clustering skipped: {0}".format(cluster["skipped"]))
            else:
                manifest.add_step("cluster", status=STATUS_COMPLETED)
                for f in cluster.get("outputFiles", []):
                    manifest.add_output("summaryTables", f)
                manifest.set_parameter("nClusters", cluster.get("nClusters"))
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

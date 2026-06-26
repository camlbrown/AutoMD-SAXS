"""Top-level OpenMM workflow orchestration.

:func:`plan_stages` is pure and testable -- it describes the pipeline for a given
config without importing OpenMM. :class:`Workflow` runs the stages; it imports the
MD/prep layers lazily so a real run requires OpenMM but planning does not.
"""

from typing import List, Tuple

from ..manifest import Manifest, STATUS_PLANNED, STATUS_COMPLETED, STATUS_RUNNING
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

    def __init__(self, config: OpenMMConfig, work_dir: str, dry_run: bool = True):
        self.config = config
        self.work_dir = work_dir
        self.dry_run = dry_run

    def plan(self) -> List[Stage]:
        return plan_stages(self.config)

    def run(self):
        """Execute the pipeline (requires OpenMM/FoXS). Not yet implemented end-to-end.

        Planning is available now via :func:`plan_stages`; the executing stages
        (:mod:`automd_saxs.openmm.prepare`, :mod:`automd_saxs.openmm.md`) import
        OpenMM lazily and raise a clear error if it is absent.
        """
        from . import prepare, md  # noqa: F401  (lazy: imports OpenMM)

        raise NotImplementedError(
            "End-to-end OpenMM execution is under development. Use plan_stages() "
            "for the dry-run plan; prepare/md provide the OpenMM-backed stages.")

"""Machine-readable run manifest.

BilboMD (and any worker) should not have to guess output filenames from a loose
directory tree. This module produces the manifest described in the project brief
-- ``pipeline``, ``status``, ``inputs``, ``outputs``, ``metrics`` -- extended
with the effective configuration and parameters so a run is self-describing.

The manifest is built incrementally (``add_output``, ``set_metric``,
``set_status``) and serialised with :meth:`Manifest.to_json`.
"""

import json
from typing import Any, Dict, List, Optional

from .config import JobConfig

PIPELINE_NAME = "automd-saxs-gromacs"

#: Output categories, matching the brief's schema (camelCase for BilboMD).
_OUTPUT_CATEGORIES = (
    "trajectories",
    "structures",
    "saxsFits",
    "summaryTables",
    "plots",
    "logs",
)

# Recognised status values.
STATUS_PLANNED = "planned"
STATUS_RUNNING = "running"
STATUS_COMPLETED = "completed"
STATUS_FAILED = "failed"


class Manifest:
    """Builder for the AutoMD-SAXS result manifest."""

    def __init__(self, config: Optional[JobConfig] = None, status: str = STATUS_PLANNED):
        self.pipeline = PIPELINE_NAME
        self.status = status
        self.config = config
        self.inputs = {"pdb": None, "saxs": None}  # type: Dict[str, Optional[str]]
        self.parameters = {}  # type: Dict[str, Any]
        self.steps = []  # type: List[Dict[str, Any]]
        self.outputs = {cat: [] for cat in _OUTPUT_CATEGORIES}  # type: Dict[str, List[str]]
        self.metrics = {"bestChi2": None, "bestFrame": None, "rgMean": None}  # type: Dict[str, Any]
        self.notes = []  # type: List[str]

        if config is not None:
            self.inputs["pdb"] = config.protein_file
            self.inputs["saxs"] = config.saxs_file
            self.parameters = self._derive_parameters(config)

    @staticmethod
    def _derive_parameters(config: JobConfig) -> Dict[str, Any]:
        return {
            "system": config.system.value,
            "forceField": config.force_field_name,
            "waterModel": config.water_model(),
            "boxShape": config.box_shape.value,
            "ionicConcentrationM": config.ionic_concentration_M,
            "simulationTimeNs": config.simulation_time_ns,
            "timestepFs": config.timestep_fs,
            "productionNsteps": config.number_of_steps(),
            "nRepeats": config.n_repeats,
            "disulfide": config.disulfide,
            "usesSaxs": config.uses_saxs,
        }

    # --- mutation --------------------------------------------------------
    def set_status(self, status: str) -> "Manifest":
        self.status = status
        return self

    def set_parameter(self, key: str, value: Any) -> "Manifest":
        self.parameters[key] = value
        return self

    def add_step(self, name: str, command=None, status: str = STATUS_PLANNED) -> "Manifest":
        step = {"name": name, "status": status}  # type: Dict[str, Any]
        if command is not None:
            step["command"] = list(command)
        self.steps.append(step)
        return self

    def add_output(self, category: str, path: str) -> "Manifest":
        if category not in self.outputs:
            raise ValueError(
                "Unknown output category {0!r}; expected one of {1}".format(
                    category, list(self.outputs)
                )
            )
        self.outputs[category].append(path)
        return self

    def set_metric(self, key: str, value: Any) -> "Manifest":
        if key not in self.metrics:
            raise ValueError(
                "Unknown metric {0!r}; expected one of {1}".format(key, list(self.metrics))
            )
        self.metrics[key] = value
        return self

    def add_note(self, note: str) -> "Manifest":
        self.notes.append(note)
        return self

    # --- serialisation ---------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        data = {
            "pipeline": self.pipeline,
            "status": self.status,
            "inputs": dict(self.inputs),
            "parameters": dict(self.parameters),
            "steps": [dict(s) for s in self.steps],
            "outputs": {k: list(v) for k, v in self.outputs.items()},
            "metrics": dict(self.metrics),
            "notes": list(self.notes),
        }
        if self.config is not None:
            data["config"] = self.config.to_dict()
        return data

    def to_json(self, **kwargs: Any) -> str:
        kwargs.setdefault("indent", 2)
        kwargs.setdefault("sort_keys", True)
        return json.dumps(self.to_dict(), **kwargs)

    def write(self, path: str) -> str:
        with open(path, "w") as handle:
            handle.write(self.to_json())
        return path


def build_plan_manifest(
    config: JobConfig,
    box_padding_nm: Optional[float] = None,
    steps=None,
) -> Manifest:
    """Build a ``planned`` manifest capturing the effective parameters of a dry run."""
    manifest = Manifest(config, status=STATUS_PLANNED)
    if box_padding_nm is not None:
        manifest.set_parameter("boxPaddingNm", box_padding_nm)
    manifest.set_parameter("genionConcentrationM", config.genion_concentration())
    for item in (steps or []):
        if isinstance(item, tuple):
            name, command = item
            manifest.add_step(name, command=command)
        else:
            manifest.add_step(str(item))
    return manifest

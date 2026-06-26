"""Typed OpenMM job configuration (JSON/YAML).

The Phase 2 job contract. Pure standard library (no pydantic), 3.6-compatible,
so it validates in a bare interpreter and in tests. Reuses
:class:`automd_saxs.config.SystemType` for protein vs protein-ligand.

Fields cover the UI targets from the project brief (job name, inputs, system
type, force field, simulation length, repeats, ionic concentration, pH,
disulfide, box padding, with advanced MD settings defaulted) plus reproducibility
(seed) and platform selection.
"""

import json
from enum import Enum
from typing import Any, Dict, Optional

from ..config import SystemType

_FEMTOSECONDS_PER_NS = 1_000_000.0


class OpenMMForceField(Enum):
    """Supported OpenMM force-field families (explicit solvent)."""

    AMBER14 = "amber14"
    CHARMM36 = "charmm36"

    @classmethod
    def from_str(cls, value: str) -> "OpenMMForceField":
        norm = value.strip().lower()
        for member in cls:
            if member.value == norm:
                return member
        raise ValueError(
            "Unknown force field {0!r}; expected one of {1}".format(
                value, [m.value for m in cls]))


class WaterModel(Enum):
    TIP3P = "tip3p"
    TIP3PFB = "tip3pfb"
    SPCE = "spce"

    @classmethod
    def from_str(cls, value: str) -> "WaterModel":
        norm = value.strip().lower()
        for member in cls:
            if member.value == norm:
                return member
        raise ValueError(
            "Unknown water model {0!r}; expected one of {1}".format(
                value, [m.value for m in cls]))


# OpenMM force-field XML bundles. The water XML is chosen by (ff, water).
_FF_PROTEIN_XML = {
    OpenMMForceField.AMBER14: "amber14-all.xml",
    OpenMMForceField.CHARMM36: "charmm36.xml",
}
_FF_WATER_XML = {
    (OpenMMForceField.AMBER14, WaterModel.TIP3P): "amber14/tip3p.xml",
    (OpenMMForceField.AMBER14, WaterModel.TIP3PFB): "amber14/tip3pfb.xml",
    (OpenMMForceField.AMBER14, WaterModel.SPCE): "amber14/spce.xml",
    (OpenMMForceField.CHARMM36, WaterModel.TIP3P): "charmm36/water.xml",
    (OpenMMForceField.CHARMM36, WaterModel.SPCE): "charmm36/spce.xml",
}

_FIELDS = (
    "job_name", "pdb", "saxs", "system", "force_field", "water_model",
    "simulation_time_ns", "timestep_fs", "n_repeats", "temperature_K",
    "ionic_concentration_M", "ph", "disulfide", "box_padding_nm",
    "equilibration_ns", "minimize_max_iterations", "nonbonded_cutoff_nm",
    "friction_per_ps", "report_interval_steps", "frame_stride", "seed",
    "platform", "extra",
)


class OpenMMConfig:
    """Validated OpenMM job description (validates on construction)."""

    def __init__(
        self,
        pdb,
        job_name="automd_saxs_openmm",
        saxs=None,
        system=SystemType.PROTEIN,
        force_field=OpenMMForceField.AMBER14,
        water_model=WaterModel.TIP3P,
        simulation_time_ns=100,
        timestep_fs=2.0,
        n_repeats=3,
        temperature_K=300.0,
        ionic_concentration_M=0.15,
        ph=7.0,
        disulfide=False,
        box_padding_nm=None,         # None -> derive from model Dmax at runtime
        equilibration_ns=0.2,
        minimize_max_iterations=0,   # 0 -> OpenMM runs until converged
        nonbonded_cutoff_nm=1.0,
        friction_per_ps=1.0,
        report_interval_steps=5000,
        frame_stride=2,
        seed=None,
        platform=None,               # None/"auto" -> fastest available
        extra=None,
    ):
        self.pdb = pdb
        self.job_name = job_name
        self.saxs = saxs
        self.system = system
        self.force_field = force_field
        self.water_model = water_model
        self.simulation_time_ns = simulation_time_ns
        self.timestep_fs = timestep_fs
        self.n_repeats = n_repeats
        self.temperature_K = temperature_K
        self.ionic_concentration_M = ionic_concentration_M
        self.ph = ph
        self.disulfide = disulfide
        self.box_padding_nm = box_padding_nm
        self.equilibration_ns = equilibration_ns
        self.minimize_max_iterations = minimize_max_iterations
        self.nonbonded_cutoff_nm = nonbonded_cutoff_nm
        self.friction_per_ps = friction_per_ps
        self.report_interval_steps = report_interval_steps
        self.frame_stride = frame_stride
        self.seed = seed
        self.platform = platform
        self.extra = {} if extra is None else dict(extra)
        self.validate()

    def __repr__(self):
        return "OpenMMConfig(job_name={0!r}, pdb={1!r}, system={2})".format(
            self.job_name, self.pdb, self.system)

    # ------------------------------------------------------------------ #
    def validate(self) -> "OpenMMConfig":
        if not self.pdb or not str(self.pdb).lower().endswith(".pdb"):
            raise ValueError("pdb must be a .pdb path, got {0!r}".format(self.pdb))
        if self.saxs is not None and not str(self.saxs).lower().endswith(".dat"):
            raise ValueError("saxs must be a .dat path or None, got {0!r}".format(self.saxs))
        if not isinstance(self.system, SystemType):
            raise TypeError("system must be a SystemType")
        if not isinstance(self.force_field, OpenMMForceField):
            raise TypeError("force_field must be an OpenMMForceField")
        if not isinstance(self.water_model, WaterModel):
            raise TypeError("water_model must be a WaterModel")
        if (self.force_field, self.water_model) not in _FF_WATER_XML:
            raise ValueError("unsupported force-field/water combination: {0}/{1}".format(
                self.force_field.value, self.water_model.value))
        if not self.simulation_time_ns > 0 or int(self.simulation_time_ns) != self.simulation_time_ns:
            raise ValueError("simulation_time_ns must be a positive integer")
        if not self.timestep_fs > 0:
            raise ValueError("timestep_fs must be > 0")
        if int(self.n_repeats) != self.n_repeats or self.n_repeats < 1:
            raise ValueError("n_repeats must be a positive integer")
        if not self.temperature_K > 0:
            raise ValueError("temperature_K must be > 0")
        if not self.ionic_concentration_M >= 0:
            raise ValueError("ionic_concentration_M must be >= 0")
        if not 0 < self.ph < 14:
            raise ValueError("ph must be between 0 and 14")
        if self.box_padding_nm is not None and not self.box_padding_nm > 0:
            raise ValueError("box_padding_nm must be > 0 or None")
        if not self.equilibration_ns >= 0:
            raise ValueError("equilibration_ns must be >= 0")
        if self.frame_stride < 1:
            raise ValueError("frame_stride must be >= 1")
        return self

    # ------------------------------------------------------------------ #
    @property
    def uses_saxs(self) -> bool:
        return self.saxs is not None

    def forcefield_files(self):
        """OpenMM ``ForceField(...)`` XML arguments for this ff + water model."""
        return [
            _FF_PROTEIN_XML[self.force_field],
            _FF_WATER_XML[(self.force_field, self.water_model)],
        ]

    def _steps(self, ns) -> int:
        raw = ns * _FEMTOSECONDS_PER_NS / self.timestep_fs
        rounded = round(raw)
        if abs(raw - rounded) > 1e-6:
            raise ValueError(
                "{0} ns is not an integer number of {1} fs steps".format(ns, self.timestep_fs))
        return int(rounded)

    def production_steps(self) -> int:
        return self._steps(self.simulation_time_ns)

    def equilibration_steps(self) -> int:
        return self._steps(self.equilibration_ns)

    # ------------------------------------------------------------------ #
    def to_dict(self) -> Dict[str, Any]:
        data = {name: getattr(self, name) for name in _FIELDS}
        data["system"] = self.system.value
        data["force_field"] = self.force_field.value
        data["water_model"] = self.water_model.value
        data["extra"] = dict(self.extra)
        return data

    def to_json(self, **kwargs: Any) -> str:
        kwargs.setdefault("indent", 2)
        kwargs.setdefault("sort_keys", True)
        return json.dumps(self.to_dict(), **kwargs)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "OpenMMConfig":
        data = dict(data)
        if "system" in data and not isinstance(data["system"], SystemType):
            data["system"] = SystemType.from_str(str(data["system"]))
        if "force_field" in data and not isinstance(data["force_field"], OpenMMForceField):
            data["force_field"] = OpenMMForceField.from_str(str(data["force_field"]))
        if "water_model" in data and not isinstance(data["water_model"], WaterModel):
            data["water_model"] = WaterModel.from_str(str(data["water_model"]))
        known = set(_FIELDS)
        kwargs = {k: v for k, v in data.items() if k in known}
        unknown = {k: v for k, v in data.items() if k not in known}
        if unknown:
            kwargs.setdefault("extra", {}).update(unknown)
        return cls(**kwargs)

    @classmethod
    def from_json(cls, text: str) -> "OpenMMConfig":
        return cls.from_dict(json.loads(text))

    @classmethod
    def from_file(cls, path: str) -> "OpenMMConfig":
        """Load from a ``.json`` or ``.yaml``/``.yml`` file (YAML needs PyYAML)."""
        with open(path, "r") as handle:
            text = handle.read()
        if path.lower().endswith((".yaml", ".yml")):
            try:
                import yaml
            except ImportError:
                raise ValueError("YAML config requires PyYAML; use JSON or install pyyaml")
            return cls.from_dict(yaml.safe_load(text))
        return cls.from_json(text)

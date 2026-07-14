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
    "equilibration_ns", "equilibration_nvt_fraction",
    "equilibration_restraint_k",
    "minimize_max_iterations", "nonbonded_cutoff_nm",
    "friction_per_ps", "report_interval_steps", "frame_stride",
    "frame_interval_ns", "seed", "platform", "hmr", "use_propka",
    "ligand_resnames", "ligand_smiles", "keep_crystallisation_agents",
    "keep_ions", "ion_resnames", "keep_waters", "protonation_overrides", "extra",
)


class OpenMMConfig:
    """Validated OpenMM job description (validates on construction).

    Relationship to BilboMD's OpenMM pipeline (auto/pdb/crd): that pipeline does
    *implicit-solvent, Rg-restrained conformational sampling* (heat to ~600 K, NVE
    Verlet, RadiusOfGyration restraint) to generate an ensemble. AutoMD-SAXS is
    the distinct *higher-accuracy, explicit-solvent, unbiased equilibrium*
    refinement path, so several settings intentionally differ:

    * explicit solvent (PME + water box + ions) and NPT, vs implicit/NVE;
    * temperature 300 K (physiological equilibrium) vs BilboMD's 600 K sampling;
    * friction 1 ps^-1 (standard equilibrium Langevin) vs BilboMD's 0.1 ps^-1;
    * no Rg restraint (unbiased).

    Peripheral engine conventions are kept consistent with BilboMD where it does
    not affect the science: HBonds constraints, and a 1.2 nm nonbonded cutoff.
    """

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
        # Total equilibration time, split NVT then NPT. 0.5 ns is a safer default
        # than 0.2 ns: unrestrained equilibration needs time for solvent density
        # and side chains to relax off the minimised structure before production.
        equilibration_ns=0.5,
        # Fraction of equilibration run as NVT (constant volume) before switching
        # to NPT (MonteCarloBarostat): NVT relaxes temperature/velocities, then
        # NPT relaxes the box density.
        equilibration_nvt_fraction=0.4,
        # Harmonic position-restraint force constant (kJ/mol/nm^2) applied to
        # SOLUTE heavy atoms (protein/nucleic + ligand; water, ions and H left
        # free) during equilibration, so solvent relaxes around a held structure.
        # Released for production. 0 -> no restraints. ~1000 kJ/mol/nm^2 ~= 2.4
        # kcal/mol/A^2 is a moderate hold.
        equilibration_restraint_k=1000.0,
        minimize_max_iterations=0,   # 0 -> OpenMM runs until converged
        # 1.0 nm: standard PME real-space cutoff; ~1.2-1.4x faster than 1.2 nm
        # (direct-space cost ~ cutoff^3) at equivalent accuracy with PME.
        nonbonded_cutoff_nm=1.0,
        friction_per_ps=1.0,
        report_interval_steps=5000,
        frame_stride=2,
        # Extract one frame per this many ns of simulation (default 0.5 ns), so a
        # 3 ns repeat yields ~6 frames in the trajectory / analysis plots / viewer.
        # Set to None/0 to fall back to the raw frame_stride.
        frame_interval_ns=0.5,
        seed=None,
        platform=None,               # None/"auto" -> fastest available
        hmr=False,                   # Hydrogen Mass Repartitioning (4 fs, ~2x)
        use_propka=True,             # structure-based pH protonation (propka)
        # Protein-ligand: explicit ligand residue names (None -> auto-detect any
        # non-standard, non-water, non-ion, non-crystallisation-agent HETATM);
        # optional SMILES hints {resname: smiles} for robust bond perception.
        ligand_resnames=None,
        ligand_smiles=None,
        keep_crystallisation_agents=False,
        # Keep bound/structural ions (Ca2+, Zn2+, Pb2+, Fe, Mg, ...) in the
        # simulation, parameterised by amber14 (default True). Set False to strip.
        keep_ions=True,
        # Optional subset of ion residue names (PDB names, e.g. ["ZN"]) to keep
        # when keep_ions is True; other ion species are stripped. None/empty ->
        # keep ALL ions (so a structure with Ca2+ and Zn2+ can keep one, not the
        # other). Ignored when keep_ions is False.
        ion_resnames=None,
        # Keep crystallographic (structural) waters from the input in the system;
        # bulk explicit solvent is still added around them at solvation (default
        # False -> crystal waters stripped, all water re-added as bulk solvent).
        keep_waters=False,
        # Manual protonation overrides from the review UI, keyed "chain:resSeq:resname"
        # -> amber14 variant (ASH/GLH/HID/HIE/HIP/LYN) or "default". Takes
        # precedence over the propka-predicted state for that residue.
        protonation_overrides=None,
        extra=None,
    ):
        self.pdb = pdb
        self.job_name = job_name
        self.saxs = saxs
        self.system = system
        self.force_field = force_field
        self.water_model = water_model
        self.simulation_time_ns = simulation_time_ns
        # Hydrogen Mass Repartitioning enables a larger stable timestep. When the
        # caller leaves the default 2 fs, HMR bumps it to 4 fs (~2x throughput);
        # an explicit non-default timestep is respected.
        self.hmr = bool(hmr)
        self.timestep_fs = 4.0 if (self.hmr and timestep_fs == 2.0) else timestep_fs
        self.n_repeats = n_repeats
        self.temperature_K = temperature_K
        self.ionic_concentration_M = ionic_concentration_M
        self.ph = ph
        self.disulfide = disulfide
        self.box_padding_nm = box_padding_nm
        self.equilibration_ns = equilibration_ns
        self.equilibration_nvt_fraction = min(1.0, max(0.0, equilibration_nvt_fraction))
        self.equilibration_restraint_k = max(0.0, float(equilibration_restraint_k))
        self.minimize_max_iterations = minimize_max_iterations
        self.nonbonded_cutoff_nm = nonbonded_cutoff_nm
        self.friction_per_ps = friction_per_ps
        self.report_interval_steps = report_interval_steps
        self.frame_stride = frame_stride
        self.frame_interval_ns = frame_interval_ns
        self.seed = seed
        self.platform = platform
        self.use_propka = bool(use_propka)
        self.ligand_resnames = list(ligand_resnames) if ligand_resnames else None
        self.ligand_smiles = dict(ligand_smiles) if ligand_smiles else None
        self.keep_crystallisation_agents = bool(keep_crystallisation_agents)
        self.keep_ions = bool(keep_ions)
        self.ion_resnames = ([str(r).strip().upper() for r in ion_resnames]
                             if ion_resnames else None)
        self.keep_waters = bool(keep_waters)
        self.protonation_overrides = dict(protonation_overrides) if protonation_overrides else None
        # Runtime-only paths (set by the workflow after prepare; not serialised):
        # the ligand SDF (perceived bonds/charges) and the GAFF AM1-BCC charge
        # cache, so createSystem/addSolvent can parameterise ligands.
        self.ligand_sdf = None
        self.gaff_cache = None
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
        # Positive number of ns. Fractional values are allowed (e.g. 0.5 ns
        # refinements, sub-ns smoke tests); production_steps() enforces that the
        # value yields a whole number of integration steps for the timestep.
        if not self.simulation_time_ns > 0:
            raise ValueError("simulation_time_ns must be > 0")
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

    def equilibration_nvt_steps(self) -> int:
        return int(self.equilibration_steps() * self.equilibration_nvt_fraction)

    def equilibration_npt_steps(self) -> int:
        return self.equilibration_steps() - self.equilibration_nvt_steps()

    def effective_frame_stride(self) -> int:
        """DCD-frame stride for extraction, derived from ``frame_interval_ns``.

        When ``frame_interval_ns`` is set (default 0.5 ns), frames are sampled at
        that simulation-time spacing regardless of timestep: a 3 ns repeat yields
        ~6 frames. Falls back to the explicit ``frame_stride`` when
        ``frame_interval_ns`` is unset/zero.
        """
        if self.frame_interval_ns and self.frame_interval_ns > 0:
            per_dcd_ns = self.report_interval_steps * self.timestep_fs / 1.0e6
            if per_dcd_ns > 0:
                return max(1, int(round(self.frame_interval_ns / per_dcd_ns)))
        return self.frame_stride

    def frame_interval_effective_ns(self) -> float:
        """Actual ns between extracted frames (for labelling plot x-axes)."""
        return self.effective_frame_stride() * self.report_interval_steps \
            * self.timestep_fs / 1.0e6

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

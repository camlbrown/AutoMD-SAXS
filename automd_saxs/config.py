"""Typed job configuration for the AutoMD-SAXS workflow.

This replaces the legacy ``configurations.txt`` -- a shell-sourceable file that
mixed environment state, derived directory paths, and (dangerously) executable
bash function definitions -- with a validated, serialisable Python model that
captures the *scientific and execution choices* a job actually needs.

Design notes
------------
* Pure standard library (``dataclasses`` + ``enum``). No pydantic/numpy
  dependency, so this loads and validates in a bare interpreter and in tests.
* Models *choices*, not the sprawling derived directory layout the old script
  wrote (``MINIM1_DIR``, ``REPEAT_DIR1`` ...). Those are reproducible from a job
  directory + protein name and belong in a future ``paths.py``.
* Reading an old ``configurations.txt`` is delegated to ``automd_saxs.legacy``
  so existing jobs keep loading.

Bug fixed at the model level
----------------------------
The legacy ``run_MD.sh`` collected an ionic concentration from the user but then
hardcoded ``genion ... -conc 0.15``, silently ignoring the choice. Here the
chosen concentration is a first-class field (:attr:`JobConfig.ionic_concentration_M`)
and :meth:`JobConfig.genion_concentration` is what command planning must use.
"""

import json
from enum import Enum
from typing import Any, Dict, Optional


# --------------------------------------------------------------------------- #
# Enumerations
# --------------------------------------------------------------------------- #

class SystemType(Enum):
    """System composition. Matches the legacy ``SYSTEM`` values."""

    PROTEIN = "Protein"
    PROTEIN_LIGAND = "Protein-ligand"

    @classmethod
    def from_str(cls, value: str) -> "SystemType":
        norm = value.strip().lower().replace("_", "-")
        for member in cls:
            if member.value.lower() == norm:
                return member
        raise ValueError(
            "Unknown system type {0!r}; expected one of {1}".format(
                value, [m.value for m in cls]
            )
        )


class ForceField(Enum):
    """Supported force fields (must match a directory in ``ff_files/``)."""

    AMBER14SB = "amber14sb"
    CHARMM36M = "charmm36m"

    @classmethod
    def from_str(cls, value: str) -> "ForceField":
        # Accept a bare name ("amber14sb") or a path whose basename is the name
        # (the legacy ``FORCE_FIELD`` was an absolute path ending in the name).
        norm = value.strip().replace("\\", "/").rstrip("/").split("/")[-1].lower()
        for member in cls:
            if member.value == norm:
                return member
        raise ValueError(
            "Unknown force field {0!r}; expected one of {1}".format(
                value, [m.value for m in cls]
            )
        )

    @property
    def water_model(self) -> str:
        """Default water model the legacy ``pdb2gmx`` used for this force field.

        The legacy scripts used ``tip3p`` on the protein-only path and ``spce``
        on the protein-ligand path. The water model therefore depends on the
        system as well; this property gives the force-field default and
        :meth:`JobConfig.water_model` applies the system-specific override.
        """
        return "tip3p"


class BoxShape(Enum):
    """Solvation box shape, mapped to ``gmx editconf -bt`` names.

    ``RECTANGULAR`` is the label the legacy front-end offered; GROMACS calls it
    ``triclinic`` so they share a value. ``AUTO`` defers the choice to a
    radius-of-gyration heuristic at planning time (legacy ``BOX_SHAPE=auto``).
    """

    DODECAHEDRON = "dodecahedron"
    TRICLINIC = "triclinic"
    OCTAHEDRON = "octahedron"
    CUBIC = "cubic"
    AUTO = "auto"

    @classmethod
    def from_str(cls, value: str) -> "BoxShape":
        norm = value.strip().lower()
        aliases = {"rectangular": "triclinic"}
        norm = aliases.get(norm, norm)
        for member in cls:
            if member.value == norm:
                return member
        raise ValueError(
            "Unknown box shape {0!r}; expected one of {1} (or 'rectangular')".format(
                value, [m.value for m in cls]
            )
        )


_FEMTOSECONDS_PER_NS = 1_000_000.0
_DMAX_MODEL_SENTINEL = "Model"


# --------------------------------------------------------------------------- #
# Job configuration
# --------------------------------------------------------------------------- #

#: Field names in construction order. Used for serialisation and for filtering
#: unknown keys in :meth:`JobConfig.from_dict`. Kept explicit (rather than via
#: ``dataclasses``) so this module runs on Python 3.6, which lacks that stdlib
#: module, as well as on the project's 3.9 conda environment.
_FIELDS = (
    "protein_file",
    "system",
    "force_field",
    "saxs_file",
    "box_shape",
    "ionic_concentration_M",
    "simulation_time_ns",
    "timestep_fs",
    "n_repeats",
    "disulfide",
    "dmax_nm",
    "python_cmd",
    "gmx_module",
    "email",
    "partition",
    "extra",
)


class JobConfig:
    """Validated description of a single AutoMD-SAXS job.

    Required scientific inputs have no defaults; everything with a sensible
    default is optional. :meth:`validate` runs automatically on construction.

    Fields
    ------
    protein_file : str (``.pdb``)
    system : SystemType (default PROTEIN)
    force_field : ForceField (default AMBER14SB)
    saxs_file : Optional[str] (``.dat``; ``None`` -> no SAXS)
    box_shape : BoxShape (default DODECAHEDRON)
    ionic_concentration_M : float > 0 (default 0.15 -- physiological)
    simulation_time_ns : positive int (default 100)
    timestep_fs : float > 0 (default 2.0)
    n_repeats : positive int (default 3)
    disulfide : bool
    dmax_nm : Optional[float] > 0 (``None`` -> use model Dmax, legacy "Model")
    python_cmd, gmx_module, email, partition : execution/environment settings
    extra : dict for forward-compatible / provenance keys
    """

    def __init__(
        self,
        protein_file,
        system=SystemType.PROTEIN,
        force_field=ForceField.AMBER14SB,
        saxs_file=None,
        box_shape=BoxShape.DODECAHEDRON,
        ionic_concentration_M=0.15,
        simulation_time_ns=100,
        timestep_fs=2.0,
        n_repeats=3,
        disulfide=False,
        dmax_nm=None,
        python_cmd="python",
        gmx_module=None,
        email=None,
        partition=None,
        extra=None,
    ):
        self.protein_file = protein_file
        self.system = system
        self.force_field = force_field
        self.saxs_file = saxs_file
        self.box_shape = box_shape
        self.ionic_concentration_M = ionic_concentration_M
        self.simulation_time_ns = simulation_time_ns
        self.timestep_fs = timestep_fs
        self.n_repeats = n_repeats
        self.disulfide = disulfide
        self.dmax_nm = dmax_nm
        self.python_cmd = python_cmd
        self.gmx_module = gmx_module
        self.email = email
        self.partition = partition
        self.extra = {} if extra is None else dict(extra)
        self.validate()

    def __repr__(self):
        return "JobConfig(protein_file={0!r}, system={1}, force_field={2})".format(
            self.protein_file, self.system, self.force_field
        )

    # ------------------------------------------------------------------ #
    # Validation
    # ------------------------------------------------------------------ #
    def validate(self) -> "JobConfig":
        if not self.protein_file or not self.protein_file.lower().endswith(".pdb"):
            raise ValueError(
                "protein_file must be a .pdb path, got {0!r}".format(self.protein_file)
            )
        if self.saxs_file is not None and not self.saxs_file.lower().endswith(".dat"):
            raise ValueError(
                "saxs_file must be a .dat path or None, got {0!r}".format(self.saxs_file)
            )
        if not isinstance(self.system, SystemType):
            raise TypeError("system must be a SystemType")
        if not isinstance(self.force_field, ForceField):
            raise TypeError("force_field must be a ForceField")
        if not isinstance(self.box_shape, BoxShape):
            raise TypeError("box_shape must be a BoxShape")

        # Legacy invariant: the protein-ligand path always used amber14sb
        # (acpype/gaff2 parameters are spliced into an AMBER topology).
        if (
            self.system is SystemType.PROTEIN_LIGAND
            and self.force_field is not ForceField.AMBER14SB
        ):
            raise ValueError(
                "Protein-ligand systems require the amber14sb force field "
                "(got {0}).".format(self.force_field.value)
            )

        if not self.ionic_concentration_M > 0:
            raise ValueError(
                "ionic_concentration_M must be > 0, got {0}".format(
                    self.ionic_concentration_M
                )
            )
        if int(self.simulation_time_ns) != self.simulation_time_ns or self.simulation_time_ns <= 0:
            raise ValueError(
                "simulation_time_ns must be a positive integer, got {0}".format(
                    self.simulation_time_ns
                )
            )
        if not self.timestep_fs > 0:
            raise ValueError(
                "timestep_fs must be > 0, got {0}".format(self.timestep_fs)
            )
        if int(self.n_repeats) != self.n_repeats or self.n_repeats < 1:
            raise ValueError(
                "n_repeats must be a positive integer, got {0}".format(self.n_repeats)
            )
        if self.dmax_nm is not None and not self.dmax_nm > 0:
            raise ValueError(
                "dmax_nm must be > 0 or None, got {0}".format(self.dmax_nm)
            )
        return self

    # ------------------------------------------------------------------ #
    # Derived scientific quantities
    # ------------------------------------------------------------------ #
    @property
    def uses_saxs(self) -> bool:
        return self.saxs_file is not None

    @property
    def force_field_name(self) -> str:
        return self.force_field.value

    def water_model(self) -> str:
        """Water model fed to ``pdb2gmx``.

        Preserves the legacy split: the protein-ligand path used ``spce`` while
        the protein-only path used ``tip3p``.
        """
        if self.system is SystemType.PROTEIN_LIGAND:
            return "spce"
        return self.force_field.water_model

    def number_of_steps(self) -> int:
        """Production ``nsteps`` from simulation time and timestep.

        Generalises the legacy ``SIMULATION_TIME * 500000`` (which assumed a
        2 fs step): ``steps = ns * 1e6 / timestep_fs``. Raises if the result is
        not a whole number of steps for the chosen timestep.
        """
        raw = self.simulation_time_ns * _FEMTOSECONDS_PER_NS / self.timestep_fs
        rounded = round(raw)
        if abs(raw - rounded) > 1e-6:
            raise ValueError(
                "simulation_time_ns={0} is not an integer number of {1} fs steps".format(
                    self.simulation_time_ns, self.timestep_fs
                )
            )
        return int(rounded)

    def genion_concentration(self) -> float:
        """Ionic concentration (M) to pass to ``gmx genion -conc``.

        This is the field the legacy ``run_MD.sh`` *should* have used; it
        hardcoded ``-conc 0.15`` instead. Command planning must read this.
        """
        return self.ionic_concentration_M

    def partition_flag(self) -> str:
        """sbatch partition flag, or empty string (matches legacy ``$PARTITION``)."""
        return "--partition={0}".format(self.partition) if self.partition else ""

    # ------------------------------------------------------------------ #
    # Serialisation
    # ------------------------------------------------------------------ #
    def to_dict(self) -> Dict[str, Any]:
        data = {name: getattr(self, name) for name in _FIELDS}
        data["system"] = self.system.value
        data["force_field"] = self.force_field.value
        data["box_shape"] = self.box_shape.value
        data["extra"] = dict(self.extra)
        return data

    def to_json(self, **kwargs: Any) -> str:
        kwargs.setdefault("indent", 2)
        kwargs.setdefault("sort_keys", True)
        return json.dumps(self.to_dict(), **kwargs)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "JobConfig":
        data = dict(data)  # shallow copy; do not mutate caller's dict
        if "system" in data and not isinstance(data["system"], SystemType):
            data["system"] = SystemType.from_str(str(data["system"]))
        if "force_field" in data and not isinstance(data["force_field"], ForceField):
            data["force_field"] = ForceField.from_str(str(data["force_field"]))
        if "box_shape" in data and not isinstance(data["box_shape"], BoxShape):
            data["box_shape"] = BoxShape.from_str(str(data["box_shape"]))

        known = set(_FIELDS)
        kwargs = {k: v for k, v in data.items() if k in known}
        unknown = {k: v for k, v in data.items() if k not in known}
        if unknown:
            kwargs.setdefault("extra", {}).update(unknown)
        return cls(**kwargs)

    @classmethod
    def from_json(cls, text: str) -> "JobConfig":
        return cls.from_dict(json.loads(text))

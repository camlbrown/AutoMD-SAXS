"""Compatibility layer for the legacy ``configurations.txt``.

The old ``simulation_setup.sh`` wrote a shell-sourceable ``configurations.txt``
containing three kinds of content:

1. shell *function definitions* (``load_gmx() { module load ...; }``) -- a
   code-execution hazard when ``source``-d;
2. environment/path variables (``PYTHON_CMD=...``, ``MINIM1_DIR=...``);
3. the scientific choices (``SYSTEM=``, ``FORCE_FIELD=``, ``IONIC_CONCENTRATION=``,
   ``SIMULATION_TIME=``, ``DMAX=``, ...).

This module parses that file *as data* (no ``source``, no code execution) and
maps the scientific choices onto a typed :class:`~automd_saxs.config.JobConfig`,
so existing jobs keep loading after the refactor.
"""

import re
from typing import Dict, Optional

from .config import BoxShape, ForceField, JobConfig, SystemType

# Matches ``KEY=VALUE`` at the start of a line. Shell function definitions such
# as ``load_gmx()   { ... }`` do NOT match (the ``()`` breaks the identifier),
# so they are skipped without us needing to understand them.
_ASSIGNMENT = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)=(.*)$")

# Legacy "no value" tokens.
_NONE_TOKENS = {"", "none"}


def _strip_quotes(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        return value[1:-1]
    return value


def parse_configurations(path: str) -> Dict[str, str]:
    """Read a legacy ``configurations.txt`` into a ``{KEY: value}`` dict.

    Comments, blank lines, and shell function definitions are ignored. Values
    have surrounding quotes stripped. The file is never executed.
    """
    result: Dict[str, str] = {}
    with open(path, "r") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            match = _ASSIGNMENT.match(stripped)
            if not match:
                continue  # function bodies, braces, etc.
            key, value = match.group(1), match.group(2)
            result[key] = _strip_quotes(value)
    return result


def _opt(value: Optional[str]) -> Optional[str]:
    """Map legacy 'None'/'' tokens to a real ``None``."""
    if value is None:
        return None
    if value.strip().lower() in _NONE_TOKENS:
        return None
    return value.strip()


def _disulfide_to_bool(value: Optional[str]) -> bool:
    return (value or "").strip().lower() in ("y", "yes")


def _partition_name(value: Optional[str]) -> Optional[str]:
    """Extract the bare partition name from a legacy ``PARTITION`` value.

    The old script stored either an empty string or the full sbatch flag
    ``--partition=<name>``.
    """
    value = _opt(value)
    if value is None:
        return None
    if value.startswith("--partition="):
        value = value[len("--partition="):]
    return value or None


def config_from_legacy_dict(values: Dict[str, str]) -> JobConfig:
    """Build a :class:`JobConfig` from parsed ``configurations.txt`` values."""
    protein_file = values.get("PROTEIN_FILE")
    if not protein_file:
        raise ValueError("configurations.txt is missing PROTEIN_FILE")

    system = SystemType.from_str(values["SYSTEM"]) if values.get("SYSTEM") else SystemType.PROTEIN
    force_field = (
        ForceField.from_str(values["FORCE_FIELD"])
        if values.get("FORCE_FIELD")
        else ForceField.AMBER14SB
    )
    box_shape = (
        BoxShape.from_str(values["BOX_SHAPE"]) if values.get("BOX_SHAPE") else BoxShape.DODECAHEDRON
    )

    dmax_raw = _opt(values.get("DMAX"))
    if dmax_raw is None or dmax_raw.lower() == "model":
        dmax_nm = None
    else:
        dmax_nm = float(dmax_raw)

    kwargs = dict(
        protein_file=protein_file,
        system=system,
        force_field=force_field,
        saxs_file=_opt(values.get("SAXS_FILE")),
        box_shape=box_shape,
        disulfide=_disulfide_to_bool(values.get("DISULFIDE")),
        dmax_nm=dmax_nm,
        python_cmd=_opt(values.get("PYTHON_CMD")) or "python",
        gmx_module=_opt(values.get("GMX_MODULE")),
        email=_opt(values.get("EMAIL_ADDR")),
        partition=_partition_name(values.get("PARTITION")),
    )
    if values.get("IONIC_CONCENTRATION"):
        kwargs["ionic_concentration_M"] = float(values["IONIC_CONCENTRATION"])
    if values.get("SIMULATION_TIME"):
        kwargs["simulation_time_ns"] = int(values["SIMULATION_TIME"])

    return JobConfig(**kwargs)


def config_from_legacy_file(path: str) -> JobConfig:
    """Parse a legacy ``configurations.txt`` and return a :class:`JobConfig`."""
    return config_from_legacy_dict(parse_configurations(path))

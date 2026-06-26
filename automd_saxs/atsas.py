"""ATSAS SAXS analysis: command planning + output parsing.

The legacy ``saxs_analysis.sh`` shelled out to ATSAS binaries (``shanum``,
``crysol``, ``gajoe``) and then post-processed their text output with several
ad-hoc Python scripts. This module:

* provides **argv builders** for the ATSAS tools (no execution); and
* ports the output **parsers** to tested functions with fixtures, so chi^2/Rg
  ranking, fit curves, and ensemble Rg distributions can be read without ATSAS.

ATSAS is *not* required to use the parsers. In Phase 2 the OpenMM branch will
replace these tools with FoXS/MultiFoXS; the parser/metric shapes here are kept
deliberately tool-agnostic so that swap is localised.
"""

import re
from typing import Dict, List, Optional, Tuple

_STRUCTURE_FRAME = re.compile(r"structure_(\d+)")


# --------------------------------------------------------------------------- #
# Command builders (no execution)
# --------------------------------------------------------------------------- #

def shanum_command(dat_file: str) -> List[str]:
    """``shanum`` estimates the useful Smax from experimental data."""
    return ["shanum", dat_file]


def crysol_command(
    dat_file: str,
    structures: str = "structure_*.pdb",
    smax: Optional[float] = None,
    lm: int = 30,
    constant_subtraction: bool = True,
) -> List[str]:
    """Theoretical SAXS profiles / fit (legacy ``crysol <dat> structure* -lm 30 -cst -sm <smax>``)."""
    argv = ["crysol", dat_file, structures, "-lm", str(lm)]
    if constant_subtraction:
        argv.append("-cst")
    if smax is not None:
        argv += ["-sm", ("%g" % smax)]
    return argv


def gajoe_command(profiles_file: str) -> List[str]:
    """Ensemble (genetic algorithm) selection over a profile list (legacy ``gajoe -p``).

    Only the input-list flag is modelled here; the legacy script also passed
    numerous GA tuning parameters whose exact flag spellings are version
    dependent, so they are intentionally left to the caller/config rather than
    guessed.
    """
    return ["gajoe", "-p", profiles_file]


# --------------------------------------------------------------------------- #
# Parsers
# --------------------------------------------------------------------------- #

def parse_crysol_summary(text: str) -> List[Dict[str, object]]:
    """Parse ``crysol_summary.txt`` into per-model records.

    Faithful to ``chi2vsRg.py``: whitespace columns where column 1 is the model
    filename (``structure_<N>.pdb``), column 3 is Rg, column 7 is chi^2. Lines
    without enough columns are skipped. Returns dicts with keys
    ``model``, ``frame`` (int or None), ``rg``, ``chi2``.
    """
    records = []  # type: List[Dict[str, object]]
    for line in text.splitlines():
        parts = line.split()
        if len(parts) < 8:
            continue
        try:
            rg = float(parts[3])
            chi2 = float(parts[7])
        except ValueError:
            continue
        model = parts[1]
        match = _STRUCTURE_FRAME.search(model)
        frame = int(match.group(1)) if match else None
        records.append({"model": model, "frame": frame, "rg": rg, "chi2": chi2})
    return records


def parse_fit(text: str) -> List[Tuple[float, ...]]:
    """Parse a CRYSOL/GAJOE ``.fit`` file (skip 1 header line; numeric columns).

    Faithful to ``ensemble_scattering.py`` (``genfromtxt(skip_header=1)``):
    typically ``q, I(q), error, ensemble_fit``. Returns a list of float tuples;
    non-numeric / short lines are skipped.
    """
    rows = []  # type: List[Tuple[float, ...]]
    for line in text.splitlines()[1:]:
        parts = line.split()
        if not parts:
            continue
        try:
            rows.append(tuple(float(p) for p in parts))
        except ValueError:
            continue
    return rows


def parse_rg_distribution(text: str) -> List[Tuple[float, float, float]]:
    """Parse ``Rg_distr_*.txt`` (skip 5 header lines; ``Rg, pool_freq, sel_freq``).

    Faithful to ``plot_rg_ensemble.py`` (``skiprows=5``).
    """
    rows = []  # type: List[Tuple[float, float, float]]
    for line in text.splitlines()[5:]:
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
        except ValueError:
            continue
    return rows


def parse_xvg(text: str) -> Tuple[Dict[str, str], List[Tuple[float, ...]]]:
    """Parse a GROMACS ``.xvg`` file into ``(metadata, data_rows)``.

    ``@`` (grace) and ``#`` (comment) lines are metadata; ``@ title "..."`` style
    entries are collected into the metadata dict. Remaining lines are numeric rows.
    """
    metadata = {}  # type: Dict[str, str]
    rows = []  # type: List[Tuple[float, ...]]
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            continue
        if stripped.startswith("@"):
            m = re.match(r'@\s+(\w+)\s+"([^"]*)"', stripped)
            if m:
                metadata[m.group(1)] = m.group(2)
            continue
        parts = stripped.split()
        try:
            rows.append(tuple(float(p) for p in parts))
        except ValueError:
            continue
    return metadata, rows


# --------------------------------------------------------------------------- #
# Derived metrics
# --------------------------------------------------------------------------- #

def best_fit(records: List[Dict[str, object]]) -> Optional[Dict[str, object]]:
    """Record with the lowest chi^2, or ``None`` if empty."""
    if not records:
        return None
    return min(records, key=lambda r: r["chi2"])


def summarize_fits(records: List[Dict[str, object]]) -> Dict[str, object]:
    """Manifest-ready metrics from crysol records: bestChi2/bestFrame/rgMean."""
    if not records:
        return {"bestChi2": None, "bestFrame": None, "rgMean": None}
    best = best_fit(records)
    rgs = [r["rg"] for r in records]
    return {
        "bestChi2": best["chi2"],
        "bestFrame": best["frame"],
        "rgMean": sum(rgs) / len(rgs),
    }

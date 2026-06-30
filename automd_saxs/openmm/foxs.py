"""FoXS / MultiFoXS: command planning + output parsing (no execution).

The OpenMM-branch replacement for ATSAS ``crysol``/``gajoe``. Provides argv
builders for FoXS (per-structure theoretical profile / experimental fit) and
MultiFoXS (ensemble selection), plus parsers for their ``.fit``/``.dat`` outputs
and a chi^2/Rg summary that feeds the manifest.

Invocation mirrors the BilboMD worker scripts:

* ``foxs -p <model.pdb>``                       -> ``<model.pdb>.dat`` (profile)
* ``foxs <exp.dat> <model.pdb>``                -> fit with chi^2 in the header
* ``multi_foxs <exp.dat> <profiles...>``        -> ensemble models

Parsers need no FoXS binary; they read files a real run produced and are tested
against fixtures.
"""

import re
from typing import Dict, List, Optional, Tuple

FOXS = "foxs"
MULTIFOXS = "multi_foxs"

# chi^2 appears in a FoXS header/comment line, e.g. "Chi^2 = 1.23" or "chi2=1.23".
_CHI2 = re.compile(r"chi\^?2\s*=\s*([0-9.eE+-]+)", re.IGNORECASE)
_C1 = re.compile(r"c1\s*=\s*([0-9.eE+-]+)", re.IGNORECASE)
_C2 = re.compile(r"c2\s*=\s*([0-9.eE+-]+)", re.IGNORECASE)
_MODEL_FRAME = re.compile(r"(?:structure|md|frame)[_-]?(\d+)", re.IGNORECASE)


# --------------------------------------------------------------------------- #
# Command builders
# --------------------------------------------------------------------------- #

def foxs_profile_command(pdb: str) -> List[str]:
    """Compute a theoretical SAXS profile for one structure (``foxs -p``)."""
    return [FOXS, "-p", pdb]


def foxs_fit_command(experimental_dat: str, pdb: str, extra=None) -> List[str]:
    """Fit a structure to experimental data (``foxs <exp.dat> <model.pdb>``)."""
    argv = [FOXS, experimental_dat, pdb]
    if extra:
        argv.extend(extra)
    return argv


def multifoxs_command(experimental_dat: str, profiles, output: Optional[str] = None) -> List[str]:
    """Ensemble selection over profile ``.dat`` files (``multi_foxs <exp> <profiles>``)."""
    argv = [MULTIFOXS, experimental_dat]
    if output:
        argv += ["-o", output]
    argv.extend(list(profiles))
    return argv


# --------------------------------------------------------------------------- #
# Parsers
# --------------------------------------------------------------------------- #

def parse_foxs_fit(text: str) -> Dict[str, object]:
    """Parse a FoXS ``.fit``/``.dat`` file.

    Returns ``{chi2, c1, c2, data}`` where ``data`` is a list of float tuples
    (typically ``q, I_exp, error, I_model``). chi^2/c1/c2 are read from comment
    lines if present (``None`` otherwise). Comment lines start with ``#``.
    """
    chi2 = c1 = c2 = None
    data = []  # type: List[Tuple[float, ...]]
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#") or _CHI2.search(stripped) and not stripped[0].isdigit():
            m = _CHI2.search(stripped)
            if m and chi2 is None:
                chi2 = float(m.group(1))
            mc1 = _C1.search(stripped)
            if mc1 and c1 is None:
                c1 = float(mc1.group(1))
            mc2 = _C2.search(stripped)
            if mc2 and c2 is None:
                c2 = float(mc2.group(1))
            continue
        parts = stripped.split()
        try:
            data.append(tuple(float(p) for p in parts))
        except ValueError:
            continue
    return {"chi2": chi2, "c1": c1, "c2": c2, "data": data}


def parse_foxs_log(text: str) -> Dict[str, Optional[float]]:
    """Extract ``c1``/``c2`` (and chi^2 if present) from a FoXS log."""
    out = {"chi2": None, "c1": None, "c2": None}  # type: Dict[str, Optional[float]]
    for pat, key in ((_CHI2, "chi2"), (_C1, "c1"), (_C2, "c2")):
        m = pat.search(text)
        if m:
            out[key] = float(m.group(1))
    return out


def frame_index(model_name: str):
    """Extract a frame number from a model filename, or ``None``."""
    m = _MODEL_FRAME.search(model_name)
    return int(m.group(1)) if m else None


# --------------------------------------------------------------------------- #
# Summary
# --------------------------------------------------------------------------- #

def run_foxs_fits(frame_pdbs, experimental_dat, runner, cwd=None):
    """Fit each frame to experimental data via the runner; parse chi^2 (real runs).

    In dry-run the commands are recorded and an empty record list is returned. In
    a real run, after each ``foxs`` call the produced fit file (``<frame>.dat`` /
    ``<frame>.fit`` in ``cwd``) is parsed for chi^2. Returns per-frame records
    ``{model, frame, chi2, fitFile}``.
    """
    import glob
    import os

    records = []
    for pdb in frame_pdbs:
        argv = foxs_fit_command(experimental_dat, pdb)
        runner.run(argv, label="foxs:{0}".format(os.path.basename(pdb)), cwd=cwd)
        if runner.dry_run:
            continue
        base = os.path.basename(pdb)
        # FoXS writes its outputs next to the input PDB regardless of cwd, naming
        # the experimental fit "<frame_stem>_<exp_stem>.fit" whose header carries
        # Chi^2. The "<base>.pdb.dat" is only the theoretical profile (no chi^2),
        # so match the .fit, not the .dat.
        out_dir = os.path.dirname(pdb) or "."
        stem = os.path.splitext(base)[0]  # e.g. structure_0
        fits = sorted(glob.glob(os.path.join(out_dir, stem + "*.fit")))
        chi2 = None
        fit_file = fits[0] if fits else None
        if fit_file:
            with open(fit_file) as fh:
                chi2 = parse_foxs_fit(fh.read()).get("chi2")
        records.append({"model": base, "frame": frame_index(base),
                        "chi2": chi2, "fitFile": fit_file})
    return records


def run_multifoxs(profiles, experimental_dat, runner, output, cwd=None):
    """Run MultiFoXS ensemble selection over profile files via the runner."""
    import os

    argv = multifoxs_command(experimental_dat, profiles, output=output)
    runner.run(argv, label="multifoxs", cwd=cwd)
    return output


def best_fit(records: List[Dict[str, object]]) -> Optional[Dict[str, object]]:
    """Record with the lowest chi^2 (records with chi2 None are ignored)."""
    scored = [r for r in records if r.get("chi2") is not None]
    if not scored:
        return None
    return min(scored, key=lambda r: r["chi2"])


def summarize_fits(records: List[Dict[str, object]]) -> Dict[str, object]:
    """Manifest-ready metrics: bestChi2 / bestFrame / rgMean (rg if present)."""
    if not records:
        return {"bestChi2": None, "bestFrame": None, "rgMean": None}
    best = best_fit(records)
    rgs = [r["rg"] for r in records if r.get("rg") is not None]
    return {
        "bestChi2": best["chi2"] if best else None,
        "bestFrame": best.get("frame") if best else None,
        "rgMean": (sum(rgs) / len(rgs)) if rgs else None,
    }

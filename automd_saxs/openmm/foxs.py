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


def multifoxs_command(experimental_dat: str, profiles, num_states: int = 5) -> List[str]:
    """Ensemble selection over structures (``multi_foxs -s N <exp> <inputs...>``).

    ``multi_foxs`` writes its outputs (``ensembles_size_<k>.txt`` and
    ``multi_state_model_<k>_1_1.dat``) to the CURRENT WORKING DIRECTORY, so the
    caller must run it in the desired output dir. It has NO output-file flag
    (``-o`` is ``--offset``). ``inputs`` may be PDB files (partial profiles
    computed internally, matching the BilboMD/Carbonara invocation) or profile
    ``.dat`` files. ``-s`` caps the maximal ensemble size.
    """
    return [MULTIFOXS, "-s", str(num_states), experimental_dat] + list(profiles)


def parse_multifoxs_ensembles(text: str) -> Optional[Dict[str, object]]:
    """Parse a ``ensembles_size_<k>.txt`` file into the BEST ensemble.

    Format (port of the Carbonara ``parseMultiFoxsEnsembles``): ensemble header
    lines start with the rank digit in column 0 and carry chi^2 as the 2nd
    ``|``-field; indented species lines carry ``<weight> | <profile/pdb>``.
    Returns ``{chi2, members:[{pdb, weight}]}`` for the first (best) ensemble, or
    ``None``.
    """
    best = None
    cur = None
    for raw in text.split("\n"):
        if not raw.strip():
            continue
        if raw[:1].isdigit():
            if best is not None:
                break  # only the first (best) ensemble is needed
            parts = raw.split("|")
            chi2 = None
            if len(parts) > 1:
                try:
                    chi2 = float(parts[1].strip())
                except ValueError:
                    chi2 = None
            cur = {"chi2": chi2, "members": []}
            best = cur
        elif cur is not None:
            parts = raw.split("|")
            if len(parts) >= 3:
                try:
                    weight = float(parts[1].strip().split()[0])
                except (ValueError, IndexError):
                    continue
                tokens = parts[2].strip().split()
                if tokens:
                    cur["members"].append({"pdb": tokens[0], "weight": weight})
    if best and best["members"]:
        return best
    return None


def parse_multifoxs_fit(text: str) -> List[Dict[str, float]]:
    """Parse a ``multi_state_model_*.dat``/``.fit`` ensemble fit curve.

    Returns rows ``{q, exp, error, model}`` (port of Carbonara ``parseMultiFoxsFit``).
    """
    rows = []
    for raw in text.split("\n"):
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        p = s.split()
        if len(p) < 4:
            continue
        try:
            q, exp, error, model = float(p[0]), float(p[1]), float(p[2]), float(p[3])
        except ValueError:
            continue
        rows.append({"q": q, "exp": exp, "error": error, "model": model})
    return rows


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


def run_multifoxs(inputs, experimental_dat, runner, out_dir, num_states=5,
                  frame_of=None):
    """Run MultiFoXS ensemble selection and parse the result.

    ``inputs`` are the structures/profiles to combine (PDB frames, matching the
    Carbonara/BilboMD invocation). multi_foxs runs in ``out_dir`` (its outputs
    land in cwd) and, for each ensemble size ``k`` from 1..num_states, writes
    ``ensembles_size_<k>.txt`` + ``multi_state_model_<k>_1_1.dat``. Returns a
    structured summary::

        {"ensembles": [{"size", "chi2", "members": [{"frame", "weight"}]}],
         "best": {"size", "chi2", "curve": [{"q","exp","error","model"}]}}

    ``frame_of`` maps a member filename back to a frame index (defaults to the
    ``structure_<n>`` pattern). Returns ``None`` in dry-run or if nothing parsed.
    """
    import glob
    import os

    argv = multifoxs_command(experimental_dat, inputs, num_states=num_states)
    runner.run(argv, label="multifoxs", cwd=out_dir)
    if runner.dry_run:
        return None

    if frame_of is None:
        def frame_of(name):
            return frame_index(os.path.basename(name))

    def _read_fit_curve(size):
        for ext in (".dat", ".fit"):
            fit_path = os.path.join(
                out_dir, "multi_state_model_{0}_1_1{1}".format(size, ext))
            if os.path.exists(fit_path):
                with open(fit_path) as fh:
                    curve = parse_multifoxs_fit(fh.read())
                if curve:
                    return curve
        return []

    ensembles = []
    for ens_path in sorted(glob.glob(os.path.join(out_dir, "ensembles_size_*.txt"))):
        m = re.search(r"ensembles_size_(\d+)\.txt$", ens_path)
        if not m:
            continue
        size = int(m.group(1))
        with open(ens_path) as fh:
            parsed = parse_multifoxs_ensembles(fh.read())
        if not parsed:
            continue
        members = [{"frame": frame_of(mem["pdb"]), "weight": mem["weight"]}
                   for mem in parsed["members"]]
        ensembles.append({"size": size, "chi2": parsed["chi2"],
                          "members": members, "curve": _read_fit_curve(size)})

    if not ensembles:
        return None

    ensembles.sort(key=lambda e: e["size"])
    # Best = lowest chi^2 across sizes (guard None chi2).
    scored = [e for e in ensembles if isinstance(e.get("chi2"), (int, float))]
    best_ens = min(scored, key=lambda e: e["chi2"]) if scored else ensembles[0]
    best = {"size": best_ens["size"], "chi2": best_ens.get("chi2"),
            "curve": best_ens.get("curve", [])}
    return {"ensembles": ensembles, "best": best}


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

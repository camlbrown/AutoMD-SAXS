"""Structure-based protonation at a target pH (propka + OpenMM variants).

OpenMM's ``Modeller.addHydrogens(pH=...)`` only meaningfully varies HIS (and
detects disulfides) with pH -- it does NOT protonate ASP/GLU at low pH or
deprotonate LYS at high pH. So on its own the user's pH barely affects
protonation. This module closes that gap: it runs `propka3` to predict
structure-specific pKa values for each ionizable group, then compares each to the
requested pH to choose the protonation state, applied via ``addHydrogens``'s
per-residue ``variants`` argument.

Only the variants amber14 actually supports are applied (probed):
    ASH, GLH  (protonated Asp/Glu, below their pKa)
    HIP       (protonated His, below its pKa; the neutral HID/HIE tautomer is
               left to OpenMM to pick from the local H-bond environment)
    LYN       (neutral Lys, above its pKa)
amber14's addHydrogens has no deprotonated-CYS/TYR variant, so those stay
protonated (rare in the pH 3-11 range anyway). Disulfides are detected by OpenMM.

Everything is best-effort: if propka is unavailable or fails, callers fall back
to plain ``addHydrogens(pH=...)`` and note it -- protonation is advisory and a
structure should always be checked visually (see the pipeline disclaimer).
"""

import os
import string
import subprocess

# amber14 addHydrogens variants we can safely force (empirically supported).
SUPPORTED_VARIANTS = {"ASH", "GLH", "HID", "HIE", "HIP", "LYN", "CYX"}

# Ionizable residues we titrate (those with a usable amber14 variant).
TITRATABLE = ("ASP", "GLU", "HIS", "LYS")


def _propka_bin():
    return os.environ.get("PROPKA_BIN", "propka3")


def chain_labels():
    """A, B, ... Z, then AA-style is not needed (26 chains is plenty)."""
    return list(string.ascii_uppercase)


def relabel_chains(topology):
    """Give every chain a distinct A/B/C... id so propka's residue identifiers
    line up with the topology (PDBFixer leaves blank chain ids, which propka and
    OpenMM then represent inconsistently)."""
    for chain, cid in zip(topology.chains(), chain_labels()):
        chain.id = cid


def parse_propka(text):
    """Parse a .pka file's SUMMARY into ``{(chain, resSeq, resName): pKa}``."""
    pkas = {}
    in_summary = False
    for line in text.splitlines():
        if "SUMMARY OF THIS PREDICTION" in line:
            in_summary = True
            continue
        if not in_summary:
            continue
        s = line.strip()
        if s.startswith("Group"):
            continue
        if not s or s.startswith("-"):
            break
        parts = s.split()
        if len(parts) >= 4:
            name, num, chain = parts[0], parts[1], parts[2]
            try:
                pka = float(parts[3])
            except ValueError:
                continue
            pkas[(chain, num, name)] = pka
    return pkas


def run_propka(pdb_path, work_dir, timeout=600):
    """Run ``propka3`` on a heavy-atom PDB; return the parsed pKa map.

    Returns ``{}`` on any failure so the caller can fall back to plain pH
    protonation. ``pdb_path`` should live in (or be reachable from) ``work_dir``;
    propka writes ``<stem>.pka`` beside its input in the cwd.
    """
    try:
        result = subprocess.run(
            [_propka_bin(), os.path.basename(pdb_path)],
            cwd=work_dir, capture_output=True, text=True, timeout=timeout, check=False)
        stem = os.path.splitext(os.path.basename(pdb_path))[0]
        pka_path = os.path.join(work_dir, stem + ".pka")
        if not os.path.exists(pka_path):
            return {}
        with open(pka_path) as fh:
            return parse_propka(fh.read())
    except Exception:  # noqa: BLE001 - protonation prediction is advisory
        return {}


def protonation_variant(resname, pka, ph):
    """OpenMM variant to force for ``resname`` at ``ph`` given predicted ``pka``,
    or ``None`` to keep OpenMM's default. Acids protonate below their pKa; Lys/His
    (bases) are protonated (charged) below their pKa."""
    if resname == "ASP":
        var = "ASH" if ph < pka else None
    elif resname == "GLU":
        var = "GLH" if ph < pka else None
    elif resname == "HIS":
        var = "HIP" if ph < pka else None  # neutral tautomer left to OpenMM
    elif resname == "LYS":
        var = None if ph < pka else "LYN"
    else:
        var = None
    return var if var in SUPPORTED_VARIANTS else None


# Selectable protonation variants per residue type for the review-UI dropdowns
# ("default" = let OpenMM choose from pH / H-bond environment).
VARIANT_CHOICES = {
    "ASP": ["default", "ASH"],       # deprotonated (default) / protonated
    "GLU": ["default", "GLH"],
    "HIS": ["default", "HID", "HIE", "HIP"],
    "LYS": ["default", "LYN"],       # protonated (default) / neutral
}


def build_variants(topology, pkas, ph, overrides=None):
    """Return ``(variants, audit, table)`` for ``Modeller.addHydrogens``.

    ``variants`` is a list of length ``n_residues`` (``None`` = OpenMM default,
    else a variant name). ``overrides`` (from the review UI) is keyed
    ``"chain:resSeq:resname"`` -> amber14 variant (or ``"default"``) and takes
    precedence over the propka-predicted state. ``audit`` records the residues
    whose state differs from the default; ``table`` records EVERY ionizable
    residue (chain/resSeq/residue/pKa/state/source/choices) for the review-page
    editor. Chain ids on ``topology`` must match the propka run (call
    :func:`relabel_chains` first).
    """
    overrides = overrides or {}
    variants = [None] * topology.getNumResidues()
    audit = []
    table = []
    for i, res in enumerate(topology.residues()):
        if res.name not in TITRATABLE:
            continue
        okey = "{0}:{1}:{2}".format(res.chain.id, res.id, res.name)
        pka = pkas.get((res.chain.id, res.id, res.name))
        source = "propka"
        state = "default"
        if okey in overrides:
            raw = overrides[okey]
            var = None if raw in (None, "default", "") else raw
            if var is None or var in SUPPORTED_VARIANTS:
                variants[i] = var
                source = "override"
                state = var or "default"
                audit.append({
                    "chain": res.chain.id, "resSeq": res.id, "residue": res.name,
                    "pKa": round(pka, 2) if pka is not None else None,
                    "variant": var or "default", "source": "override",
                })
        elif pka is not None:
            var = protonation_variant(res.name, pka, ph)
            if var:
                variants[i] = var
                state = var
                audit.append({
                    "chain": res.chain.id, "resSeq": res.id, "residue": res.name,
                    "pKa": round(pka, 2), "variant": var, "source": "propka",
                })
        table.append({
            "key": okey,
            "chain": res.chain.id, "resSeq": res.id, "residue": res.name,
            "pKa": round(pka, 2) if pka is not None else None,
            "state": state, "source": source,
            "choices": VARIANT_CHOICES.get(res.name, ["default"]),
        })
    return variants, audit, table

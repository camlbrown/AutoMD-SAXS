"""Model Dmax and box-padding calculations.

This module is a faithful port of the scientific logic that the legacy workflow
buried in shell. Two pieces are reproduced here so they can be unit-tested
without GROMACS, ATSAS, or Slurm:

1. **Model Dmax** -- the maximum interatomic distance of the input structure,
   used as a fallback when no experimental SAXS Dmax is supplied and as the
   reference for box padding. Ported from the inline ``python3`` heredoc in
   ``run_MD.sh`` (the block computing ``MODEL_DMAX`` via
   ``scipy.spatial.distance.pdist``).

2. **Box padding** -- the solvation buffer fed to ``gmx editconf -d``. Ported
   from the ``bc`` expression in ``run_MD.sh``::

       BOX_PADDING = sqrt((DMAX - MODEL_DMAX)^2) + 2     # == |Δ| + 2 nm
       if BOX_PADDING < 2: BOX_PADDING = 2               # 2 nm minimum

Legacy behaviour deliberately preserved
---------------------------------------
* Only ``ATOM`` records contribute to model Dmax. ``HETATM`` records (ligands,
  ions, crystallographic waters) are ignored, matching the legacy
  ``line.startswith("ATOM")`` test.
* Coordinates are read from the fixed PDB columns 31-54 (0-based slices
  ``[30:38] [38:46] [46:54]``), not whitespace-split, matching the legacy parser.
* Distances are computed in Angstrom (PDB units) then divided by 10 to report
  nanometres, matching the legacy ``/10.0``.
* A structure with fewer than two parseable atoms yields a model Dmax of 0.0.

``numpy`` is used when importable (fast, vectorised) and is otherwise replaced
by an equivalent pure-Python computation, so this module has no hard runtime
dependency and can be tested in a bare interpreter.
"""

import math
from typing import List, Optional, Sequence, Tuple

Coordinate = Tuple[float, float, float]

# Sentinel used by the legacy ``simulation_setup.sh`` when no experimental SAXS
# Dmax is provided: ``DMAX=Model``. ``run_MD.sh`` then substitutes the computed
# model Dmax. We accept it (case-insensitively) so old config files still load.
MODEL_SENTINEL = "Model"

ANGSTROM_PER_NM = 10.0

DEFAULT_BUFFER_NM = 2.0
DEFAULT_MIN_PADDING_NM = 2.0


def parse_atom_coordinates(pdb_path: str) -> List[Coordinate]:
    """Return ``(x, y, z)`` triples (Angstrom) for every ``ATOM`` record.

    Mirrors the legacy parser exactly: only lines beginning with ``ATOM`` are
    considered, coordinates are sliced from fixed columns, and any line whose
    coordinate field cannot be parsed as a float is silently skipped (the legacy
    code wrapped the conversion in ``try/except ValueError: continue``).
    """
    coords: List[Coordinate] = []
    with open(pdb_path, "r") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            coords.append((x, y, z))
    return coords


def _max_pairwise_distance(coords: Sequence[Coordinate]) -> float:
    """Largest Euclidean distance between any two points (same units as input).

    Returns ``0.0`` for fewer than two points, matching the legacy guard
    ``np.max(pdist(coords)) if coords.shape[0] > 1 else 0.0``.
    """
    n = len(coords)
    if n < 2:
        return 0.0

    try:
        import numpy as np  # type: ignore
    except ImportError:
        return _max_pairwise_distance_py(coords)

    arr = np.asarray(coords, dtype=float)
    max_sq = 0.0
    # Row-by-row to keep memory O(n) rather than materialising the full O(n^2)
    # distance matrix that ``pdist`` would build for large structures.
    for i in range(n - 1):
        diffs = arr[i + 1:] - arr[i]
        sq = (diffs * diffs).sum(axis=1)
        local = float(sq.max())
        if local > max_sq:
            max_sq = local
    return math.sqrt(max_sq)


def _max_pairwise_distance_py(coords: Sequence[Coordinate]) -> float:
    """Pure-Python fallback for :func:`_max_pairwise_distance`."""
    n = len(coords)
    max_sq = 0.0
    for i in range(n - 1):
        xi, yi, zi = coords[i]
        for j in range(i + 1, n):
            xj, yj, zj = coords[j]
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            d_sq = dx * dx + dy * dy + dz * dz
            if d_sq > max_sq:
                max_sq = d_sq
    return math.sqrt(max_sq)


def model_dmax_nm(pdb_path: str) -> float:
    """Maximum interatomic distance of ``ATOM`` records, in nanometres.

    Faithful port of the ``MODEL_DMAX`` heredoc in ``run_MD.sh``.
    """
    coords = parse_atom_coordinates(pdb_path)
    return _max_pairwise_distance(coords) / ANGSTROM_PER_NM


def resolve_dmax_nm(
    experimental_dmax: Optional[object],
    model_dmax: float,
) -> float:
    """Resolve the effective Dmax used for box padding.

    ``experimental_dmax`` may be:

    * ``None`` or the ``"Model"`` sentinel -> use ``model_dmax`` (the no-SAXS
      path, where ``simulation_setup.sh`` wrote ``DMAX=Model``);
    * a number or numeric string -> use that experimental value.

    Raises ``ValueError`` for a non-empty, non-sentinel value that is not a
    valid positive number, replacing the legacy ``exit 1`` regex guards with an
    actionable Python error.
    """
    if experimental_dmax is None:
        return model_dmax
    if isinstance(experimental_dmax, str):
        text = experimental_dmax.strip()
        if text == "" or text.lower() == MODEL_SENTINEL.lower():
            return model_dmax
        try:
            value = float(text)
        except ValueError:
            raise ValueError(
                "DMAX is not a valid number: {0!r}".format(experimental_dmax)
            )
    else:
        value = float(experimental_dmax)

    if value <= 0:
        raise ValueError(
            "DMAX must be greater than zero, got {0}".format(value)
        )
    return value


def box_padding_nm(
    experimental_dmax_nm: float,
    model_dmax_nm: float,
    buffer_nm: float = DEFAULT_BUFFER_NM,
    minimum_nm: float = DEFAULT_MIN_PADDING_NM,
) -> float:
    """Solvation box padding (``gmx editconf -d``) in nanometres.

    ``padding = |experimental - model| + buffer``, floored at ``minimum``.
    Faithful port of the ``bc`` block in ``run_MD.sh`` (default buffer 2 nm,
    default floor 2 nm).
    """
    padding = abs(experimental_dmax_nm - model_dmax_nm) + buffer_nm
    if padding < minimum_nm:
        padding = minimum_nm
    return padding

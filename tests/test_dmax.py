"""Unit tests for ``automd_saxs.dmax``.

These tests require no GROMACS, ATSAS, Slurm, numpy, scipy, or pytest. They run
under pytest (``pytest tests/test_dmax.py``) and also standalone on a bare
interpreter (``python tests/test_dmax.py``) via the ``__main__`` runner at the
bottom -- useful because the project's scientific conda env is not present on
every development machine.
"""

import os
import sys
import tempfile

# Make the package importable when run as a plain script from the repo root.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs import dmax  # noqa: E402

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "tiny.pdb")


def _almost(a, b, tol=1e-6):
    assert abs(a - b) <= tol, "expected {0}, got {1}".format(b, a)


# --- parse_atom_coordinates -------------------------------------------------

def test_parses_only_atom_records():
    coords = dmax.parse_atom_coordinates(FIXTURE)
    # 3 ATOM records; the HETATM ZN at (1000,0,0) must be ignored.
    assert coords == [(0.0, 0.0, 0.0), (30.0, 0.0, 0.0), (0.0, 40.0, 0.0)]


def test_parser_skips_unparseable_coordinate_lines():
    # An ATOM line whose coordinate columns are not floats must be skipped,
    # matching the legacy try/except ValueError: continue.
    line = "ATOM      5  CA  ALA A   5     XXX.XXX   0.000   0.000  1.00  0.00           C\n"
    good = "ATOM      6  CA  ALA A   6       0.000   0.000   0.000  1.00  0.00           C\n"
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as fh:
        fh.write(line)
        fh.write(good)
        path = fh.name
    try:
        coords = dmax.parse_atom_coordinates(path)
        assert coords == [(0.0, 0.0, 0.0)]
    finally:
        os.unlink(path)


# --- model_dmax_nm ----------------------------------------------------------

def test_model_dmax_known_geometry():
    # max pairwise distance is B-C = 50 Angstrom = 5.000 nm.
    _almost(dmax.model_dmax_nm(FIXTURE), 5.0)


def test_model_dmax_single_atom_is_zero():
    one = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\n"
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as fh:
        fh.write(one)
        path = fh.name
    try:
        _almost(dmax.model_dmax_nm(path), 0.0)
    finally:
        os.unlink(path)


def test_max_pairwise_pure_python_matches_default():
    coords = [(0.0, 0.0, 0.0), (30.0, 0.0, 0.0), (0.0, 40.0, 0.0)]
    _almost(dmax._max_pairwise_distance_py(coords), 50.0)
    _almost(dmax._max_pairwise_distance(coords), 50.0)


# --- resolve_dmax_nm --------------------------------------------------------

def test_resolve_dmax_model_sentinel_uses_model():
    _almost(dmax.resolve_dmax_nm("Model", 4.2), 4.2)
    _almost(dmax.resolve_dmax_nm("model", 4.2), 4.2)
    _almost(dmax.resolve_dmax_nm(None, 4.2), 4.2)
    _almost(dmax.resolve_dmax_nm("", 4.2), 4.2)


def test_resolve_dmax_numeric_string_and_number():
    _almost(dmax.resolve_dmax_nm("7.5", 4.2), 7.5)
    _almost(dmax.resolve_dmax_nm(7.5, 4.2), 7.5)


def test_resolve_dmax_rejects_garbage():
    for bad in ("abc", "-1", "0"):
        try:
            dmax.resolve_dmax_nm(bad, 4.2)
        except ValueError:
            continue
        raise AssertionError("expected ValueError for {0!r}".format(bad))


# --- box_padding_nm ---------------------------------------------------------

def test_box_padding_difference_plus_buffer():
    # |8 - 5| + 2 = 5 nm
    _almost(dmax.box_padding_nm(8.0, 5.0), 5.0)


def test_box_padding_is_symmetric_in_difference():
    _almost(dmax.box_padding_nm(5.0, 8.0), 5.0)


def test_box_padding_floored_at_minimum():
    # equal Dmax -> |0| + 2 = 2, exactly the 2 nm floor
    _almost(dmax.box_padding_nm(5.0, 5.0), 2.0)


def test_box_padding_custom_buffer_can_hit_floor():
    # tiny buffer + tiny difference would drop below the floor -> clamped to 2 nm
    _almost(dmax.box_padding_nm(5.0, 5.05, buffer_nm=0.1), 2.0)


def _run_standalone():
    """Discover and run every ``test_*`` function in this module."""
    funcs = sorted(
        (name, obj)
        for name, obj in globals().items()
        if name.startswith("test_") and callable(obj)
    )
    failures = []
    for name, fn in funcs:
        try:
            fn()
            print("PASS", name)
        except Exception as exc:  # noqa: BLE001 - report and continue
            failures.append((name, exc))
            print("FAIL", name, "->", exc)
    print("\n{0}/{1} passed".format(len(funcs) - len(failures), len(funcs)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_run_standalone())

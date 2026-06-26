"""Unit tests for ``automd_saxs.ligand``. No external deps; runs standalone."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs import ligand  # noqa: E402


def _atom(serial, name, res, chain, resnum, x=0.0, y=0.0, z=0.0, record="ATOM"):
    fmt = "%-6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n"
    return fmt % (record, serial, name, res, chain, resnum, x, y, z)


COMPLEX = (
    _atom(1, " N  ", "ALA", "A", 1) +
    _atom(2, " CA ", "ALA", "A", 1) +
    _atom(3, " C1 ", "LIG", "B", 500, record="HETATM") +
    _atom(4, " C2 ", "LIG", "B", 500, record="HETATM") +
    _atom(5, " O1 ", "DRG", "B", 501, record="HETATM")
)


# --- structure splitting ----------------------------------------------------

def test_split_complex_separates_protein_and_ligands():
    protein, ligands = ligand.split_complex(COMPLEX)
    assert "ALA" in protein
    assert "LIG" not in protein and "DRG" not in protein
    names = [(rn, num) for rn, num, _ in ligands]
    assert names == [("LIG", 500), ("DRG", 501)]
    # ligand groups carry their own atoms
    lig_lines = ligands[0][2]
    assert len(lig_lines) == 2


def test_atom_to_hetatm_relabels_only_nonstandard():
    out = ligand.atom_to_hetatm(_atom(1, " CA ", "ALA", "A", 1) +
                                _atom(2, " C1 ", "LIG", "A", 9))
    lines = out.splitlines()
    assert lines[0].startswith("ATOM")     # ALA kept
    assert lines[1].startswith("HETATM")   # LIG relabelled


def test_set_residue_number_and_chain():
    line = _atom(1, " C1 ", "LIG", "B", 500)
    assert ligand.set_residue_number(line, 7).splitlines()[0][22:26] == "   7"
    assert ligand.set_chain_identifier(line, "X").splitlines()[0][21] == "X"


def test_first_residue_number_and_chain():
    assert ligand.first_residue_number(COMPLEX) == 1
    assert ligand.first_chain_identifier(COMPLEX) == "A"


# --- topology splicing ------------------------------------------------------

TOP = (
    "; topology\n"
    '#include "amber14sb.ff/forcefield.itp"\n'
    "\n"
    "[ molecules ]\n"
    "; name  count\n"
    "Protein   1\n"
)


def test_add_include_after_forcefield():
    out = ligand.add_include_after_forcefield(TOP, "LIG.itp")
    lines = out.splitlines()
    i = next(k for k, l in enumerate(lines) if "forcefield.itp" in l)
    assert lines[i + 1] == '#include "LIG.itp"'


def test_add_include_without_forcefield_raises():
    try:
        ligand.add_include_after_forcefield("[ molecules ]\nProtein 1\n", "x.itp")
    except ValueError:
        return
    raise AssertionError("expected ValueError")


def test_add_molecule_appends_entry():
    out = ligand.add_molecule(TOP, "LIG", 1)
    assert out.rstrip().endswith("LIG         1")


def test_extract_atomtypes_splits_block():
    itp = (
        "; ligand itp\n"
        "[ atomtypes ]\n"
        " ca  ca  0.0  0.0  A  1.0  1.0\n"
        " c3  c3  0.0  0.0  A  1.0  1.0\n"
        "[ moleculetype ]\n"
        "LIG  3\n"
    )
    atomtypes, remaining = ligand.extract_atomtypes(itp)
    assert "[ atomtypes ]" in atomtypes
    assert "ca  ca" in atomtypes
    assert "[ atomtypes ]" not in remaining
    assert "[ moleculetype ]" in remaining


def test_extract_atomtypes_absent_returns_original():
    itp = "[ moleculetype ]\nLIG 3\n"
    atomtypes, remaining = ligand.extract_atomtypes(itp)
    assert atomtypes == ""
    assert remaining == itp


# --- command builders / plan ------------------------------------------------

def test_acpype_uses_gaff2():
    argv = ligand.acpype_command("LIG.pdb", "LIG")
    assert argv == ["acpype", "-a", "gaff2", "-i", "LIG.pdb", "-b", "LIG"]


def test_protein_pdb2gmx_ligand_path_flags():
    argv = ligand.protein_pdb2gmx_command(disulfide=True)
    assert "-ignh" in argv
    assert argv[argv.index("-water") + 1] == "spce"
    assert "-ss" in argv


def test_plan_ligand_setup_order():
    steps = ligand.plan_ligand_setup(["ligand_500", "ligand_501"], disulfide=False)
    labels = [s[0] for s in steps]
    assert labels[0] == "protein_pdb2gmx"
    assert labels.count("acpype") == 2
    assert labels[-2:] == ["pdb_tidy", "editconf_complex"]
    # reduce step records the stdout redirect target
    reduce_step = next(s for s in steps if s[0] == "reduce")
    assert reduce_step[2] == "ligand_500_H.pdb"


def _run_standalone():
    funcs = sorted((n, o) for n, o in globals().items()
                   if n.startswith("test_") and callable(o))
    failures = []
    for name, fn in funcs:
        try:
            fn()
            print("PASS", name)
        except Exception as exc:  # noqa: BLE001
            failures.append((name, exc))
            print("FAIL", name, "->", exc)
    print("\n{0}/{1} passed".format(len(funcs) - len(failures), len(funcs)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_run_standalone())

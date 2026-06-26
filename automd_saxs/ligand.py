"""Protein-ligand setup: structure splitting, topology splicing, tool planning.

The legacy protein-ligand path in ``run_MD.sh`` interleaved external tools
(``reduce``, ``pdb4amber``, ``acpype``, ``pdb_tidy``, ``pdb2gmx``) with brittle
``sed``/``awk`` topology edits and several one-off Python scripts under
``slurms/ligand_setup/``. This module:

* ports the structure/topology manipulations to **pure, testable functions**
  (no ``sed``, no in-place file surgery); and
* provides **argv builders** + an ordered per-ligand plan for the external tools,
  so the protein-ligand pipeline can be planned without running anything.

Faithful to the legacy behaviour (same residue classification, same tool flags,
same topology include placement); the implementation is just testable.
"""

from typing import Dict, List, Tuple

# Residues treated as protein/solvent (everything else is a ligand). Union of the
# lists used by the legacy separate_multi_lig.py / ATOM_to_HETATM.py (deduped).
STANDARD_RESIDUES = frozenset({
    "THR", "ACE", "GLU", "VAL", "GLN", "LEU", "SOL", "SER", "ARG", "HOH", "CYS",
    "TYR", "ILE", "TRP", "PRO", "LYS", "ASN", "ASP", "MET", "NME", "GLY", "PHE",
    "HIS", "HIE", "HID", "HIP", "ALA", "NHE", "CTHR", "CACE", "CGLU", "CVAL",
    "CGLN", "CLEU", "CSER", "CARG", "CCYS", "CCYX", "CYX", "NCYX", "CTYR", "CILE",
    "CTRP", "CPRO", "CLYS", "CASN", "CASP", "CMET", "CGLY", "CPHE", "CHIS", "CHIE",
    "HcID", "CHIP", "CALA", "GLH", "CYM", "CCYM", "NCYM", "HYP", "CHYP", "NHYP",
    "LYN", "NLYN", "ORN", "DAB", "ASH", "NGLH", "CASH", "NASH", "NTHR", "NGLU",
    "NVAL", "NGLN", "NLEU", "NSER", "NARG", "NCYS", "NTYR", "NILE", "NH2", "URE",
    "HO4", "NTRP", "NPRO", "NLYS", "NASN", "NASP", "NMET", "NGLY", "NPHE", "NHIS",
    "NHIE", "NHID", "NHIP", "NALA", "MSE",
})

DEFAULT_GMX = "gmx_mpi"


# --------------------------------------------------------------------------- #
# Structure manipulation (pure)
# --------------------------------------------------------------------------- #

def split_complex(pdb_text: str):
    """Split a complex PDB into protein text and per-ligand records.

    Port of ``separate_multi_lig.py``: ``ATOM``/``HETATM`` records whose residue
    name is not in :data:`STANDARD_RESIDUES` are ligand atoms, grouped by
    ``(residue_name, residue_number)`` in first-seen order; everything else
    (including non-coordinate lines) is protein.

    Returns ``(protein_text, ligands)`` where ``ligands`` is a list of
    ``(residue_name, residue_number, lines)`` tuples.
    """
    protein_lines = []  # type: List[str]
    ligands = []  # type: List[Tuple[str, int, List[str]]]
    index = {}  # type: Dict[Tuple[str, int], int]

    for line in pdb_text.splitlines(keepends=True):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            resname = line[17:20].strip()
            try:
                resnum = int(line[22:26])
            except ValueError:
                protein_lines.append(line)
                continue
            if resname not in STANDARD_RESIDUES:
                key = (resname, resnum)
                if key not in index:
                    index[key] = len(ligands)
                    ligands.append((resname, resnum, []))
                ligands[index[key]][2].append(line)
            else:
                protein_lines.append(line)
        else:
            protein_lines.append(line)

    return "".join(protein_lines), ligands


def atom_to_hetatm(pdb_text: str) -> str:
    """Relabel non-standard ``ATOM`` residues as ``HETATM`` (port of ATOM_to_HETATM.py)."""
    out = []
    for line in pdb_text.splitlines(keepends=True):
        if line.startswith("ATOM"):
            resname = line[17:20].strip()
            if resname not in STANDARD_RESIDUES:
                line = line.replace("ATOM  ", "HETATM", 1)
        out.append(line)
    return "".join(out)


def set_residue_number(pdb_text: str, residue_number: int) -> str:
    """Force every coordinate record to ``residue_number`` (port of replace_residue_numbers)."""
    out = []
    rj = str(residue_number).rjust(4)
    for line in pdb_text.splitlines(keepends=True):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            line = line[:22] + rj + line[26:]
        out.append(line)
    return "".join(out)


def set_chain_identifier(pdb_text: str, chain: str) -> str:
    """Force every coordinate record to ``chain`` (port of replace_chain_identifier)."""
    out = []
    for line in pdb_text.splitlines(keepends=True):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            line = line[:21] + chain + line[22:]
        out.append(line)
    return "".join(out)


def first_residue_number(pdb_text: str):
    """Residue number of the first coordinate record, or ``None``."""
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                return int(line[22:26].strip())
            except ValueError:
                return None
    return None


def first_chain_identifier(pdb_text: str):
    """Chain identifier of the first coordinate record, or ``None``."""
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            return line[21]
    return None


# --------------------------------------------------------------------------- #
# Topology manipulation (pure) -- replaces the legacy sed/awk splicing
# --------------------------------------------------------------------------- #

_FORCEFIELD_INCLUDE = "forcefield.itp"


def add_include_after_forcefield(top_text: str, include_target: str) -> str:
    """Insert ``#include "<target>"`` immediately after the forcefield include.

    Replaces the legacy ``sed`` insert that placed ligand ``.itp`` /
    ``_atomtypes.txt`` includes after the ``#include ".../forcefield.itp"`` line.
    """
    line_to_add = '#include "{0}"'.format(include_target)
    out = []
    inserted = False
    for line in top_text.splitlines():
        out.append(line)
        if not inserted and _FORCEFIELD_INCLUDE in line and line.lstrip().startswith("#include"):
            out.append(line_to_add)
            inserted = True
    if not inserted:
        raise ValueError("no forcefield.itp #include found in topology")
    return "\n".join(out) + ("\n" if top_text.endswith("\n") else "")


def add_molecule(top_text: str, name: str, count: int = 1) -> str:
    """Append a molecule entry to the ``[ molecules ]`` section.

    Replaces the legacy ``echo "<name> 1" >> top``. The entry is appended after
    the last existing molecule line (or at end of file if the section is last).
    """
    entry = "{0}{1}{2}".format(name, " " * 9, count)
    text = top_text if top_text.endswith("\n") else top_text + "\n"
    return text + entry + "\n"


def extract_atomtypes(itp_text: str):
    """Split an ``.itp`` into its ``[ atomtypes ]`` block and the remainder.

    Port of the legacy ``sed -n '<start>,<end>p'`` / ``sed -i '<start>,<end>d'``
    that lifted the ``[ atomtypes ]`` section (up to the following
    ``[ moleculetype ]``) out of each ligand ``.itp`` so it could be included
    once near the top of the topology. Returns ``(atomtypes_block, remaining_itp)``.
    If there is no ``[ atomtypes ]`` section, returns ``("", itp_text)``.
    """
    lines = itp_text.splitlines(keepends=True)
    start = None
    end = None
    for i, line in enumerate(lines):
        stripped = line.strip().replace(" ", "")
        if start is None and stripped.startswith("[atomtypes]"):
            start = i
        elif start is not None and stripped.startswith("[moleculetype]"):
            end = i
            break
    if start is None:
        return "", itp_text
    if end is None:
        end = len(lines)
    atomtypes = "".join(lines[start:end])
    remaining = "".join(lines[:start] + lines[end:])
    return atomtypes, remaining


# --------------------------------------------------------------------------- #
# External-tool argv builders + per-ligand plan
# --------------------------------------------------------------------------- #

def reduce_command(input_pdb: str) -> List[str]:
    """``reduce`` adds hydrogens; legacy redirected stdout to ``<name>_H.pdb``."""
    return ["reduce", input_pdb]


def pdb4amber_command(input_pdb: str, output_pdb: str) -> List[str]:
    return ["pdb4amber", "-i", input_pdb, "-o", output_pdb]


def acpype_command(input_pdb: str, basename: str, atom_type: str = "gaff2") -> List[str]:
    return ["acpype", "-a", atom_type, "-i", input_pdb, "-b", basename]


def pdb_tidy_command(input_pdb: str) -> List[str]:
    """``pdb_tidy`` cleans the merged complex; legacy redirected stdout."""
    return ["pdb_tidy", input_pdb]


def protein_pdb2gmx_command(
    disulfide: bool,
    input_pdb: str = "ProteinAmber.pdb",
    output_pdb: str = "Protein_pdb2gmx.pdb",
    topology: str = "Protein.top",
    force_field: str = "amber14sb",
    water: str = "spce",
    gmx: str = DEFAULT_GMX,
) -> List[str]:
    """Protein ``pdb2gmx`` for the ligand path (legacy used ``-ignh -water spce``)."""
    argv = [
        gmx, "pdb2gmx",
        "-ff", force_field,
        "-f", input_pdb,
        "-o", output_pdb,
        "-p", topology,
        "-ter",
        "-water", water,
        "-ignh",
    ]
    if disulfide:
        argv.append("-ss")
    return argv


def plan_ligand_parameterisation(ligand_name: str, atom_type: str = "gaff2"):
    """Ordered ``(label, argv, stdout)`` steps to parameterise one ligand.

    Mirrors the legacy per-ligand loop: add H (``reduce``) -> clean
    (``pdb4amber``) -> parameterise (``acpype``). ``stdout`` names the file the
    legacy redirected ``reduce`` output to (``None`` when not redirected).
    """
    h_pdb = "{0}_H.pdb".format(ligand_name)
    clean_pdb = "{0}.pdb".format(ligand_name)
    src_pdb = "ligands/{0}.pdb".format(ligand_name)
    return [
        ("reduce", reduce_command(src_pdb), h_pdb),
        ("pdb4amber", pdb4amber_command(h_pdb, clean_pdb), None),
        ("acpype", acpype_command(clean_pdb, ligand_name, atom_type), None),
    ]


def plan_ligand_setup(ligand_names, disulfide: bool = False):
    """Full protein-ligand setup plan as ``(label, argv, stdout)`` steps.

    Order matches ``run_MD.sh``: split the complex (handled separately via
    :func:`split_complex`), build the protein topology, then per-ligand
    parameterise, then tidy/relabel and build the complex ``.gro``.
    """
    steps = [
        ("protein_pdb2gmx", protein_pdb2gmx_command(disulfide), None),
    ]
    for name in ligand_names:
        steps.extend(plan_ligand_parameterisation(name))
    steps.append(("pdb_tidy", pdb_tidy_command("Complex.pdb"), "Complex_tidy.pdb"))
    steps.append(
        ("editconf_complex",
         [DEFAULT_GMX, "editconf", "-f", "GMX.pdb", "-o", "GMX.gro"],
         None)
    )
    return steps

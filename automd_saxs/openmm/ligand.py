"""Small-molecule (ligand) handling for protein-ligand simulations.

Separates a PDB into protein / ligand(s) / water / ions, parameterises each
ligand with GAFF (openmmforcefields, AM1-BCC charges via AmberTools), and builds
a :class:`ForceField` that knows about both the amber14 protein and the GAFF
ligands. A ligand's bond orders + hydrogens come from RDKit perception on the
PROTONATED ligand (a net-charge search), or from a user-supplied SMILES
(``Molecule.from_pdb_and_smiles``) for robustness.

Ligand parameterisation is genuinely hard for arbitrary PDBs (no bond orders,
sometimes no hydrogens). Best practice, and the reliable path here: give the
ligand already protonated, or supply a SMILES hint. Results should be checked.

All heavy deps (openff-toolkit, openmmforcefields, rdkit, AmberTools) ship in the
runtime image and are imported lazily.
"""

import os

from ..command_runner import MissingDependencyError

# Residues we do NOT treat as GAFF ligands. Standard amino acids + the Task-1
# protonation variants + terminal caps, nucleic acids, and water.
_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # protonation / disulfide variants
    "ASH", "GLH", "HID", "HIE", "HIP", "LYN", "CYX", "CYM", "HYP",
    # terminal caps
    "ACE", "NME", "NMA", "NH2", "FOR",
}
_NUCLEIC = {
    "DA", "DT", "DG", "DC", "DU", "DI", "A", "U", "G", "C", "I", "T",
    "RA", "RU", "RG", "RC", "5MC", "PSU",
}
WATER = {"HOH", "WAT", "TIP", "TIP3", "TIP4", "TIP5", "SOL", "H2O", "DOD", "SPC"}

# Monatomic ions (kept for Task 3 ion support; excluded from the GAFF path).
IONS = {
    "NA", "CL", "K", "MG", "CA", "ZN", "FE", "FE2", "MN", "CU", "CU1", "CO",
    "NI", "LI", "RB", "CS", "SR", "BA", "F", "BR", "IOD", "CD", "HG", "PB",
    "AL", "AU", "AG", "PT", "PD", "CR", "MO", "V", "W", "SB", "SN", "NA+",
    "CL-", "K+", "MG2", "CA2", "ZN2", "SOD", "CLA", "POT", "MG+2", "CA+2",
}

# Common crystallisation / buffer additives stripped by default (not usually
# part of the biological assembly the user wants to simulate).
CRYSTALLISATION_AGENTS = {
    "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "P6G", "MPD", "DMS", "DMSO",
    "SO4", "PO4", "ACT", "FMT", "MES", "EPE", "TRS", "IMD", "BME", "MRD",
    "BOG", "NAG", "MAN", "BMA", "FUC", "GAL", "CIT", "TAR", "MLI", "ACY",
    "TBU", "MOH", "IPA", "ETX", "BU3", "PGO", "12P", "15P", "2PE", "XPE",
    "SIN", "SCN", "AZI", "NO3", "NH4", "PEO", "OLC", "BME", "DTT", "GSH",
}

STANDARD_RESIDUES = _AA | _NUCLEIC

# Map common PDB ion residue names to the amber14 ion template name so a bound
# ion is kept and parameterised (Task 3). Most match directly; these are aliases.
ION_ALIASES = {
    "CAL": "CA", "CA2": "CA", "CA+2": "CA", "ZN2": "ZN", "MG2": "MG",
    "MG+2": "MG", "MN2": "MN", "FE3": "FE", "NA+": "NA", "CL-": "CL", "K+": "K",
    "POT": "K", "SOD": "NA", "CLA": "CL", "IOD": "IOD", "I": "IOD",
}


def amber_ion_resname(pdb_resname):
    """amber14 ion template name for a PDB ion residue name (alias or identity)."""
    return ION_ALIASES.get(pdb_resname, pdb_resname)


_MODIFIED_AA_CACHE = None


def modified_residues():
    """Residue names PDBFixer can convert to standard amino acids (MSE->MET, ...).

    These are kept WITH the protein so PDBFixer's replaceNonstandardResidues can
    convert them -- they must not be split off as GAFF ligands. Sourced from
    PDBFixer's own substitution table (authoritative), with a common-case
    fallback if it can't be imported.
    """
    global _MODIFIED_AA_CACHE
    if _MODIFIED_AA_CACHE is None:
        try:
            from pdbfixer.pdbfixer import substitutions
            _MODIFIED_AA_CACHE = set(substitutions.keys())
        except Exception:  # noqa: BLE001 - fall back to common modified residues
            _MODIFIED_AA_CACHE = {
                "MSE", "SEP", "TPO", "PTR", "CSO", "CME", "MLY", "KCX", "PCA",
                "HYP", "FME", "CSD", "OCS", "M3L", "CAS", "CSS",
            }
    return _MODIFIED_AA_CACHE


def _require():
    try:
        from openff.toolkit import Molecule
        from openmmforcefields.generators import GAFFTemplateGenerator
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except ImportError as exc:
        raise MissingDependencyError(
            "Protein-ligand support requires openff-toolkit, openmmforcefields, "
            "rdkit and AmberTools (all in the runtime image). Original import "
            "error: {0}".format(exc))
    return Molecule, GAFFTemplateGenerator, Chem, rdDetermineBonds


def _resname(line):
    return line[17:20].strip()


def classify_residues(pdb_path, strip_agents=True, ligand_resnames=None):
    """Split a PDB's ATOM/HETATM lines into categories.

    Returns ``{"protein": [lines], "ligand": {resname: [lines]}, "water": [lines],
    "ion": {resname: [lines]}, "stripped": {resname: n}}``.

    - Standard amino/nucleic residues -> protein (kept, prepped).
    - Water -> water (stripped; re-added at solvation).
    - Monatomic ions -> ion (kept aside for Task 3; excluded from GAFF).
    - Crystallisation agents -> stripped by default (``strip_agents``).
    - Everything else (or ``ligand_resnames`` if given) -> ligand (GAFF).
    """
    want = set(ligand_resnames) if ligand_resnames else None
    modified = modified_residues()
    out = {"protein": [], "ligand": {}, "water": [], "ion": {}, "stripped": {}}
    for line in open(pdb_path):
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        name = _resname(line)
        # Standard + modified amino acids stay with the protein (PDBFixer converts
        # modified ones, e.g. MSE->MET); never split them off as GAFF ligands.
        if name in STANDARD_RESIDUES or name in modified:
            out["protein"].append(line)
        elif name in WATER:
            out["water"].append(line)
        elif name in IONS:
            out["ion"].setdefault(name, []).append(line)
        elif want is not None:
            # explicit ligand list: only those are ligands, others stripped
            if name in want:
                out["ligand"].setdefault(name, []).append(line)
            else:
                out["stripped"][name] = out["stripped"].get(name, 0) + 1
        elif strip_agents and name in CRYSTALLISATION_AGENTS:
            out["stripped"][name] = out["stripped"].get(name, 0) + 1
        else:
            out["ligand"].setdefault(name, []).append(line)
    return out


def write_lines(lines, out_path):
    with open(out_path, "w") as fh:
        fh.writelines(lines)
        fh.write("END\n")
    return out_path


def _set_uniform_residue(mol, resname):
    """Force every atom of an OpenFF molecule into one residue.

    Hydrogens added by RDKit ``AddHs`` (the SMILES path) carry no residue
    metadata, so ``to_openmm()`` would split the ligand into two residues and
    break GAFF matching. Setting uniform metadata makes it a single residue.
    """
    name = (resname or "LIG")[:3].upper()
    for atom in mol.atoms:
        atom.metadata["residue_name"] = name
        atom.metadata["residue_number"] = 1
        atom.metadata["insertion_code"] = " "
        atom.metadata["chain_id"] = "L"
    return mol


def build_ligand_molecule(ligand_pdb, resname="LIG", smiles=None):
    """Build an OpenFF ``Molecule`` for one ligand residue.

    With ``smiles`` (recommended for robustness) the PDB coordinates are mapped
    onto the SMILES graph (correct bonds + hydrogens, even for an un-protonated
    ligand). Without it, RDKit perceives bond orders from the PROTONATED 3D
    structure using a net-charge search -- so the ligand must carry its
    hydrogens. Raises ``ValueError`` with a clear message if neither works. The
    returned molecule is normalised to a single ``resname`` residue.
    """
    Molecule, _, Chem, rdDetermineBonds = _require()
    if smiles:
        # Assign bond orders from the SMILES template onto the PDB heavy atoms
        # (matched by graph), then add hydrogens at 3D positions. This works even
        # when the ligand is UN-protonated (heavy atoms only) -- the SMILES
        # defines the bonds, charge and protonation.
        from rdkit.Chem import AllChem
        template = Chem.RemoveHs(Chem.MolFromSmiles(smiles))
        pdbmol = Chem.MolFromPDBFile(ligand_pdb, removeHs=True, sanitize=False)
        if template is None or pdbmol is None:
            raise ValueError("could not read SMILES/PDB for ligand {0}".format(
                os.path.basename(ligand_pdb)))
        if template.GetNumAtoms() != pdbmol.GetNumAtoms():
            raise ValueError(
                "ligand {0}: SMILES heavy-atom count ({1}) != PDB heavy atoms "
                "({2}); check the SMILES hint.".format(
                    os.path.basename(ligand_pdb), template.GetNumAtoms(),
                    pdbmol.GetNumAtoms()))
        mol = AllChem.AssignBondOrdersFromTemplate(template, pdbmol)
        mol = Chem.AddHs(mol, addCoords=True)
        return _set_uniform_residue(
            Molecule.from_rdkit(mol, allow_undefined_stereo=True), resname)

    rd = Chem.MolFromPDBFile(ligand_pdb, removeHs=False, sanitize=False)
    if rd is None:
        raise ValueError("RDKit could not read ligand PDB {0}".format(ligand_pdb))
    has_h = any(a.GetAtomicNum() == 1 for a in rd.GetAtoms())
    if not has_h:
        raise ValueError(
            "ligand {0} has no hydrogens; bond perception needs a protonated "
            "ligand or a SMILES hint (ligand_smiles).".format(os.path.basename(ligand_pdb)))
    last = None
    for charge in (0, 1, -1, 2, -2, 3, -3):
        try:
            m = Chem.Mol(rd)
            rdDetermineBonds.DetermineBonds(m, charge=charge)
            Chem.SanitizeMol(m)
            return _set_uniform_residue(
                Molecule.from_rdkit(m, allow_undefined_stereo=True), resname)
        except Exception as exc:  # noqa: BLE001 - try the next candidate charge
            last = exc
    raise ValueError(
        "could not perceive bond orders for ligand {0} (tried net charges "
        "0, +-1, +-2, +-3); supply a SMILES hint. Last error: {1}".format(
            os.path.basename(ligand_pdb), last))


def load_ligand_molecules(sdf_path):
    """Load one or many ligand ``Molecule`` objects from an SDF file."""
    Molecule, _, _, _ = _require()
    mols = Molecule.from_file(sdf_path)
    return mols if isinstance(mols, list) else [mols]


def write_ligand_sdf(molecules, sdf_path):
    """Write all ligand molecules (with perceived bonds/charges) to one SDF."""
    if not molecules:
        return None
    # openff writes a single molecule per to_file; append for multiples.
    import tempfile
    with open(sdf_path, "w") as out:
        for mol in molecules:
            tmp = sdf_path + ".one"
            mol.to_file(tmp, "SDF")
            with open(tmp) as fh:
                out.write(fh.read())
            os.remove(tmp)
    return sdf_path


def make_gaff_generator(molecules, cache=None, forcefield="gaff-2.11"):
    """A GAFF template generator for the given ligand molecules (cached charges)."""
    _, GAFFTemplateGenerator, _, _ = _require()
    return GAFFTemplateGenerator(molecules=molecules, forcefield=forcefield, cache=cache)


def build_forcefield(config, app):
    """``app.ForceField`` for the protein (amber14) plus, if a ligand SDF is set
    on the config (``config.ligand_sdf``), a registered GAFF template generator so
    ``createSystem`` / ``addSolvent`` can parameterise the ligand(s). AM1-BCC
    charges are cached (``config.gaff_cache``) so they are computed only once."""
    ff = app.ForceField(*config.forcefield_files())
    sdf = getattr(config, "ligand_sdf", None)
    if sdf and os.path.exists(sdf):
        mols = load_ligand_molecules(sdf)
        cache = getattr(config, "gaff_cache", None)
        generator = make_gaff_generator(mols, cache=cache)
        ff.registerTemplateGenerator(generator.generator)
    return ff

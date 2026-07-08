"""Structure preparation and solvation (OpenMM, imported lazily).

Follows the BilboMD worker pattern: PDBFixer to add missing atoms/residues and
``Modeller.addHydrogens(forcefield, pH=...)`` for protonation, then explicit
solvation with ``Modeller.addSolvent`` using a padding derived from the model
Dmax and the requested ionic strength.

OpenMM/PDBFixer are imported only when these functions run, so importing this
module never requires them. Audit trail: callers should record what PDBFixer
changed (missing residues/atoms/terminals) in the manifest rather than silently
renaming atoms.
"""

import os

from ..command_runner import MissingDependencyError
from ..dmax import box_padding_nm, model_dmax_nm, resolve_dmax_nm
from . import ligand, protonation
from .schema import OpenMMConfig


def _require_openmm():
    try:
        import openmm
        from openmm import app, unit
        import pdbfixer
    except ImportError as exc:
        raise MissingDependencyError(
            "OpenMM structure preparation requires openmm and pdbfixer "
            "(conda-forge). Original import error: {0}".format(exc))
    return openmm, app, unit, pdbfixer


def resolve_box_padding(config: OpenMMConfig, pdb_path: str):
    """Box padding (nm): explicit config value, else |Dmax_expt - Dmax_model| + 2.

    Reuses the Phase 1 :mod:`automd_saxs.dmax` calculators so the OpenMM and
    GROMACS branches size boxes identically.
    """
    if config.box_padding_nm is not None:
        return config.box_padding_nm
    model = model_dmax_nm(pdb_path)
    effective = resolve_dmax_nm(None, model)  # no experimental Dmax field yet
    return box_padding_nm(effective, model)


def prepare_structure(config: OpenMMConfig, pdb_path: str, out_pdb: str):
    """Clean the input PDB and add hydrogens at the requested pH.

    Protonation is structure-based: `propka3` predicts per-residue pKa values and
    the requested pH decides each ionizable group's state (see
    :mod:`automd_saxs.openmm.protonation`), applied via ``addHydrogens`` variants.
    Falls back to plain ``addHydrogens(pH=...)`` if propka is unavailable/fails.

    Returns an audit dict describing what PDBFixer changed and the protonation
    decisions (method + per-residue overrides), for the manifest / UI. Protonation
    is a best-effort prediction and should be checked visually.
    """
    openmm, app, unit, pdbfixer = _require_openmm()
    work_dir = os.path.dirname(os.path.abspath(out_pdb))

    # Split the input into protein / ligand(s) / water / ions. Water and (for now)
    # ions and crystallisation agents are set aside; the protein is prepped +
    # protonated (propka) and each ligand is parameterised with GAFF.
    classes = ligand.classify_residues(
        pdb_path,
        strip_agents=not getattr(config, "keep_crystallisation_agents", False),
        ligand_resnames=getattr(config, "ligand_resnames", None))
    protein_pdb = ligand.write_lines(
        classes["protein"], os.path.join(work_dir, "protein_only.pdb"))

    fixer = pdbfixer.PDBFixer(filename=protein_pdb)
    fixer.findMissingResidues()
    # Convert modified residues (e.g. MSE selenomethionine -> MET) to their
    # standard forms so amber14 has templates for them (robustness on raw PDBs).
    fixer.findNonstandardResidues()
    nonstandard = [(str(r), std) for r, std in fixer.nonstandardResidues]
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    missing_residues = dict(fixer.missingResidues)
    missing_atoms = {str(k): [a.name for a in v] for k, v in fixer.missingAtoms.items()}
    fixer.addMissingAtoms()
    fixer.findMissingAtoms()

    # Consistent chain ids so propka residue identifiers line up with the topology.
    protonation.relabel_chains(fixer.topology)

    forcefield = app.ForceField(*config.forcefield_files())
    modeller = app.Modeller(fixer.topology, fixer.positions)

    # Structure-based pH protonation (propka) with graceful fallback.
    proton_method = "openmm-default"
    proton_changes = []
    pkas = {}
    if getattr(config, "use_propka", True):
        heavy_pdb = os.path.join(work_dir, "prepared_heavy.pdb")
        with open(heavy_pdb, "w") as handle:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
        pkas = protonation.run_propka(heavy_pdb, work_dir)

    if pkas:
        variants, proton_changes = protonation.build_variants(
            modeller.topology, pkas, config.ph)
        modeller.addHydrogens(forcefield, pH=config.ph, variants=variants)
        proton_method = "propka"
    else:
        modeller.addHydrogens(forcefield, pH=config.ph)

    # Ligand parameterisation: build an OpenFF molecule per ligand residue
    # (RDKit bond perception on the protonated ligand, or a SMILES hint), persist
    # to ligand.sdf for downstream GAFF, and merge the ligand atoms (with their
    # own hydrogens) into the prepared complex.
    ligands_audit = []
    lig_mols = []
    ligand_warnings = []
    smiles_hints = getattr(config, "ligand_smiles", None) or {}
    explicit = set(getattr(config, "ligand_resnames", None) or []) | set(smiles_hints)
    for resname, lines in classes["ligand"].items():
        lig_pdb = ligand.write_lines(
            lines, os.path.join(work_dir, "ligand_{0}.pdb".format(resname)))
        try:
            mol = ligand.build_ligand_molecule(
                lig_pdb, resname=resname, smiles=smiles_hints.get(resname))
        except Exception as exc:  # noqa: BLE001
            # A residue the user explicitly asked to keep must not be silently
            # dropped; otherwise strip the unparameterisable HETATM (e.g. a
            # crystallisation agent with no hydrogens) and warn, so the job still
            # runs. The Task 4 UI lets the user force-keep with a SMILES hint.
            if resname in explicit:
                raise
            classes["stripped"][resname] = classes["stripped"].get(resname, 0) + len(lines)
            ligand_warnings.append(
                "stripped ligand candidate {0} (could not parameterise: {1})".format(
                    resname, str(exc)[:120]))
            continue
        lig_mols.append(mol)
        modeller.add(mol.to_topology().to_openmm(), mol.conformers[0].to_openmm())
        ligands_audit.append({
            "resname": resname,
            "nAtoms": mol.n_atoms,
            "formalCharge": float(mol.total_charge.magnitude)
            if hasattr(mol.total_charge, "magnitude") else float(mol.total_charge / mol.total_charge.unit),
            "smiles": mol.to_smiles(explicit_hydrogens=False),
            "source": "smiles" if resname in smiles_hints else "perceived",
        })
    if lig_mols and getattr(config, "ligand_sdf", None):
        ligand.write_ligand_sdf(lig_mols, config.ligand_sdf)

    # Keep bound/structural ions (Ca2+, Zn2+, Pb2+, Fe, Mg, ...): merge them back
    # as their own residues so amber14 parameterises them and solvation
    # neutralises accounting for their charge. Merged after the protein so they
    # never confuse chain-terminus detection.
    ions_audit = {}
    if getattr(config, "keep_ions", True) and classes["ion"]:
        ions_audit = _merge_ions(modeller, classes["ion"], app, openmm, unit)

    with open(out_pdb, "w") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle)

    return {
        "missingResidues": {str(k): v for k, v in missing_residues.items()},
        "nonstandardResidues": nonstandard,
        "missingAtoms": missing_atoms,
        "pH": config.ph,
        "protonationMethod": proton_method,
        "protonationChanges": proton_changes,
        "ligands": ligands_audit,
        "ligandWarnings": ligand_warnings,
        "strippedResidues": classes["stripped"],
        "ions": ions_audit if ions_audit else {},
        "ionsSetAside": ({} if getattr(config, "keep_ions", True)
                         else {k: len(v) for k, v in classes["ion"].items()}),
        "output": out_pdb,
    }


def _merge_ions(modeller, ion_classes, app, openmm, unit):
    """Merge bound ions into the prepared complex as amber14-named residues.

    ``ion_classes`` is ``{pdb_resname: [pdb lines]}``. Each ion becomes its own
    single-atom residue named for the amber14 ion template (see
    :func:`ligand.amber_ion_resname`) with the element taken from the PDB element
    column (falling back to the residue name). Returns ``{amber_resname: count}``.
    """
    topology = app.Topology()
    chain = topology.addChain("I")
    positions = []
    kept = {}
    for pdb_resname, lines in ion_classes.items():
        amber = ligand.amber_ion_resname(pdb_resname)
        for line in lines:
            element_symbol = (line[76:78].strip() or pdb_resname).capitalize()
            try:
                element = app.Element.getBySymbol(element_symbol)
            except Exception:  # noqa: BLE001 - fall back to the residue name
                element = app.Element.getBySymbol(pdb_resname.capitalize())
            residue = topology.addResidue(amber, chain)
            topology.addAtom(line[12:16].strip() or amber, element, residue)
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            positions.append(openmm.Vec3(x, y, z) * 0.1)  # Angstrom -> nm
            kept[amber] = kept.get(amber, 0) + 1
    modeller.add(topology, positions * unit.nanometer)
    return kept


def solvate(config: OpenMMConfig, prepared_pdb: str, out_pdb: str, padding_nm: float):
    """Add explicit water and ions to a neutral box (Modeller.addSolvent)."""
    openmm, app, unit, pdbfixer = _require_openmm()

    pdb = app.PDBFile(prepared_pdb)
    forcefield = ligand.build_forcefield(config, app)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(
        forcefield,
        model=config.water_model.value,
        padding=padding_nm * unit.nanometer,
        ionicStrength=config.ionic_concentration_M * unit.molar,
        neutralize=True,
    )
    with open(out_pdb, "w") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle)
    return out_pdb

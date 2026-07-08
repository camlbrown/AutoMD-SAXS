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
from . import protonation
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

    fixer = pdbfixer.PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
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

    with open(out_pdb, "w") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle)

    return {
        "missingResidues": {str(k): v for k, v in missing_residues.items()},
        "missingAtoms": missing_atoms,
        "pH": config.ph,
        "protonationMethod": proton_method,
        "protonationChanges": proton_changes,
        "output": out_pdb,
    }


def solvate(config: OpenMMConfig, prepared_pdb: str, out_pdb: str, padding_nm: float):
    """Add explicit water and ions to a neutral box (Modeller.addSolvent)."""
    openmm, app, unit, pdbfixer = _require_openmm()

    pdb = app.PDBFile(prepared_pdb)
    forcefield = app.ForceField(*config.forcefield_files())
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

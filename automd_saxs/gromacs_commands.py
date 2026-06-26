"""GROMACS command construction (no execution).

Each function returns an ``argv`` list -- never a shell string -- for one step
of the setup pipeline that the legacy ``run_MD.sh`` ran on the login node before
submitting jobs:

    pdb2gmx -> editconf (box) -> editconf (center) -> solvate
            -> grompp (ions) -> genion

These builders are pure and unit-testable. They deliberately differ from the
legacy script in two safe ways:

* **No wildcards.** The legacy used ``-c *.gro`` / ``-p *.top``; here filenames
  are explicit arguments, removing the "exactly one match" fragility.
* **Ionic concentration is honoured.** :func:`genion_command` takes the
  concentration from :class:`~automd_saxs.config.JobConfig`, fixing the legacy
  hardcoded ``-conc 0.15``.

The mdrun-bearing stages (minim/nvt/npt/production) run inside Slurm scripts and
are planned in :mod:`automd_saxs.slurm`, not here.
"""

from typing import List

from .config import BoxShape, JobConfig

#: Default GROMACS executable (the legacy scripts used the MPI build).
DEFAULT_GMX = "gmx_mpi"

#: Standard solvent coordinate file shipped with GROMACS (legacy ``solvate -cs``).
SOLVENT_BOX = "spc216.gro"

#: Index group piped to ``genion`` on stdin to choose which molecules to replace
#: with ions (legacy ``echo "SOL" | gmx genion ...``).
GENION_REPLACE_GROUP = "SOL"

# GROMACS ``editconf -bt`` accepts these box types. ``BoxShape.AUTO`` is not a
# real box type -- it signals a runtime radius-of-gyration heuristic and must be
# resolved to a concrete shape before a command can be built.
_BT_NAMES = {
    BoxShape.DODECAHEDRON: "dodecahedron",
    BoxShape.TRICLINIC: "triclinic",
    BoxShape.OCTAHEDRON: "octahedron",
    BoxShape.CUBIC: "cubic",
}


def pdb2gmx_command(
    config: JobConfig,
    input_pdb: str = "GMX.pdb",
    output_gro: str = "GMX.gro",
    topology: str = "topol.top",
    gmx: str = DEFAULT_GMX,
) -> List[str]:
    """Build the topology with ``pdb2gmx`` (protein-only path).

    Mirrors the legacy protein branch::

        gmx pdb2gmx -f GMX.pdb -o GMX.gro -p topol.top -chainsep ter \\
            -ff <ff> -water <water> -ter -merge all [-ss]
    """
    argv = [
        gmx, "pdb2gmx",
        "-f", input_pdb,
        "-o", output_gro,
        "-p", topology,
        "-chainsep", "ter",
        "-ff", config.force_field_name,
        "-water", config.water_model(),
        "-ter",
        "-merge", "all",
    ]
    if config.disulfide:
        argv.append("-ss")
    return argv


def editconf_box_command(
    config: JobConfig,
    box_padding_nm: float,
    input_gro: str = "GMX.gro",
    output_gro: str = "1.gro",
    gmx: str = DEFAULT_GMX,
) -> List[str]:
    """Define the periodic box with ``editconf -bt <shape> -d <padding>``."""
    if config.box_shape is BoxShape.AUTO:
        raise ValueError(
            "BoxShape.AUTO cannot be planned statically; resolve it to a "
            "concrete shape (the legacy code used a runtime gyrate heuristic)."
        )
    bt = _BT_NAMES[config.box_shape]
    return [
        gmx, "editconf",
        "-f", input_gro,
        "-o", output_gro,
        "-bt", bt,
        "-d", _fmt(box_padding_nm),
    ]


def editconf_center_command(
    input_gro: str = "1.gro",
    output_gro: str = "centered.gro",
    gmx: str = DEFAULT_GMX,
) -> List[str]:
    """Centre the solute in the box (legacy ``editconf -c``)."""
    return [gmx, "editconf", "-f", input_gro, "-o", output_gro, "-c"]


def solvate_command(
    topology: str = "topol.top",
    centered_gro: str = "centered.gro",
    output_gro: str = "solvate.gro",
    solvent: str = SOLVENT_BOX,
    gmx: str = DEFAULT_GMX,
) -> List[str]:
    """Add explicit water (legacy ``solvate -cp centered.gro -cs spc216.gro``)."""
    return [
        gmx, "solvate",
        "-cp", centered_gro,
        "-cs", solvent,
        "-o", output_gro,
        "-p", topology,
    ]


def ions_grompp_command(
    ions_mdp: str,
    structure_gro: str = "solvate.gro",
    topology: str = "topol.top",
    output_tpr: str = "ions.tpr",
    maxwarn: int = 100,
    gmx: str = DEFAULT_GMX,
) -> List[str]:
    """Assemble the ions ``.tpr`` (legacy ``grompp -f ions.mdp ... -maxwarn 100``)."""
    return [
        gmx, "grompp",
        "-f", ions_mdp,
        "-c", structure_gro,
        "-p", topology,
        "-o", output_tpr,
        "-maxwarn", str(maxwarn),
    ]


def genion_command(
    config: JobConfig,
    input_tpr: str = "ions.tpr",
    output_gro: str = "solvate.gro",
    topology: str = "topol.top",
    pname: str = "Na",
    nname: str = "Cl",
    gmx: str = DEFAULT_GMX,
) -> List[str]:
    """Add neutralising ions at the chosen concentration.

    Mirrors the legacy command but reads the concentration from the config
    instead of hardcoding ``-conc 0.15``. The replace-group selection
    (:data:`GENION_REPLACE_GROUP`) is fed on stdin by the runner, matching
    ``echo "SOL" | gmx genion``.
    """
    return [
        gmx, "genion",
        "-s", input_tpr,
        "-o", output_gro,
        "-p", topology,
        "-pname", pname,
        "-nname", nname,
        "-neutral",
        "-conc", _fmt(config.genion_concentration()),
    ]


def setup_pipeline(config: JobConfig, box_padding_nm: float, ions_mdp: str):
    """Return the ordered list of ``(label, argv)`` setup steps.

    This is the protein-only login-node pipeline run by ``run_MD.sh`` before
    Slurm submission. The protein-ligand path (acpype/itp splicing) is not yet
    ported and is reported as unsupported by the planner.
    """
    return [
        ("pdb2gmx", pdb2gmx_command(config)),
        ("editconf_box", editconf_box_command(config, box_padding_nm)),
        ("editconf_center", editconf_center_command()),
        ("solvate", solvate_command()),
        ("ions_grompp", ions_grompp_command(ions_mdp)),
        ("genion", genion_command(config)),
    ]


def _fmt(value: float) -> str:
    """Format a float for the command line without trailing noise."""
    return ("%g" % value)

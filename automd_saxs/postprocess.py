"""Trajectory post-processing: per-repeat cleanup + combine (trjcat/align).

Ports the GROMACS trajectory-processing pipeline from the legacy
``post_processing.sh`` into argv builders and ordered command plans. For each
production repeat it: centres the final frame, builds a solvent-free index,
strips solvent from the gro/tpr/xtc, removes PBC jumps, fits to remove
rotation/translation, builds a reference ``final.pdb``, and extracts frames for
SAXS. Then it concatenates the three ``final.xtc`` and aligns them into
``combined_aligned.xtc`` -- the input the ``cluster`` stage consumes.

Each step is ``(label, argv, stdin)``: GROMACS index-group selections that the
legacy script fed via here-docs are carried as the ``stdin`` string. Pure and
testable; execution happens through :class:`~automd_saxs.command_runner.CommandRunner`.
"""

from typing import List, Optional, Tuple

DEFAULT_GMX = "gmx_mpi"

# GROMACS default index groups used by the legacy here-docs.
_GROUP_SYSTEM = "0"
_GROUP_PROTEIN = "1"
_GROUP_CALPHA = "3"

Step = Tuple[str, List[str], Optional[str]]

# Working filenames produced inside each processed/<rep> directory (legacy names).
MODEL_PDB = "model.pdb"
NOSOLV_NDX = "nosolv.ndx"
NOSOLV_GRO = "nosolv.gro"
NOSOLV_TPR = "nosolv.tpr"
NOSOLV_XTC = "nosolv.xtc"
NOJUMP_XTC = "nojump.xtc"
FINAL_XTC = "final.xtc"
FINAL_PDB = "final.pdb"
COMBINED_XTC = "combined.xtc"
COMBINED_ALIGNED_XTC = "combined_aligned.xtc"


def plan_repeat_processing(
    src_gro: str,
    src_tpr: str,
    src_xtc: str,
    reference_gro: str,
    skip: int = 2,
    gmx: str = DEFAULT_GMX,
) -> List[Step]:
    """Ordered steps to clean one repeat's trajectory (run inside its processed dir).

    ``reference_gro`` is the topology-stage ``GMX.gro`` used to build ``final.pdb``.
    """
    return [
        ("editconf_center",
         [gmx, "editconf", "-f", src_gro, "-o", MODEL_PDB, "-c"], None),
        ("make_ndx",
         [gmx, "make_ndx", "-f", MODEL_PDB, "-o", NOSOLV_NDX],
         _GROUP_PROTEIN + "\nq\n"),
        ("editconf_nosolv",
         [gmx, "editconf", "-f", MODEL_PDB, "-n", NOSOLV_NDX, "-o", NOSOLV_GRO],
         _GROUP_PROTEIN + "\n"),
        ("convert_tpr",
         [gmx, "convert-tpr", "-s", src_tpr, "-n", NOSOLV_NDX, "-o", NOSOLV_TPR],
         _GROUP_PROTEIN + "\n"),
        ("trjconv_nosolv",
         [gmx, "trjconv", "-f", src_xtc, "-s", NOSOLV_TPR, "-n", NOSOLV_NDX,
          "-o", NOSOLV_XTC], _GROUP_PROTEIN + "\n"),
        ("trjconv_nojump",
         [gmx, "trjconv", "-f", NOSOLV_XTC, "-s", NOSOLV_TPR, "-pbc", "nojump",
          "-o", NOJUMP_XTC], _GROUP_SYSTEM + "\n"),
        ("trjconv_fit",
         [gmx, "trjconv", "-f", NOJUMP_XTC, "-s", NOSOLV_GRO, "-o", FINAL_XTC,
          "-fit", "rot+trans"], _GROUP_SYSTEM + "\n" + _GROUP_SYSTEM + "\n"),
        ("editconf_final_pdb",
         [gmx, "editconf", "-f", reference_gro, "-n", NOSOLV_NDX, "-o", FINAL_PDB],
         _GROUP_PROTEIN + "\n"),
        ("extract_frames",
         [gmx, "trjconv", "-f", FINAL_XTC, "-s", NOSOLV_TPR, "-o", "structure_.pdb",
          "-sep", "-skip", str(skip)], _GROUP_PROTEIN + "\n"),
    ]


def plan_combine(
    final_xtcs: List[str],
    reference_pdb: str,
    skip: int = 2,
    gmx: str = DEFAULT_GMX,
) -> List[Step]:
    """Steps to concatenate per-repeat ``final.xtc`` and align them.

    Produces ``combined.xtc`` then ``combined_aligned.xtc`` (fit on C-alpha,
    output Protein), matching the legacy ``trjcat -cat`` + ``trjconv -fit``.
    """
    trjcat = [gmx, "trjcat", "-f"] + list(final_xtcs) + ["-o", COMBINED_XTC, "-cat"]
    align = [gmx, "trjconv", "-s", reference_pdb, "-f", COMBINED_XTC,
             "-o", COMBINED_ALIGNED_XTC, "-skip", str(skip), "-fit", "rot+trans"]
    return [
        ("trjcat", trjcat, None),
        ("align_combined", align, _GROUP_CALPHA + "\n" + _GROUP_PROTEIN + "\n"),
    ]

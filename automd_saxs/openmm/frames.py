"""Trajectory frame extraction and combination (mdtraj, imported lazily).

Replaces the GROMACS ``trjconv``/``trjcat`` post-processing with mdtraj
operations on OpenMM DCD output: extract per-frame PDBs for FoXS, and combine
the per-repeat trajectories for clustering. mdtraj is imported only when these
run, so importing this module never requires it.

The OpenMM trajectories are explicit-solvent (protein + water + ions). For SAXS
fitting (FoXS models the hydration layer implicitly via c1/c2) and for structural
clustering we want the SOLUTE only, so frames/combined trajectories are stripped
to the protein selection by default. This matches the GROMACS branch (which
removed solvent before SAXS) and avoids writing huge ~400k-atom water-box PDBs.
"""

import os

from ..command_runner import MissingDependencyError

# mdtraj selection for the macromolecule (drops water, ions). Falls back to all
# atoms if a topology has no matching atoms (e.g. unusual residue naming).
DEFAULT_SELECTION = "protein"


def _require_mdtraj():
    try:
        import mdtraj
    except ImportError as exc:
        raise MissingDependencyError(
            "Frame extraction requires mdtraj (conda-forge). "
            "Original import error: {0}".format(exc))
    return mdtraj


def _solute(traj, selection):
    """Return the solute-only slice of a trajectory, or the original if empty."""
    idx = traj.topology.select(selection)
    if len(idx) == 0:
        return traj
    return traj.atom_slice(idx)


def extract_frames(trajectory, topology, out_dir, stride=2, selection=DEFAULT_SELECTION):
    """Write stride-sampled SOLUTE frames as ``structure_<i>.pdb`` into ``out_dir``.

    Water/ions are stripped (``selection``) so each frame is the protein only —
    correct for FoXS and far smaller on disk. Returns the written PDB paths.
    """
    mdtraj = _require_mdtraj()
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    traj = _solute(mdtraj.load(trajectory, top=topology, stride=stride), selection)
    written = []
    for i in range(traj.n_frames):
        path = os.path.join(out_dir, "structure_{0}.pdb".format(i))
        traj[i].save_pdb(path)
        written.append(path)
    return written


def structural_timeseries(trajectory, topology, stride=1, selection=DEFAULT_SELECTION):
    """Per-frame structural metrics for one repeat's trajectory.

    Returns a list of ``{frame, rg, rmsd, sasa}`` (Rg and Cα-RMSD-to-first-frame in
    Angstrom; total SASA in nm^2) computed on the solute. Mirrors the GROMACS
    branch's gyrate/rms/sasa analysis but via mdtraj so it runs in the OpenMM
    pipeline. Best-effort: returns [] if it cannot be computed.
    """
    mdtraj = _require_mdtraj()
    try:
        traj = _solute(mdtraj.load(trajectory, top=topology, stride=stride), selection)
        rg = mdtraj.compute_rg(traj) * 10.0  # nm -> Angstrom
        ca = traj.topology.select("name CA")
        rmsd = (mdtraj.rmsd(traj, traj, 0, atom_indices=ca) * 10.0
                if len(ca) else [0.0] * traj.n_frames)
        try:
            sasa = mdtraj.shrake_rupley(traj, mode="atom").sum(axis=1)
        except Exception:  # noqa: BLE001 - SASA is advisory
            sasa = [None] * traj.n_frames
        out = []
        for i in range(traj.n_frames):
            out.append({
                "frame": i,
                "rg": float(rg[i]),
                "rmsd": float(rmsd[i]),
                "sasa": (float(sasa[i]) if sasa[i] is not None else None),
            })
        return out
    except Exception:  # noqa: BLE001 - time-series is advisory; never fail the run
        return []


def read_pca_coords(path):
    """Read a clustering PCA_coords.txt into a list of [pc1, pc2, ...] float rows.

    Header line (``PC1(0.43) PC2(0.21)``) is skipped. Returns [] on any error.
    """
    try:
        rows = []
        with open(path) as fh:
            lines = fh.read().splitlines()
        for line in lines[1:]:
            parts = line.split()
            if not parts:
                continue
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
        return rows
    except Exception:  # noqa: BLE001
        return []


def rg_of_pdb(pdb_path):
    """Radius of gyration (Angstrom) of a single-structure PDB, or None on error."""
    mdtraj = _require_mdtraj()
    try:
        t = mdtraj.load(pdb_path)
        # mdtraj returns Rg in nm; report Angstrom to match FoXS/SAXS conventions.
        return float(mdtraj.compute_rg(t)[0] * 10.0)
    except Exception:  # noqa: BLE001 - Rg is advisory; never fail the run
        return None


def combine_trajectories(trajectories, topology, out_path, stride=1,
                         selection=DEFAULT_SELECTION):
    """Concatenate per-repeat trajectories into one solute-only DCD. Returns out_path.

    Also writes a solute-only reference topology (``<out>.pdb``) alongside the DCD
    so the combined trajectory can be loaded without the full solvated topology.
    """
    mdtraj = _require_mdtraj()
    frames = [_solute(mdtraj.load(t, top=topology, stride=stride), selection)
              for t in trajectories]
    if not frames:
        raise ValueError("no trajectories to combine")
    combined = frames[0]
    for extra in frames[1:]:
        combined = combined.join(extra)
    combined.save_dcd(out_path)
    # Solute-only topology for downstream loaders (clustering).
    top_pdb = os.path.splitext(out_path)[0] + ".pdb"
    combined[0].save_pdb(top_pdb)
    return out_path

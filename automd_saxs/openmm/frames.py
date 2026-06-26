"""Trajectory frame extraction and combination (mdtraj, imported lazily).

Replaces the GROMACS ``trjconv``/``trjcat`` post-processing with mdtraj
operations on OpenMM DCD output: extract per-frame PDBs for FoXS, and combine
the per-repeat trajectories for clustering. mdtraj is imported only when these
run, so importing this module never requires it.
"""

import os

from ..command_runner import MissingDependencyError


def _require_mdtraj():
    try:
        import mdtraj
    except ImportError as exc:
        raise MissingDependencyError(
            "Frame extraction requires mdtraj (conda-forge). "
            "Original import error: {0}".format(exc))
    return mdtraj


def extract_frames(trajectory: str, topology: str, out_dir: str, stride: int = 2):
    """Write stride-sampled frames as ``structure_<i>.pdb`` into ``out_dir``.

    Returns the list of written PDB paths.
    """
    mdtraj = _require_mdtraj()
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    traj = mdtraj.load(trajectory, top=topology, stride=stride)
    written = []
    for i in range(traj.n_frames):
        path = os.path.join(out_dir, "structure_{0}.pdb".format(i))
        traj[i].save_pdb(path)
        written.append(path)
    return written


def combine_trajectories(trajectories, topology: str, out_path: str, stride: int = 1):
    """Concatenate per-repeat trajectories into one DCD (mdtraj). Returns out_path."""
    mdtraj = _require_mdtraj()
    frames = [mdtraj.load(t, top=topology, stride=stride) for t in trajectories]
    if not frames:
        raise ValueError("no trajectories to combine")
    combined = frames[0]
    for extra in frames[1:]:
        combined = combined.join(extra)
    combined.save_dcd(out_path)
    return out_path

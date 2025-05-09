import MDAnalysis as mda
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--coords')
parser.add_argument('--traj')
parser.add_argument('--skip', type=int)
parser.add_argument('--out')
args = parser.parse_args()

u = mda.Universe(args.coords, args.traj)

atoms = u.select_atoms('all')

with mda.Writer(f'{args.out}.xtc', atoms.n_atoms) as f:
    for ts in u.trajectory[::args.skip]:
        f.write(atoms)

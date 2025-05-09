from MDAnalysis.analysis import align
import MDAnalysis as mda
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--coords')
parser.add_argument('--traj')
parser.add_argument('--out')
args = parser.parse_args()

u = mda.Universe(args.coords, args.traj)
u_ref = mda.Universe(args.coords)

alignment = align.AlignTraj(u, u_ref, select="protein and name CA", filename=f'{args.out}.dcd')
alignment.run()

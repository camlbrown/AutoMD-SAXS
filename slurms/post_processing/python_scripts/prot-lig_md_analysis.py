import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.ticker import ScalarFormatter


# 1) Lig RMSD
fig, ax = plt.subplots(figsize=(6,4))
data = np.loadtxt('rmsd_ligand.xvg', comments=('@','#'))
time_ns = data[:, 0] /1000
rmsd    = data[:, 1]
ax.plot(time_ns, rmsd, linewidth=1)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Ligand RMSD (Å)")
fig.tight_layout()
fig.savefig('ligand_rmsd.svg')
fig.savefig('ligand_rmsd.pdf')
plt.close(fig)

# 2) Lig RMSF
fig, ax = plt.subplots(figsize=(6,4))
data = np.loadtxt('rmsf_ligand.xvg', comments=('@','#'))
time_ns = data[:, 0]
rmsf    = data[:, 1]
ax.plot(time_ns, rmsf, linewidth=1)
ax.set_xlabel("Ligand Atom Number")
ax.set_ylabel("RMSF (Å)")
fig.tight_layout()
fig.savefig('ligand_rmsf.svg')
fig.savefig('ligand_rmsf.pdf')
plt.close(fig)

# 4) Prot-lig H-bond count
fig, ax = plt.subplots(figsize=(6,4))
data = np.loadtxt('hb_lig.xvg', comments=('@','#'))
time_ns = data[:, 0] / 1000.0 
hbond   = data[:, 1]
ax.plot(time_ns, hbond, linewidth=1)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Protein-Ligand H-Bond Count")
fig.tight_layout()
fig.savefig('prot-lig_hbonds.svg')
fig.savefig('prot-lig_hbonds.pdf')
plt.close(fig)


fig = plt.figure(figsize=(6, 12))
gs = gridspec.GridSpec(
    nrows=3, ncols=1,
    height_ratios=[1, 1, 1],
    hspace=0.4, wspace=0.3
)

ax_rmsd   = fig.add_subplot(gs[0, 0])
ax_rmsf   = fig.add_subplot(gs[1, 0])
ax_hbond  = fig.add_subplot(gs[2, 0])

# 1) RMSD
data = np.loadtxt('rmsd_ligand.xvg', comments=('@','#'))
ax_rmsd.plot(data[:,0]/1000, data[:,1], linewidth=1)
ax_rmsd.set_xlabel("Time (ns)")
ax_rmsd.set_ylabel("RMSD (Å)")

# 1) RMSF
data = np.loadtxt('rmsf_ligand.xvg', comments=('@','#'))
ax_rmsf.plot(data[:,0], data[:,1], linewidth=1)
ax_rmsf.set_xlabel("Ligand Atom Number)")
ax_rmsf.set_ylabel("RMSF (Å)")

# 4) System H-bonds
data = np.loadtxt('hb_lig.xvg', comments=('@','#'))
ax_hbond.plot(data[:,0]/1000, data[:,1], linewidth=1)
ax_hbond.set_xlabel("Time (ns)")
ax_hbond.set_ylabel("Protein-Ligand H-Bond Count")

plt.tight_layout()
fig.savefig('prot-lig_md_analysis.svg')
fig.savefig('prot-lig_md_analysis.pdf')


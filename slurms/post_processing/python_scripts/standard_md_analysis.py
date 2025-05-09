import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.ticker import ScalarFormatter


# 1) RMSD
fig, ax = plt.subplots(figsize=(6,4))
data = np.loadtxt('rmsd_backbone.xvg', comments=('@','#'))
time_ns, rmsd = data[:,0], data[:,1]
ax.plot(time_ns, rmsd, linewidth=1)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("C-\u03B1 RMSD (Å)")
fig.tight_layout()
fig.savefig('rmsd_calpha.svg')
fig.savefig('rmsd_calpha.pdf')
plt.close(fig)

# 2) Gyration
fig, ax = plt.subplots(figsize=(6,4))
data = np.loadtxt('gyrate.xvg', comments=('@','#'))
time_ns, gyr = data[:,0]/1e3, data[:,1] * 10
ax.plot(time_ns, gyr, linewidth=1)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Rg (Å)")
fig.tight_layout()
fig.savefig('gyration.svg')
fig.savefig('gyration.pdf')
plt.close(fig)

# 3) SASA
fig, ax = plt.subplots(figsize=(6,4))
data = np.loadtxt('sasa_total.xvg', comments=('@','#'))
time_ns, sasa = data[:,0]/1e3, data[:,1]
ax.plot(time_ns, sasa, linewidth=1)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("SASA (Å²)")
fig.tight_layout()
fig.savefig('sasa.svg')
fig.savefig('sasa.pdf')
plt.close(fig)

# 4) System H-bond count
fig, ax = plt.subplots(figsize=(6,4))
data = np.loadtxt('hb_system.xvg', comments=('@','#'))
time_ns, hbond = data[:,0]/1e3, data[:,1]
ax.plot(time_ns, hbond, linewidth=1)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("System H-Bond Count")
fig.tight_layout()
fig.savefig('hbonds.svg')
fig.savefig('hbonds.pdf')
plt.close(fig)

# 5) Energy 
fig, ax = plt.subplots(figsize=(6,4))
data      = np.loadtxt('energy.xvg', comments=('@','#'))
times_ns  = data[:,0]/1e3
dt        = times_ns[1] - times_ns[0]
stride    = int(round(1.0/dt))
t_sample  = times_ns[::stride]
e_sample  = data[:,1][::stride] * 1e-6
ax.plot(t_sample, e_sample, linewidth=1)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Total Energy (×10⁶ kJ/mol)")
formatter = ScalarFormatter(useOffset=False)
formatter.set_scientific(False)
ax.yaxis.set_major_formatter(formatter)
fig.tight_layout()
fig.savefig('energy.svg')
fig.savefig('energy.pdf')
plt.close(fig)

fig = plt.figure(figsize=(12, 14))
gs = gridspec.GridSpec(
    nrows=3, ncols=2,
    height_ratios=[1, 1, 0.8],
    hspace=0.4, wspace=0.3
)

ax_sasa  = fig.add_subplot(gs[0, 0])
ax_gyr   = fig.add_subplot(gs[0, 1])
ax_rmsd  = fig.add_subplot(gs[1, 0])
ax_hbond = fig.add_subplot(gs[1, 1])

# 1) RMSD
data = np.loadtxt('rmsd_backbone.xvg', comments=('@','#'))
ax_rmsd.plot(data[:,0], data[:,1], linewidth=1)
ax_rmsd.set_xlabel("Time (ns)")
ax_rmsd.set_ylabel("RMSD (Å)")

# 2) Gyration
data = np.loadtxt('gyrate.xvg', comments=('@','#'))
ax_gyr.plot(data[:,0]/1e3, data[:,1]*10, linewidth=1)
ax_gyr.set_xlabel("Time (ns)")
ax_gyr.set_ylabel("Rg (Å)")

# 3) SASA
data = np.loadtxt('sasa_total.xvg', comments=('@','#'))
ax_sasa.plot(data[:,0]/1e3, data[:,1], linewidth=1)
ax_sasa.set_xlabel("Time (ns)")
ax_sasa.set_ylabel("SASA (Å²)")

# 4) System H-bonds
data = np.loadtxt('hb_system.xvg', comments=('@','#'))
ax_hbond.plot(data[:,0]/1e3, data[:,1], linewidth=1)
ax_hbond.set_xlabel("Time (ns)")
ax_hbond.set_ylabel("System H-Bond Count")

# 5) Energy (custom axes)
data     = np.loadtxt('energy.xvg', comments=('@','#'))
times_ns = data[:,0]/1e3
stride   = int(round(1.0/ (times_ns[1] - times_ns[0])))
t_samp   = times_ns[::stride]
e_samp   = data[:,1][::stride] * 1e-6

bb     = ax_rmsd.get_position()
w      = bb.width
h      = bb.height * 1.0
left   = 0.5 - w/2
bottom = 0.08

ax_energy = fig.add_axes([left, bottom, w, h])
ax_energy.plot(t_samp, e_samp, linewidth=1)
ax_energy.set_xlabel("Time (ns)")
ax_energy.set_ylabel("Total Energy (×10⁶ kJ/mol)")
fmt = ScalarFormatter(useOffset=False)
fmt.set_scientific(False)
ax_energy.yaxis.set_major_formatter(fmt)

plt.tight_layout()
fig.savefig('standard_md_analysis.svg')
fig.savefig('standard_md_analysis.pdf')


#!/bin/bash

import matplotlib.pyplot as plt
import seaborn as sns

# read data from file
with open('crysol_summary.txt', 'r') as f:
    lines = f.readlines()

# extract data from lines
pdb_files = []
chi2_values = []
rg_values = []

min_rg = float('inf')
max_rg = -float('inf')
min_chi2 = float('inf')
max_chi2 = -float('inf')

for line in lines:
    line = line.strip()
    if line:
        parts = line.split()
        pdb_file = parts[1].replace('structure_*', '').replace('.pdb', '').split('structure_')[1]
        pdb_files.append(int(pdb_file))
        chi2 = float(parts[7])
        rg = float(parts[3])
        chi2_values.append(chi2)
        rg_values.append(rg)
        
        # Update min and max values for Rg and Chi^2
        min_rg = min(min_rg, rg)
        max_rg = max(max_rg, rg)
        min_chi2 = min(min_chi2, chi2)
        max_chi2 = max(max_chi2, chi2)

# sort data by pdb_files
data = sorted(zip(pdb_files, chi2_values, rg_values))

pdb_files, chi2_values, rg_values = zip(*data)

sns.set_palette('deep')
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

x = range(len(pdb_files))
ax1.set_xlabel('Time (nanoseconds)')

ax1.set_ylim([min_chi2, max_chi2])
ax1.set_ylabel('\u03C7\u00B2')
line1, = ax1.plot(x, chi2_values, color='orange', label='\u03C7\u00B2')
ax2.set_ylim([min_rg, max_rg])
ax2.set_ylabel('Rg')
line2, = ax2.plot(x, rg_values, color='royalblue', label='Rg')

ax1.legend(handles=[line1, line2], loc='upper left')

fig.savefig('chi2_vs_rg.svg', bbox_inches='tight')
fig.savefig('chi2_vs_rg.pdf')


import re
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
from sklearn.decomposition import PCA

# Step 1: Read chi2 values from the summary file
with open('crysol_summary.txt', 'r') as f:
    chi2_dict = {}
    for line in f:
        match = re.search(r'Model: (.+)\s+Rg: .+ Chi\^2:\s+([\d.]+)', line)
        if match:
            model_name = match.group(1)
            chi2 = float(match.group(2))
            chi2_dict[model_name] = chi2

# Step 2: Sort the chi2 values and write them to a file
sorted_keys = sorted(chi2_dict.keys(), key=lambda x: int(re.search(r'\d+', x).group()))
chi2_values = []
with open('chi2_scores.txt', 'w') as f:
    for key in sorted_keys:
        chi2_value = chi2_dict[key]
        chi2_values.append(chi2_value)  # Store chi2 values for later use
        f.write(f'{key} {chi2_value:.2f}\n')

# Convert chi2_values to numpy array
chi2_values = np.array(chi2_values)

# Step 3: Load and prepare the MD trajectory
traj = md.load('../../combined_aligned.xtc', top='../../final.pdb')
ca_atoms = traj.topology.select('name CA')
ca_traj = traj.xyz[:, ca_atoms, :]
ca_traj_reshaped = ca_traj.reshape(ca_traj.shape[0], ca_traj.shape[1] * 3)


#PC1 vs PC2

pca = PCA(n_components=3)
pca_traj = pca.fit_transform(ca_traj_reshaped)
fig, ax = plt.subplots(figsize=(12, 12))

# Scatter plot using reversed chi2 values as color (mapped to a reversed gradient)
sc = ax.scatter(pca_traj[:, 0], pca_traj[:, 1], c=chi2_values, cmap='viridis_r', alpha=0.75)

# color bar
cbar = plt.colorbar(sc)
cbar.set_label('$\chi^2$ Value', fontsize=22)
cbar.ax.tick_params(labelsize=18)

ax.set_xlabel('PC1 ({:.2f}%)'.format(pca.explained_variance_ratio_[0] * 100), fontsize=22)
ax.set_ylabel('PC2 ({:.2f}%)'.format(pca.explained_variance_ratio_[1] * 100), fontsize=22)

# grid and ticks formatting
ax.grid(which='major', color='grey', linewidth=0.5, alpha=0.2)
ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True, top=False, labeltop=False)
ax.tick_params(axis='y', which='major', left=True, labelleft=True, right=False, labelright=False)
x_ticks = ax.get_xticks()
y_ticks = ax.get_yticks()
ax.set_xticklabels(['{:,.1f}'.format(x) for x in x_ticks], fontsize=18)
ax.set_yticklabels(['{:,.1f}'.format(y) for y in y_ticks], fontsize=18)

plt.subplots_adjust(left=0.12, hspace=0.05, wspace=0.05)
plt.savefig('PC1vsPC2.svg', bbox_inches='tight')
plt.savefig('PC1vsPC2.pdf')


#PC1 vs PC3

pca = PCA(n_components=3)
pca_traj = pca.fit_transform(ca_traj_reshaped)
fig, ax = plt.subplots(figsize=(12, 12))

sc = ax.scatter(pca_traj[:, 0], pca_traj[:, 2], c=chi2_values, cmap='viridis_r', alpha=0.75)

cbar = plt.colorbar(sc)
cbar.set_label('$\chi^2$ Value', fontsize=22)
cbar.ax.tick_params(labelsize=18)

ax.set_xlabel('PC1 ({:.2f}%)'.format(pca.explained_variance_ratio_[0] * 100), fontsize=22)
ax.set_ylabel('PC3 ({:.2f}%)'.format(pca.explained_variance_ratio_[2] * 100), fontsize=22)

ax.grid(which='major', color='grey', linewidth=0.5, alpha=0.2)
ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True, top=False, labeltop=False)
ax.tick_params(axis='y', which='major', left=True, labelleft=True, right=False, labelright=False)
x_ticks = ax.get_xticks()
y_ticks = ax.get_yticks()
ax.set_xticklabels(['{:,.1f}'.format(x) for x in x_ticks], fontsize=18)
ax.set_yticklabels(['{:,.1f}'.format(y) for y in y_ticks], fontsize=18)

plt.subplots_adjust(left=0.12, hspace=0.05, wspace=0.05)
plt.savefig('PC1vsPC3.svg', bbox_inches='tight')
plt.savefig('PC1vsPC3.pdf')


#PC2 vs PC3

pca = PCA(n_components=3)
pca_traj = pca.fit_transform(ca_traj_reshaped)
fig, ax = plt.subplots(figsize=(12, 12))

sc = ax.scatter(pca_traj[:, 1], pca_traj[:, 2], c=chi2_values, cmap='viridis_r', alpha=0.75)

cbar = plt.colorbar(sc)
cbar.set_label('$\chi^2$ Value', fontsize=22)
cbar.ax.tick_params(labelsize=18)

ax.set_xlabel('PC2 ({:.2f}%)'.format(pca.explained_variance_ratio_[1] * 100), fontsize=22)
ax.set_ylabel('PC3 ({:.2f}%)'.format(pca.explained_variance_ratio_[2] * 100), fontsize=22)

ax.grid(which='major', color='grey', linewidth=0.5, alpha=0.2)
ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True, top=False, labeltop=False)
ax.tick_params(axis='y', which='major', left=True, labelleft=True, right=False, labelright=False)
x_ticks = ax.get_xticks()
y_ticks = ax.get_yticks()
ax.set_xticklabels(['{:,.1f}'.format(x) for x in x_ticks], fontsize=18)
ax.set_yticklabels(['{:,.1f}'.format(y) for y in y_ticks], fontsize=18)

plt.subplots_adjust(left=0.12, hspace=0.05, wspace=0.05)
plt.savefig('PC2vsPC3.svg', bbox_inches='tight')
plt.savefig('PC2vsPC3.pdf')

#output motions captured by each PC
#project traj onto each PC to visualize the motion along PC
pc_vectors = pca.components_.reshape(3, ca_traj.shape[1], 3)

#save PC motion in PDB format 
def save_pc_motion(pc_vector, pc_label, filename, n_frames=20, scale_factor=10):
    mean_structure = np.mean(ca_traj, axis=0)

    #interpolate between the mean and displaced structures
    frames = []
    ca_topology = traj.top.subset(ca_atoms)

    #forward motion frames
    for i in range(n_frames):
        displacement = scale_factor * ((i / (n_frames - 1)) - 0.5) * pc_vector
        interpolated_structure = mean_structure + displacement
        frames.append(interpolated_structure)

    #backward motion frames 
    for i in range(n_frames - 1, -1, -1):  # Reverse iteration
        frames.append(frames[i])  # Append the reverse frame

    motion_traj = md.Trajectory(np.array(frames), ca_topology)

    motion_traj.save(f'{filename}_motion.pdb')

save_pc_motion(pc_vectors[0], 'PC1', 'PC1_motion')
save_pc_motion(pc_vectors[1], 'PC2', 'PC2_motion')
save_pc_motion(pc_vectors[2], 'PC3', 'PC3_motion')




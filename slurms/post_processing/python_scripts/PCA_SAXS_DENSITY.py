#!/bin/bash

	#####PCA + SAXS COLOURING. EXTRACTS CRYSOL DATA, CONVERTS CHI^2 TO COLOURS AND RUNS PCA WITH FRAMES COLOURES BY SAXS FIT + DENSITY DISTRIBITONS AND SCREE#####

#####To do#####
#Expand plotting to three dimensions 
#Incorporate clustering either into this script (e.g. with DBSCAN) or incorporate this script into CLoNe
#Visualise the principal motions

import os
import re
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import mdtraj as md
from sklearn.decomposition import PCA
import warnings 
from scipy.signal import savgol_filter

#temporary housekeeping
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

									#####extract crysol data#####

with open('crysol_summary.txt', 'r') as f:
    #dict to stor Chi^2 & model 
    chi2_dict = {}
    for line in f:
        #search for model & Chi^2 and add to dict 
        match = re.search(r'Model: (.+)\s+Rg: .+ Chi\^2:\s+([\d.]+)', line)
        if match:
            model_name = match.group(1)
            chi2 = float(match.group(2))
            chi2_dict[model_name] = chi2

#sort dict in numerical order, loop through and add to chi2_scores.txt
sorted_keys = sorted(chi2_dict.keys(), key=lambda x: int(re.search(r'\d+', x).group()))

with open('chi2_scores.txt', 'w') as f:
    for key in sorted_keys:
        f.write(f'{key} {chi2_dict[key]:.2f}\n')

									######chi^2 to colour######

#using american spelling (colors) for seaborn compatability
colors = []
with open('chi2_scores.txt', 'r') as f:
    for line in f:
        filename, chi2 = line.strip().split()
        chi2 = float(chi2)
        if 0.8 > chi2:
            colors.append('pink')
        elif 0.8 < chi2 <=2:
            colors.append('green')
        elif 2.0 < chi2 <= 3:
            colors.append('blue')
        elif 3 < chi2 <= 4:
            colors.append('purple')
        elif 4 < chi2 <= 5:
            colors.append('red')
        elif chi2 > 5:
            colors.append('black')

with open('colours.txt', 'w') as f:
    for color in colors:
        f.write(color + '\n')


#convert to RBG colouring
with open('colours.txt', 'r') as f:
    colors = [line.strip() for line in f.readlines()]

colors_rgb = sns.color_palette(colors, n_colors=len(colors))

										#####PCA#####

#PC VISUALISATION
# Load traj + topo

traj = md.load('combined_aligned.xtc', top='final.pdb')


ca_atoms = traj.topology.select('name CA')

ca_traj = traj.xyz[:, ca_atoms, :]

ca_traj_reshaped = ca_traj.reshape(ca_traj.shape[0], ca_traj.shape[1] * 3)



# Reshape trajectory for PCA
#traj_reshaped = traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3)

# Perform PCA
pca = PCA(n_components=3)
pca.fit(ca_traj_reshaped)

# Transform the trajectory into the PCA space
pca_traj = pca.transform(ca_traj_reshaped)

# Define the number of frames to select for each principal component
num_frames = 50

# Save the motions of the first three principal components as separate PDB files with a looping effect
for i in range(3):
    principal_motion = pca.components_.T[:, i]
    projection = np.sum(pca_traj[:, i:i+1] * principal_motion[:i+1], axis=1)
    principal_frames = np.argsort(projection)[:num_frames]

    # Duplicate frames in reverse order
    principal_frames_loop = np.concatenate((principal_frames, principal_frames[::-1]))

    selected_frames = traj[principal_frames_loop]
    selected_frames_alpha = selected_frames.atom_slice(selected_frames.topology.select('name CA'))
    selected_frames_alpha.save(f'pc{i+1}_motion_loop.pdb')

# calculate the cosine content for the first principal componentis
first_principal_motion = pca.components_.T[0]
cosine_content_1 = np.dot(first_principal_motion, first_principal_motion)
second_principal_motion = pca.components_.T[1]
cosine_content_2 = np.dot(first_principal_motion, second_principal_motion)
third_principal_motion = pca.components_.T[2]
cosine_content_3 = np.dot(first_principal_motion, third_principal_motion)

# print the cosine content for the first three principal components
print("Cosine content for PC 1 = {:.8f}".format(cosine_content_1))
print("Cosine content for PC 2 = {:.8f}".format(cosine_content_2))
print("Cosine content for PC 3 = {:.8f}".format(cosine_content_3))

cosine_values = np.array([cosine_content_1, cosine_content_2, cosine_content_3])
np.savetxt('cosine_values.txt', cosine_values)

# plot transformed components over time
time = np.arange(pca_traj.shape[0])
plt.plot(time, savgol_filter(pca_traj[:, 0], window_length=51, polyorder=3), label='PC1')
plt.plot(time, savgol_filter(pca_traj[:, 1], window_length=51, polyorder=3), label='PC2')
plt.plot(time, savgol_filter(pca_traj[:, 2], window_length=51, polyorder=3), label='PC3')
plt.xlabel('Time')
plt.ylabel('Transformed Components')
plt.legend()
plt.show()
plt.savefig('transformed_components_over_time.png')



#find dat file for legend
dat_file = next((f for f in os.listdir() if f.endswith('.dat')), None)

###PC2 vs PC1###

#scatter plot
fig, ax = plt.subplots(figsize=(16, 16))
fig.suptitle('MD Frames Projected onto PC2 vs PC1', fontsize=13)
plt.subplots_adjust(top=0.85)
sc = ax.scatter(pca.transform(ca_traj_reshaped)[:, 0],
                pca.transform(ca_traj_reshaped)[:, 1],
                c=colors_rgb, alpha=0.5)
ax.grid(which='major', color='grey', linewidth=0.5, alpha=0.5)
ax.set_xlabel('PC1 ({:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize=14)
ax.set_ylabel('PC2 ({:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize=14)

legend_elements = [
    plt.Line2D([0], [0], marker='', color='w', label=f'{dat_file}', linestyle='', markersize=0),
    plt.Line2D([0], [0], marker='o', color='w', label='0.8 > $\chi^2$', markerfacecolor='pink', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='0.8 < $\chi^2$ < 2.0', markerfacecolor='green', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='2.0 < $\chi^2$ < 3.0', markerfacecolor='blue', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='3.0 < $\chi^2$ < 4.0', markerfacecolor='purple', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='4.0 < $\chi^2$ < 5.0', markerfacecolor='red', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='$\chi^2$ > 5.0', markerfacecolor='black', markersize=10)
]
plt.legend(handles=legend_elements)

#density plots
x = pca.transform(ca_traj_reshaped)[:, 0] #0=PC1
y = pca.transform(ca_traj_reshaped)[:, 1] #1=PC2
color_order = ['green', 'blue', 'purple', 'red', 'black']
for i, pc in enumerate([x, y]):  # iterate through the two PCs and create horizontal and vertical density plots
    if i == 0:  # x-axis PC
        ax_histx = fig.add_axes([ax.get_position().x0, ax.get_position().y1 + 0.01, ax.get_position().width, 0.08], xticklabels=[], sharex=ax)
        ax_histx.spines['top'].set_visible(False)  # remove spines
        ax_histx.spines['right'].set_visible(False)  # remove spines
        for color in color_order:  # iterate over color names
            if color in colors:
                kde = sns.kdeplot(pc[np.where(np.array(colors) == color)[0]], color=colors_rgb[colors.index(color)], fill=True, alpha=0.4, ax=ax_histx)
        ax_histx.spines['left'].set_visible(False)
        ax_histx.set_yticks([])
        ax_histx.tick_params(axis='both', which='both', length=0, bottom=False, labelbottom=False, top=False, labeltop=False)
    else:  # y-axis PC
        ax_histy = fig.add_axes([ax.get_position().x1 + 0.01, ax.get_position().y0, 0.08, ax.get_position().height], yticklabels=[], sharey=ax)
        ax_histy.spines['top'].set_visible(False)
        ax_histy.spines['right'].set_visible(False)
        for color in color_order:  # iterate over color names
            if color in colors:
                kde = sns.kdeplot(pc[np.where(np.array(colors) == color)[0]], color=colors_rgb[colors.index(color)], fill=True, alpha=0.4, ax=ax_histy, vertical=True)
        ax_histy.spines['bottom'].set_visible(False)
        ax_histy.set_xticks([])
        ax_histy.tick_params(axis='both', which='both', length=0, left=False, labelleft=False, right=False, labelright=False)

for ax_hist in [ax_histx, ax_histy]:
    ax_hist.set_xlabel('')
    ax_hist.set_xlabel('')
    ax_hist.set_ylabel('')
    ax_hist.spines['bottom'].set_visible(False)
    ax_hist.spines['left'].set_visible(False)
    ax_hist.spines['right'].set_visible(False)
    ax_hist.spines['top'].set_visible(False)

ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True, top=False, labeltop=False)
ax.tick_params(axis='y', which='major', left=True, labelleft=True, right=False, labelright=False)
x_ticks = ax.get_xticks()
y_ticks = ax.get_yticks()
ax.set_xticklabels(['{:,.1f}'.format(x) for x in x_ticks])
ax.set_yticklabels(['{:,.1f}'.format(y) for y in y_ticks])

plt.subplots_adjust(left=0.05, hspace=0.05, wspace=0.05)
plt.savefig('PC2vsPC1.png', dpi=600, bbox_inches='tight')
plt.show()

###PC3 vs PC1###
fig, ax = plt.subplots(figsize=(16, 16))
fig.suptitle('MD Frames Projected onto PC3 vs PC1', fontsize=13)
plt.subplots_adjust(top=0.85)
sc = ax.scatter(pca.transform(ca_traj_reshaped)[:, 0],
                pca.transform(ca_traj_reshaped)[:, 2],
                c=colors_rgb, alpha=0.5)
ax.grid(which='major', color='grey', linewidth=0.5, alpha=0.5)
ax.set_xlabel('PC1 ({:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize=14)
ax.set_ylabel('PC3 ({:.2f}%)'.format(pca.explained_variance_ratio_[2]*100), fontsize=14)


legend_elements = [
    plt.Line2D([0], [0], marker='', color='w', label=f'{dat_file}', linestyle='', markersize=0),
    plt.Line2D([0], [0], marker='o', color='w', label='0.8 > $\chi^2$', markerfacecolor='pink', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='0.8 < $\chi^2$ < 2.0', markerfacecolor='green', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='2.0 < $\chi^2$ < 3.0', markerfacecolor='blue', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='3.0 < $\chi^2$ < 4.0', markerfacecolor='purple', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='4.0 < $\chi^2$ < 5.0', markerfacecolor='red', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='$\chi^2$ > 5.0', markerfacecolor='black', markersize=10)
]
plt.legend(handles=legend_elements)

#denity plots
x = pca.transform(ca_traj_reshaped)[:, 0] #PC1
y = pca.transform(ca_traj_reshaped)[:, 2] #PC3
color_order = ['green', 'blue', 'purple', 'red', 'black']
for i, pc in enumerate([x, y]):  # iterate through the two PCs and create horizontal and vertical density plots
    if i == 0:  # x-axis PC
        ax_histx = fig.add_axes([ax.get_position().x0, ax.get_position().y1 + 0.01, ax.get_position().width, 0.08], xticklabels=[], sharex=ax)
        ax_histx.spines['top'].set_visible(False)  # remove spines
        ax_histx.spines['right'].set_visible(False)  # remove spines
        for color in color_order:  # iterate over color names
            if color in colors:
                kde = sns.kdeplot(pc[np.where(np.array(colors) == color)[0]], color=colors_rgb[colors.index(color)], fill=True, alpha=0.4, ax=ax_histx)
        ax_histx.spines['left'].set_visible(False)
        ax_histx.set_yticks([])
        ax_histx.tick_params(axis='both', which='both', length=0, bottom=False, labelbottom=False, top=False, labeltop=False)
    else:  # y-axis PC
        ax_histy = fig.add_axes([ax.get_position().x1 + 0.01, ax.get_position().y0, 0.08, ax.get_position().height], yticklabels=[], sharey=ax)
        ax_histy.spines['top'].set_visible(False)
        ax_histy.spines['right'].set_visible(False)
        for color in color_order:  # iterate over color names
            if color in colors:
                kde = sns.kdeplot(pc[np.where(np.array(colors) == color)[0]], color=colors_rgb[colors.index(color)], fill=True, alpha=0.4, ax=ax_histy, vertical=True)
        ax_histy.spines['bottom'].set_visible(False)
        ax_histy.set_xticks([])
        ax_histy.tick_params(axis='both', which='both', length=0, left=False, labelleft=False, right=False, labelright=False)

for ax_hist in [ax_histx, ax_histy]:
    ax_hist.set_xlabel('')
    ax_hist.set_xlabel('')
    ax_hist.set_ylabel('')
    ax_hist.spines['bottom'].set_visible(False)
    ax_hist.spines['left'].set_visible(False)
    ax_hist.spines['right'].set_visible(False)
    ax_hist.spines['top'].set_visible(False)

ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True, top=False, labeltop=False)
ax.tick_params(axis='y', which='major', left=True, labelleft=True, right=False, labelright=False)
x_ticks = ax.get_xticks()
y_ticks = ax.get_yticks()
ax.set_xticklabels(['{:,.1f}'.format(x) for x in x_ticks])
ax.set_yticklabels(['{:,.1f}'.format(y) for y in y_ticks])

plt.subplots_adjust(left=0.05, hspace=0.05, wspace=0.05)
plt.savefig('PC3vsPC1.png', dpi=600, bbox_inches='tight')
plt.show()


###PC3 vs PC2###
fig, ax = plt.subplots(figsize=(16, 16))
fig.suptitle('MD Frames Projected onto PC3 vs PC2', fontsize=13)
plt.subplots_adjust(top=0.85)
sc = ax.scatter(pca.transform(ca_traj_reshaped)[:, 1],
                pca.transform(ca_traj_reshaped)[:, 2],
                c=colors_rgb, alpha=0.5)
ax.grid(which='major', color='grey', linewidth=0.5, alpha=0.5)
ax.set_xlabel('PC2 ({:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize=14)
ax.set_ylabel('PC3 ({:.2f}%)'.format(pca.explained_variance_ratio_[2]*100), fontsize=14)

legend_elements = [
    plt.Line2D([0], [0], marker='', color='w', label=f'{dat_file}', linestyle='', markersize=0),
    plt.Line2D([0], [0], marker='o', color='w', label='0.8 > $\chi^2$', markerfacecolor='pink', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='0.8 < $\chi^2$ < 2.0', markerfacecolor='green', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='2.0 < $\chi^2$ < 3.0', markerfacecolor='blue', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='3.0 < $\chi^2$ < 4.0', markerfacecolor='purple', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='4.0 < $\chi^2$ < 5.0', markerfacecolor='red', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='$\chi^2$ > 5.0', markerfacecolor='black', markersize=10)
]
plt.legend(handles=legend_elements)

#density plots
x = pca.transform(ca_traj_reshaped)[:, 1] #PC2
y = pca.transform(ca_traj_reshaped)[:, 2] #PC3
color_order = ['green', 'blue', 'purple', 'red', 'black']
for i, pc in enumerate([x, y]):  # iterate through the two PCs and create horizontal and vertical density plots
    if i == 0:  # x-axis PC
        ax_histx = fig.add_axes([ax.get_position().x0, ax.get_position().y1 + 0.01, ax.get_position().width, 0.08], xticklabels=[], sharex=ax)
        ax_histx.spines['top'].set_visible(False)  # remove spines
        ax_histx.spines['right'].set_visible(False)  # remove spines
        for color in color_order:  # iterate over color names
            if color in colors:
                kde = sns.kdeplot(pc[np.where(np.array(colors) == color)[0]], color=colors_rgb[colors.index(color)], fill=True, alpha=0.4, ax=ax_histx)
        ax_histx.spines['left'].set_visible(False)
        ax_histx.set_yticks([])
        ax_histx.tick_params(axis='both', which='both', length=0, bottom=False, labelbottom=False, top=False, labeltop=False)
    else:  # y-axis PC
        ax_histy = fig.add_axes([ax.get_position().x1 + 0.01, ax.get_position().y0, 0.08, ax.get_position().height], yticklabels=[], sharey=ax)
        ax_histy.spines['top'].set_visible(False)
        ax_histy.spines['right'].set_visible(False)
        for color in color_order:  # iterate over color names
            if color in colors:
                kde = sns.kdeplot(pc[np.where(np.array(colors) == color)[0]], color=colors_rgb[colors.index(color)], fill=True, alpha=0.4, ax=ax_histy, vertical=True)
        ax_histy.spines['bottom'].set_visible(False)
        ax_histy.set_xticks([])
        ax_histy.tick_params(axis='both', which='both', length=0, left=False, labelleft=False, right=False, labelright=False)

for ax_hist in [ax_histx, ax_histy]:
    ax_hist.set_xlabel('')
    ax_hist.set_xlabel('')
    ax_hist.set_ylabel('')
    ax_hist.spines['bottom'].set_visible(False)
    ax_hist.spines['left'].set_visible(False)
    ax_hist.spines['right'].set_visible(False)
    ax_hist.spines['top'].set_visible(False)

ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True, top=False, labeltop=False)
ax.tick_params(axis='y', which='major', left=True, labelleft=True, right=False, labelright=False)
x_ticks = ax.get_xticks()
y_ticks = ax.get_yticks()
ax.set_xticklabels(['{:,.1f}'.format(x) for x in x_ticks])
ax.set_yticklabels(['{:,.1f}'.format(y) for y in y_ticks])

plt.subplots_adjust(left=0.05, hspace=0.05, wspace=0.05)
plt.savefig('PC3vsPC2.png', dpi=600, bbox_inches='tight')
plt.show()


# Calculate cumulative explained variance ratio
cumulative_variance_ratio = np.cumsum(pca.explained_variance_ratio_)

# Plot cumulative explained variance ratio
fig, ax = plt.subplots()
ax.plot(range(1, pca.n_components_+1), cumulative_variance_ratio, '-o')
ax.set_xlabel('Principal Components')
ax.set_ylabel('Cumulative Variance Explained (%)')

ax.set_xticks(range(1, pca.n_components_+1))

plt.savefig('cumulative_variance.png', dpi=600, bbox_inches='tight')
plt.show()




#plot scree
#fig, ax = plt.subplots()
#ax.plot(range(1, pca.n_components_+1), pca.explained_variance_ratio_, '-o')
#ax.set_xlabel('Principal Component')
#ax.set_ylabel('Variance Explained (%)')
#plt.savefig('scree.png', dpi=600, bbox_inches='tight')
#plt.show()
#


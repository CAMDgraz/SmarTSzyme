#!/bin/python

"""
@autors: Daniel Platero Rochart [daniel.platero-rochart@medunigraz.at]
"""

from evcouplings.couplings import CouplingsModel
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('./images.mplstyle')

positions = [355, 373, 446]
model = CouplingsModel('mhet.model')

single_mut = model.smm(positions)

# Plot single_mut
min_val = single_mut.min()
max_val = single_mut.max()
cbar_limit = abs(min_val) if abs(min_val) > abs(max_val) else abs(max_val)    

fig, ax = plt.subplots()
sns.heatmap(single_mut.T, ax=ax, center=0, cmap='coolwarm_r', vmin=-cbar_limit,
            vmax=cbar_limit)
ax.set_ylabel(r'Amino Acids')
ax.set_yticks(np.arange(0.5, len(model.alphabet) + 0.5, 1))
ax.set_yticklabels(model.alphabet, rotation=0)

ax.set_xticks(np.arange(0.5, len(positions) + 0.5, 1))
ax.set_xticklabels(positions)
ax.set_xlabel(r'Positions')

# Get original amino acids
for pos_idx, pos in enumerate(positions):
    original_aa = np.where(model.alphabet == model.seq(pos))[0]
    y_pos = np.arange(0.5, len(model.alphabet) + 0.5, 1)[original_aa]
    x_pos = np.arange(0.5, len(positions) + 0.5, 1)[pos_idx]
    ax.scatter(x_pos, y_pos, color='black', marker='o')
fig.savefig('sm_mhet_improving_papers.png')


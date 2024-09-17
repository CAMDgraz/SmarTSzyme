#!/bin/python

"""
@autors: Daniel Platero Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sánchez-Murcia [pedro.murcia@medunigraz.at]
"""

"""
Reduce:
Selection of key residues based on the structure. The method implemented is 
based on the work of Amor et. al. (doi.org/10.1038/ncomms12477). 
"""

import pickle
import heapq
import numpy as np
import matplotlib.pyplot as plt

import reduce_functions as rf
import extra

print("""
********************************************************************************
* SmarTSzyme-reduce:                                                           *
*      Selection of important residues for enzyme engineering                  *
********************************************************************************
""")

# Reading fluxes  and calculating overall contribution =========================
flux_stab = rf.load_pickle('./test/test_dir/flux_stabilizing.pickle')
flux_destab = rf.load_pickle('./test/test_dir/flux_destabilizing.pickle')

flux_total = flux_destab - flux_stab
# ==============================================================================

# Plot results =================================================================
flux_positive = np.copy(flux_total)
flux_positive[np.where(flux_positive < 0)[0]] = 0

flux_negative = np.copy(flux_total)
flux_negative[np.where(flux_negative > 0)[0]] = 0

extra.mplstyle()

# Bar plot divided
fig, ax = plt.subplots()
ax.bar(np.arange(len(flux_total)), flux_positive, color='red',
                 label='Destabilizing')
ax.bar(np.arange(len(flux_total)), flux_negative, color='blue',
                 label='Stabilizing')
ax.set_xticks(np.arange(1, len(flux_total) + 1, 1))
ax.set_xlabel(r'Residues')
ax.set_ylabel(r'Normalized Flux')
fig.savefig('flux_bar.png')
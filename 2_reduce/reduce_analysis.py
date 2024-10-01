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

import heapq
import numpy as np
import argparse
import matplotlib.pyplot as plt

import reduce_functions as rf
import extra

def parse_arguments():
    """
    Parse arguments of the cli
    """
    desc = '''\nReduce: Identification of key residues in order to reduce the 
              the mutational landscape'''
    parser = argparse.ArgumentParser(prog='Reduce',
                                     description = desc,
                                     add_help=True,
                                     allow_abbrev = False)
    inputs = parser.add_argument_group(title='Input options')
    inputs.add_argument('-maxres', dest = 'max_residues', action = 'store', 
                      help = 'Top n residues to select.',
                      type=int, required = True)
    inputs.add_argument('-offset', dest = 'offset', action = 'store',
                        help = 'Offset for the first residue', default = 0,
                        type=int)
    user_inputs = parser.parse_args()
    return user_inputs

args = parse_arguments()

# Reading fluxes  and calculating overall contribution =========================
flux_stab = rf.load_pickle('./flux_stabilizing.pickle')
flux_destab = rf.load_pickle('./flux_destabilizing.pickle')

flux_total = flux_destab - flux_stab

del flux_destab
del flux_stab

flux_positive = np.copy(flux_total)
flux_positive[np.where(flux_positive < 0)[0]] = 0

flux_negative = np.copy(flux_total)
flux_negative[np.where(flux_negative > 0)[0]] = 0

# ==============================================================================

# Write results to csv =========================================================
with open('results.csv', 'w') as f:
    f.write('residue,flux\n')
    for resid, flux in enumerate(flux_total):
        f.write(f'{resid + 1},{flux:.4f}\n')

del flux_total

# Plot results =================================================================
extra.mplstyle()

# Bar plot divided
fig, ax = plt.subplots()
ax.bar(np.arange(1, len(flux_positive) + 1), flux_positive, color='red',
                 label='Destabilizing')
ax.bar(np.arange(1, len(flux_negative) + 1), flux_negative, color='blue',
                 label='Stabilizing')
ax.set_xlabel(r'Residues')
ax.set_ylabel(r'Normalized Flux')
fig.savefig('flux_bar.png')

# Select residues for next step =================================================
ordered_positive = []
ordered_negative = []

heapq.heapify(ordered_negative)
heapq.heapify(ordered_positive)

for res_id, flux in enumerate(flux_negative):
    heapq.heappush(ordered_negative, (flux, res_id + 1))
nnegative = heapq.nsmallest(args.max_residues, ordered_negative)
nnegative_res = np.asarray([neg[1] for neg in nnegative]) + args.offset
nnegative_flux = np.asarray([neg[0] for neg in nnegative])

for res_id, flux in enumerate(flux_positive):
    heapq.heappush(ordered_positive, (flux, res_id + 1))
npositive = heapq.nlargest(args.max_residues, ordered_positive)
npositive_res = np.asarray([pos[1] for pos in npositive]) + args.offset
npositive_flux = np.asarray([pos[0] for pos in npositive])

fig, ax = plt.subplots()
ax.bar(np.arange(len(nnegative_res)), -nnegative_flux, color='blue')
ax.set_xlabel(r'Residues')
ax.set_xticks(np.arange(len(nnegative_res)))
ax.set_xticklabels(nnegative_res, rotation=45, size=10)
ax.set_ylabel(r'|Normalized flux|')
fig.savefig('top_negative_flux.png')

fig, ax = plt.subplots()
ax.bar(np.arange(len(npositive_res)), npositive_flux, color='red')
ax.set_xlabel(r'Residues')
ax.set_xticks(np.arange(len(npositive_res)))
ax.set_xticklabels(npositive_res, rotation=45, size=10)
ax.set_ylabel(r'|Normalized flux|')
fig.savefig('top_positive_flux.png')
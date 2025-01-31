#!/bin/python
# -*- coding: utf-8 -*-

"""
@authors: Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sánchez-Murcia [pedro.murcia@medunigraz.at]
"""

"""
Reduce:
Selection of key residues based on the structure. The method implemented is 
based on the work of Amor et. al. (doi.org/10.1038/ncomms12477). 
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import reduce_functions as rf
import mdtraj as md

print("""
********************************************************************************
* SmarTSzyme-reduce:                                                           *
*      Selection of important residues for enzyme engineering                  *
********************************************************************************
""")

# Parse arguments ==============================================================
args = rf.parse_arguments()
catalytic_residues = np.asarray(args.catalytic_residues, dtype=int) - 1

import os
import sys
import shutil
if args.output == './':
    print('Output path can not be ./')
    sys.exit()
try:
    os.mkdir(args.output)
    os.mkdir(f'{args.output}/matrices/')
except FileExistsError:
    if args.force:
        shutil.rmtree(args.output, ignore_errors=True)
        os.mkdir(args.output)
        os.mkdir(f'{args.output}/matrices/')
    else:
        raise FileExistsError
# ==============================================================================

# Read list of qmmm jobs =======================================================
print(f'\nReading QMMM jobs from: {args.qmmm_list}')
jobs_df = pd.read_csv(args.qmmm_list, header=None)
jobs = np.asarray(jobs_df.iloc[:, 0])
print(f'Number of trajectories: {len(jobs)}')
# ==============================================================================

# Get TS indices ===============================================================
qmmm_type = 'smd'
traj = md.load(f'{jobs[0]}/traj_{args.sufix}.nc', top=f'{jobs[0]}/top_{args.sufix}.parm7')
nframes = traj.n_frames
del traj

if qmmm_type == 'smd':
    ts_indices = np.zeros(len(jobs), dtype=np.int16)
    for job_idx, job in enumerate(jobs):
        ts_index, work_lines = rf.identify_smd_TS(f'{job}/smd_{args.sufix}.txt') # args
        factor = int(work_lines/nframes)
        if ts_index % factor != 0:
            ts_index = int(ts_index//factor + 1)
        else:
            ts_index = int(ts_index//factor)
        ts_indices[job_idx] = ts_index
    with open(f'{args.output}/ts_indexes.dat', 'w') as f:
        for job, ts in zip(jobs, ts_indices):
            f.write(f'{job},{ts}\n')
# ==============================================================================

# Calculate matrices ===========================================================
interactions = ['vdw', 'hbonds', 'coulomb']

# Parallel calculation
print(f'Calculating interactions for Enzyme-Substrate complex:')
for interaction in interactions:
    print(f'    * calculating {interaction}')
    matrix_int = np.zeros((args.nresidues, args.nresidues))
    arg1 = [interaction for i in range(len(jobs))]
    arg2 = [f'{job}/traj_{args.sufix}.nc' for job in jobs]
    arg3 = [f'{job}/top_{args.sufix}.parm7' for job in jobs]
    arg4 = [args.cutoff/10 for i in range(len(jobs))]
    with mp.Pool(processes=args.ncpus) as pool:
        results = pool.starmap(rf.calculate_matrix, zip(arg1, arg2, arg3,
                                                              arg4))
    for result in results:
        matrix_int += result
    del results # clean memory
    rf.write_pickle(matrix_int,
                    f'{args.output}/matrices/{interaction}_es.pickle')
del matrix_int

print(f'Calculating interactions for (pseudo) Transition State complex:')
for interaction in interactions:
    print(f'    * calculating {interaction}')
    matrix_int = np.zeros((args.nresidues, args.nresidues))
    arg1 = [interaction for i in range(len(jobs))]
    arg2 = [f'{job}/traj_{args.sufix}.nc' for job in jobs]
    arg3 = [f'{job}/top_{args.sufix}.parm7' for job in jobs]
    arg4 = [args.cutoff/10 for i in range(len(jobs))]
    with mp.Pool(processes=args.ncpus) as pool:
        results = pool.starmap(rf.calculate_matrix, zip(arg1, arg2, arg3,
                                                         arg4, ts_indices-1))
        # -1 to math the 0-based index
    for result in results:
        matrix_int += result
    del results # clean memory
    rf.write_pickle(matrix_int,
                    f'{args.output}/matrices/{interaction}_pts.pickle')
del matrix_int
# ==============================================================================

# Read matrices one by one =====================================================
print(f'\nReading matrices ...')
matrix_es = np.zeros((args.nresidues, args.nresidues))
for interaction in ['vdw', 'hbonds', 'coulomb']:
    matrix_es += rf.load_pickle(
                            f'{args.output}/matrices/{interaction}_es.pickle')
matrix_es /= len(jobs)

matrix_pts = np.zeros((args.nresidues, args.nresidues))
for interaction in ['vdw', 'hbonds', 'coulomb']:
    matrix_pts += rf.load_pickle(
                            f'{args.output}/matrices/{interaction}_pts.pickle')
matrix_pts /= len(jobs)
# ==============================================================================

# Separate contributions =======================================================
matrix_diff = matrix_pts - matrix_es
del matrix_pts
del matrix_es

stabilizing = np.copy(matrix_diff)
stabilizing[np.where(stabilizing > 0)] = 0
stabilizing[catalytic_residues] = matrix_diff[catalytic_residues]
stabilizing.T[catalytic_residues] = matrix_diff.T[catalytic_residues]

destabilizing = np.copy(matrix_diff)
destabilizing[np.where(destabilizing < 0)] = 0
destabilizing[catalytic_residues] = matrix_diff[catalytic_residues]
destabilizing.T[catalytic_residues] = matrix_diff.T[catalytic_residues]
# ==============================================================================

# Calculate edge-to-edge propensity ============================================
print('\nCalculating edge_to_edge for stabilizing residues ...')
edges_catal = np.asarray([])
for catalytic_residue in catalytic_residues:
    edges_catal = np.concatenate((edges_catal,
                                  np.where(stabilizing.nonzero()[0] ==
                                          catalytic_residue)[0]))
edge_to_edge = rf.edge_transfer_matrix(stabilizing, args.nresidues)
for diag in range(edge_to_edge.shape[0]):
    edge_to_edge[diag][diag] = 0

flux_stabilizing = np.zeros(args.nresidues)
for residue in range(args.nresidues):
    if residue in catalytic_residues:
        continue
    edges_residue = np.where(stabilizing.nonzero()[0] == residue)[0]
    edges_residue_repeat = np.repeat(edges_residue, len(edges_catal))
    edges_residue_repeat = edges_residue_repeat.astype(np.int16)
    edges_catal_tile = np.tile(edges_catal, len(edges_residue))
    edges_catal_tile = edges_catal_tile.astype(np.int16)
    flux_stabilizing[residue] = np.sum(abs(edge_to_edge[edges_residue_repeat,
                                                        edges_catal_tile]))

# Normalize fluxes and write results
flux_stabilizing /= np.sum(flux_stabilizing)
rf.write_pickle(flux_stabilizing, f'{args.output}/flux_stabilizing.pickle')

print('Calculating edge_to_edge for destabilizing residues ...')
edges_catal = np.asarray([])
for catalytic_residue in catalytic_residues:
    edges_catal = np.concatenate((edges_catal,
                                  np.where(destabilizing.nonzero()[0] ==
                                          catalytic_residue)[0]))
edge_to_edge = rf.edge_transfer_matrix(destabilizing, args.nresidues)
for diag in range(edge_to_edge.shape[0]):
    edge_to_edge[diag][diag] = 0


flux_destabilizing = np.zeros(args.nresidues)
for residue in range(args.nresidues):
    if residue in catalytic_residues:
        continue
    edges_residue = np.where(destabilizing.nonzero()[0] == residue)[0]
    edges_catal_tile = np.tile(edges_catal, len(edges_residue))
    edges_catal_tile = edges_catal_tile.astype(np.int16)
    edges_residue_repeat = np.repeat(edges_residue, len(edges_catal))
    edges_residue_repeat = edges_residue_repeat.astype(np.int16)
    flux_destabilizing[residue] = np.sum(abs(edge_to_edge[edges_residue_repeat,
                                                          edges_catal_tile]))

# Normalize fluxes and write results
flux_destabilizing /= np.sum(flux_destabilizing)
rf.write_pickle(flux_destabilizing, f'{args.output}/flux_destabilizing.pickle')
# ==============================================================================

# Plot and save results ========================================================
flux_total = flux_destabilizing - flux_stabilizing

del flux_destabilizing
del flux_stabilizing

flux_positive = np.copy(flux_total)
flux_positive[np.where(flux_positive < 0)[0]] = 0

flux_negative = np.copy(flux_total)
flux_negative[np.where(flux_negative > 0)[0]] = 0

# Write results to csv
with open(f'{args.output}/results.csv', 'w') as f:
    f.write('residue,flux\n')
    for resid, flux in enumerate(flux_total):
        f.write(f'{resid + 1},{flux:.4f}\n')

del flux_total

rf.mplstyle()

# Bar plot divided
fig, ax = plt.subplots()
ax.bar(np.arange(1, len(flux_positive) + 1), flux_positive, color='red',
                 label='Destabilizing')
ax.bar(np.arange(1, len(flux_negative) + 1), flux_negative, color='blue',
                 label='Stabilizing')
ax.set_xlabel(r'Residues')
ax.set_ylabel(r'Normalized Flux')
fig.savefig(f'{args.output}/flux_bar.png')

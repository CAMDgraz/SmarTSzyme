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
import os
import sys
import shutil

print("""
********************************************************************************
* SmarTSzyme-reduce:                                                           *
*      Selection of key residues for enzyme engineering                  *
********************************************************************************
""")

# Handle arguments =============================================================
args = rf.parse_arguments()

# Check if qmmm_list or interactions_list provided
if args.qmmm_list and args.interactions_list:
    print('\nError!!! Provide either qmmm_list or interactions_list.')
    sys.exit()
elif not args.qmmm_list and not args.interactions_list:
    print('\nError!!! None qmmm_list or interactions_list provided.')
    sys.exit()

# Check the output argument and create the file tree
if args.output == './':
    print('Error!!! Output path can not be ./')
    sys.exit()
try:
    os.mkdir(args.output)
    os.mkdir(f'{args.output}/matrices/')
    os.mkdir(f'{args.output}/coupling/')
except FileExistsError:
    if args.force: # Delete output folder if -f
        shutil.rmtree(args.output, ignore_errors=True)
        os.mkdir(args.output)
        os.mkdir(f'{args.output}/matrices/')
        os.mkdir(f'{args.output}/coupling/')
    else:
        print('Error!!! Output file exist. Use -f True to overwrite.')
        raise FileExistsError

catalytic_residues = np.asarray(args.catalytic_residues, dtype=int) - 1
interactions = ['vdw', 'hbonds', 'coulomb']
jobs_df = pd.read_csv(args.qmmm_list, header=None)
jobs = np.asarray(jobs_df.iloc[:, 0])
# ==============================================================================

# Input options ================================================================
print(f'\n*** User inputs ***')
if args.qmmm_list:
    only_flux = False
    print(f'\nReading QMMM jobs from:          {args.qmmm_list}')
elif args.interactions_list:
    only_flux = True
    print(f'\nReading interactions paths from: {args.interactions_list}')
print(f'\nNumber of residues:              {args.nres}')
print(f'\nCatalytic residues:              {args.catalytic_residues}')
print(f'\nNumber of trajectories:          {len(jobs)}')
print(f'\nNumber of CPUs:                  {args.ncpus}')
# ==============================================================================

# Get TS indices ===============================================================
if not only_flux:
    # Get number of frames
    traj = md.load(f'{jobs[0]}/traj_{args.sufix}.nc',
                top=f'{jobs[0]}/top_{args.sufix}.parm7')
    nframes = traj.n_frames
    del traj

    # Get frame ID of the pseudo TS
    print('\n*** Finding frame ID of the Transition State structures ***')
    ts_indices = np.zeros(len(jobs), dtype=np.int16)
    for job_idx, job in enumerate(jobs):
        ts_index, work_lines = rf.identify_smd_TS(f'{job}/smd_{args.sufix}.txt')
        factor = int(work_lines/nframes)
        if ts_index % factor != 0:
            ts_index = int(ts_index//factor + 1)
        else:
            ts_index = int(ts_index//factor)
        ts_indices[job_idx] = ts_index

    # Write frame ID of the pseudo TS
    print('\nWritting frames ID')
    with open(f'{args.output}/ts_indexes.dat', 'w') as f:
        for job, ts in zip(jobs, ts_indices):
            f.write(f'{job},{ts}\n')
    print('\nDone')
# ==============================================================================

# Calculate matrices ===========================================================
    # Non-parallelized calculation
    if args.ncpus == 1:
        print(f'\n*** Calculating interactions for Enzyme-Substrate complex ***')
        for interaction in interactions:
            print(f'    - calculating {interaction}')
            for jobid, job in enumerate(jobs):
                rf.calculate_matrix(interaction,
                                    f'{job}/traj_{args.sufix}.nc',
                                    f'{job}/top_{args.sufix}.parm7',
                                    jobid, args.output)
        print('\nDone')
        print(f'\n*** Calculating interactions for (pseudo) TS complex ***')
        for interaction in interactions:
            print(f'    - calculating {interaction}')
            for jobid, job in enumerate(jobs):
                rf.calculate_matrix(interaction,
                                    f'{job}/traj_{args.sufix}.nc',
                                    f'{job}/top_{args.sufix}.parm7',
                                    jobid, args.output, ts_indices[jobid] - 1)
        print('\nDone')
    # Parallelized calculation
    elif args.ncpus != 1:
        ncpus = args.ncpus
        if args.ncpus > mp.cpu_count():
            print(f'\nWarning!!! Not enough CPUs')
            print(f'\nUsing the CPUs availables: {mp.cpu_count()}')
            ncpus = mp.cpu_count()
        print(f'\n*** Calculating interactions for Enzyme-Substrate complex ***')
        for interaction in interactions:
            print(f'    - calculating {interaction}')
            arg1 = [interaction for _ in range(len(jobs))]
            arg2 = [f'{job}/traj_{args.sufix}.nc' for job in jobs]
            arg3 = [f'{job}/top_{args.sufix}.parm7' for job in jobs]
            arg5 = [args.output for _ in jobs]

            with mp.Pool(processes=ncpus) as pool:
                _ = pool.starmap_async(rf.calculate_matrix,
                                    zip(arg1, arg2, arg3, range(jobs),
                                        arg5))
        print('\nDone')
        print(f'\n*** Calculating interactions for (pseudo) TS complex ***')
        for interaction in interactions:
            print(f'    - calculating {interaction}')
            arg1 = [interaction for _ in range(len(jobs))]
            arg2 = [f'{job}/traj_{args.sufix}.nc' for job in jobs]
            arg3 = [f'{job}/top_{args.sufix}.parm7' for job in jobs]
            arg5 = [args.output for _ in jobs]

            with mp.Pool(processes=ncpus) as pool:
                _ = pool.starmap_async(rf.calculate_matrix,
                                    zip(arg1, arg2, arg3, range(jobs),
                                        arg5, ts_indices - 1))
        print('\nDone')
    with open(f'{args.output}/matrices.dat', 'w') as f:
        for interaction in interactions:
            for jobid in range(job):
                path = f'{args.output}/matrices'
                f.write(f'{path}/es_{interaction}.{jobid}.pickle')
                f.write(f'{path}/pts_{interaction}.{jobid}.pickle')
# ==============================================================================

# Read interaction matrices ====================================================

print('\n*** Loading interaction matrices ***')
if 

matrix_es = np.zeros((args.nresidues, args.nresidues))
for interaction in interactions:
    matrix_es += rf.load_pickle(
                            f'{args.output}/matrices/{interaction}_es.pickle')

matrix_pts = np.zeros((args.nresidues, args.nresidues))
for interaction in interactions:
    matrix_pts += rf.load_pickle(
                            f'{args.output}/matrices/{interaction}_pts.pickle')
# ==============================================================================

# Separate contributions =======================================================
matrix_diff = matrix_pts - matrix_es
del matrix_pts
del matrix_es

# Calculate edge to edge propensity ============================================
print(f'*** Calculating edge to edge ***')

edges_catal = np.asarray([])
for catalytic_residue in catalytic_residues:
    edges_catal = np.concatenate((edges_catal,
                                  np.where(matrix_diff.nonzero()[0] ==
                                          catalytic_residue)[0]))
edge_to_edge = rf.edge_transfer_matrix(matrix_diff, args.nresidues)
for diag in range(edge_to_edge.shape[0]):
    edge_to_edge[diag][diag] = 0

flux = np.zeros(args.nresidues)
for residue in range(args.nresidues):
    if residue in catalytic_residues:
        continue
    edges_residue = np.where(matrix_diff.nonzero()[0] == residue)[0]
    edges_catal_tile = np.tile(edges_catal, len(edges_residue))
    edges_catal_tile = edges_catal_tile.astype(np.int16)
    edges_residue_repeat = np.repeat(edges_residue, len(edges_catal))
    edges_residue_repeat = edges_residue_repeat.astype(np.int16)
    flux[residue] = np.sum(abs(edge_to_edge[edges_residue_repeat,
                                                          edges_catal_tile]))
del edge_to_edge

# separate stab from destab ====================================================
total_contribution = matrix_diff.sum(axis=1)
stab_index = np.where(total_contribution < 0)[0]
destab_index = np.where(total_contribution > 0)[0]

flux_stab = flux[stab_index]
flux_destab = flux[destab_index]

print(f'*** Writing output ***') 
with open(f'{args.output}/results.csv', 'w') as f:
    f.write('residue,flux\n')
    for res in range(args.nresidues):
        if res in stab_index:
            flux_res = -flux_stab[np.where(stab_index == res)[0]]
        elif res in destab_index:
            flux_res = flux_destab[np.where(destab_index == res)[0]]
        f.write(f'{res+1},{flux_res[0]}\n')

with open(f'{args.output}/stab.dat', 'w') as f:
    for res in stab_index:
        f.write(f'{res+1},{total_contribution[res]}\n')

with open(f'{args.output}/destab.dat', 'w') as f:
    for res in destab_index:
        f.write(f'{res+1},{total_contribution[res]}\n')
print('*** Normal termination ***')

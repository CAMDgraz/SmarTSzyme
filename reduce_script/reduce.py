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

print("""
********************************************************************************
*                              SmarTSzyme-reduce                               *       
*               Selection of key residues for enzyme engineering               *
********************************************************************************
""")

# Parse arguments from the CLI
from cli import parse_reduce
import sys

if len(sys.argv) == 1:
    parse_reduce().print_help()
    sys.exit(1)
args = parse_reduce().parse_args()

# Remaining imports
import pandas as pd
import numpy as np
import multiprocessing as mp
import mdtraj as md
import os
import glob
import shutil
import utils as ut

# Check if qmmm_list or interactions_list provided
if args.qmmm_list and args.int_matrices:
    print('\nError!!! Provide either qmmm_list or interaction_matrices.')
    sys.exit()
elif not args.qmmm_list and not args.int_matrices:
    print('\nError!!! Neither qmmm_list or interaction_matrices provided.')
    sys.exit()

# Check the output argument and create the file tree
if args.output == './':
    print('\nError!!! Output path can not be ./')
    sys.exit()
try:
    os.mkdir(args.output)
    if args.qmmm_list:
        os.mkdir(f'{args.output}/matrices/')
    os.mkdir(f'{args.output}/coupling/')
except FileExistsError:
    if args.force: # Delete output folder if -f
        shutil.rmtree(args.output, ignore_errors=True)
        os.mkdir(args.output)
        if args.qmmm_list:
            os.mkdir(f'{args.output}/matrices/')
        os.mkdir(f'{args.output}/coupling/')
    else:
        print('\nError!!! Output file exist. Use -f True to overwrite.')
        raise FileExistsError

catalytic_residues = np.asarray(args.catalytic_residues, dtype=int) - 1
interactions = ['vdw', 'hbonds', 'coulomb']

# ==============================================================================

# Input options ================================================================
print(f'*** User inputs ***')
if args.qmmm_list:
    only_flux = False
    jobs_df = pd.read_csv(args.qmmm_list, header=None)
    jobs = np.asarray(jobs_df.iloc[:, 0])
    print(f'Reading QMMM jobs from:          {args.qmmm_list}')
    print(f'Number of trajectories:          {len(jobs)}')
elif args.int_matrices:
    only_flux = True
    print(f'Reading interactions from:       {args.int_matrices}')
print(f'Number of residues:              {args.nresidues}')
print(f'Catalytic residues:              {args.catalytic_residues}')
print(f'Number of CPUs:                  {args.ncpus}')
# ==============================================================================

# Get TS indices ===============================================================
if not only_flux:
    # Get number of frames
    traj = md.load(f'{jobs[0]}/traj_{args.suffix}.nc',
                top=f'{jobs[0]}/top_{args.suffix}.parm7')
    nframes = traj.n_frames
    del traj

    # Get frame ID of the pseudo TS
    print('\n*** Finding frame ID of the Transition State structures ***')
    ts_indices = np.zeros(len(jobs), dtype=np.int16)
    for job_idx, job in enumerate(jobs):
        ts_index, work_lines = ut.identify_smd_TS(f'{job}/smd_{args.suffix}.txt')
        factor = int(work_lines/nframes)
        if ts_index % factor != 0:
            ts_index = int(ts_index//factor + 1)
        else:
            ts_index = int(ts_index//factor)
        ts_indices[job_idx] = ts_index

    # Write frame ID of the pseudo TS
    print('Writting frames ID')
    with open(f'{args.output}/ts_indexes.dat', 'w') as f:
        for job, ts in zip(jobs, ts_indices):
            f.write(f'{job},{ts}\n')
    print('Done')
# ==============================================================================

# Calculate matrices ===========================================================
    # Non-parallelized calculation
    if args.ncpus == 1:
        print(f'\n*** Calculating Enzyme-Substrate complex ***')
        for interaction in interactions:
            print(f'- calculating {interaction}')
            for jobid, job in enumerate(jobs):
                ut.calculate_matrix(interaction,
                                    job, args.suffix,
                                    jobid, args.output)
        print('Done')
        print(f'\n*** Calculating (pseudo) TS complex ***')
        for interaction in interactions:
            print(f'- calculating {interaction}')
            for jobid, job in enumerate(jobs):
                ut.calculate_matrix(interaction,
                                    job, args.suffix,
                                    jobid, args.output, ts_indices[jobid] - 1)
        print('Done')
    # Parallelized calculation
    elif args.ncpus != 1:
        ncpus = args.ncpus
        if args.ncpus > mp.cpu_count():
            print(f'\nWarning!!! Not enough CPUs')
            print(f'Using the CPUs availables: {mp.cpu_count()}')
            ncpus = mp.cpu_count()
        print(f'\n*** Calculating Enzyme-Substrate complex ***')
        for interaction in interactions:
            print(f'- calculating {interaction}')
            arg1 = [interaction for _ in range(len(jobs))]
            arg2 = [f'{job}' for job in jobs]
            arg3 = [f'{args.suffix}' for _ in jobs]
            arg4 = np.arange(0, len(jobs))
            arg5 = [args.output for _ in jobs]

            with mp.Pool(processes=ncpus) as pool:
                _ = pool.starmap(ut.calculate_matrix,
                                 zip(arg1, arg2, arg3, arg4, arg5))
        print('Done')
        print(f'\n*** Calculating (pseudo) TS complex ***')
        for interaction in interactions:
            print(f'- calculating {interaction}')
            arg1 = [interaction for _ in range(len(jobs))]
            arg2 = [f'{job}' for job in jobs]
            arg3 = [f'{args.suffix}' for _ in jobs]
            arg4 = np.arange(0, len(jobs))
            arg5 = [args.output for _ in jobs]

            with mp.Pool(processes=ncpus) as pool:
                _ = pool.starmap(ut.calculate_matrix,
                                 zip(arg1, arg2, arg3, arg4, arg5,
                                     ts_indices - 1))
        print('Done')

# ==============================================================================

# Calculate coupling ===========================================================

print('\n*** Loading interaction matrices and calculating coupling ***')
if only_flux:
    path_matrices = f'{args.int_matrices}'
elif not only_flux:
    path_matrices = f'{args.output}/matrices/'

matrices_id = ut.find_matrices(path_matrices)

# Calculate coupling for each reaction mechanism
for job in matrices_id:
    # es matrix
    matrix_es = np.zeros((args.nresidues, args.nresidues))
    int_es = glob.glob(f'{path_matrices}/es_*.{job}.pickle')
    for matrix in int_es:
        matrix_es += ut.load_pickle(matrix)

    # pts matrix
    matrix_pts = np.zeros((args.nresidues, args.nresidues))
    int_pts = glob.glob(f'{path_matrices}/pts_*.{job}.pickle')
    for matrix in int_pts:
        matrix_pts += ut.load_pickle(matrix)
    
    matrix_diff = matrix_pts - matrix_es

    del matrix_pts
    del matrix_es

    # calculate edge-to-edge
    edge_to_edge, edges = ut.edge_transfer_matrix(matrix_diff)
    dweights = ut.get_weights(matrix_diff, edges)

    flux_stab = np.zeros(args.nresidues)
    flux_dest = np.zeros(args.nresidues)

    # get edges involved in the catalytic residues
    edges_catal = []
    for cr in catalytic_residues:
        for edge_id, vertices in enumerate(edges):
            if cr == vertices[0] or cr == vertices[1]:
                edges_catal.append(edge_id)

    # do not consider the diagonal elements (or embedness)
    for diag in range(edge_to_edge.shape[0]):
        edge_to_edge[diag][diag] = 0
    
    for residue in range(args.nresidues):
        if residue in catalytic_residues:
            continue
        edges_residue_stab = []
        edges_residue_dest = []
        for edge_id, vertices in enumerate(edges):
            if residue == vertices[0] or residue == vertices[1]:
                if dweights[edge_id][edge_id] < 0:
                    edges_residue_stab.append(edge_id)
                elif dweights[edge_id][edge_id] > 0:
                    edges_residue_dest.append(edge_id)
    
        # stab fluxes
        edges_catal_tile = np.tile(edges_catal, len(edges_residue_stab))
        edges_catal_tile = edges_catal_tile.astype(np.int16)
        edges_residue_repeat = np.repeat(edges_residue_stab, len(edges_catal))
        edges_residue_repeat = edges_residue_repeat.astype(np.int16)
        flux_stab[residue] = np.sum(abs(edge_to_edge[edges_residue_repeat,
                                                    edges_catal_tile]))

        # dest fluxes
        edges_catal_tile = np.tile(edges_catal, len(edges_residue_dest))
        edges_catal_tile = edges_catal_tile.astype(np.int16)
        edges_residue_repeat = np.repeat(edges_residue_dest, len(edges_catal))
        edges_residue_repeat = edges_residue_repeat.astype(np.int16)
        flux_dest[residue] = np.sum(abs(edge_to_edge[edges_residue_repeat,
                                                    edges_catal_tile]))
    # normalize fluxes
    total_flux = abs(edge_to_edge[:, edges_catal]).sum()

    flux_stab /= total_flux
    flux_dest /= total_flux

    # Print stab and destab ====================================================
    with open(f'{args.output}/coupling/stab_flux.{job}.csv', 'w') as f:
        f.write('residue,flux\n')
        for res, flux in enumerate(flux_stab):
            f.write(f'{res},{flux}\n')

    with open(f'{args.output}/coupling/dest_flux.{job}.csv', 'w') as f:
        f.write('residue,flux\n')
        for res, flux in enumerate(flux_dest):
            f.write(f'{res},{flux}\n')
print('Done')
# ==============================================================================
print('\n*** Normal termination of SmarTSzyme ***')

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

# General imports
import numpy as np
import pandas as pd
import multiprocessing as mp
import argparse

# SmarTSyzme imports
import topology_loaders
import interactions
import extra

# Functions ====================================================================
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
    inputs.add_argument('-qmmm_list', dest = 'qmmm_list', action = 'store', 
                      help = 'List of QMMM jobs to analyze.',
                      type=str, required = True)
    inputs.add_argument('-sufix', dest='sufix', action='store', type=str,
                        required=True,
                        help='Sufix for the top_, traj_ and smd_ files')
    inputs.add_argument('-nres', dest='nresidues', action='store', type=int,
                      required=True, help='Number of residues')
    inputs.add_argument('-cr', dest='catalytic_residues', nargs='+',
                        action='store', required=True,
                        help='Catalytic residues')
    inputs.add_argument('-cutoff', dest='cutoff', type=float, required=True,
                        help='''Maximum distance between residues to be consider
                          in the pairwise interactions (in A)''')
    inputs.add_argument('-ncpus', dest='ncpus', required=False, default=1,
                        type=int, help='Number of CPUs to use [default: 1]')
    outputs = parser.add_argument_group(title='Output options')
    outputs.add_argument('-out', dest='output', action='store', type=str,
                      required=True, help='prefix for the outputs')
    user_inputs = parser.parse_args()
    return user_inputs


def get_matrix(trajectory, topology, ts_index, cutoff):
    """
    Calculate pairwise interactions used to determine the flux.
    
    Parameters
    ----------
    trajectory: mdtraj.trajectory
        Steered MD trajectory
    topology: mdtraj.topology
        System topology
    ts_index: int
        Index of the (pseudo) Transition State structure
    cutoff: float
        Cutoff (A) for the non-bonded interactions.

    Return
    ------
        matrix_es: numpy.array
            2D array with the interactions in the Enzyme-Substrate
            complex.
        matrix_ts: numpy.array
            2D array with the interactions in the Transition State complex.
    """
    traj = extra.load_traj(trajectory, topology=topology)
    top_info = topology_loaders.load_top(topology)
    matrix_es = np.zeros((traj.n_residues, traj.n_residues))
    matrix_es += interactions.compute_vdw(traj[0], top_info, cutoff)
    matrix_es += interactions.compute_hbonds(traj[0], top_info, cutoff)
    matrix_es += interactions.compute_coulomb(traj[0], top_info, cutoff)

    matrix_ts = np.zeros((traj.n_residues, traj.n_residues))
    matrix_ts += interactions.compute_vdw(traj[ts_index], top_info, cutoff)
    matrix_ts += interactions.compute_hbonds(traj[ts_index], top_info, cutoff)
    matrix_ts += interactions.compute_coulomb(traj[ts_index], top_info, cutoff)
    return matrix_es, matrix_ts


def find_TS(work_file):
    """
    Finds maximum in the pulling work of a smd job.

    Parameters
    ----------
        work_file: str
            Path to the steered MD output file.
    
    Return
    ------
        maximum: float
            Maximum found for the reaction.
    """
    works = pd.read_csv(work_file, skiprows=3, skipfooter=3, engine='python',
                        header=None, sep=r'\s+')
    works = np.asarray(works.iloc[:, -1])
    works_diff = np.diff(works)
    maximum = 0
    maximum_value = 0
    for frame, diff in enumerate(works_diff):
        if diff < 0 and maximum_value < works[frame]:
            maximum = frame
            maximum_value = works[frame]
    return maximum

def edge_transfer_matrix(matrix, N):
    """
    Calculate the edge to edge transfer matrix used to select the most important
    residues based of the flux.

    Parameters
    ----------
        matrix: numpy.array
            2D array with the interactions of the system.
        N: int
            Numeber of residues of the system.
    
    Return
    ------
        edge_to_edge: numpy.array
            2D array with the edge_to_edge propensity values.
    """
    edges = matrix.nonzero()
    incidence = np.zeros((N, len(edges[0])))
    for edge_id in range(incidence.shape[1]):
        incidence[edges[0][edge_id]][edge_id] = 1
        incidence[edges[1][edge_id]][edge_id] = -1
    dweights = np.diag(matrix[edges[0], edges[1]])

    # Calculate the laplacian matrix
    laplacian_matrix = np.matmul(incidence, dweights)
    laplacian_matrix = np.matmul(laplacian_matrix, incidence.T)

    # Monroe-Penrose pseudoinverse
    pseudoinverse = np.linalg.pinv(laplacian_matrix)

    # Edge-To-Edge
    edge_to_edge = np.matmul(dweights, incidence.T)
    edge_to_edge = np.matmul(edge_to_edge, pseudoinverse)
    edge_to_edge = np.matmul(edge_to_edge, incidence)

    return edge_to_edge

def get_difference_matrix(job_path):
    """
    Return the difference (matrix) between the (pseudo) Transition State  and
    the Enzyme-Substrate complexes.

    Parameters
    ----------
        job_path: str
            Path to the folder containing the Steered MD job.
        
    Return
    ------
        matrix_dif: numpy.array
            2D array with the difference between the (pseudo) Transition State
            and the Enzyme-Substrate complexes.
    """
    ts_index = find_TS(f'{job_path}/smd_{args.sufix}.txt')
    if ts_index % 2 == 1:
        ts_index = int(ts_index/2 + 0.5)
    else:
        ts_index = int(ts_index/2)
    matrix_es, matrix_ts = get_matrix(f'{job_path}/traj_{args.sufix}.nc',
                                      f'{job_path}/top_{args.sufix}.parm7',
                                        ts_index, cutoff)
    matrix_diff = matrix_ts - matrix_es
    return matrix_diff


# Main script ==================================================================
# Get options from the CLI
args = parse_arguments()
cutoff = args.cutoff / 10
catalytic_residues = args.catalytic_residues
catalytic_residues = np.asarray(catalytic_residues, dtype=int)
ncpus = args.ncpus
sufix = args.sufix

# Read list of qmmm jobs
jobs_df = pd.read_csv(args.qmmm_list, header=None)
jobs = np.asarray(jobs_df.iloc[:, 0])

# Split jobs to calculate in parallel
job_iterable = [(jobs[chunk*ncpus:(chunk+1)*ncpus]) for
                chunk in range(len(jobs)//ncpus + 1) if len(jobs) > chunk*ncpus]

# TODO This is parallelization part which has to be improved
matrix_diff_total = np.zeros((args.nresidues, args.nresidues))

# Iterate over the job_iterable
for batch, job_list in enumerate(job_iterable):
    print(f'Calculating interaction matrix for batch {batch + 1} of' +
          f' {len(job_iterable)} with {len(job_list)} reaction(s) ...')
    chunksize = int(1 if len(job_list)//ncpus == 0 else len(job_list)/ncpus)
    with mp.Pool(processes=10) as pool:
       results = pool.map(get_difference_matrix, job_list, chunksize=1)
    for result in results:
        matrix_diff_total += result

# Separation of the difference in stabilizing and destabilizing
matrix_diff_total /= len(jobs)

stabilizing = np.copy(matrix_diff_total)
stabilizing[np.where(stabilizing > 0)] = 0
stabilizing[catalytic_residues] = matrix_diff_total[catalytic_residues]
stabilizing.T[catalytic_residues] = matrix_diff_total.T[catalytic_residues]

destabilizing = np.copy(matrix_diff_total)
destabilizing[np.where(destabilizing < 0)] = 0
destabilizing[catalytic_residues] = matrix_diff_total[catalytic_residues]
destabilizing.T[catalytic_residues] = matrix_diff_total.T[catalytic_residues]

# Calculate stabilizing part
print('Calculating edge_to_edge for stabilizing residues')
edges_catal = np.asarray([])
for catalytic_residue in catalytic_residues:
    edges_catal = np.concatenate((edges_catal,
                                  np.where(stabilizing.nonzero()[0] ==
                                          catalytic_residue)[0]))
edge_to_edge = edge_transfer_matrix(stabilizing, args.nresidues)
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

# Calculate destabilizing part
print('Calculating edge_to_edge for destabilizing residues')
edges_catal = np.asarray([])
for catalytic_residue in catalytic_residues:
    edges_catal = np.concatenate((edges_catal,
                                  np.where(destabilizing.nonzero()[0] ==
                                          catalytic_residue)[0]))
edge_to_edge = edge_transfer_matrix(destabilizing, args.nresidues)
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

# Output results ===============================================================
flux_diff = flux_destabilizing - flux_stabilizing
extra.write_results(flux_diff, args.output)
extra.plot_results(flux_stabilizing, flux_destabilizing, args.output,
                   args.nresidues, nres_out=10)

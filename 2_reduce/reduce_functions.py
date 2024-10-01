#!/bin/python

"""
@autors: Daniel Platero Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sánchez-Murcia [pedro.murcia@medunigraz.at]
"""

"""
Functions used in reduce.py
"""

# General imports
import numpy as np
import pandas as pd
import multiprocessing as mp
import argparse
import pickle

# SmarTSyzme imports
import topology_loaders
import interactions
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


def indentify_smd_TS(work_file):
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

def calculate_matrix(interaction, trajectory, topology, cutoff, ts_index=None):
    """
    Calculate the Van der Waals interactions for the ES and TS.
    """

    traj = extra.load_traj(trajectory, topology=topology)
    top_info = topology_loaders.load_top(topology)

    if not ts_index:
        if interaction == 'vdw':
            matrix = interactions.compute_vdw(traj[0], top_info, cutoff)
        elif interaction == 'coulomb':
            matrix = interactions.compute_coulomb(traj[0], top_info, cutoff)
        elif interaction == 'hbonds':
            matrix = interactions.compute_hbonds(traj[0], top_info, cutoff)
    else:
        if interaction == 'vdw':
            matrix = interactions.compute_vdw(traj[ts_index], top_info, cutoff)
        elif interaction == 'coulomb':
            matrix = interactions.compute_coulomb(traj[ts_index], top_info, cutoff)
        elif interaction == 'hbonds':
            matrix = interactions.compute_hbonds(traj[ts_index], top_info, cutoff)
    return matrix

def write_pickle(matrix, file):
    with open(file, 'wb') as f:
        pickle.dump(matrix, f)
    return None

def load_pickle(file):
    with open(file, 'rb') as f:
        data = pickle.load(f)
    return data
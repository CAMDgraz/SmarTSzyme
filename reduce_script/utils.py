#!/bin/python

"""
@autors: Daniel Platero Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sánchez-Murcia [pedro.murcia@medunigraz.at]
"""

"""
Functions used along the reduce.py script
"""

# General imports
import numpy as np
import pandas as pd
import multiprocessing as mp
import argparse
import pickle
import os
import mdtraj as md
import matplotlib as mpl
import glob
import sys

# SmarTSyzme imports
import topology_loaders
import interactions

# Arguments parsing
def parse_arguments_reduce():
    """
    Parse arguments of the cli
    """
    desc = '''\nReduce: Reduction of the mutational landscape'''
    parser = argparse.ArgumentParser(prog='Reduce',
                                     usage='reduce.py [OPTIONS]',
                                     description = desc,
                                     add_help=True,
                                     allow_abbrev = False)
    inputs = parser.add_argument_group(title='Input options')
    inputs.add_argument('-qmmm_list', dest='qmmm_list', action='store', 
                        help='File with the paths to the QMMM results.', 
                        type=str, required = False)
    
    inputs.add_argument('-int_matrices', dest='int_matrices', 
                        help='Path to the interaction matrices',
                        action='store', required=False, type=str)
    
    inputs.add_argument('-suffix', dest='suffix', action='store', type=str,
                        help='Suffix for the top_, traj_ and smd_ files',
                        required=True)
    
    inputs.add_argument('-nres', dest='nresidues', action='store', type=int,
                      required=True, help='Number of residues')
    
    inputs.add_argument('-catres', dest='catalytic_residues', nargs='+',
                        action='store', required=True,
                        help='ResID of the catalytic residues')
    
    inputs.add_argument('-ncpus', dest='ncpus', required=False, default=1,
                        type=int, help='Number of CPUs to use [default: 1]')
    
    inputs.add_argument('-f', dest='force', action='store',
                        help='Delete output folder if it exist',
                        required=False, default=False, type=bool)
    
    outputs = parser.add_argument_group(title='Output options')
    outputs.add_argument('-out', dest='output', action='store', type=str,
                      required=False, help='Output path (different from ./)',
                      default='./out')
    user_inputs = parser.parse_args()
    return user_inputs

def parse_arguments_analysis():
    """
    Parse arguments of the cli
    """
    parser = argparse.ArgumentParser(prog='Analysis',
                                     usage='analysis.py [OPTIONS]',
                                     add_help=True,
                                     allow_abbrev = False)
    inputs = parser.add_argument_group(title='Input')
    inputs.add_argument('-coupling_path', dest='coupling_path', action='store',
                        required=True, type=str,
                        help='Path to the coupling folder of SmarTSzyme')
    
    inputs.add_argument('-int_path', dest='int_path', action='store',
                        required=True, type=str,
                        help='Path to the matrices folder of SmarTSzyme')
    inputs.add_argument('-nres', dest='nres', type=int, action='store',
                        required=True, help='Number of residues')
    
    inputs.add_argument('-batch', dest='batch', required=False, type=int,
                        default=None, action='store',
                        help='Increments in the number of trajectories to analyze')
    
    inputs.add_argument('-topn', dest='topn', action='store', type=int, 
                        default=10, required=False, 
                        help='Top n residues to output')
    
    inputs.add_argument('-exp_res', dest='exp_res', action='store', nargs='+',
                        required=False, default=None,
                        help='ID of experimental residues to check (1-based)')
    
    output = parser.add_argument_group(title='Output')
    output.add_argument('-out', dest='output', action='store', required=False,
                        default='./', type=str,
                        help='Output directory')
    user_inputs = parser.parse_args()
    return user_inputs

# Topology and trajectory
def load_traj(trajectory, topology):
    """
    Load trajectory using mdtraj

    Parameters
    ----------
    trajectory: str
        Path to the trajectory file
    topology: str
        Path to the topology file
    
    Return
    ------
    traj: md.Trajectory
        Trajectory loaded with mdtraj
    """
    name, extension = os.path.splitext(trajectory)
    # Load trajectory
    traj = md.load(trajectory, top=topology)
    return traj

# Main functionalities of reduce
def pairwise_distance(atom1: np.ndarray,
                      atom2: np.ndarray,
                      coord: np.ndarray) -> np.ndarray:
    """
    Compute the pairwise distance between two arrays of atoms.

    Parameters
    ----------
    atom1 : np.ndarray
        Numpy array with the atom indices
    atom2 : np.ndarray
        Numpy array with the second atom indices
    
    Return
    ------
    distance : np.ndarray
        Numpy array with the pairwise distance
    """
    points_pair = np.vstack((atom1, atom2))
    distance_ = np.diff(coord[points_pair.T], axis=1)**2
    distance = (distance_.sum(axis=2)**(1/2)).flatten()
    return distance

def identify_smd_TS(work_file: str) -> tuple[float, int]:
    """
    Finds maximum in the pulling work of a smd job.

    Parameters
    ----------
        work_file: str
            Path to the steered MD output file.
    
    Return
    ------
        maximum: float
            Maximum found in the sMD job.
        num_works: int
            Number of pulling works.
    """
    works = pd.read_csv(work_file, skiprows=3, skipfooter=3, engine='python',
                        header=None, sep=r'\s+')
    works = np.asarray(works.iloc[:, -1])
    works_diff = np.diff(works)
    maximum = 0
    maximum_value = 0
    for frame, diff in enumerate(works_diff):
        if diff < 0 and maximum_value < works[frame]:
            maximum = frame + 1
            maximum_value = works[frame]
    return maximum, len(works)

def get_incidence_matrix(matrix: np.ndarray) -> tuple[list, np.ndarray]:
    """
    Construct the incidence matrix of the interaction graph

    Parameters
    ----------
        matrix : numpy.array
            2D numpy array matrix of the graph.
    
    Return
    ------
        edges : list
            List of tuples (node1,node2) containing the nodes involved in the
            each edge.
        incidence_matrix : np.ndarray
            2D Numpy array with the incidence matrix. The direction is assign 
            from the node with the lower index to the highest.
    """
    # Get edges from the matrix
    edges = []
    connections = matrix.nonzero()
    for node1, node2 in zip(connections[0], connections[1]):
        if (node1,node2) not in edges and (node2,node1) not in edges:
            edges.append((node1,node2))
    
    incidence_matrix = np.zeros((len(matrix), len(edges)))
    for edge, vertices in enumerate(edges):
        incidence_matrix[vertices[0]][edge] = 1
        incidence_matrix[vertices[1]][edge] = -1
    return edges, incidence_matrix

def get_weights(matrix: np.ndarray, edges: list) -> np.ndarray:
    """
    Get diagonal matrix of weights

    Parameters
    ----------
        matrix : numpy.array
            2D interaction matrix
        edges : list
            List with the tuples containing the nodes in each edge.

    Return
    ------
        dweights : numpy.array
            Diagonal matrix with the weights of each edge.
    """
    weights = np.zeros(len(edges))
    for edge, vertices in enumerate(edges):
        weights[edge] = matrix[vertices[0]][vertices[1]]
    dweights = np.diag(weights)
    return dweights

def edge_transfer_matrix(matrix: np.ndarray) -> tuple[np.ndarray, list]:
    """
    Calculate the edge to edge transfer matrix used to select the most important
    residues based of the flux.

    Parameters
    ----------
        matrix: numpy.array
            2D array with the interactions of the system.
    
    Return
    ------
        edge_to_edge: numpy.array
            2D array with the edge_to_edge propensity values.
        esges : list
            List with the tuples containing the nodes in each edge.
    """
    edges, incidence = get_incidence_matrix(matrix)
    dweights = get_weights(matrix, edges)

    # Calculate the laplacian matrix
    laplacian_matrix = incidence @ dweights
    laplacian_matrix = laplacian_matrix @ incidence.T

    # Monroe-Penrose pseudoinverse
    pseudoinverse = np.linalg.pinv(laplacian_matrix)

    # Edge-To-Edge
    edge_to_edge = dweights @ incidence.T
    edge_to_edge = edge_to_edge @ pseudoinverse
    edge_to_edge = edge_to_edge @ incidence

    return edge_to_edge, edges

def calculate_matrix(interaction: str, trajectory: str, topology: str,
                     jobid: int, output: str, ts_index=None) -> None:
    """
    Calculate interactions for the ES and TS.

    Parameters
    ----------
        interaction: str
            Interaction to calculate.
        trajectory: str
            Path to the trajectory file.
        topology: str
            Path to the topology file.
        jobid: int
            Jobid value
        output: str
            Path to output
        ts_index: int
            Frame id of the pseudo TS structure
    
    Return
    ------
        matrix: np.ndarray
            Numpy 2D matrix (nres x nres) with the residues pairwise
            interaction.
    """

    traj = load_traj(trajectory, topology=topology)
    top_info = topology_loaders.load_top(topology)
    cutoff = 1

    if not ts_index:
        if interaction == 'vdw':
            matrix = interactions.compute_vdw(traj[0], top_info, cutoff)
        elif interaction == 'coulomb':
            matrix = interactions.compute_coulomb(traj[0], top_info, cutoff)
        elif interaction == 'hbonds':
            matrix = interactions.compute_hbonds(traj[0], top_info, cutoff)
        outfile = f'{output}/matrices/es_{interaction}.{jobid}.pickle'
    else:
        if interaction == 'vdw':
            matrix = interactions.compute_vdw(traj[ts_index], top_info, cutoff)
        elif interaction == 'coulomb':
            matrix = interactions.compute_coulomb(traj[ts_index], top_info,
                                                  cutoff)
        elif interaction == 'hbonds':
            matrix = interactions.compute_hbonds(traj[ts_index], top_info,
                                                 cutoff)
        outfile = f'{output}/matrices/pts_{interaction}.{jobid}.pickle'
    
    write_pickle(matrix, outfile)
    
    return None

def find_matrices(path: str) -> np.ndarray:
    """
    Return a list with the matrices id in the path

    Parameters
    ----------
        path: str
            Path to the interaction matrices
    
    Return
    ------
        matrices_id: np .ndarray
            Numpy array with the matrices id in path
    """
    es_matrices = glob.glob(f'{path}/es_*.pickle')
    pts_matrices = glob.glob(f'{path}/pts_*.pickle')

    if len(es_matrices) != len(pts_matrices):
        print('\nError!!! Number of matrices for the es different than for the pts')
        sys.exit()
    
    matrices_id = np.unique([int(matrix.split('.')[-2]) for matrix in
                             es_matrices])
    return matrices_id


# Output functions
def write_pickle(matrix: np.ndarray, file: str) -> None:
    """
    Write pickle file

    Parameters
    ----------
        matrix: np.ndarray
            Numpy matrix to save
        file: str
            Filename for storing the matrix
    """
    with open(file, 'wb') as f:
        pickle.dump(matrix, f)
    return None

def load_pickle(file: str) -> np.ndarray:
    """
    Load data from a pickle file

    Parameters
    ----------
        file: str
            Path the pickle file
    
    Return
    ------
        data: np.ndarray
            Data contained in the pickle file
    """
    with open(file, 'rb') as f:
        data = pickle.load(f)
    return data
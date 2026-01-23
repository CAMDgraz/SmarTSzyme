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
import pickle
import os
import mdtraj as md
import glob
import sys

# SmarTSyzme imports
import topology_loaders
import interactions

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

def calculate_matrix(interaction: str, job: str, suffix: str,
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

    traj_file = f'{job}/traj_{suffix}.nc'
    top_file = f'{job}/top_{suffix}.parm7'
    qmfile = f'{job}/qmmm.out'
    traj = load_traj(traj_file, topology=top_file)
    top_info = topology_loaders.load_top(top_file)
    cutoff = 1
    charges = get_charges(top_info, qmfile, ts_index)

    if not ts_index:
        if interaction == 'vdw':
            matrix = interactions.compute_vdw(traj[0], top_info, cutoff)
        elif interaction == 'coulomb':
            matrix = interactions.compute_coulomb(traj[0], top_info, cutoff,
                                                  charges)
        elif interaction == 'hbonds':
            matrix = interactions.compute_hbonds(traj[0], top_info, cutoff)
        outfile = f'{output}/matrices/es_{interaction}.{jobid}.pickle'
    else:
        if interaction == 'vdw':
            matrix = interactions.compute_vdw(traj[ts_index], top_info, cutoff)
        elif interaction == 'coulomb':
            matrix = interactions.compute_coulomb(traj[ts_index], top_info,
                                                  cutoff, charges)
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

# Find pattern in file
def find_pattern(file:str, pattern:str) -> list:
    """
    Docstring for find_pattern_line
    
    :param file: Text file to find the pattern
    :type file: str
    :param pattern: Pattern to find in <file>
    :type pattern: str
    :return: List of tuples with the line number and the line matching the 
             pattern
    :rtype: list
    """
    matches = []
    with open(file, 'r') as f:
        for lno, line in enumerate(f):
            if pattern in line:
                matches.append((lno, line.strip()))
    return matches

# Get QM point charges
def quantum_charges(file:str) -> tuple:
    """
    Extract quantum charges from smd (amber+DFTB3)
    
    :param file: Output file of the smd job
    :type file: str
    :return: Tuple with the mapping of the atom ID from the QM to the MM region
    and the point charges for each atom in each structure.
    :rtype: tuple
    """
    # Patterns in the out file to look for
    qatoms_pattern = 'nquant'
    qmapping_pattern = 'QM Region Cartesian Coordinates (*=link atom)'
    qcharges_pattern = 'Atomic Charges for Step'

    # Number of atoms in the QM region (without the link dummy atoms)
    _, qatoms_match = find_pattern(file, qatoms_pattern)[0]
    qatoms = int(qatoms_match.split()[-1])

    # Mapping of the QM atoms in the MM region
    qmapping_start, _ = find_pattern(file, qmapping_pattern)[0]
    qmapping = pd.read_csv(file, engine='python', header=None, sep=r'\s+',
                           skiprows=qmapping_start + 2, nrows=qatoms)
    qmapping = np.asarray(qmapping.iloc[:, 2])

    # Quantum point charges
    qcharges_starts = find_pattern(file, qcharges_pattern)
    qcharges = []
    for lno, _ in qcharges_starts:
        charges = pd.read_csv(file, engine='python', header=None, sep=r'\s+',
                              skiprows=lno+2, nrows=qatoms)
        qcharges.append(np.asarray(charges.iloc[:, 2]))
    
    return (qmapping, qcharges)

def get_charges(top_info, qmfile, ts_index=None):
    """
    Mix the MM charges and the QM charges
    """

    # MM charges
    charges = np.asarray(top_info['CHARGE'])
    charges /= 18.2223 # Factor for the charges in amber

    # QM charges
    qmmapping, qmcharges = quantum_charges(qmfile)

    # Replace the MM by QM in the active site
    if not ts_index:
        charges[qmmapping - 1] = qmcharges[0]
    else:
        charges[qmmapping - 1] = qmcharges[ts_index]
    
    return charges
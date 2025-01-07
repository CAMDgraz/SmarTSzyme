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

# SmarTSyzme imports
import topology_loaders
import interactions

# Arguments parsing
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
    inputs.add_argument('-f', dest='force', action='store',
                        help='Delete output folder if it exist',
                        required=False, default=False, type=bool)
    outputs = parser.add_argument_group(title='Output options')
    outputs.add_argument('-out', dest='output', action='store', type=str,
                      required=False, help='Output path (different from ./)',
                      default='./out')
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
def pairwise_distance(atom1, atom2, coord):
    """
    Compute the pairwise distance between two arrays of atoms.

    Parameters
    ----------
    atom1 : np.array
        Numpy array with the atom indices
    atom2 : np.array
        Numpy array with the second atom indices
    
    Return
    ------
    distance : np.array
        Numpy array with the pairwise distance
    """
    points_pair = np.vstack((atom1, atom2))
    distance_ = np.diff(coord[points_pair.T], axis=1)**2
    distance = (distance_.sum(axis=2)**(1/2)).flatten()
    return distance

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
            maximum = frame + 1
            maximum_value = works[frame]
    return maximum, len(works)

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

def calculate_matrix(interaction, trajectory, topology, cutoff, ts_index=None):
    """
    Calculate the Van der Waals interactions for the ES and TS.
    """

    traj = load_traj(trajectory, topology=topology)
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

# Output
def write_pickle(matrix, file):
    with open(file, 'wb') as f:
        pickle.dump(matrix, f)
    return None

def load_pickle(file):
    with open(file, 'rb') as f:
        data = pickle.load(f)
    return data

def mplstyle():
    mpl.rcParams['axes.titlesize'] = 30
    mpl.rcParams['axes.labelsize'] = 30
    mpl.rcParams['axes.spines.top'] = True
    mpl.rcParams['axes.spines.bottom'] = True
    mpl.rcParams['axes.spines.left'] = True
    mpl.rcParams['axes.spines.right'] = True
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['axes.titlepad'] = 10

    # Latex configuration
    mpl.rcParams['text.usetex'] = False

    # Legend configuration
    mpl.rcParams['legend.fancybox'] = True
    mpl.rcParams['legend.loc'] = 'lower right'
    mpl.rcParams['legend.fontsize'] = 20
    mpl.rcParams['legend.handletextpad'] = 0.1

    # Ticks configuration
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.color'] = (0.2, 0.2, 0.2)
    mpl.rcParams['ytick.color'] = (0.2, 0.2, 0.2)

    # Figure configuration
    mpl.rcParams['figure.figsize'] = 10, 10
    mpl.rcParams['figure.dpi'] = 300

    # Layout
    mpl.rcParams['figure.constrained_layout.use'] = True

    # Saving
    mpl.rcParams['savefig.dpi'] = 300
    mpl.rcParams['savefig.transparent'] = False
    return None

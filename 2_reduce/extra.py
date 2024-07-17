#!/bin/python
# -*- coding: utf-8 -*-

"""
@authors: Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sánchez-Murcia [pedro.murcia@medunigraz.at]
"""

"""
Functions used along the reduce script 
"""

# Imports
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import heapq
import mdtraj as md

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

def write_results(flux, filename):
    """
    Write resulting flux to a .dat file

    Parameters
    ----------
    flux: np.array
        Numpy array containing the flux for each residue
    filename: str
        Path to .dat file
    
    Return
    ------
    None
    """
    with open(filename + 'flux.dat', 'w') as outfile:
        outfile.write(f'residue,flux\n')
        for residue, flux in enumerate(flux):
            outfile.write(f'{residue + 1},{flux}\n')
    return None


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

def plot_results(flux_stabilizing, flux_destabilizing, filename_prefix, nres,
                 nres_out):
    """
    Plot the results from SmarTSzyme

    Parameters
    ----------
    flux_stabilizing: np.array
        Numpy array with the fluxes of the stabilizing interactions
    flux_destabilizing: np.array
        Numpy array with the fluxes of the destabilizing interactions
    filename_prefix: str
        Prefix of the plots
    nres: int
        Numeber of residues

    Return
    ------
        None
    """

    mplstyle()

    # Bar plot
    nresidues = np.arange(1, nres + 1, 1)
    per_res_contribution = flux_destabilizing - flux_stabilizing
    stab_contribution = np.where(per_res_contribution < 0)[0]
    destab_contribution = np.where(per_res_contribution > 0)[0]

    fig, ax = plt.subplots()
    ax.bar(nresidues[stab_contribution], per_res_contribution[stab_contribution]
           ,color='blue', label='Stabilizing interactions')
    ax.bar(nresidues[destab_contribution],
           per_res_contribution[destab_contribution], color='red',
           label='Destabilizing interactions')
    ax.set_xticks(nresidues[::5])
    ax.set_xticklabels([f'{label}' for label in nresidues[::5]], rotation=90,
                       size=10)

    ax.set_xlabel('Residues')
    ax.set_ylabel('Flux')
    ax.legend(loc='upper right')
    fig.savefig(filename_prefix  + '_fluxes.png')
    plt.close(fig)

    # Bar plot higher fluxes stabilizing and destabilizing
    ordered_flux = []
    heapq.heapify(ordered_flux)
    for resid, flux in enumerate(per_res_contribution):
        heapq.heappush(ordered_flux, (flux, resid + 1))

    nstabilizing = heapq.nsmallest(nres_out, ordered_flux)
    ndestabilizing = heapq.nlargest(nres_out, ordered_flux)    
    
    plot_flux_stabilizing = []
    plot_resid_stabilizing = []
    plot_flux_destabilizing = []
    plot_resid_destabilizing = []
    for i_stab, i_destab in zip(nstabilizing, ndestabilizing):
        plot_flux_stabilizing.append(-i_stab[0])
        plot_resid_stabilizing.append(i_stab[1])
        plot_flux_destabilizing.append(-i_destab[0])
        plot_resid_destabilizing.append(i_destab[1])
    

    fig, ax = plt.subplots()
    ax.bar(plot_resid_stabilizing, plot_flux_stabilizing, color='blue',
           label='Stabilizing')
    ax.bar(plot_resid_destabilizing, plot_flux_destabilizing, color='red',
           label='Destabilizing')
    ax.set_xlabel(r'Residues')
    ax.set_ylabel(r'Flux')
    fig.savefig(filename_prefix  + '_nres.png')
    plt.close(fig)
    return

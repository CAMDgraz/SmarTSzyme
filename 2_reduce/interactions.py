#!/bin/python
# -*- coding: utf-8 -*-

"""
@authors: Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sánchez-Murcia [pedro.murcia@medunigraz.at]
"""

"""
Calculate interactions between residues 
"""
import numpy as np
import reduce_functions as rf
from sklearn.neighbors import KDTree as kdtree

def compute_vdw(traj, top_info, cutoff):
    """
    Calculate van der Waals interaction interresidue using a kdtree structure 
    for faster computation. Only half of the matrix is computed.
    """

    n_excluded = np.asarray(top_info['NUMBER_EXCLUDED_ATOMS'])
    excluded_atoms = np.asarray(top_info['EXCLUDED_ATOMS_LIST'])
    type_index = np.asarray(top_info['ATOM_TYPE_INDEX'])
    acoef = np.asarray(top_info['LENNARD_JONES_ACOEF'])
    bcoef = np.asarray(top_info['LENNARD_JONES_BCOEF'])
    nonbonded_index = np.asarray(top_info['NONBONDED_PARM_INDEX'])
    ntypes = top_info['POINTERS'][1]
    residue_pointer = np.asarray(top_info['RESIDUE_POINTER']) - 1
    names = np.asarray(top_info['ATOM_NAME'])
    coord = traj.xyz[0]
    van_der_Waals = np.zeros((traj.n_residues, traj.n_residues))

    for residue1 in range(traj.n_residues):
        atoms1 = traj.topology.select(f'resid {residue1}')
        to_delete = []
        for atom_id, atom in enumerate(atoms1):
            if names[atom] in ['C', 'N', 'O', 'CA']:
                to_delete.append(atom_id)
        atoms1 = np.delete(atoms1, to_delete)

        frame_tree = kdtree(coord)
        in_cutoff = np.concatenate(frame_tree.query_radius(coord[atoms1],
                                                           r=cutoff))
        neighbors = [np.where(residue_pointer <= atom)[0][-1]
                     for atom in in_cutoff]
        neighbors = np.unique(np.asarray(neighbors))
        for residue2 in neighbors[np.where(neighbors > residue1)[0]]:
            atoms2 = traj.topology.select(f'resid {residue2}')
            to_delete = []
            for atom_id, atom in enumerate(atoms2):
                if names[atom] in ['C', 'N', 'O', 'CA']:
                    to_delete.append(atom_id)
            atoms2 = np.delete(atoms2, to_delete)

            atoms1_repeated = np.repeat(atoms1, len(atoms2))
            atoms2_tiled = np.tile(atoms2, len(atoms1))
            vdw_atoms1 = np.zeros(len(atoms1_repeated))
            atom_id = 0
            for atom_i, atom_j in zip(atoms1_repeated, atoms2_tiled):
                r2 = np.square(coord[atom_i] - coord[atom_j])
                r6 = (r2.sum())**3 * (10**6)

                exclusion = np.asarray([])
                if atom_j > atom_i:
                    higher_atom = atom_j
                    if n_excluded[atom_i] != 1: # 1 stands for non exclusion
                        exclude_from = n_excluded[:atom_i].sum()
                        exclude_to = n_excluded[:atom_i + 1].sum()
                        exclusion = excluded_atoms[exclude_from:exclude_to]

                else: # check minimum index for excluded atoms
                    higher_atom = atom_i
                    if n_excluded[atom_j] != 1:
                        exclude_from = n_excluded[:atom_j].sum()
                        exclude_to = n_excluded[:atom_j + 1].sum()
                        exclusion = excluded_atoms[exclude_from : exclude_to]
                if higher_atom + 1 not in exclusion:
                    pair_index = (ntypes * (type_index[atom_i] - 1)
                                  + type_index[atom_j])
                    lj_index = nonbonded_index[pair_index - 1]

                    # Calculate vdW
                    lj_alpha = acoef[lj_index - 1]
                    lj_beta = bcoef[lj_index - 1]
                    vdw_ = (lj_alpha/(r6)**2) - (lj_beta/(r6))
                    vdw_atoms1[atom_id] = vdw_
                atom_id += 1

            van_der_Waals[residue1][residue2] = np.sum(vdw_atoms1)
            van_der_Waals[residue2][residue1] = np.sum(vdw_atoms1)
    return van_der_Waals

def compute_hbonds(traj, top_info, cutoff):
    """
    From MM-ISMSA
    """
    # r-values
    r_min = 1.5
    r_min_ideal = 1.8
    r_max_ideal = 2.4
    r_max = 2.7
    # alpha-values
    a_min = 100
    a_min_ideal = 130
    a_max_ideal = 165
    a_max = 180
    # beta-values
    b_min = 90
    b_min_ideal = 115
    b_max_ideal = 145
    b_max = 180

    atomic_numbers = np.asarray(top_info['ATOMIC_NUMBER'])
    residue_pointer = np.asarray(top_info['RESIDUE_POINTER']) - 1
    charges = np.asarray(top_info['CHARGE'])
    charges /= 18.2223
    charges = np.round(charges)
    names = np.asarray(top_info['ATOM_NAME'])
    # Get Hydrogen atoms and donors
    bonds_w_H = np.asarray(top_info['BONDS_INC_HYDROGEN'])
    bonds_w_H = np.delete(bonds_w_H, np.arange(2, len(bonds_w_H), 3))
    bonds_w_H = (bonds_w_H/3 + 1).astype(np.int32)
    no_donor_index = np.concatenate((np.where(atomic_numbers[bonds_w_H - 1]
                                              == 7)[0], 
                                     np.where(atomic_numbers[bonds_w_H - 1]
                                              == 8)[0]))
    donors = bonds_w_H[no_donor_index]
    hydrogens_index = no_donor_index + (-1)**(no_donor_index%2)
    hydrogens = bonds_w_H[hydrogens_index]
    # Get Acceptors and X atoms
    bonds_wo_H = np.asarray(top_info['BONDS_WITHOUT_HYDROGEN'])
    bonds_wo_H = np.delete(bonds_wo_H, np.arange(2, len(bonds_wo_H), 3))
    bonds_wo_H = (bonds_wo_H/3 + 1).astype(np.int32)
    no_acceptor_index = np.concatenate((np.where(atomic_numbers[bonds_wo_H - 1]
                                                 == 7)[0],
                                        np.where(atomic_numbers[bonds_wo_H - 1]
                                                 == 8)[0]))
    acceptors = np.concatenate((bonds_wo_H[no_acceptor_index], donors))
    x_index = no_acceptor_index + (-1)**(no_acceptor_index%2)
    x_atoms = bonds_wo_H[x_index]    
    # Get coordinates
    coord = traj.xyz[0]


    hbonds_matrix = np.zeros((traj.n_residues, traj.n_residues))
    for residue1 in range(traj.n_residues):
        atoms1 = traj.topology.select(f'resid {residue1}')
        to_delete = []
        for atom_id, atom in enumerate(atoms1):
            if names[atom] in ['C', 'N', 'O', 'CA']:
                to_delete.append(atom_id)
        atoms1 = np.delete(atoms1, to_delete)

        frame_tree = kdtree(coord)
        in_cutoff = np.concatenate(frame_tree.query_radius(coord[atoms1], r=1))
        neighbors = [np.where(residue_pointer <= atom)[0][-1]
                     for atom in in_cutoff]
        neighbors = np.unique(np.asarray(neighbors))
        # Hydrogens and donors in res1
        h1_info = np.intersect1d(atoms1 + 1, hydrogens, return_indices=True)
        hydrogens_res1 = h1_info[0].astype(np.int32)
        donors_res1 = donors[h1_info[2]].astype(np.int32)
        # Acceptors and X atoms in res1
        a1_info = np.intersect1d(atoms1 + 1, acceptors, return_indices=True)
        acceptors_res1 = a1_info[0].astype(np.int32)
        x_atoms_res1 = x_atoms[a1_info[2]].astype(np.int32) 

        for residue2 in neighbors[np.where(neighbors > residue1)[0]]:
            atoms2 = traj.topology.select(f'resid {residue2}')
            to_delete = []
            for atom_id, atom in enumerate(atoms2):
                if names[atom] in ['C', 'N', 'O', 'CA']:
                    to_delete.append(atom_id)
            atoms2 = np.delete(atoms2, to_delete)

            # Hydrogens and donors in res1
            h2_info = np.intersect1d(atoms2 + 1, hydrogens, return_indices=True)
            hydrogens_res2 = h2_info[0].astype(np.int32)
            donors_res2 = donors[h2_info[2]].astype(np.int32)
            # Acceptors and X atoms in res1
            a2_info = np.intersect1d(atoms2 + 1, acceptors, return_indices=True)
            acceptors_res2 = a2_info[0].astype(np.int32)
            x_atoms_res2 = x_atoms[a2_info[2]].astype(np.int32)
            # Join all type of atoms for faster calculation
            hydrogens_res12 = np.concatenate((hydrogens_res1, hydrogens_res2))
            donors_res12 = np.concatenate((donors_res1, donors_res2))
            acceptors_res12 = np.concatenate((acceptors_res1, acceptors_res2))
            x_atoms_res12 = np.concatenate((x_atoms_res1, x_atoms_res2))
            # Repeating and tiling atoms
            hydrogens_repeated = np.repeat(hydrogens_res12,
                                           len(acceptors_res12))
            donors_repeated = np.repeat(donors_res12, len(acceptors_res12))
            acceptors_tiled = np.tile(acceptors_res12, len(hydrogens_res12))
            x_repeated = np.repeat(x_atoms_res12, len(hydrogens_res12))
            # Delete H and acceptors in the same residue
            to_delete = []
            for hyd_idx in range(len(hydrogens_repeated)):
                hyd, don = hydrogens_repeated[hyd_idx], acceptors_tiled[hyd_idx] 
                if hyd - 1 in atoms1 and don - 1 in atoms1:
                    to_delete.append(hyd_idx)
                elif hyd - 1 in atoms2 and don - 1 in atoms2:
                    to_delete.append(hyd_idx)
            hydrogens_repeated = np.delete(hydrogens_repeated, to_delete)
            donors_repeated = np.delete(donors_repeated, to_delete)
            acceptors_tiled = np.delete(acceptors_tiled, to_delete)
            x_repeated = np.delete(x_repeated, to_delete)
            # Eliminate possibility of same atom for donor and acceptor
            #to_delete = np.where((acceptors_tiled == donors_repeated)
            #                     == True)[0]
            #hydrogens_repeated = np.delete(hydrogens_repeated, to_delete)
            #donors_repeated = np.delete(donors_repeated, to_delete)
            #acceptors_tiled = np.delete(acceptors_tiled, to_delete)
            #x_repeated = np.delete(x_repeated, to_delete)
            # Calculate relevant distances 
            acceptor_hydrogen = rf.pairwise_distance(hydrogens_repeated - 1,
                                          acceptors_tiled - 1, coord)
            donor_acceptor = rf.pairwise_distance(donors_repeated - 1,
                                               acceptors_tiled - 1, coord)
            donor_hydrogen = rf.pairwise_distance(hydrogens_repeated - 1,
                                               donors_repeated - 1, coord)
            acceptor_x = rf.pairwise_distance(acceptors_tiled - 1, x_repeated - 1,
                                           coord)
            x_hydrogen = rf.pairwise_distance(hydrogens_repeated - 1,
                                           x_repeated - 1, coord)
            # Calculate angles
            cos_alpha = np.clip((donor_hydrogen**2 + acceptor_hydrogen**2 - 
                                 donor_acceptor**2)/(2*acceptor_hydrogen*
                                                     donor_hydrogen), -1, 1)
            cos_beta = np.clip((acceptor_x**2 + acceptor_hydrogen**2 - 
                                 x_hydrogen**2)/(2*acceptor_hydrogen*
                                                     acceptor_x), -1, 1)
            alpha = np.arccos(cos_alpha)*180/np.pi
            beta = np.arccos(cos_beta)*180/np.pi
            # Calculate scores
            r_score = np.zeros(len(hydrogens_repeated), dtype=np.float32)
            a_score = np.zeros(len(hydrogens_repeated), dtype=np.float32)
            b_score = np.zeros(len(hydrogens_repeated), dtype=np.float32)

            for r, a, b, ind in zip(acceptor_hydrogen*10, alpha, beta,
                                    range(len(hydrogens_repeated))):
                if r >= r_min_ideal and r <= r_max_ideal:
                    r_score[ind] = 1
                elif r >= r_min and r <= r_min_ideal:
                    r_score[ind] = 1 - (r_min_ideal - r)/(r_min_ideal - r_min)
                elif r >= r_max_ideal and r <= r_max:
                    r_score[ind] = 1 - (r - r_max_ideal)/(r_max - r_max_ideal)
                else:
                    r_score[ind] = 0

                if a >= a_min_ideal and a <= a_max_ideal:
                    a_score[ind] = 1
                elif a >= a_min and a <= a_min_ideal:
                    a_score[ind] = 1 - (a_min_ideal - a)/(a_min_ideal - a_min)
                elif a >= a_max_ideal and a <= a_max:
                    a_score[ind] = 1 - (a - a_max_ideal)/(a_max - a_max_ideal)
                else:
                    a_score[ind] = 0
                
                if b >= b_min_ideal and b <= b_max_ideal:
                    b_score[ind] = 1
                elif b >= b_min and b <= b_min_ideal:
                    b_score[ind] = 1 - (b_min_ideal - b)/(b_min_ideal - b_min)
                elif b >= b_max_ideal and b <= b_max:
                    b_score[ind] = 1 - (b - b_max_ideal)/(b_max - b_max_ideal)
                else:
                    b_score[ind] = 0
            # Checking charge
            charge_factor = np.zeros(len(hydrogens_repeated), dtype=np.float32)
            for donor, acceptor, ind in zip(donors_repeated-1, acceptors_tiled-1,
                                            range(len(donors_repeated))):
                
                if charges[donor] != 0 and charges[acceptor] != 0:
                    charge_factor[ind] = -3
                elif charges[donor] == 0 and charges[acceptor] == 0:
                    charge_factor[ind] = -1
                else:
                    charge_factor[ind] = -2
            hbonds_matrix[residue1][residue2] = np.sum(r_score*a_score*b_score*charge_factor)
            hbonds_matrix[residue2][residue1] = np.sum(r_score*a_score*b_score*charge_factor)
    return hbonds_matrix

def compute_coulomb(traj, top_info, cutoff):
    """
    Compute electrostatic interaction with taking into account the dielectric
    constant of the media (from 10.1021/acs.jctc.7b00125)
    """
    n_excluded = np.asarray(top_info['NUMBER_EXCLUDED_ATOMS'])
    excluded_atoms = np.asarray(top_info['EXCLUDED_ATOMS_LIST'])
    residue_pointer = np.asarray(top_info['RESIDUE_POINTER']) - 1
    charges = np.asarray(top_info['CHARGE'])
    names = np.asarray(top_info['ATOM_NAME'])
    charges /= 18.2223
    coulomb_constant = 322
    max_dist = 5.5

    coord = traj.xyz[0]
    coulomb_matrix = np.zeros((traj.n_residues, traj.n_residues))
    for residue1 in range(traj.n_residues):
        atoms1 = traj.topology.select(f'resid {residue1}')
        to_delete = []
        for atom_id, atom in enumerate(atoms1):
            if names[atom] in ['C', 'N', 'O', 'CA']:
                to_delete.append(atom_id)
        atoms1 = np.delete(atoms1, to_delete)
        frame_tree = kdtree(coord)
        in_cutoff = np.concatenate(frame_tree.query_radius(coord[atoms1],
                                                           r=cutoff))
        neighbors = [np.where(residue_pointer <= atom)[0][-1]
                     for atom in in_cutoff]
        neighbors = np.unique(np.asarray(neighbors))
        for residue2 in neighbors[np.where(neighbors > residue1)[0]]:
            atoms2 = traj.topology.select(f'resid {residue2}')
            to_delete = []
            for atom_id, atom in enumerate(atoms2):
                if names[atom] in ['C', 'N', 'O', 'CA']:
                    to_delete.append(atom_id)
            atoms2 = np.delete(atoms2, to_delete)
            atoms1_repeated = np.repeat(atoms1, len(atoms2))
            atoms2_tiled = np.tile(atoms2, len(atoms1))
            coulomb_atoms1 = np.zeros(len(atoms1_repeated))
            atom_id = 0
            for atom_i, atom_j in zip(atoms1_repeated, atoms2_tiled):
                r2 = np.square(coord[atom_i] - coord[atom_j])
                r = np.sqrt(r2.sum())
                exclusion = np.asarray([])
                if atom_j > atom_i:
                    higher_atom = atom_j
                    if n_excluded[atom_i] != 1: # 1 stands for non exclusion
                        exclude_from = n_excluded[:atom_i].sum()
                        exclude_to = n_excluded[:atom_i + 1].sum()
                        exclusion = excluded_atoms[exclude_from:exclude_to]

                else: # check minimum index for excluded atoms
                    higher_atom = atom_i
                    if n_excluded[atom_j] != 1:
                        exclude_from = n_excluded[:atom_j].sum()
                        exclude_to = n_excluded[:atom_j + 1].sum()
                        exclusion = excluded_atoms[exclude_from : exclude_to]
                if higher_atom + 1 not in exclusion:
                    # Calculate vdW
                    if r > max_dist:
                        coulomb_atoms1[atom_id] = 0
                    elif r < 1.45:
                        coulomb_ = (coulomb_constant*charges[atom_i]*charges[atom_j])/(1.45*4)
                    else:
                        coulomb_ = (coulomb_constant*charges[atom_i]*charges[atom_j])/(r*4)
                    coulomb_atoms1[atom_id] = coulomb_
                atom_id += 1

            coulomb_matrix[residue1][residue2] = np.sum(coulomb_atoms1)
            coulomb_matrix[residue2][residue1] = np.sum(coulomb_atoms1)
    return coulomb_matrix

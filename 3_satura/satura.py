#!/bin/python

"""
@authors: Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sánchez-Murcia [pedro.murcia@medunigraz.at]
"""

"""
Satura:
Scoring of mutations based on MSA model (evcouplings) and sequence model (ESM)
"""

# General imports
import os
import heapq
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# satura_functions.py import
import satura_functions as sf

print("""
********************************************************************************
* SmarTSzyme-satura:                                                           *
*     Scoring single point mutations for the selected residues                 *
********************************************************************************
""")

args = sf.parse_arguments()
sf.mplstyle()
aa_list = np.asarray(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                     'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])

os.mkdir(args.output)

if not args.positions:
    # Read reduce results ==================================================
    flux_pd = pd.read_csv(args.flux_file)
    flux_total = np.asarray(flux_pd.iloc[:, 1])

    flux_positive = np.copy(flux_total)
    flux_positive[np.where(flux_positive < 0)[0]] = 0

    flux_negative = np.copy(flux_total)
    flux_negative[np.where(flux_negative > 0)[0]] = 0

    del flux_total

        
    # Select residues for saturation =======================================
    ordered_positive = []
    ordered_negative = []

    heapq.heapify(ordered_negative)
    heapq.heapify(ordered_positive)

    for res_id, flux in enumerate(flux_negative):
        heapq.heappush(ordered_negative, (flux, res_id + 1))
    n_negative = heapq.nsmallest(args.max_residues, ordered_negative)
    n_negative_res = np.asarray([neg[1] for neg in n_negative])
    n_negative_flux = np.asarray([neg[0] for neg in n_negative])

    for res_id, flux in enumerate(flux_positive):
        heapq.heappush(ordered_positive, (flux, res_id + 1))
    n_positive = heapq.nlargest(args.max_residues, ordered_positive)
    n_positive_res = np.asarray([pos[1] for pos in n_positive])
    n_positive_flux = np.asarray([pos[0] for pos in n_positive])

    fig, ax = plt.subplots()
    ax.bar(np.arange(len(n_negative_res)), -n_negative_flux, color='blue')
    ax.set_xlabel(r'Residues')
    ax.set_xticks(np.arange(len(n_negative_res)))
    ax.set_xticklabels(n_negative_res, rotation=45, size=10)
    ax.set_ylabel(r'|Normalized flux|')
    fig.savefig(f'{args.output}/top_negative_flux.png')

    fig, ax = plt.subplots()
    ax.bar(np.arange(len(n_positive_res)), n_positive_flux, color='red')
    ax.set_xlabel(r'Residues')
    ax.set_xticks(np.arange(len(n_positive_res)))
    ax.set_xticklabels(n_positive_res, rotation=45, size=10)
    ax.set_ylabel(r'|Normalized flux|')
    fig.savefig(f'{args.output}/top_positive_flux.png')

    positions = n_positive_res

else:
    positions = args.positions

# Predicting models
total_scores = np.zeros((len(positions), 20), dtype=float)

if 'esm' in args.model and 'None' not in args.model:
    print('Running ESM model(s) ...')
    # Read sequence from fasta
    with open(args.fasta, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('>'):
                pass
            else:
                sequence = list(str(line.strip()))

    column = ['mutation']
    data = []
    positions_esm = []
    original_aas = []

    for position in positions:
        if position > len(sequence):
            position_first_chain = (position - len(sequence)) + args.offset_chain
            positions_esm.append(position_first_chain)
            original_aas.append(sequence[position_first_chain - 1])
            print(f"Position {position} equivalent to {position_first_chain} in the provided chain.")
            for aa in aa_list:
                data.append(f'{sequence[position_first_chain - 1]}{position_first_chain}{aa}')
        else:
            positions_esm.append(position)
            original_aas.append(sequence[position - 1])
            for aa in aa_list:
                data.append(f'{sequence[position - 1]}{position}{aa}')
    df = pd.DataFrame(data=data, columns=column)
    sequence = ''.join(sequence)


    esm_df = sf.run_esm(sequence, df, column, arguments=args)

    # Plot ESM results
    scores = np.asarray(esm_df.iloc[:, -1]).reshape((-1, 20))
    np.savetxt(f'{args.output}/esm.csv', scores, delimiter=',')
    total_scores += scores
    del esm_df

    min_val = scores.min()
    max_val = scores.max()
    cbar_limit = abs(min_val) if abs(min_val) > abs(max_val) else abs(max_val)
    fig, ax = plt.subplots()
    sns.heatmap(scores.T, ax=ax, cmap='coolwarm_r', center=0, vmin=-cbar_limit,
                vmax=cbar_limit)
    ax.set_yticks(np.arange(0.5, len(aa_list) + 0.5, 1))
    ax.set_yticklabels(aa_list, rotation=0)
    ax.set_ylabel(r'Mutation')
    ax.set_xticks(np.arange(0.5, len(positions_esm) + 0.5, 1))
    ax.set_xticklabels(positions_esm, rotation=45)
    ax.set_xlabel(r'Residues')

    for pos_idx, pos in enumerate(positions_esm):
        original_aa = np.where(aa_list == original_aas[pos_idx])[0]
        y_pos = np.arange(0.5, len(aa_list) + 0.5, 1)[original_aa]
        x_pos = np.arange(0.5, len(positions_esm) + 0.5, 1)[pos_idx]
        ax.scatter(x_pos, y_pos, color='black', marker='o')
    fig.savefig(f'{args.output}/single_mutations_esm.png')

if 'evcouplings' in args.model and 'None' not in args.model:
    print('Running EVcouplings model ...')
    from evcouplings.couplings import CouplingsModel

    positions_evc = []

    for position in positions:
        if position > len(sequence):
            position_first_chain = (position - len(sequence)) + args.offset_chain
            positions_evc.append(position_first_chain)
            print(f"Position {position} equivalent to {position_first_chain} in the provided chain.")
        else:
            positions_evc.append(position)

    model = CouplingsModel(args.ev_model)

    # Single mutations
    scores = np.zeros((len(positions_evc), 20))
    for pos_idx, pos in enumerate(positions_evc):
        try:
            score_ = model.smm(pos)
            scores[pos_idx] = score_
        except KeyError:
            print(f'Position {pos} not present in the model... skipping')
            pass
    np.savetxt(f'{args.output}/evmut.csv', scores, delimiter=',')
    #scores = model.smm(positions)
    total_scores += scores

    # Plot results
    min_val = scores.min()
    max_val = scores.max()
    cbar_limit = abs(min_val) if abs(min_val) > abs(max_val) else abs(max_val)
    fig, ax = plt.subplots()
    sns.heatmap(scores.T, ax=ax, center=0, cmap='coolwarm_r',
                vmin=-cbar_limit, vmax=cbar_limit)
    ax.set_ylabel(r'Mutations')
    ax.set_yticks(np.arange(0.5, len(aa_list) + 0.5, 1))
    ax.set_yticklabels(aa_list, rotation=0)

    ax.set_xlabel(r'Residues')
    ax.set_xticks(np.arange(0.5, len(positions) + 0.5, 1))
    ax.set_xticklabels(positions, rotation=45)

    # Get original Amino Acids
    for pos_idx, pos in enumerate(positions):
        try:
            original_aa = np.where(aa_list == model.seq(pos))[0]
            y_pos = np.arange(0.5, len(aa_list) + 0.5, 1)[original_aa]
            x_pos = np.arange(0.5, len(positions) + 0.5, 1)[pos_idx]
            ax.scatter(x_pos, y_pos, color='black', marker='o')
        except KeyError:
            pass
    fig.savefig(f'{args.output}/single_mutations_evc.png')

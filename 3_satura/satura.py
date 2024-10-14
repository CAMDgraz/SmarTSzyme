#!/bin/python

"""
@autors: Daniel Platero Rochart [daniel.platero-rochart@medunigraz.at]
         Pedro A. Sánchez-Murcia [pedro.murcia@medunigraz.at]
"""

"""
Satura:
Scoring of mutations based on MSA model (evcoupling) and sequence model (ESM)
"""

# General imports
import os
import heapq
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Satura imports
import satura_functions as sf

# ESM imports
import torch
from esm import pretrained
from tqdm import tqdm

print("""
********************************************************************************
* SmarTSzyme-satura:                                                           *
*     Scoring single point mutations based on the selected residues            *
********************************************************************************
""")

args = sf.parse_arguments()
sf.mplstyle()
aa_list = np.asarray(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                     'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])

if not args.positions:
    # Read reduce results ==========================================================
    flux_pd = pd.read_csv(args.flux_file)
    flux_total = np.asarray(flux_pd.iloc[:, 1])

    flux_positive = np.copy(flux_total)
    flux_positive[np.where(flux_positive < 0)[0]] = 0

    flux_negative = np.copy(flux_total)
    flux_negative[np.where(flux_negative > 0)[0]] = 0

    del flux_total

    os.mkdir(args.output)

    # Select residues for next step =================================================
    ordered_positive = []
    ordered_negative = []

    heapq.heapify(ordered_negative)
    heapq.heapify(ordered_positive)

    for res_id, flux in enumerate(flux_negative):
        heapq.heappush(ordered_negative, (flux, res_id + 1))
    nnegative = heapq.nsmallest(args.max_residues, ordered_negative)
    nnegative_res = np.asarray([neg[1] for neg in nnegative])
    nnegative_flux = np.asarray([neg[0] for neg in nnegative])

    for res_id, flux in enumerate(flux_positive):
        heapq.heappush(ordered_positive, (flux, res_id + 1))
    npositive = heapq.nlargest(args.max_residues, ordered_positive)
    npositive_res = np.asarray([pos[1] for pos in npositive])
    npositive_flux = np.asarray([pos[0] for pos in npositive])

    fig, ax = plt.subplots()
    ax.bar(np.arange(len(nnegative_res)), -nnegative_flux, color='blue')
    ax.set_xlabel(r'Residues')
    ax.set_xticks(np.arange(len(nnegative_res)))
    ax.set_xticklabels(nnegative_res, rotation=45, size=10)
    ax.set_ylabel(r'|Normalized flux|')
    fig.savefig(f'{args.output}/top_negative_flux.png')

    fig, ax = plt.subplots()
    ax.bar(np.arange(len(npositive_res)), npositive_flux, color='red')
    ax.set_xlabel(r'Residues')
    ax.set_xticks(np.arange(len(npositive_res)))
    ax.set_xticklabels(npositive_res, rotation=45, size=10)
    ax.set_ylabel(r'|Normalized flux|')
    fig.savefig(f'{args.output}/top_positive_flux.png')

    positions = npositive_res

else:
    positions = args.positions

if 'esm' in args.model:
    print('Running ESM model(s)')
    # Read sequence from fasta
    with open(args.fasta, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('>'):
                pass
            else:
                sequence = list(str(line))

    column = ['mutation']
    data = []
    positions_esm = []
    original_aas = []

    for position in positions:
        if position > 155:
            continue
        else:
            positions_esm.append(position)
            original_aas.append(sequence[position - 1])
        for aa in aa_list:
            data.append(f'{sequence[position - 1]}{position}{aa}')
    df = pd.DataFrame(data=data, columns=column)
    sequence = ''.join(sequence)

    esm_df = sf.run_esm(sequence, df, args=args)

    # Plot ESM results
    scores = np.asarray(df.iloc[:, -1]).reshape((-1, 20))
    del df

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
    ax.set_xticklabels(positions_esm, rotation=0)
    ax.set_xlabel(r'Residues')

    for pos_idx, pos in enumerate(positions_esm):
        original_aa = np.where(aa_list == original_aas[pos_idx])[0]
        y_pos = np.arange(0.5, len(aa_list) + 0.5, 1)[original_aa]
        x_pos = np.arange(0.5, len(positions_esm) + 0.5, 1)[pos_idx]
        ax.scatter(x_pos, y_pos, color='black', marker='o')
    fig.savefig(f'{args.output}/single_mutations.png')

if 'evcouplings' in args.model:
    print('Under development!!!!!!')
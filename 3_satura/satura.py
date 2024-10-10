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
import matplotlib.pyplot as plt

# Satura imports
import satura_functions as sf

# ESM imports
import torch
from esm import Alphabet, pretrained, MSATransformer
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import itertools

print("""
********************************************************************************
* SmarTSzyme-satura:                                                           *
*     Scoring single point mutations based on the selected residues            *
********************************************************************************
""")

aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
           'R', 'S', 'T', 'V', 'W', 'Y']
args = sf.parse_arguments()

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

# Run ESM ======================================================================
with open(args.fasta, 'r') as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith('>'):
            pass
        else:
            sequence = list(str(line))

# Provisional, remove the second chain residues #####
column = ['mutation']
data = []
for position in npositive_res:
    if position > 155:
        continue
    for aa in aa_list:
        data.append(f'{sequence[position - 1]}{position}{aa}')
df = pd.DataFrame(data=data, columns=column)

device = torch.device("cuda" if torch.cuda.is_available() and
                      not args.nogpu else "cpu")

# load models
for model_location in args.model_location:
    model, alphabet = pretrained.load_model_and_alphabet(model_location)
    model = model.to(device)
    model.eval()
    batch_converter = alphabet.get_batch_converter()


    data = [("protein1", sequence)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(device)

    if args.scoring_strategy == "wt-marginals":
        with torch.no_grad():
            token_probs = torch.log_softmax(model(batch_tokens)["logits"],
                                            dim=-1)
        df[model_location] = df.apply(
            lambda row: sf.label_row(row[args.mutation_col], args.sequence,
                                     token_probs, alphabet, args.offset_idx),
                                    axis=1)
    elif args.scoring_strategy == "masked-marginals":
        all_token_probs = []
        for i in tqdm(range(batch_tokens.size(1))):
            batch_tokens_masked = batch_tokens.clone()
            batch_tokens_masked[0, i] = alphabet.mask_idx
            with torch.no_grad():
                token_probs = torch.log_softmax(model(batch_tokens_masked)["logits"],
                                                dim=-1)
            all_token_probs.append(token_probs[:, i])
        token_probs = torch.cat(all_token_probs, dim=0).unsqueeze(0)
        df[model_location] = df.apply(
            lambda row: sf.label_row(row[args.mutation_col], args.sequence,
                                     token_probs, alphabet, args.offset_idx),
            axis=1)

df.to_csv('out_test_esm.csv')
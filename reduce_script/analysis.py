"""
@autors: Daniel Platero Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sánchez-Murcia [pedro.murcia@medunigraz.at]
"""

"""
Analaysis of SmarTSzyme's results
"""
print("""
********************************************************************************
*                            SmarTSzyme-analysis                               *       
*                      Analysis of SmarTSzyme's results                        *
********************************************************************************
""")

# Parser arguments from the CLI
from cli import parse_analysis
import sys

if len(sys.argv) == 1:
    parse_analysis().print_help()
    sys.exit(1)
args = parse_analysis().parse_args()

# Remaining imports
import glob
import numpy as np
import pandas as pd
import heapq
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import utils as ut

# Input options ================================================================
print(f'*** User inputs ***')
print(f'Reading matrices from:          {args.int_path}')
print(f'Reading coupling from:          {args.coupling_path}')
if args.exp_res:
    print(f'Checking experimental residues: {args.exp_res}')
print(f'Output top:                     {args.topn}')
print(f'Batch size:                     {args.batch}')
# ==============================================================================

# Check interaction matrices and coupling files ================================
# Get number of couplings
coupling_stab = glob.glob(f'{args.coupling_path}/stab_flux.*.csv')
coupling_dest = glob.glob(f'{args.coupling_path}/dest_flux.*.csv')

# Get interaction matrices 
matrices_id = ut.find_matrices(args.int_path)

if len(coupling_dest) != len(coupling_stab):
    print('\nError!!! The number of stabilizing and destabilizing files differ')
    sys.exit()
elif len(coupling_stab) != len(matrices_id):
    print('\nError!!! The number of interaction matrices and couplings differ')
    sys.exit()

batches = np.arange(args.batch, len(coupling_stab) + args.batch, args.batch)
print(f'No. of couplings file:          {len(coupling_stab)}')
print(f'No. Batches:                    {len(batches)}')

# ==============================================================================

# Main analysis ================================================================
score_batches = np.zeros((len(batches), args.nres))

print(f'\n*** Analysing results ***')
for batch_id, batch in enumerate(batches):
    # Variables for stab and destab contributions
    stab_dest = np.zeros((batch, args.nres))
    flux_stab_batch = np.zeros((batch, args.nres))
    flux_dest_batch = np.zeros((batch, args.nres))

    for job_id, job in enumerate(matrices_id[:batch]):
        results_stab = pd.read_csv(f'{args.coupling_path}/stab_flux.{job}.csv',
                                   sep=',')
        results_dest = pd.read_csv(f'{args.coupling_path}/dest_flux.{job}.csv',
                                   sep=',')
        flux_stab_batch[job] = np.asarray(results_stab.iloc[:, 1])
        flux_dest_batch[job] = np.asarray(results_dest.iloc[:, 1])

        # get es and pts interactions
        es_int = np.zeros((args.nres, args.nres))
        pts_int = np.zeros((args.nres, args.nres))

        es_files = glob.glob(f'{args.int_path}/es_*.{job}.pickle')
        pts_files = glob.glob(f'{args.int_path}/pts_*.{job}.pickle')

        for es_file in es_files:
            es_int += ut.load_pickle(es_file)
        for pts_file in pts_files:
            pts_int += ut.load_pickle(pts_file)
        

        diff_matrix = pts_int - es_int
        stab_dest[job_id] = diff_matrix.sum(axis=0)/abs(diff_matrix.sum(axis=0))

    dscore = (stab_dest > 0).sum(axis=0)/batch
    sscore = (stab_dest < 0).sum(axis=0)/batch
    dest_score = (dscore*flux_dest_batch).sum(axis=0) 
    stab_score = (sscore*flux_stab_batch).sum(axis=0)
    score_batches[batch_id] = dest_score - stab_score
    print(f'Batch {batch_id + 1}: Done')

# Get topn
print(f'\n*** Getting top {args.topn} residues ***')
top_stab_residues = np.zeros(score_batches.shape)
top_dest_residues = np.zeros(score_batches.shape)

for batch_id, score in enumerate(score_batches):
    ordered = []
    heapq.heapify(ordered)
    for res in range(args.nres):
        heapq.heappush(ordered, (score[res], res))
    
    nlargest = [i[1] for i in heapq.nlargest(args.topn, ordered)]
    nsmallest = [i[1] for i in heapq.nsmallest(args.topn, ordered)]

    top_dest_residues[batch_id, nlargest] = 1
    top_stab_residues[batch_id, nsmallest] = 1

print(f'Maximum mutatibility score:')
for batch_id, res_batch in enumerate(top_stab_residues):
    print(f'Batch {batch_id}: {np.where(res_batch == 1)[0] + 1}')

print(f'Minimum mutatibility score:')
for batch_id, res_batch in enumerate(top_dest_residues):
    print(f'Batch {batch_id}: {np.where(res_batch == 1)[0] + 1}')

print('Done')            

# Plotting and writting results ================================================
print('\n*** Saving results ***')

# Write score
for batch_id, score_batch in enumerate(score_batches):
    with open(f'{args.output}/score_{batch_id}.csv', 'w') as f:
        f.write('Residue,Score\n')
        for res, score in enumerate(score_batch):
            f.write(f'{res + 1},{score}\n')

# Plot top score
mpl.rcParams['text.usetex'] = False
mpl.rcParams['svg.fonttype'] = 'none'

top_stab = np.unique(np.where(top_stab_residues == 1)[1])
top_dest = np.unique(np.where(top_dest_residues == 1)[1])

# Maximum mutability
fig, ax = plt.subplots()
sns.heatmap(top_dest_residues[:, top_dest], cmap='cubehelix_r', vmin=0, vmax=1,
            center=0.5, cbar=False, linewidths=0.5, linecolor='grey', ax=ax,
            square=True)
ax.set_yticks(np.arange(0.5, len(batches) + 0.5, 1))
ax.set_yticklabels(batches, rotation=0, size=10)
ax.set_xticks(np.arange(0.5, len(top_dest) + 0.5, 1))
ax.set_xticklabels(top_dest + 1, rotation=45, size=10)

ax.set_ylabel(r'Number of Trajectories')
ax.set_xlabel(r'Residue id')
fig.savefig(f'{args.output}/top_max.svg')
plt.close(fig)

# Minimum mutability
fig, ax = plt.subplots()
sns.heatmap(top_stab_residues[:, top_stab], cmap='cubehelix_r', vmin=0, vmax=1,
            center=0.5, cbar=False, linewidths=0.5, linecolor='grey', ax=ax,
            square=True)
ax.set_yticks(np.arange(0.5, len(batches) + 0.5, 1))
ax.set_yticklabels(batches, rotation=0, size=10)
ax.set_xticks(np.arange(0.5, len(top_stab) + 0.5, 1))
ax.set_xticklabels(top_stab + 1, rotation=45, size=10)

ax.set_ylabel(r'Number of Trajectories')
ax.set_xlabel(r'Residue id')
fig.savefig(f'{args.output}/top_min.svg')
plt.close(fig)

# Experimental residues
if args.exp_res:        
    exp_res = np.asarray(args.exp_res, dtype=int) - 1
    exp_scores = score_batches[:, exp_res]
    with open(f'{args.output}/experimental.csv', 'w') as f:
        for res, scores in zip(exp_res, exp_scores.T):
            f.write(f'{res+1}')
            for score in scores:
                f.write(f',{score}')
            f.write('\n')

    max_val = score_batches[:, exp_res].max()
    min_val = score_batches[:, exp_res].min()

    vmax = max_val if max_val > abs(min_val) else abs(min_val)

    fig, ax = plt.subplots()
    sns.heatmap(score_batches[:, exp_res], cmap='vlag', center=0,
                ax=ax, vmin=-vmax, vmax=vmax, square=True)
    ax.set_yticks(np.arange(0.5, len(batches) + 0.5, 1))
    ax.set_yticklabels(batches, rotation=0, size=10)
    ax.set_xticks(np.arange(0.5, len(args.exp_res) + 0.5, 1))
    ax.set_xticklabels(exp_res + 1, rotation=45, size=10)

    ax.set_ylabel(r'Number of Trajectories')
    ax.set_xlabel(r'Residue id')
    fig.savefig(f'{args.output}/experimental.svg')
    plt.close(fig)
print('Done')

print('\n*** Normal termination of SmarTSzyme ***')
#!/bin/python

"""
@autors: Daniel Platero Rochart [daniel.platero-rochart@medunigraz.at]
         Pedro A. Sánchez-Murcia [pedro.murcia@medunigraz.at]
"""

"""
Functions used along satura.py script. The esm functions are taken from 
(https://huggingface.co/blog/AmelieSchreiber/mutation-scoring)
"""

# General imports
import matplotlib as mpl
import argparse

# ESM imports
import pathlib
import string

# Arguments parsering
def parse_arguments():
    """
    Parse arguments of the cli
    """
    desc = '''\nSatura: Scoring of single point mutations based on the position
                        identified by reduce'''
    parser = argparse.ArgumentParser(prog='Reduce',
                                     description = desc,
                                     add_help=True,
                                     allow_abbrev = False)
    inputs = parser.add_argument_group(title='Input options')
    inputs.add_argument('-file', dest='flux_file', action = 'store',
                        type=str, required = True,
                        help = 'Path to the csv file output of reduce.')
    inputs.add_argument('-maxres', dest = 'max_residues', action = 'store', 
                      help = 'Top n residues to select.',
                      type=int, required = True)
    inputs_esm = parser.add_argument_group(title='ESM inputs')
    inputs_esm.add_argument("-model-location", type=str,
                            help="Pretrained model name", nargs="+",
                            default=['esm1v_t33_650M_UR90S_1',
                                     'esm1v_t33_650M_UR90S_2',
                                     'esm1v_t33_650M_UR90S_3',
                                     'esm1v_t33_650M_UR90S_4',
                                     'esm1v_t33_650M_UR90S_5'])
    inputs_esm.add_argument("-fasta", type=str, 
                            help="fasta file containing the sequence")
    inputs_esm.add_argument("-scoring-strategy", type=str,
                            default="wt-marginals",
                            choices=["wt-marginals", "masked-marginals"])
    inputs_esm.add_argument("-nogpu", action="store_true",
                            help="Do not use GPU even if available")
    outputs = parser.add_argument_group(title='Output options')
    outputs.add_argument('-out', dest = 'output', type = str,
                         default = './', help = 'Output path.')
    user_inputs = parser.parse_args()
    return user_inputs

# ESM functions
def remove_insertions(sequence: str) -> str:
    deletekeys = dict.fromkeys(string.ascii_lowercase)
    deletekeys["."] = None
    deletekeys["*"] = None
    translation = str.maketrans(deletekeys)
    return sequence.translate(translation)

def label_row(row, sequence, token_probs, alphabet, offset_idx):
    wt, idx, mt = row[0], int(row[1:-1]) - offset_idx, row[-1]
    assert sequence[idx] == wt, "The listed wildtype does not match the provided sequence"
    wt_encoded, mt_encoded = alphabet.get_idx(wt), alphabet.get_idx(mt)
    score = token_probs[0, 1 + idx, mt_encoded] - token_probs[0, 1 + idx, wt_encoded]
    return score.item()

# Output
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

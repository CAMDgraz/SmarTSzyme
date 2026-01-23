"""
@autors: Daniel Platero Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sánchez-Murcia [pedro.murcia@medunigraz.at]
"""

"""
Parsing of the command line interface
"""

import argparse

def parse_reduce():
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
    
    return parser

def parse_analysis():
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
    return parser
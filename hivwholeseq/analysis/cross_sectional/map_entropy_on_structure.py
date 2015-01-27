# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/15
content:    Map site entropy onto protein structure for visualization.
'''
# Modules
import argparse

from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment_entropy)
from hivwholeseq.structure.get_PDB import get_PDB
from hivwholeseq.structure_utils import get_chainseq


# Functions



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Map subtype entropy onto structure',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--region', required=True,
                        help='Region to analyze (e.g. p17)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot results')

    args = parser.parse_args()
    region = args.region
    VERBOSE = args.verbose
    use_plot = args.plot


    if VERBOSE >= 1:
        print 'Get entropy'
    S = get_subtype_reference_alignment_entropy(region, subtype='B',
                                                refname='HXB2',
                                                type='aa')

    if VERBOSE >= 1:
        print 'Get structure PDB'
    pdb = get_PDB(region, VERBOSE=VERBOSE)

    # Note: PR is a dimer
    ch = list(pdb.get_chains())[0]

    # Note: they often put drugs after the proteins in the same PDB chain
    seq = get_chainseq(ch).rstrip('X')


    # Check for length
    if len(seq) != len(S):
        raise ValueError('Entropy and structure not aligned')

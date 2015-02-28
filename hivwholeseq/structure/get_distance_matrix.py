# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/15
content:    Get PDB of a structure.
'''
# Modules
import argparse
from Bio.PDB import PDBParser
import matplotlib.pyplot as plt

from hivwholeseq.utils.structure import get_chainseq, get_distance_matrix
from hivwholeseq.structure.get_PDB import get_PDB



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
        print 'Get structure PDB and index of chains'
    pdb, chinds = get_PDB(region, VERBOSE=VERBOSE)

    ch = list(pdb.get_chains())[chinds[0]]

    # Note: they often put drugs after the proteins in the same PDB chain
    seq = get_chainseq(ch).rstrip('X')

    ds = get_distance_matrix(ch)

    if use_plot:
        fig, ax = plt.subplots()
        h = ax.imshow(ds, interpolation='nearest')
        cb = plt.colorbar(h)

        ax.set_title(region)
        cb.set_label('Distance [Angstrom]', rotation=270, labelpad=30)

        plt.ion()
        plt.show()

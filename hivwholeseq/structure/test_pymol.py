# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/15
content:    Get PDB of a structure.
'''
# Modules
import os
import argparse
from Bio.PDB import PDBParser

from hivwholeseq.structure.filenames import get_PDB_filename



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
        print 'Get structure PDB'
    fn_pdb = get_PDB_filename(region, VERBOSE=VERBOSE)


    # NOTE: pymol_manager is to be used as a script, not a module
    # (there are all kinds of race conditions)

    import ipymol
    mol=ipymol.MolViewer()
    mol.server.do('fetch 3odu; as cartoon; bg white; png /home/fabio/Desktop/test.png')

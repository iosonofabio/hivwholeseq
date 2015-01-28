# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/15
content:    Get PDB of a structure.
'''
# Modules
import argparse
from Bio.PDB import PDBParser

from hivwholeseq.utils.structure import get_chainseq



# Functions
def get_PDB(region, seqid=None, VERBOSE=0):
    '''Get a PDB structure from file'''
    from hivwholeseq.structure.filenames import get_PDB_filename
    fn = get_PDB_filename(region, seqid=seqid, VERBOSE=VERBOSE)
    pdbp = PDBParser()
    pdb = pdbp.get_structure(region, fn)
    return pdb



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
    pdb = get_PDB(region, VERBOSE=VERBOSE)

    # Note: PR is a dimer
    ch = list(pdb.get_chains())[0]

    # Note: they often put drugs after the proteins in the same PDB chain
    seq = get_chainseq(ch).rstrip('X')

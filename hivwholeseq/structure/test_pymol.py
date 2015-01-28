# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/15
content:    Get PDB of a structure.
'''
# Modules
import os
import argparse
import numpy as np
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

    # Get subtype entropy
    from hivwholeseq.cross_sectional.get_subtype_entropy import (
        get_subtype_reference_alignment_entropy)
    S = get_subtype_reference_alignment_entropy(region, subtype='B',
                                                refname='HXB2',
                                                type='aa')

    #S = np.log(S + 1e-2)
    vdws = 0.1 + 3.0 * (S - S.min()) / (S.max() - S.min())

    # NOTE: pymol_manager is to be used as a script, not a module
    # (there are all kinds of race conditions)

    import ipymol
    mol=ipymol.MolViewer()
    mol.server.do('load '+fn_pdb+';')
    mol.server.do('zoom; as cartoon; show spheres, chain A; hide spheres, resn HOH')
    for pos, vdw in enumerate(vdws):
        #mol.server.do('alter resi '+str(pos+1)+', vdw='+str(vdw)+';')
        #mol.server.do('color '+'blue'+', resi '+str(pos+1)+';')
    mol.server.do('rebuild;')
    mol.server.do('bg white; png /home/fabio/Desktop/'+region+'.png')

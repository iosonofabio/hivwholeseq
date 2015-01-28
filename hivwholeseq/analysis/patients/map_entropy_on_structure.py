# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/15
content:    Map site entropy onto protein structure for visualization.
'''
# Modules
import os
import sys
import argparse
import numpy as np

from hivwholeseq.patients.patients import load_patient
from hivwholeseq.structure.get_PDB import get_PDB
from hivwholeseq.structure.filenames import get_PDB_filename
from hivwholeseq.structure_utils import get_chainseq
from hivwholeseq.one_site_statistics import get_entropy



# Globals
times = [('early', 1),
         ('midway', -5),
         ('late', -1)]


# Functions



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Map subtype entropy onto structure',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--region', required=True,
                        help='Region to analyze (e.g. p17)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot results')

    args = parser.parse_args()
    pname = args.patient
    region = args.region
    VERBOSE = args.verbose
    use_plot = args.plot
    
    patient = load_patient(pname)
    
    if VERBOSE >= 1:
        print 'Get allele frequencies'
    aft, ind = patient.get_allele_frequency_trajectories(region, cov_min=100)

    for ii, (tlabel, i) in enumerate(times, 1):
        if VERBOSE >= 1:
            print 'Get entropy'
        # FIXME: take only last time point for now
        af = aft[i]
        # FIXME: We have nucleotides and want amino acids, but ok for now
        Snuc = get_entropy(af)
        S = Snuc.reshape((Snuc.shape[0] // 3, 3)).mean(axis=1)


        if VERBOSE >= 1:
            print 'Get structure PDB filename'
        fn_pdb = get_PDB_filename(region, VERBOSE=VERBOSE)

        # Get subtype entropy
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
        

        # NOTE: pymol_manager is to be used as a script, not a module
        # (there are all kinds of race conditions)

        ## Plot with different radii based on entropy
        #vdws = 0.1 + 3.0 * (S - S.min()) / (S.max() - S.min())
        #import ipymol
        #mol=ipymol.MolViewer()
        #mol.server.do('load '+fn_pdb+';')
        #mol.server.do('zoom; as cartoon; show spheres, chain A; hide spheres, resn HOH')
        #for pos, vdw in enumerate(vdws):
        #    mol.server.do('alter resi '+str(pos+1)+', vdw='+str(vdw)+';')
        #    #mol.server.do('color '+'blue'+', resi '+str(pos+1)+';')
        #mol.server.do('rebuild;')
        #mol.server.do('bg white; png /home/fabio/Desktop/'+region+'.png')

        # Plot with different colors based on entropy
        x = np.log10(S + 1e-5)
        xmax = np.log(2)
        xmin = -5
        color_levels = 100 * (x - xmin) / (xmax - xmin)

        # NOTE: pymol_manager is to be used as a script, not a module
        # (there are all kinds of race conditions)

        import ipymol
        mol=ipymol.MolViewer()
        mol.server.do('load '+fn_pdb+';')
        mol.server.do('zoom center, 20; as cartoon;')
        for pos, Spos in enumerate(S):
            mol.server.do('alter resi '+str(pos+1)+', b='+str(color_levels[pos])+';')
        mol.server.do('rebuild;')
        mol.server.do('spectrum b, rainbow, minimum=0, maximum=100;')
        mol.server.do('bg white; png /home/fabio/Desktop/'+region+'_'+pname+'_'+str(ii)+tlabel+'.png')

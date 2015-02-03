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
from hivwholeseq.utils.structure import get_chainseq
from hivwholeseq.one_site_statistics import get_entropy
from hivwholeseq.utils import ipymol



# Globals



# Functions



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Map subtype entropy onto structure',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--protein', required=True,
                        help='Protein to analyze (e.g. p17)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot results')

    args = parser.parse_args()
    pname = args.patient
    protein = args.protein
    VERBOSE = args.verbose
    use_plot = args.plot
    
    # NOTE: pymol_manager is to be used as a script, not a module
    # (there are all kinds of race conditions)

    if VERBOSE >= 1:
        print 'Get structure PDB filename'
    fn_pdb = get_PDB_filename(protein, VERBOSE=VERBOSE)

    if VERBOSE >= 1:
        print 'Get structure PDB'
    pdb = get_PDB(protein, VERBOSE=VERBOSE)

    # Note: PR is a dimer
    ch = list(pdb.get_chains())[0]

    # Note: they often put drugs after the proteins in the same PDB chain
    seq = get_chainseq(ch).rstrip('X')

    if VERBOSE >= 1:
        print 'Load protein into pymol'
    mol=ipymol.MolViewer()
    cmd = ('delete all; ' +
           'load '+fn_pdb+'; ' +
           'center; ' +
           'zoom center, 20; as cartoon; ')
    if VERBOSE >= 2:
        print cmd
    mol.server.do(cmd)

    # Go to the patient
    patient = load_patient(pname)
    
    if VERBOSE >= 1:
        print 'Get allele frequencies of amino acids'
    aft, ind = patient.get_allele_frequency_trajectories_aa(protein, cov_min=100,
                                                            depth_min=50)

    Ss = []

    for it, af in enumerate(aft):
        if VERBOSE >= 1:
            print 'Time point n.'+str(it+1)

        if VERBOSE >= 1:
            print 'Get entropy'
        S = get_entropy(af)
        Ss.append(S)

        # Check for length
        if len(seq) != len(S):
            raise ValueError('Entropy and structure not aligned')
        
        # Output parameters
        #folder_out = '/home/fabio/Desktop/'
        folder_out = '/home/fabio/university/phd/talks/IGIM_2015/figures/animations/'
        fn_out = folder_out+protein+'/'+pname+'_'+'{:d}'.format(it)+'.png'

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
        #mol.server.do('bg white; png /home/fabio/Desktop/'+protein+'.png')

        # Plot with different colors based on entropy
        x = np.log10(S + 1e-5)
        xmax = np.log(2)
        xmin = -5
        color_levels = 100 * (x - xmin) / (xmax - xmin)

        for pos, Spos in enumerate(S):
            cmd = 'alter resi '+str(pos+1)+', b='+str(color_levels[pos])+'; '
            if VERBOSE >= 2:
                print cmd
            mol.server.do(cmd)
        cmd = ('rebuild; ' +
                'spectrum b, rainbow, minimum=0, maximum=100; ' +
                'bg white; png '+fn_out+'; ')
        if VERBOSE >= 2:
            print cmd
        mol.server.do(cmd)

    Ss = np.array(Ss)


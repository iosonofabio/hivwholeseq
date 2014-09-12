# vim: fdm=marker
'''
author:     Fabio Zanini
date:       12/09/14
content:    Plot the phylogenetic tree of the consensus sequences.
'''
# Modules
import sys
import os
import argparse
import numpy as np
from Bio import SeqIO, AlignIO, Phylo
import matplotlib.pyplot as plt

from hivwholeseq.sequencing.samples import SampleSeq
from hivwholeseq.patients.patients import load_patients, Patient, SamplePat
from hivwholeseq.patients.filenames import get_tree_consensi_figure_filename

# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Align consensi',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=['all'],
                        help='Patients to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Store tree plot to file')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    use_save = args.save

    patients = load_patients()
    if pnames != ['all']:
        patients = patients.iloc[patients.index.isin(pnames)]

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        for fragment in fragments:
            if VERBOSE >= 1:
                print pname, fragment

            fn = patient.get_consensi_alignment_filename(fragment)
            patient.ali = AlignIO.read(fn, 'fasta')
            tree = Phylo.read(patient.get_consensi_tree_filename(fragment),
                              'newick')

            fig, ax = plt.subplots(figsize=(15, 12))
            Phylo.draw(tree, do_show=False, axes=ax)
            ax.set_title(pname+', '+fragment)

            x_max = max(tree.depths().itervalues())
            ax.set_xlim(0.995, 0.995 + (x_max - 0.995) * 1.4)
            ax.grid(True)

            if use_save:
                fn = get_tree_consensi_figure_filename(pname, fragment)
                from hivwholeseq.generic_utils import mkdirs
                mkdirs(os.path.dirname(fn))
                fig.savefig(fn)
                plt.close(fig)
            
    if not use_save:
        plt.ion()
        plt.show()






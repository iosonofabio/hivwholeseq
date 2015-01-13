# vim: fdm=marker
'''
author:     Fabio Zanini
date:       13/01/15
content:    Get the tree of consensi from a patient.
'''
# Modules
import sys
import os
import argparse
from operator import attrgetter
import numpy as np
from Bio import SeqIO, AlignIO

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.sequencing.samples import SampleSeq
from hivwholeseq.patients.patients import load_patients, Patient, SamplePat
from hivwholeseq.tree_utils import build_tree_fasttree



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Align consensi',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=['all'],
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='*',
                        help='Regions to analyze (e.g. V3 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save alignment to file')
    parser.add_argument('--plot', action='store_true',
                        help='Plot phylogenetic tree. Requires --save and --tree')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    use_save = args.save
    use_plot = args.plot

    patients = load_patients()
    if pnames != ['all']:
        patients = patients.iloc[patients.index.isin(pnames)]

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        if regions is None:
            refseq_gw = patient.get_reference('genomewide', 'gb')
            regionspat = map(attrgetter('id'), refseq_gw.features) + ['genomewide']
        else:
            regionspat = regions

        for region in regionspat:
            if VERBOSE >= 1:
                print pname, region
                if VERBOSE == 1:
                    print ''

            if VERBOSE >= 2:
                print 'Get alignment'
            ali = patient.get_consensi_alignment(region)

            if VERBOSE >= 2:
                print 'Build tree'
                sys.stdout.flush()
            tree = build_tree_fasttree(ali, rootname=ali[0].id,
                                       VERBOSE=VERBOSE)
            tree.ladderize()

            if use_save:
                if VERBOSE >= 2:
                    print 'Save tree'
                from Bio import Phylo
                fn_out = patient.get_consensi_tree_filename(region)
                mkdirs(os.path.dirname(fn_out))
                Phylo.write(tree, fn_out, 'newick')

            if use_plot:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(15, 12))
                Phylo.draw(tree, do_show=False, axes=ax)
                ax.set_title(pname+', '+region)

                x_max = max(tree.depths().itervalues())
                ax.set_xlim(0.995, 0.995 + (x_max - 0.995) * 1.4)
                ax.grid(True)
                
                plt.ion()
                plt.show()



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
from hivwholeseq.patients.patients import load_patients, Patient



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get phylogenetic tree of consensi',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. V3 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    use_plot = args.plot

    patients = load_patients()
    if pnames:
        patients = patients.loc[pnames]

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        for region in regions:
            if VERBOSE >= 1:
                print pname, region

            ali = patient.get_consensi_alignment(region)
            tree = patient.get_consensi_tree(region, format='json')

            if use_plot:
                fig, ax = plt.subplots(figsize=(15, 12))
                Phylo.draw(tree, do_show=False, axes=ax)
                ax.set_title(pname+', '+region)

                x_max = max(tree.depths().itervalues())
                ax.set_xlim(0.998, 0.998 + (x_max - 0.998) * 1.5)
                ax.grid(True)
            
    if use_plot:
        plt.ion()
        plt.show()






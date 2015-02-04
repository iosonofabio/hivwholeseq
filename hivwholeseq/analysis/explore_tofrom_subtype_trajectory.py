# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/14
content:    Study whether HIV in a patient tends to explore the mutational space
            in a way that comes closer to a subtype average (entropy profile).
'''
# Modules
import os, sys
import argparse
from collections import defaultdict
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alphal, alpha
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.argparse import RoiAction
from hivwholeseq.reference import load_custom_reference, load_custom_alignment
from hivwholeseq.cross_sectional.get_subtype_allele_frequencies import (
    get_subtype_reference_alignment_allele_frequencies)


# Globals
refname = 'HXB2'



# Functions
def trim_comap(comap, ref1, ref2, VERBOSE=0):
    '''Trim comap to relevant region'''
    from seqanpy import align_overlap

    (score, ali1, ali2) = align_overlap(ref1, ref2, score_gapopen=-20)
    start = len(ali2) - len(ali2.lstrip('-'))
    end = len(ali2.rstrip('-'))

    if VERBOSE >= 3:
        from hivwholeseq.utils.sequence import pretty_print_pairwise_ali
        pretty_print_pairwise_ali((ali1, ali2), width=100, name1='pat', name2='ali')

    ind = (comap[:, 1] >= start) & (comap[:, 1] < end)
    comap = comap[ind]
    comap[:, 0] -= start

    return comap



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Mutations away from/towards subtype',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Genomic regions (e.g. V3 IN)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot results')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    use_plot = args.plot

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    

    for region in regions:
        if VERBOSE >= 1:
            print region

        # Load HXB2
        refseq = load_custom_reference(refname, 'gb')
        refseq.name = refname
        refm = np.array(refseq)

        # Load subtype alignment
        afs = get_subtype_reference_alignment_allele_frequencies(region)


        for pname, patient in patients.iterrows():
            patient = Patient(patient)

            if VERBOSE >= 1:
                print region, pname

            comap = patient.get_map_coordinates_reference(region,
                                                          refname=(refname, region))

            aft, ind = patient.get_allele_frequency_trajectories(region)
            
            # Sometimes the alignment is trimmed
            comap = trim_comap(comap, alpha[aft[0].argmax(axis=0)], alpha[afs.argmax(axis=0)],
                               VERBOSE=VERBOSE)

            afsp = afs[:, comap[:, 0]]
            aft = aft[:, :, comap[:, 1]]

            # First, check the fraction of sites for which initial consensus = subtype consensus
            consi = aft[0].argmax(axis=0)
            conss = afsp.argmax(axis=0)
            print 'Fraction of sites at which patient initial consensus agrees with subtype:',
            print '{:2.0%}'.format((consi == conss).mean())

            # Check at the end
            consf = aft[-1].argmax(axis=0)
            print 'Fraction of sites at which patient final consensus agrees with subtype:',
            print '{:2.0%}'.format((consf == conss).mean())

            # FIXME
            sys.exit()





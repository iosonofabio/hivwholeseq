# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/01/15
content:    Do mutations happen more towards the consensus or away from it?
'''
# Modules
import os
import argparse
from itertools import izip
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.sequence import translate_with_gaps
from hivwholeseq.multipatient.explore_codon_usage_patsubtype import get_degenerate_dict
import hivwholeseq.utils.plot
from hivwholeseq.multipatient.explore_entropy_patsubtype import (
    get_subtype_reference_alignment, get_ali_entropy)
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment_entropy)
from hivwholeseq.patients.one_site_statistics import get_codons_n_polymorphic
from hivwholeseq.utils.one_site_statistics import get_allele_frequencies_alignment



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Explore conservation levels across patients and subtype',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot results')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot

    data = defaultdict(list)

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for region in regions:
        if VERBOSE >= 1:
            print region

        ali = get_subtype_reference_alignment(region)
        alim = np.array(ali, 'S1')
        afs = get_allele_frequencies_alignment(alim, VERBOSE=VERBOSE)

        for ipat, (pname, patient) in enumerate(patients.iterrows()):

            patient = Patient(patient)
            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=1000,
                                                                 depth_min=300,
                                                                 VERBOSE=VERBOSE)
            if len(ind) == 0:
                continue

            # Get the coordinate map to HXB2 (and ref alignment)
            coomap = patient.get_map_coordinates_reference(region)
            coomap[:, 0] -= coomap[coomap[:, 1] == 0, 0][0]

            # NOTE: Pavel has sometimes cut the second half of the alignment (e.g. RT)
            if len(alim.T) < (coomap[:, 0]).max() + 1:
                coomap = coomap[coomap[:, 0] < len(Ssub)]
                aft = aft[:, :, :coomap[:, 1].max() + 1]

            coomap = dict(coomap[:, ::-1])

            times = patient.times[ind]

            icons = patient.get_initial_consensus_noinsertions(aft, VERBOSE=VERBOSE)
            cons = alpha[icons]

            for pospat, posali in coomap.iteritems():
                for ia, a in enumerate(alpha[:4]):
                    af = aft[pospat, ia]
                    data['muts'].append(())


                    # FIXME

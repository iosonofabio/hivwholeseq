# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/01/15
content:    Look at how long are polymorphisms staying around.
'''
# Modules
import os
import sys
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
from hivwholeseq.one_site_statistics import get_allele_frequencies_alignment



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

        #ali = get_subtype_reference_alignment(region)
        #alim = np.array(ali, 'S1')
        #afs = get_allele_frequencies_alignment(alim, VERBOSE=VERBOSE)
        S = get_subtype_reference_alignment_entropy(region)

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
            if len(S) < (coomap[:, 0]).max() + 1:
                coomap = coomap[coomap[:, 0] < len(Ssub)]
                aft = aft[:, :, :coomap[:, 1].max() + 1]

            coomap = dict(coomap[:, ::-1])

            times = patient.times[ind]

            icons = patient.get_initial_consensus_noinsertions(aft, VERBOSE=VERBOSE)
            cons = alpha[icons]

            # Take only polymorphic sites and derived alleles
            pos_poly, ia_poly = (((aft > 0.01) & (aft < 0.99)).any(axis=0) & 
                                 (aft[0] < 0.1)).T.nonzero()
            for pospat, ia in izip(pos_poly, ia_poly):
                if pospat not in coomap:
                    continue

                posali = coomap[pospat]

                af = aft[:, ia, pospat]

                it0 = (af > 0.01).nonzero()[0][0]

                dit = ((af[it0:] < 0.01) | (af[it0:] > 0.99)).nonzero()[0]

                # Stays polymorphic forever
                if len(dit) == 0:
                    data['t'].append((pname, region,
                                      posali,
                                      cons[pospat], alpha[ia],
                                      af.max(),
                                      np.inf,
                                      times[it0], times[-1],
                                      len(times) - it0,
                                      'poly',
                                      S[posali]))

                else:
                    it1 = dit[0] + it0
                    if af[it1] < 0.01:
                        fate = 'loss'
                    else:
                        fate = 'fix'
                    data['t'].append((pname, region,
                                      posali,
                                      cons[pospat], alpha[ia],
                                      af.max(),
                                      times[it1] - times[it0],
                                      times[it0], times[it1],
                                      it1 - it0,
                                      fate,
                                      S[posali]))


    data['t'] = pd.DataFrame(data=data['t'],
                             columns=['Patient', 'Region',
                                      'Pos',
                                      'Anc', 'Mut',
                                      'afmax',
                                      'dt',
                                      't0', 't1',
                                      'dind',
                                      'fate',
                                      'Ssub'])

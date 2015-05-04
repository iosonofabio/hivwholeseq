# vim: fdm=indent
'''
author:     Fabio Zanini
date:       04/05/15
content:    Describe substitutions away from and to subtype consensus and
            whether or not they are within CTL epitopes.
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
from Bio.Seq import translate

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.patients import load_patients, iterpatient
from hivwholeseq.utils.sequence import translate_with_gaps
import hivwholeseq.utils.plot
from hivwholeseq.analysis.explore_entropy_patsubtype import (
    get_subtype_reference_alignment, get_ali_entropy)
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment_entropy)
from hivwholeseq.cross_sectional.get_subtype_consensus import (
    get_subtype_reference_alignment_consensus)
from hivwholeseq.patients.one_site_statistics import get_codons_n_polymorphic
from hivwholeseq.analysis.mutation_rate.explore_divergence_synonymous import translate_masked
from hivwholeseq.utils.argparse import PatientsAction


# Globals
regions = ['p17', 'p24',
           'PR', 'RT', 'IN',
           'vif', 'vpu',
           'gp41',
           'nef']



# Functions
def collect_data(pnames, regions, VERBOSE=0):
    '''Collect data for sweep call'''
    data = []
    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    from hivwholeseq.reference import load_custom_reference
    from hivwholeseq.utils.sequence import find_annotation
    ref = load_custom_reference('HXB2', 'gb')

    for region in regions:
        if VERBOSE >= 1:
            print region

        fea = find_annotation(ref, region)
        region_start = fea.location.nofuzzy_start

        if VERBOSE >= 2:
            print 'Get subtype consensus'
        conssub = get_subtype_reference_alignment_consensus(region, VERBOSE=VERBOSE)

        for ipat, (pname, patient) in enumerate(iterpatient(patients)):
            pcode = patient.code
            if VERBOSE >= 2:
                print pcode, region

            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=100,
                                                                 depth_min=10,
                                                                 VERBOSE=VERBOSE)
            if len(ind) == 0:
                if VERBOSE >= 2:
                    print 'No time points: skip'
                continue

            times = patient.times[ind]

            if VERBOSE >= 2:
                print 'Get CTL epitopes'
            ctl_table = patient.get_ctl_epitopes(regions=[region])

            if VERBOSE >= 2:
                print 'Get coordinate map'
            coomap = patient.get_map_coordinates_reference(region, refname=('HXB2', region))

            icons = patient.get_initial_consensus_noinsertions(aft, VERBOSE=VERBOSE,
                                                               return_ind=True)
            consm = alpha[icons]

            # Get the map as a dictionary from patient to subtype
            coomapd = {'pat_to_subtype': dict(coomap[:, ::-1]),
                       'subtype_to_pat': dict(coomap)}

            # Condition on fixation
            ind_sweep = zip(*((aft[0] < 0.05) & (aft[-2:] > 0.95).any(axis=0)).T.nonzero())

            for posdna, inuc in ind_sweep:
                # Get the position in reference coordinates
                if posdna not in coomapd['pat_to_subtype']:
                    continue
                pos_sub = coomapd['pat_to_subtype'][posdna]
                if (pos_sub >= len(conssub)):
                    continue
                conspos_sub = conssub[pos_sub]
                pos_ref_gw = pos_sub + region_start

                # Get allele frequency trajectory
                aftpos = aft[:, :, posdna].T

                # Get only non-masked time points
                indpost = -aftpos[0].mask
                if indpost.sum() == 0:
                    continue
                timespos = times[indpost]
                aftpos = aftpos[:, indpost]

                anc = consm[posdna]
                ianc = icons[posdna]
                nuc = alpha[inuc]
                mut = anc+'->'+nuc

                # Ignore indels
                if (inuc >= 4) or (ianc >= 4):
                    continue

                # Skip if the site is already polymorphic at the start
                if aftpos[ianc, 0] < 0.95:
                    continue

                # Define transition/transversion
                if frozenset(nuc+anc) in (frozenset('CT'), frozenset('AG')):
                    trclass = 'ts'
                else:
                    trclass = 'tv'

                # Away/to subtype consensus
                if (anc == conspos_sub):
                    away_conssub = 'away'
                elif (nuc == conspos_sub):
                    away_conssub = 'to'
                else:
                    away_conssub = 'neither'
                    # NOTE: we could keep those
                    continue

                # Find whether is within an epitope
                if len(ctl_table):
                    is_epitope = ((pos_ref_gw >= np.array(ctl_table['start_HXB2'])) &
                                  (pos_ref_gw < np.array(ctl_table['end_HXB2']))).any()
                else:
                    is_epitope = False

                datum = {'pcode': patient.code,
                         'region': region,
                         'pos_patient': posdna,
                         'pos_ref': pos_ref_gw,
                         'mut': mut,
                         'trclass': trclass,
                         'epitope': is_epitope,
                         'awayto': away_conssub,
                        }

                data.append(datum)

    data = pd.DataFrame(data)
    return data



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Study accumulation of minor alleles for different kinds of mutations',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', action=PatientsAction,
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+', default=regions,
                        help='Region to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot results')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot


    data = collect_data(pnames, regions, VERBOSE=VERBOSE)


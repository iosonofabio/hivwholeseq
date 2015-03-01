# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/01/15
content:    Get the fixation probability of alleles in different categories.
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
from hivwholeseq.patients.patients import load_patients, Patient
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


# Globals



# Functions
def calculate_probabilities(data, binsc_af, VERBOSE=0):
    '''Calculate fixation probabilities'''
    datap = (data
             .loc[:, ['afbin', 'fix']]
             .groupby('afbin')
             .mean())
    datap['# alleles'] = (data
                          .loc[:, ['afbin', 'fix']]
                          .groupby('afbin')
                          .size())
    datap['af'] = binsc_af
    return datap



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Study accumulation of minor alleles for different kinds of mutations',
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

    regdata = []
    data = []

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for region in regions:
        if VERBOSE >= 1:
            print region

        if VERBOSE >= 2:
            print 'Get subtype consensus (for checks only)'
        conssub = get_subtype_reference_alignment_consensus(region, VERBOSE=VERBOSE)

        if VERBOSE >= 2:
            print 'Get subtype entropy'
        Ssub = get_subtype_reference_alignment_entropy(region, VERBOSE=VERBOSE)

        regdata.append({'name': region, 'consensus': conssub, 'S': Ssub, 'L': len(conssub)})

        for ipat, (pname, patient) in enumerate(patients.iterrows()):
            pcode = patient.code
            if VERBOSE >= 2:
                print pname, pcode

            patient = Patient(patient)
            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=1000,
                                                                 depth_min=300,
                                                                 VERBOSE=VERBOSE)
            if len(ind) == 0:
                if VERBOSE >= 2:
                    print 'No time points: skip'
                continue

            times = patient.times[ind]

            if VERBOSE >= 2:
                print 'Get coordinate map'
            coomap = patient.get_map_coordinates_reference(region, refname=('HXB2', region))

            icons = patient.get_initial_consensus_noinsertions(aft, VERBOSE=VERBOSE,
                                                               return_ind=True)
            consm = alpha[icons]
            protm = translate_masked(consm)
            
            # Premature stops in the initial consensus???
            if '*' in protm:
                # Trim the stop codon if still there (some proteins are also end of translation)
                if protm[-1] == '*':
                    if VERBOSE >= 2:
                        print 'Ends with a stop, trim it'
                    icons = icons[:-3]
                    consm = consm[:-3]
                    protm = protm[:-1]
                    aft = aft[:, :, :-3]
                    coomap = coomap[coomap[:, 1] < len(consm)]

                else:
                    continue

            # Get the map as a dictionary from patient to subtype
            coomapd = {'pat_to_subtype': dict(coomap[:, ::-1]),
                       'subtype_to_pat': dict(coomap)}

            # FIXME: deal better with depth (this should be already there?)
            aft[aft < 2e-3] = 0

            for posdna, pos_sub in coomapd['pat_to_subtype'].iteritems():
                # Get the entropy
                if pos_sub >= len(Ssub):
                    continue
                Ssubpos = Ssub[pos_sub]

                # Get allele frequency trajectory
                aftpos = aft[:, :, posdna].T

                # Get only non-masked time points
                indpost = -aftpos[0].mask
                if indpost.sum() < 2:
                    continue
                timespos = times[indpost]
                aftpos = aftpos[:, indpost]

                anc = consm[posdna]
                ianc = icons[posdna]

                for inuc, nuc in enumerate(alpha[:4]):
                    mut = anc+'->'+nuc

                    # Ignore ancestral allele
                    if inuc == ianc:
                        continue

                    # Transition/transversion
                    if frozenset(nuc+anc) in (frozenset('CT'), frozenset('AG')):
                        trclass = 'ts'
                    else:
                        trclass = 'tv'

                    # Syn/nonsyn
                    codanc = consm[posdna - posdna % 3: posdna - posdna % 3 + 3]
                    codder = codanc.copy()
                    codder[posdna % 3] = nuc
                    if ('-' in codanc) or ('-' in codder):
                        continue

                    if translate(''.join(codanc)) == translate(''.join(codder)):
                        mutclass = 'syn'
                    else:
                        mutclass = 'nonsyn'

                    # Analyze the trajectory
                    for it, (af, time) in enumerate(izip(aftpos[inuc, :-1], timespos[:-1])):
                        af_future = aftpos[inuc, it + 1:]
                        has_fixed = (af_future > 0.95).any()
                        has_lost = (af_future < 0.05).any()
                        if not (has_fixed or has_lost):
                            continue

                        if has_fixed and has_lost:
                            has_fixed = False
                            has_lost = False
                            it_fixed = (af_future > 0.95).nonzero()[0][0]
                            it_lost = (af_future < 0.05).nonzero()[0][0]
                            if it_fixed < it_lost:
                                has_fixed = True
                            else:
                                has_lost = True

                        if has_fixed:
                            fixlost = 'fixed'
                        else:
                            fixlost = 'lost'

                        data.append((region, pcode,
                                     posdna, pos_sub,
                                     anc, nuc, mut,
                                     mutclass, trclass,
                                     fixlost,
                                     Ssubpos,
                                     time, af))

    data = pd.DataFrame(data=data,
                        columns=['region', 'pcode',
                                 'posdna', 'possub',
                                 'anc', 'der', 'mut',
                                 'class', 'tr',
                                 'fixclass',
                                 'Ssub',
                                 'time', 'af'])

    regdata = pd.DataFrame(regdata).set_index('name', drop=False)

    if VERBOSE >= 1:
        print 'Analyze collected alleles'
    # Bin allele frequency
    #logistic_fun = lambda x: 1.0 / (1.0 + 10**(-x))
    #bins_af_logit = np.linspace(-2, 1.5, 10)
    #binsc_af_logit = 0.5 * (bins_af_logit[1:] + bins_af_logit[:-1])
    #bins_af = logistic_fun(bins_af_logit)
    #binsc_af = logistic_fun(binsc_af_logit)

    bins_af = np.linspace(1e-2, 1-1e-2, 8)
    binsc_af = 0.5 * (bins_af[1:] + bins_af[:-1])
    binsw_af = bins_af[1:] - bins_af[:-1]
    data['afbin'] = -1
    for b in bins_af:
        data.loc[data.loc[:, 'af'] >= b, 'afbin'] += 1

    # Set boolean fix/lost
    data['fix'] = data.loc[:, 'fixclass'] == 'fixed'

    # Exclude alleles that are outside of boundaries
    data = data.loc[(data.loc[:, 'afbin'] >= 0) & (data.loc[:, 'afbin'] < len(binsc_af))]


    if plot:

        fix, ax = plt.subplots()

        # Plot all alleles together
        datap = calculate_probabilities(data, binsc_af, VERBOSE=VERBOSE)
        x = np.array(datap.loc[:, 'af'])
        y = np.array(datap.loc[:, 'fix'])
        n = np.array(datap.loc[:, '# alleles'])
        ax.plot(x, y, lw=2, label='all')
        for i in xrange(len(x)):
            ax.text(x[i], y[i], str(n[i]),
                    horizontalalignment='center',
                    verticalalignment='center')


        # Plot syn/nonsyn
        datap = calculate_probabilities(data[data['class'] == 'syn'], binsc_af, VERBOSE=VERBOSE)
        x = np.array(datap.loc[:, 'af'])
        y = np.array(datap.loc[:, 'fix'])
        n = np.array(datap.loc[:, '# alleles'])
        ax.plot(x, y, lw=2, label='syn')
        for i in xrange(len(x)):
            ax.text(x[i], y[i], str(n[i]),
                    horizontalalignment='center',
                    verticalalignment='center')

        datap = calculate_probabilities(data[data['class'] == 'nonsyn'], binsc_af, VERBOSE=VERBOSE)
        x = np.array(datap.loc[:, 'af'])
        y = np.array(datap.loc[:, 'fix'])
        n = np.array(datap.loc[:, '# alleles'])
        ax.plot(x, y, lw=2, label='nonsyn')
        for i in xrange(len(x)):
            ax.text(x[i], y[i], str(n[i]),
                    horizontalalignment='center',
                    verticalalignment='center')

        # Neutral line
        xfit = np.linspace(1e-3, 1-1e-3, 1000)
        yfit = xfit
        ax.plot(xfit, yfit, lw=1.5, color='k', alpha=0.5)

        ax.grid(True)
        ax.set_xlim(1e-3, 1-1e-3)
        #ax.set_xscale('logit')
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel('Initial frequency')
        ax.set_ylabel('Fixation probability')
        ax.legend(loc='lower right')

        plt.ion()
        plt.show()


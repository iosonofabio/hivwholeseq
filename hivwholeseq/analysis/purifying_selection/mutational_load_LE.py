# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/01/15
content:    Quantify purifying selection on different subtype entropy classes.
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
import hivwholeseq.utils.plot
from hivwholeseq.analysis.explore_entropy_patsubtype import (
    get_subtype_reference_alignment, get_ali_entropy)
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment_entropy)
from hivwholeseq.cross_sectional.get_subtype_consensus import (
    get_subtype_reference_alignment_consensus)
from hivwholeseq.patients.one_site_statistics import get_codons_n_polymorphic
from hivwholeseq.analysis.mutation_rate.explore_divergence_synonymous import translate_masked
from hivwholeseq.analysis.purifying_selection.filenames import get_fitness_cost_entropy_filename


# Globals



# Functions



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Study load by deleterious mutations on a hypothetical pop in linkage eq',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot results')
    parser.add_argument('-N', type=int, default=100000,
                        help='Number of points in the distribution')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot
    N = args.N

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    data = []
    for region in regions:
        if VERBOSE >= 1:
            print region

        if VERBOSE >= 2:
            print 'Get subtype consensus'
        conssub = get_subtype_reference_alignment_consensus(region, VERBOSE=VERBOSE)
        conssubm = np.array(conssub)

        if VERBOSE >= 2:
            print 'Get subtype entropy'
        Ssub = get_subtype_reference_alignment_entropy(region, VERBOSE=VERBOSE)

        if VERBOSE >= 2:
            print 'Get fitness costs'
        fn = get_fitness_cost_entropy_filename(region)
        fitness_cost = pd.read_pickle(fn)
        fitness_cost.set_index('iSbin', drop=False, inplace=True)

        if VERBOSE >= 2:
            print 'Get entropy bins'
        iSbin = -np.ones(len(Ssub), int)
        for _, datum in fitness_cost.iterrows():
            iSbin[Ssub > datum['Smin']] += 1
        # Fix extremes
        iSbin = iSbin.clip(fitness_cost['iSbin'].min(), fitness_cost['iSbin'].max())

        if VERBOSE >= 2:
            print 'Get fitness array'
        fitness_cost_arr = np.array(fitness_cost.loc[iSbin, 's'])

        for ipat, (pname, patient) in enumerate(patients.iterrows()):
            pcode = patient.code
            if VERBOSE >= 2:
                print pname, pcode

            patient = Patient(patient)

            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=1000,
                                                                 depth_min=100)
            if len(ind) == 0:
                if VERBOSE >= 2:
                    print 'No time points: skip'
                continue

            times = patient.times[ind]

            if VERBOSE >= 2:
                print 'Get coordinate map'
            coomap = patient.get_map_coordinates_reference(region, refname=('HXB2', region))

            iconsm = patient.get_initial_consensus_noinsertions(aft, return_ind=True)
            consm = alpha[iconsm]
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

            # Get the map as a dictionary from patient to subtype and vice versa
            coomapd = {'pat_to_subtype': pd.Series(coomap[:, 0], index=coomap[:, 1]),
                       'subtype_to_pat': pd.Series(coomap[:, 1], index=coomap[:, 0])}

            
            # Restrict allele frequencies to positions in the map
            pos_tmp = coomap[:, 0]
            cons_sub_tmp = conssubm[pos_tmp]
            icons_sub_tmp = np.array(map(alphal.index, cons_sub_tmp), int)
            fitness_cost_tmp = fitness_cost_arr[pos_tmp]

            aft_tmp = aft[:, :, coomap[:, 1]]
            aft_cons_sub = aft_tmp[:, icons_sub_tmp, np.arange(aft_tmp.shape[-1])]


            if VERBOSE >= 2:
                print 'Get hypothetical LE populations and their fitnesses'
            ran = np.random.rand(N, aft_cons_sub.shape[0], aft_cons_sub.shape[-1])
            fitness_distr = np.dot((ran > aft_cons_sub), -fitness_cost_tmp).T

            data.append({'pcode': pcode, 'region': region,
                         'times': times, 'fitness': fitness_distr})

    data = pd.DataFrame(data)

    if plot: 
        if VERBOSE >= 1:
            print 'Plot'
        for _, datum in data.iterrows():
            region = datum['region']
            pcode = datum['pcode']
            times = datum['times']
            fitness = datum['fitness']

            # Plot cumulative
            fig, ax = plt.subplots()
            for it, time in enumerate(times):
                # Get the absolute value for plotting with log scale
                fit_t = -fitness[it]
                fit_t.sort()

                x = fit_t
                y = 1.0 - np.linspace(0, 1, len(x))
                color = cm.jet(1.0 * it / len(times))
                ax.plot(x, y, lw=2, color=color, label=str(int(time))+' days')

            ax.set_xlabel('Fitness cost')
            ax.set_ylabel('Fraction of haplotypes with fitness > x')
            ax.set_title(', '.join([region, pcode]))
            ax.grid(True)
            ax.set_ylim(-0.05, 1.05)
            ax.set_xlim(0, 0.2)

            plt.tight_layout()

            # Plot density
            fig, ax = plt.subplots()
            bins = np.linspace(0.96 * -fitness.max(), 1.04 * -fitness.min(), 100)
            binsc = 0.5 * (bins[1:] + bins[:-1])
            for it, time in enumerate(times):
                # Get the absolute value for plotting with log scale
                fit_t = -fitness[it]
                fit_t = np.histogram(fit_t, bins=bins, density=True)[0]

                x = binsc
                y = fit_t
                color = cm.jet(1.0 * it / len(times))
                ax.plot(x, y, lw=2,
                        color=color, alpha=0.7,
                        label=str(int(time))+' days')

            ax.set_xlabel('Fitness cost')
            ax.set_ylabel('Density')
            ax.set_title(', '.join([region, pcode]))
            ax.grid(True)
            ax.set_xlim(0, 0.2)
            ax.set_ylim(ymin=-0.04*ax.get_ylim()[1])
            ax.set_yticklabels([])
            ax.legend(loc='upper right', title='Time:', fontsize=14)

            plt.tight_layout()


        plt.ion()
        plt.show()


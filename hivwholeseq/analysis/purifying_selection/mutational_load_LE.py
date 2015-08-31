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

from hivwholeseq.utils.miseq import alpha, alphal
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
regions = ['p17', 'p24']



# Functions
def collect_data_mutational_load_LE(pnames, regions, VERBOSE=0, N=10000):
    '''Collect data on mutational load under linkage equilibrium'''
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
            print 'Get fitness costs (genomewide average)'
        fn = get_fitness_cost_entropy_filename('all_joint')
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
                    iconsm = iconsm[:-3]
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
            cons_tmp = consm[coomap[:, 1]]

            if VERBOSE >= 2:
                print 'Get hypothetical LE populations and their fitnesses'
            ran = np.random.rand(N, aft_cons_sub.shape[0], aft_cons_sub.shape[-1])
            pops = ran > aft_cons_sub
            fitness_distr = np.dot(pops, -fitness_cost_tmp).T

            if VERBOSE >= 2:
                print 'Subtract fitness of initial consensus'
            fitness_cons0 = np.dot(cons_tmp != cons_sub_tmp, -fitness_cost_tmp)

            fitness_distr -= fitness_cons0

            for time, fitn in izip(times, fitness_distr):
                data.append({'pcode': pcode, 'region': region,
                             'time': time, 'fitness': fitn,
                             'nsites': len(pos_tmp)})

    data = pd.DataFrame(data)
    data.set_index(['region', 'time'], drop=False, inplace=True)

    return data


def get_shared_times(data, VERBOSE=0):
    '''For genomewide estimates, take only shared time points'''
    # FIXME: this function should group by patient
    pcode = data['pcode'].iloc[0]

    regions = np.unique(data['region']).tolist()
    
    tmp = defaultdict(set)
    for _, datum in data.iterrows():
        tmp[datum['time']].add(datum['region'])
    timesgw = sorted(t for t, val in tmp.iteritems() if len(val) == len(regions))

    datagw = []
    for time in timesgw:
        datatmp = data.loc[data.loc[:, 'time'] == time]
        datumgw = {'fitness': datatmp.loc[:, 'fitness'].sum(),
                   'region': 'genomewide',
                   'pcode': pcode,
                   'time': time,
                   'nsites': datatmp.loc[:, 'nsites'].sum()}
        datagw.append(datumgw)
    datagw = pd.DataFrame(datagw)

    return datagw


def plot_mutational_load_LE(data, scale_gw=False, type='density', VERBOSE=0):
    '''Plot the distribution of mutation loads'''
    data = data.copy()

    region = data.iloc[0]['region']
    pcode = data.iloc[0]['pcode']
    nsites = data.iloc[0]['nsites']

    if scale_gw:
        L = 9000
        data['fitness'] *= L / nsites

    if type in ['both', 'cumulative']:
        if VERBOSE >= 1:
            print 'Plot cumulative distribution'
        fig, ax = plt.subplots()
        for it, (_, datum) in enumerate(data.iterrows()):
            time = np.array(datum['time'])
            fit_t = np.sort(-np.array(datum['fitness']))

            x = fit_t
            y = 1.0 - np.linspace(0, 1, len(x))
            color = cm.jet(1.0 * it / len(data['time']))
            ax.plot(x, y, lw=2, color=color, label=str(int(time))+' days')

        ax.set_xlabel('Fitness cost [compared to founder]')
        ax.set_ylabel('Fraction of haplotypes with fitness > x')
        ax.set_title(', '.join([region, pcode]))
        ax.grid(True)
        ax.set_ylim(-0.05, 1.05)

        plt.tight_layout()

    if type in ['both', 'density']:
        if VERBOSE >= 1:
            print 'Plot density'
        fig, ax = plt.subplots()
        fitbounds = np.array([min(map(np.min, -data['fitness'])),
                              max(map(np.max, -data['fitness']))])

        bins = np.linspace(0.96 * fitbounds[0], 1.04 * fitbounds[1], 30)
        binsc = 0.5 * (bins[1:] + bins[:-1])
        for it, (_, datum) in enumerate(data.iterrows()):
            time = np.array(datum['time'])
            fit_t = -np.array(datum['fitness'])

            fit_t = np.histogram(fit_t, bins=bins, density=True)[0]
            x = binsc
            y = fit_t
            color = cm.jet(1.0 * it / len(data['time']))
            ax.plot(x, y, lw=2,
                    color=color, alpha=0.7,
                    label=str(int(time))+' days')

        ax.set_xlabel('Fitness cost [compared to founder]')
        ax.set_ylabel('Density')
        ax.set_title(', '.join([region, pcode]))
        ax.grid(True)
        ax.set_ylim(-0.04*ax.get_ylim()[1], 1.04*ax.get_ylim()[1])
        ax.set_yticklabels([])
        ax.legend(loc='upper right', title='Time:', fontsize=14)

        plt.tight_layout()



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
    parser.add_argument('--plot', action='store_true',
                        help='Plot results')
    parser.add_argument('-N', type=int, default=10000,
                        help='Number of points in the distribution')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot
    N = args.N

    data = collect_data_mutational_load_LE(pnames, regions, VERBOSE=VERBOSE, N=N)

    datagw = get_shared_times(data, VERBOSE=VERBOSE)

    if plot: 
        if VERBOSE >= 1:
            print 'Plot'
        #for region, datum in data.groupby('region'):
        #    datum['region'] = region
        #    plot_mutational_load_LE(datum, scale_gw=False, VERBOSE=VERBOSE)

        plot_mutational_load_LE(datagw, scale_gw=True, VERBOSE=VERBOSE)

        plt.ion()
        plt.show()


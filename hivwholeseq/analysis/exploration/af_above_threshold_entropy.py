# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/07/14
content:    Explore sites that are very conserved within our patient set, across
            the whole infection, compared to their behaviour within the subtype.
'''
# Modules
import os
import sys
import argparse
from itertools import izip, combinations, chain
from collections import defaultdict
from operator import itemgetter
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO
from scipy.stats import pearsonr, spearmanr

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.filenames import root_patient_folder
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.cross_sectional.get_subtype_consensus import (
    get_subtype_reference_alignment_consensus)
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment_entropy)
from hivwholeseq.cross_sectional.get_subtype_allele_frequencies import (
    get_subtype_reference_alignment_allele_frequencies)



# Globals
regions = ['p17', 'p24', 'PR', 'RT', 'p15', 'IN', 'gp41']
pnames = ['20097', '15363', '15823', '9669', '20529', '15241', '15376', '15319']




# Functions
def collect_data(pnames, regions, VERBOSE=0):
    '''Collect data'''

    data = []
    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for region in regions:
        if VERBOSE >= 1:
            print region

        if VERBOSE >= 2:
            print 'Get subtype consensus'
        conssub = get_subtype_reference_alignment_consensus(region, VERBOSE=VERBOSE)

        if VERBOSE >= 2:
            print 'Get subtype allele frequencies'
        af_sub = get_subtype_reference_alignment_allele_frequencies(region)

        if VERBOSE >= 2:
            print 'Get subtype entropy'
        Ssub = get_subtype_reference_alignment_entropy(region, VERBOSE=VERBOSE)

        for ipat, (pname, patient) in enumerate(patients.iterrows()):
            pcode = patient.code
            if VERBOSE >= 2:
                print pname, pcode

            patient = Patient(patient)
            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=1000,
                                                                 depth_min=30,
                                                                 VERBOSE=VERBOSE)
            if len(ind) == 0:
                if VERBOSE >= 2:
                    print 'No time points found: skip.'
                continue

            times = patient.times[ind]

            if VERBOSE >= 2:
                print 'Get coordinate map'
            coomap = patient.get_map_coordinates_reference(region, refname=('HXB2', region))

            if VERBOSE >= 2:
                print 'Get initial consensus'
            icons = patient.get_initial_consensus_noinsertions(aft, VERBOSE=VERBOSE,
                                                               return_ind=True)
            consm = alpha[icons]

            # Get the map as a dictionary from patient to subtype
            coomapd = dict(coomap[:, ::-1])

            if VERBOSE >= 2:
                print 'Get alleles'
            for posdna in xrange(aft.shape[-1]):
                if VERBOSE >= 3:
                    print posdna

                # Look for this position in the subtype alignment
                if posdna not in coomapd:
                    continue
                pos_sub = coomapd[posdna]
                if (pos_sub >= af_sub.shape[1]):
                    continue

                # Ancestral allele
                ianc = icons[posdna]
                anc = alpha[ianc]

                # Discard if the initial time point is already polymorphic
                aft_anc0 = aft[0, ianc, posdna]
                if aft_anc0 < 0.9:
                    continue

                # KEY: Take sum of derived alleles
                aft_der = 1.0 - aft[:, ianc, posdna]

                # Get only non-masked time points
                indpost = -aft_der.mask
                if indpost.sum() == 0:
                    continue
                timespos = times[indpost]

                # Set subtype wide properties
                Spos_sub = Ssub[pos_sub]
                conspos_sub = conssub[pos_sub]
                if (anc == conspos_sub):
                    away_conssub = 'away'
                else:
                    away_conssub = 'to'

                for it, time in enumerate(timespos):
                    af = aft_der[it]
                    data.append((region, pcode,
                                 anc,
                                 Spos_sub, conspos_sub, away_conssub,
                                 time, af))

    data = pd.DataFrame(data=data, columns=('region', 'pcode',
                                            'anc',
                                            'Ssub', 'conssub', 'awayto',
                                            'time', 'af'))


    return data


def calculate_fraction_above_threshold(data):
    '''Calculate the fraction of allele above thresholds'''
    from hivwholeseq.utils.pandas import add_binned_column
    bins_S, binsc_S = add_binned_column(data, 'Sbin', 'Ssub', bins=5, clip=True)
    bins_t, binsc_t = add_binned_column(data, 'tbin', 'time',
                                        bins=300 * np.arange(8),
                                        clip=True)

    thresholds = 0.01 * np.arange(10)
    datap = []
    for (Sbin, tbin, pcode), datum in data.groupby(['Sbin', 'tbin', 'pcode']):
        datump = {'Ssub': binsc_S[Sbin],
                  'time': binsc_t[tbin],
                  'pcode': pcode,
                 }
        for thr in thresholds:
            datump['frac > '+'{:.1G}'.format(thr)] = (datum['af'] > thr).mean()
        datap.append(datump)

    datap = pd.DataFrame(datap)

    return datap


def plot_af_above_threshold(datap, threshold=0.01, singles=False):
    '''Plot the fraction of alleles above a certain threshold'''
    Ssubs = datap.Ssub.unique().tolist()

    colormap = cm.jet
    fs = 16

    fig, ax = plt.subplots()

    if singles:
        for (Ssub, pcode), datump in datap.groupby(['Ssub', 'pcode']):
            x = np.array(datump['time'])
            y = np.array(datump['frac > '+'{:.1G}'.format(threshold)])

            ax.plot(x, y,
                    lw=1,
                    ls='-', marker='x', ms=10,
                    color=colormap(1.0 * Ssubs.index(Ssub) / len(Ssubs)),
                    alpha=0.7,
                   )

    # Average over patients
    for Ssub, datump in datap.groupby('Ssub'):
        datumpt = datump.groupby('time').mean()
        ddatumpt = datump.groupby('time').std()
        x = np.array(datumpt.index)
        y = np.array(datumpt['frac > '+'{:.1G}'.format(threshold)])
        dy = np.array(ddatumpt['frac > '+'{:.1G}'.format(threshold)])
        label = '{:.1G}'.format(Ssub)
        ax.errorbar(x, y, yerr=dy,
                    lw=2,
                    ls='-', marker='o', ms=12,
                    color=colormap(1.0 * Ssubs.index(Ssub) / len(Ssubs)),
                    label=label,
                   )

    ax.set_xlabel('Time [days from infection]', fontsize=fs)
    ax.set_ylabel('frac > '+'{:.1G}'.format(threshold), fontsize=fs)
    ax.grid(True)
    ax.legend(loc=2, title='Subtype B entropy:', fontsize=fs)

    plt.tight_layout()



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Explore conservation levels across patients and subtype',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=pnames,
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+', default=regions,
                        help='Regions to analyze (e.g. F1 p17)')
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

    datap = calculate_fraction_above_threshold(data)

    if plot:
        plot_af_above_threshold(datap, threshold=0.01, singles=False)

        plt.ion()
        plt.show()

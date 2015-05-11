# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/01/15
content:    Study the site frequency spectrum to and against subtype consensus.
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
import hivwholeseq.utils.plot

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.cross_sectional.get_subtype_consensus import (
    get_subtype_reference_alignment_consensus)
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment_entropy)
from hivwholeseq.cross_sectional.get_subtype_allele_frequencies import (
    get_subtype_reference_alignment_allele_frequencies)



# Globals
regions = ['p17', 'p24', 'PR', 'RT', 'IN', 'vif', 'vpu', 'nef', 'gp41']



# Functions
def collect_data(pnames, regions, VERBOSE=0):
    '''Collect data for the SFS'''
    data = []

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    # Restrict to subtype B patients
    patients = patients.loc[patients.Subtype == 'B']

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
                                                                 depth_min=300,
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

                # Iterate over derived alleles
                for ider, der in enumerate(alpha[:4]):
                    # Skip non-mutations
                    if ider == ianc:
                        continue

                    # Get only non-masked time points
                    aft_der = aft[:, ider, posdna]
                    indpost = -aft_der.mask
                    if indpost.sum() == 0:
                        continue
                    timespos = times[indpost]
                    aft_der = aft_der[indpost]
                    
                    mut = anc + '->' + der

                    # Get the difference in subtype abundances
                    afanc_sub = af_sub[ianc, pos_sub]
                    afder_sub = af_sub[ider, pos_sub]
                    daf_sub = afder_sub - afanc_sub
                    Spos_sub = Ssub[pos_sub]
                    conspos_sub = conssub[pos_sub]

                    if (anc == conspos_sub):
                        away_conssub = 'away'
                    elif (der == conspos_sub):
                        away_conssub = 'to'
                    else:
                        away_conssub = 'neither'
                        continue


                    for it, time in enumerate(timespos):
                        af = aft_der[it]
                        data.append((region, pcode,
                                     anc, der, mut,
                                     afanc_sub, afder_sub, daf_sub,
                                     Spos_sub, conspos_sub, away_conssub,
                                     time, af))

    data = pd.DataFrame(data=data, columns=('region', 'pcode',
                                            'anc', 'der', 'mut',
                                            'afancsub', 'afdersub', 'dafsub',
                                            'Ssub', 'conssub', 'awayto',
                                            'time', 'af'))

    return data


def bin_data(data, bincats, VERBOSE=0):
    '''Add bin columns to the data'''
    from hivwholeseq.utils.pandas import add_binned_column
    bins = {}

    if 'time' in bincats:
        if VERBOSE >= 2:
            print 'Bin by time'
        bins_t = 30.5 * np.array([-1, 12, 24, 48, 72, 120])
        binsc_t = 0.5 * (bins_t[1:] + bins_t[:-1])
        add_binned_column(data, 'tbin', 'time', bins=bins_t, clip=True)
        bins['time'] = (bins_t, binsc_t)

    if 'time_rough' in bincats:
        if VERBOSE >= 2:
            print 'Bin by time, roughly'
        bins_tr = 365.24 * np.array([-0.1, 1.5, 3, 10])
        add_binned_column(data, 'trbin', 'time', bins=bins_tr, clip=True)
        bins['time_rough'] = (bins_tr, None)

    if 'entropy' in bincats:
        if VERBOSE >= 2:
            print 'Bin by subtype entropy'
        bins_S = np.array([0, 0.01, 0.05, 0.1, 0.5, 3])
        binsc_S = 0.5 * (bins_S[1:] + bins_S[:-1])
        add_binned_column(data, 'Sbin', 'Ssub', bins=bins_S, clip=True)
        bins['entropy'] = (bins_S, binsc_S)

    if 'af' in bincats:
        if VERBOSE >= 2:
            print 'Bin by allele freq, including fixed'
        logistic_fun = lambda x: 1.0 / (1.0 + 10**(-x))
        bins_af_logit = np.linspace(-2.5, 1.5, 8)
        binsc_af_logit = 0.5 * (bins_af_logit[1:] + bins_af_logit[:-1])
        bins_af = logistic_fun(bins_af_logit)
        binsc_af = logistic_fun(binsc_af_logit)
        add_binned_column(data, 'afbin', 'af', bins=bins_af, clip='upper')
        bins['af'] = (bins_af, binsc_af)

    return bins


def get_sfs(data, bins_af, attrnames=['tbin', 'awayto'], VERBOSE=0, normalize=True):
    '''Get SFS from data'''
    datah = data.loc[(data.loc[:, 'afbin'] != -1) &
                     (data.loc[:, 'afbin'] != len(bins_af) - 1)]
    datah = (datah
             .loc[:, attrnames + ['afbin']]
             .groupby(attrnames + ['afbin'])
             .size()
             .unstack('afbin'))

    datah[np.isnan(datah)] = 0

    if not normalize:
        return datah

    datah = (datah.T / datah.sum(axis=1)).T # The sum of counts
    # (the bin widths)
    binsw_af = bins_af[1:] - bins_af[:-1]
    datah /= binsw_af

    # Add pseudocounts
    vmin = 0.1 * datah[datah > 0].min().min()
    datah[datah < vmin] = vmin

    return datah


def plot_sfs_awayto_time(data, bins, VERBOSE=0):
    '''Plot SFS (time and away/to)'''
    bins_tr, _ = bins['time_rough']
    bins_af, binsc_af = bins['af']
    datah_t = get_sfs(data, bins_af, attrnames=['trbin', 'awayto'], VERBOSE=VERBOSE)

    import seaborn as sns
    sns.set_style('darkgrid')
    fs = 16
    from hivwholeseq.paper_figures.plots import HIVEVO_colormap
    colormap = HIVEVO_colormap('website')

    fig, ax = plt.subplots() 
    for ik, (keys, arr) in enumerate(datah_t.iterrows()):
        time_window = '{:.1f} - {:.1f}'.format(max(0, bins_tr[keys[0]] / 365.24),
                                                  (bins_tr[keys[0] + 1]) / 365.24)
        awayto = keys[1]

        x = binsc_af[np.array(arr.index, int)]
        y = np.array(arr)
        y *= 1e2 / y[0]

        color = colormap(1.0 * keys[0] / len(binsc_t))
        if keys[1] == 'away':
            ls = '--'
        else:
            ls = '-'

        if awayto == 'to':
            label = time_window
        else:
            label = ''

        ax.plot(x, y,
                ls=ls,
                lw=2,
                color=color,
                label=label)

    # Fix legend
    from matplotlib.lines import Line2D
    lines = [Line2D([], [], ls='-', color='k', label='to'),
             Line2D([], [], ls='--', color='k', label='away')]
    handles, _ = ax.get_legend_handles_labels()
    handles.extend(lines)

    ax.set_xlabel('SNV frequency [includes substitutions]', fontsize=fs)
    ax.set_xlim(bins_af[0], bins_af[-1])
    ax.set_xscale('logit')
    ax.set_yscale('log')
    ax.set_ylim(ymax=5e2)
    ax.grid(True)
    ax.legend(handles=handles, loc='upper right', ncol=2, fontsize=fs-2,
              title='  Time from       Subtype\ninfection [yrs]   consensus')
    ax.set_ylabel('Spectrum [normalized by bin width]', fontsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
    plt.tight_layout()


def plot_sfs_awayto_entropy(data, VERBOSE=0):
    '''Plot SFS (entropy and away/to)'''
    bins_S, binsc_S = bins['entropy']
    bins_af, binsc_af = bins['af']
    datah_S = get_sfs(data, bins_af, attrnames=['Sbin', 'awayto'], VERBOSE=VERBOSE)

    fig, ax = plt.subplots() 
    for ik, (keys, arr) in enumerate(datah_S.iterrows()):
        iSbin = keys[0]
        awayto = keys[1]

        x = binsc_af[np.array(arr.index, int)]
        y = np.array(arr)
        y *= 1e2 / y[0]

        color = cm.jet(1.0 * keys[0] / len(binsc_t))
        if keys[1] == 'away':
            ls = '--'
        else:
            ls = '-'
        
        label = ', '.join(['S e ['+str(bins_S[iSbin])+', '+str(bins_S[iSbin+1])+']', awayto])

        ax.plot(x, y,
                ls=ls,
                lw=2,
                color=color, label=label)

    ax.set_xlabel('Allele frequency', fontsize=fs)
    ax.set_xlim(bins_af[0], bins_af[-1])
    ax.set_xscale('logit')
    ax.set_yscale('log')
    ax.set_ylim(ymax=1e4)
    ax.grid(True)
    ax.legend(loc='upper right', ncol=2, fontsize=fs-4, title='Categories')
    ax.set_ylabel('Spectrum [normalized by bin width]', fontsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
    plt.tight_layout()


def plot_af_avg_awayto_entropy(data, VERBOSE=0):
    '''Plot average frequency (simpler plot)'''
    bins_S, binsc_S = bins['entropy']

    fig, ax = plt.subplots(figsize=(5, 4))
    datap = (data
             .loc[:, ['Sbin', 'awayto', 'af']]
             .groupby(['Sbin', 'awayto'], as_index=False)
             .mean()
             .groupby('awayto'))
    for awayto, datum in datap:
        if awayto == 'away':
            ls = '--'
        else:
            ls = '-'

        x = binsc_S[datum['Sbin']]
        y = datum['af']

        ax.plot(x, y, ls=ls, lw=2, label=awayto, color='k')

    ax.set_xlabel('Entropy in subtype [bits]', fontsize=fs)
    ax.set_ylabel('Average frequency', fontsize=fs)
    ax.set_ylim(1e-4, 1e-1)
    ax.set_xlim(1e-3, 5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper left', fontsize=fs)
    ax.grid(True)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
    plt.tight_layout()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Accumulation of minor alleles stratified by abundance difference in subtype',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', action=PatientsAction, default=None,
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+', default=regions,
                        help='Regions to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot results')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot


    data = collect_data(pnames, regions, VERBOSE=VERBOSE)

    bins = bin_data(data, ['time', 'time_rough', 'entropy', 'af'], VERBOSE=VERBOSE)

    if plot:

        plot_sfs_awayto_time(data, VERBOSE=VERBOSE)

        plot_sfs_awayto_entropy(data, VERBOSE=VERBOSE)

        plot_af_avg_awayto_entropy(data, VERBOSE=VERBOSE)

        plt.ion()
        plt.show()



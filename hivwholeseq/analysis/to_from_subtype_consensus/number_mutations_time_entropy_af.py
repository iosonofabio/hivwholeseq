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

from hivwholeseq.utils.miseq import alpha, alphal
from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.cross_sectional.get_subtype_consensus import (
    get_subtype_reference_alignment_consensus)
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment_entropy)
from hivwholeseq.cross_sectional.get_subtype_allele_frequencies import (
    get_subtype_reference_alignment_allele_frequencies)



# Globals
regions = ['p17', 'p24',
           'PR', 'RT', 'IN',
           'vif', 'vpu', 'nef',
           'gp41', 'gp120',
          ]



# Functions
def collect_data(pnames, regions, VERBOSE=0):
    '''Collect data for the SFS'''
    data = []

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    else:
        # Restrict to subtype B patients
        patients = patients.loc[patients.Subtype == 'B']
        patients = patients.loc[-patients.code.isin(['p4', 'p7'])]

    # Get HXB2 for coordinates
    from hivwholeseq.reference import load_custom_reference
    from hivwholeseq.utils.sequence import find_annotation
    ref = load_custom_reference('HXB2', 'gb')

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

        if VERBOSE >= 2:
            print 'Get position of region in HXB2 genomewide coordinates'
        reg_start = find_annotation(ref, region).location.nofuzzy_start

        for ipat, (pname, patient) in enumerate(patients.iterrows()):
            pcode = patient.code
            if VERBOSE >= 2:
                print pname, pcode

            patient = Patient(patient)
            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=500,
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

                # Discard if too many gaps in subtype alignment
                if af_sub[4, pos_sub] > 0.05:
                    continue

                # Ancestral allele
                ianc = icons[posdna]
                anc = alpha[ianc]

                # Discard if the initial time point is already polymorphic
                aft_anc0 = aft[0, ianc, posdna]
                #if aft_anc0 < 0.9:
                #    continue

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
                                     pos_sub, pos_sub + reg_start,
                                     Spos_sub, conspos_sub, away_conssub,
                                     time, af))

    data = pd.DataFrame(data=data, columns=('region', 'pcode',
                                            'anc', 'der', 'mut',
                                            'pos_sub', 'pos_HXB2',
                                            'Ssub', 'conssub', 'awayto',
                                            'time', 'af'))

    return data


def bin_data(data, bincats, include_substitutions=True, VERBOSE=0):
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
            if include_substitutions:
                print 'Bin by allele freq, including fixed'
            else:
                print 'Bin by allele freq, NOT including fixed'
        logistic_fun = lambda x: 1.0 / (1.0 + 10**(-x))
        bins_af_logit = np.linspace(-2.5, 1.5, 8)
        binsc_af_logit = 0.5 * (bins_af_logit[1:] + bins_af_logit[:-1])
        bins_af = logistic_fun(bins_af_logit)
        binsc_af = logistic_fun(binsc_af_logit)
        if include_substitutions:
            clip = 'upper'
        else:
            clip = False
        add_binned_column(data, 'afbin', 'af', bins=bins_af, clip=clip)
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


def get_n_mutations_avg(data, attrnames):
    '''Get the average number of mutations per genome, globally'''
    npat = len(data.pcode.unique())
    return (data
            .loc[:, list(attrnames) + ['af']]
            .groupby(attrnames)
            .sum()
            ['af'] /
            npat)


def get_n_mutations_patientwise(data, attrnames=['time', 'awayto']):
    '''Get the average number of mutations per genome, patient by patient'''
    return (data
            .loc[:, ['pcode'] + list(attrnames) + ['af']]
            .groupby(['pcode'] + list(attrnames))
            .sum()
            ['af'])


def average_n_muts_patients(datap, n_bootstrap=10, VERBOSE=0):
    '''Rebin by time and average over samples'''
    from hivwholeseq.utils.pandas import add_binned_column

    if isinstance(datap, pd.Series):
        datapf = pd.DataFrame(data=datap, dtype=float)
        datapf['time'] = datapf.index.get_level_values('time')
        datapf['awayto'] = datapf.index.get_level_values('awayto')
        datapf['pcode'] = datapf.index.get_level_values('pcode')
    else:
        datapf = datap.copy()
    datapf.index = np.arange(len(datapf))
    bins_t = 365.24 * np.array([0, 0.6, 1.3, 1.5, 4, 6, 10])
    bins_t, binsc_t = add_binned_column(datapf, 'tbin', 'time', bins=bins_t, clip=True)
    
    # Actually perform the mean
    n_muts = (datapf
              .loc[:, ['af', 'tbin', 'awayto']]
              .groupby(['tbin', 'awayto'])
              .mean())

    n_muts.rename(columns={'af': 'nmuts'}, inplace=True)
    n_muts['time'] = binsc_t[n_muts.index.get_level_values('tbin')]
    n_muts['awayto'] = n_muts.index.get_level_values('awayto')

    if n_bootstrap >= 1:
        if VERBOSE >= 2:
            print 'Bootstrapping over patients'
        datagr = datapf.groupby('pcode')
        pcodes = datapf['pcode'].unique().tolist()
        npat = len(pcodes)
        res = pd.DataFrame(index=n_muts.index)
        for i in xrange(n_bootstrap):
            if VERBOSE >= 2:
                print 'Bootstrap:', i+1
            datab = []
            for j in xrange(npat):
                ipat = np.random.randint(npat)
                datab.append(datagr.get_group(pcodes[ipat]).copy())
                datab[-1]['pcode'] = 'BS'+str(j+1)
            datab = pd.concat(datab)
            res['BS'+str(i+1)] = average_n_muts_patients(datab, n_bootstrap=0)['nmuts']
        std = res.std(axis=1)
        n_muts['nmuts_std'] = std

    return n_muts


def average_n_muts_patients_total(datap, n_bootstrap=10, VERBOSE=0):
    '''Rebin by time and average over samples'''
    from hivwholeseq.utils.pandas import add_binned_column

    if isinstance(datap, pd.Series):
        datapf = pd.DataFrame(data=datap, dtype=float)
        datapf['time'] = datapf.index.get_level_values('time')
        datapf['awayto'] = datapf.index.get_level_values('awayto')
        datapf['pcode'] = datapf.index.get_level_values('pcode')
    else:
        datapf = datap.copy()
    datapf.index = np.arange(len(datapf))
    bins_t = 365.24 * np.array([0, 0.6, 1.3, 1.5, 4, 6, 10])
    bins_t, binsc_t = add_binned_column(datapf, 'tbin', 'time', bins=bins_t, clip=True)
    n_muts = (datapf[['af', 'tbin', 'awayto']]
              .groupby(['tbin', 'awayto'])
              .mean()
              .unstack('awayto')
              .sum(axis=1))
    n_muts = pd.DataFrame(n_muts, columns=['nmuts'])
    n_muts['time'] = binsc_t[n_muts.index.get_level_values('tbin')]

    if n_bootstrap >= 1:
        if VERBOSE >= 2:
            print 'Bootstrapping over patients'
        datagr = datapf.groupby('pcode')
        pcodes = datapf['pcode'].unique().tolist()
        npat = len(pcodes)
        res = pd.DataFrame(index=n_muts.index)
        for i in xrange(n_bootstrap):
            if VERBOSE >= 2:
                print 'Bootstrap:', i+1
            datab = []
            for j in xrange(npat):
                ipat = np.random.randint(npat)
                datab.append(datagr.get_group(pcodes[ipat]).copy())
                datab[-1]['pcode'] = 'BS'+str(j+1)
            datab = pd.concat(datab)
            res['BS'+str(i+1)] = average_n_muts_patients_total(datab, n_bootstrap=0)['nmuts']
        std = res.std(axis=1)
        n_muts['nmuts_std'] = std

    return n_muts


def average_n_muts_patients_fraction(datap, n_bootstrap=10, VERBOSE=0):
    '''Rebin by time and average over samples'''
    from hivwholeseq.utils.pandas import add_binned_column

    if isinstance(datap, pd.Series):
        datapf = pd.DataFrame(data=datap, dtype=float)
        datapf['time'] = datapf.index.get_level_values('time')
        datapf['awayto'] = datapf.index.get_level_values('awayto')
        datapf['pcode'] = datapf.index.get_level_values('pcode')
    else:
        datapf = datap.copy()
    datapf.index = np.arange(len(datapf))
    bins_t = 365.24 * np.array([0, 0.6, 1.3, 1.5, 4, 6, 10])
    bins_t, binsc_t = add_binned_column(datapf, 'tbin', 'time', bins=bins_t, clip=True)
    n_muts = (datapf[['af', 'tbin', 'awayto']]
              .groupby(['tbin', 'awayto'])
              .mean()
              .unstack('awayto')
              ['af'])
    frac = n_muts['to'] / (n_muts['to'] + n_muts['away'])
    frac = pd.DataFrame(frac, columns=['fraction to/away'])
    frac['time'] = binsc_t[frac.index.get_level_values('tbin')]

    if n_bootstrap >= 1:
        if VERBOSE >= 2:
            print 'Bootstrapping over patients'
        datagr = datapf.groupby('pcode')
        pcodes = datapf['pcode'].unique().tolist()
        npat = len(pcodes)
        res = pd.DataFrame(index=n_muts.index)
        for i in xrange(n_bootstrap):
            if VERBOSE >= 2:
                print 'Bootstrap:', i+1
            datab = []
            for j in xrange(npat):
                ipat = np.random.randint(npat)
                datab.append(datagr.get_group(pcodes[ipat]).copy())
                datab[-1]['pcode'] = 'BS'+str(j+1)
            datab = pd.concat(datab)
            res['BS'+str(i+1)] = average_n_muts_patients_fraction(datab,
                                            n_bootstrap=0)['fraction to/away']
        std = res.std(axis=1)
        frac['fraction_std'] = std

    return frac


def get_fraction_to_baseline(pnames, regions):
    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    else:
        # Restrict to subtype B patients
        patients = patients.loc[patients.Subtype == 'B']
        patients = patients.loc[-patients.code.isin(['p4', 'p7'])]

    data = []
    for region in regions:
        conssub = get_subtype_reference_alignment_consensus(region, VERBOSE=VERBOSE)

        for ipat, (pname, patient) in enumerate(patients.iterrows()):
            pcode = patient.code
            patient = Patient(patient)
            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=500,
                                                                 depth_min=30,
                                                                 VERBOSE=VERBOSE)
            icons = patient.get_initial_consensus_noinsertions(aft, VERBOSE=VERBOSE,
                                                               return_ind=True)
            consm = alpha[icons]
            coomap = patient.get_map_coordinates_reference(region, refname=('HXB2', region))
            coomapd = dict(coomap[:, ::-1])

            datum = {'cons': 0,
                     'noncons': 0,
                     'pcode': pcode,
                     'region': region,
                    }
            
            for posdna in xrange(aft.shape[-1]):

                # Look for this position in the subtype alignment
                if posdna not in coomapd:
                    continue
                pos_sub = coomapd[posdna]
                if (pos_sub >= len(conssub)):
                    continue

                if conssub[pos_sub] == consm[posdna]:
                    datum['cons'] += 1
                else:
                    datum['noncons'] += 1

            data.append(datum)

    data = pd.DataFrame(data)
    return 1.0 * data['noncons'].sum() / np.array(data[['cons', 'noncons']]).sum()


def get_allele_frequency_entropy(data, n_bootstrap=10, VERBOSE=0, include_substitutions=True):
    '''Get the average allele frequency by entropy away/to'''
    datap = data.copy()
    if not include_substitutions:
        datap = datap.loc[datap['af'] < 0.95]

    datap = (datap
             .loc[:, ['Sbin', 'awayto', 'af']]
             .groupby(['Sbin', 'awayto'], as_index=False)
             .mean())

    # Bootstrap over patients
    if n_bootstrap >= 1:
        if VERBOSE >= 2:
            print 'Bootstrapping over patients'
        datagr = data.groupby('pcode')
        pcodes = data['pcode'].unique().tolist()
        npat = len(pcodes)
        res = pd.DataFrame(index=datap.index)
        for i in xrange(n_bootstrap):
            if VERBOSE >= 2:
                print 'Bootstrap:', i+1
            datab = []
            for j in xrange(npat):
                ipat = np.random.randint(npat)
                datab.append(datagr.get_group(pcodes[ipat]).copy())
                datab[-1]['pcode'] = 'BS'+str(j+1)
            datab = pd.concat(datab) 
            res['BS'+str(i+1)] = get_allele_frequency_entropy(datab, n_bootstrap=0)['af']
        std = res.std(axis=1)
        datap['af_std'] = std

    return datap


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


def plot_af_avg_awayto_entropy(datap, bins, VERBOSE=0):
    '''Plot average frequency (simpler plot)'''
    bins_S, binsc_S = bins['entropy']

    import seaborn as sns
    sns.set_style('darkgrid')
    fs = 16

    fig, ax = plt.subplots(figsize=(5, 4))
    for awayto, datum in datap.groupby('awayto'):
        if awayto == 'away':
            ls = '--'
        else:
            ls = '-'

        x = binsc_S[datum['Sbin']]
        y = datum['af']
        dy = datum['af_std']

        ax.errorbar(x, y, yerr=dy,
                    ls=ls, lw=2,
                    label=awayto,
                    color='k')

    ax.set_xlabel('Entropy in subtype [bits]', fontsize=fs)
    ax.set_ylabel('Average frequency', fontsize=fs)
    ax.set_ylim(1e-5, 1)
    ax.set_xlim(1e-3, 5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper left', fontsize=fs)
    ax.grid(True)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
    plt.tight_layout()


def plot_n_muts_avg(n_muts):
    '''Plot the averages'''
    import seaborn as sns
    sns.set_style('darkgrid')
    fs = 16
    colormap = cm.jet
    fig, ax = plt.subplots(figsize=(5, 4))
    for key, datum in n_muts.groupby('awayto'):
        x = np.array(datum['time'])
        y = np.array(datum['nmuts'])
        if key == 'away':
            label = 'Away from consensus'
            ls = '--'
            color = 'darkred'
        else:
            label = 'Back to consensus'
            ls = '-'
            color = 'steelblue'

        label = key

        ax.plot(x, y,
                lw=2,
                ls=ls,
                marker='o',
                markersize=15,
                color=color,
                label=label
               )

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Time from infection [days]', fontsize=fs)
    ax.set_ylabel('Average number of mutations\nper genome [no env]', fontsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.legend(loc='lower right', fontsize=fs, title='Mutations:')
    ax.grid(True)

    plt.tight_layout()


def plot_n_muts_avg_ratio(n_muts):
    '''Plot the averages'''

    data = []
    for time, datum in n_muts.groupby('time'):
        datum = datum.set_index('awayto')
        data.append({'time': time,
                     'ratio': 1.0 * datum.loc['to', 'af'] / datum.loc['away', 'af'],
                     'sum': datum['af'].sum(),
                    })
    data = pd.DataFrame(data)

    import seaborn as sns
    sns.set_style('dark')
    fs = 16
    colormap = cm.jet
    fig, ax = plt.subplots(figsize=(5, 4))

    # Plot sum
    x = np.array(data['time'])
    y = np.array(data['sum'])
    ax.plot(x, y,
            lw=2,
            marker='o',
            markersize=15,
            color='steelblue',
           )

    ax.set_xlabel('Time from infection [days]', fontsize=fs)
    ax.set_ylabel('Average number of mutations\nper genome [no env]', fontsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.grid(False)

    # Plot ratio
    ax = ax.twinx()
    x = np.array(data['time'])
    y = np.array(data['ratio'])
    ax.plot(x, y,
            lw=2,
            marker='o',
            markersize=15,
            color='darkred',
           )

    ax.set_ylim(0, 0.6)
    ax.set_ylabel('Ratio to/away', fontsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.grid(False)

    plt.tight_layout()


def plot_n_muts_single(datap):
    import seaborn as sns
    sns.set_style('darkgrid')
    fs = 16
    colormap = cm.jet
    fig, ax = plt.subplots()
    pcodes = datap.index.get_level_values('pcode').unique().tolist()
    pcodes.sort(key=lambda x: int(x[1:]))
    for pcode in pcodes:
        datum = datap.loc[pcode].unstack('awayto')
        for key in ('away', 'to'):
            datump = datum[key]
            x = np.array(datump.index)
            y = np.array(datump)
            if key == 'away':
                ls = '--'
                label = ''
            else:
                ls = '-'
                label = pcode

            ax.plot(x, y,
                    lw=2,
                    ls=ls,
                    marker='o',
                    markersize=15,
                    color=colormap(1.0 * pcodes.index(pcode) / len(pcodes)),
                    label=label
                   )

    ax.set_xlabel('Time from infection [days]', fontsize=fs)
    ax.set_ylabel('Average number\nof mutations', fontsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.legend(loc='upper left', fontsize=fs)
    ax.grid(True)

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


    # FIXME: why does p10 have only one time point??
    data = collect_data(pnames, regions, VERBOSE=VERBOSE)

    # Plot afs of the low-entropy to-consensus alleles, since that is the main
    # difference between Richard's and my script
    def plot_afs_to_lowS(data, Smax=0.145):
        import matplotlib.pyplot as plt
        from matplotlib import cm

        d = data.groupby('awayto').get_group('to')
        d = d.loc[d['Ssub'] < Smax]
        for pcode, dp in d.groupby('pcode'):
            fig, ax = plt.subplots()
            ax.set_xlabel('Time [DSI]')
            ax.set_ylabel('af')
            ax.set_title(pcode)
            ax.grid(True)

            for pos, datum in dp.groupby('pos_HXB2'):
                datum = datum.sort('time')
                x = np.array(datum['time'])
                y = np.array(datum['af'])
                ax.plot(x, y, label=str(pos),
                        lw=2,
                        color=cm.jet(1.0 * pos / 9800),
                       )

            plt.tight_layout()

        plt.ion()
        plt.show()

    #plot_afs_to_lowS(data, Smax=0.145)

    bins = bin_data(data, ['time', 'time_rough', 'entropy', 'af'], VERBOSE=VERBOSE,
                    include_substitutions=True)

    # Get global counts
    counts = get_sfs(data, bins['af'][0],
                     attrnames=['trbin', 'awayto'],
                     VERBOSE=VERBOSE, normalize=False)

    # Plot the avg number of mutations per genome, patient by patient
    datap = get_n_mutations_patientwise(data, attrnames=['time', 'awayto'])

    # Rebin by time and average over samples
    n_muts = average_n_muts_patients(datap)

    frac = average_n_muts_patients_fraction(datap, n_bootstrap=100, VERBOSE=3)
    frac0 = get_fraction_to_baseline(None, regions)

    ## NOTE: I tried the avg af without substitutions, tiny changes only
    #af_avg = get_allele_frequency_entropy(data, include_substitutions=False)
    #plot_af_avg_awayto_entropy(af_avg, VERBOSE=VERBOSE)

    af_avg = get_allele_frequency_entropy(data, include_substitutions=True)

    if plot:
        
        plot_n_muts_single(datap)

        plot_n_muts_avg(n_muts)

        plot_af_avg_awayto_entropy(af_avg, bins, VERBOSE=VERBOSE)

        plt.ion()
        plt.show()

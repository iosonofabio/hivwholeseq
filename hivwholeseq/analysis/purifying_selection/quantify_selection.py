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
pnames = ['20097', '15376', '15823', '15241', '9669', '15319']
regions = ['p17', 'p24', 'nef', 'PR', 'RT', 'IN', 'vif']
fun = lambda x, l, u: l * (1 - np.exp(- u/l * x))



# Functions
def collect_data_fitness_cost(regions, pnames, VERBOSE=0):
    '''Collect data for fitness cost estimate'''
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

            # Get only codons with at most one polymorphic site, to avoid obvious epistasis
            ind_poly, _ = get_codons_n_polymorphic(aft, icons, n=[0, 1], VERBOSE=VERBOSE)
            ind_poly_dna = [i * 3 + j for i in ind_poly for j in xrange(3)]

            for posdna in ind_poly_dna:
                # Get the entropy
                if posdna not in coomapd['pat_to_subtype']:
                    continue
                pos_sub = coomapd['pat_to_subtype'][posdna]
                if pos_sub >= len(Ssub):
                    continue
                Ssubpos = Ssub[pos_sub]

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

                # Skip if the site is already polymorphic at the start
                if aftpos[ianc, 0] < 0.95:
                    continue

                # Skip if the site has sweeps (we are looking at purifying selection only)
                # Obviously, it is hard to distinguish between sweeps and unconstrained positions
                if (aftpos[ianc] < 0.6).any():
                    continue

                for inuc, af in enumerate(aftpos[:4]):
                    nuc = alpha[inuc]
                    if nuc == anc:
                        continue

                    mut = anc+'->'+nuc

                    # Define transition/transversion
                    if frozenset(nuc+anc) in (frozenset('CT'), frozenset('AG')):
                        trclass = 'ts'
                    else:
                        trclass = 'tv'

                    # Get the whole trajectory for plots against time
                    for af, time in izip(aftpos[inuc], timespos):
                        data.append((region, pcode,
                                     posdna, pos_sub,
                                     anc, nuc, mut,
                                     trclass,
                                     Ssubpos,
                                     time, af))

    data = pd.DataFrame(data=data,
                        columns=['region', 'pcode',
                                 'posdna', 'possub',
                                 'anc', 'der', 'mut',
                                 'tr',
                                 'Ssub',
                                 'time', 'af'])

    return data


def fit_fitness_cost(x, y, mu, VERBOSE=0):
    '''Fit saturation curve for fitness costs
    
    NOTE: as can be seen from the plots below, the fit for l is not very
    sensitive to mu.
    '''
    from scipy.optimize import minimize_scalar

    fun_min_scalar = lambda s: ((y - fun(x, mu/s, mu))**2).sum()
    s = minimize_scalar(fun_min_scalar, bounds=[1e-5, 1e0]).x

    return s


def fit_saturation(data, bins_S, binsc_S, method='group', VERBOSE=0):
    '''Fit saturation curves to the data
    
    Parameters:
       method (str): whether to fit single allele trajectories ('single') or
       time averages ('group'). The fit result is similar, but noise statistics
       differ.
    '''
    # Every single observation counts one
    if method == 'single':
        dataf = (data
                 .loc[:, ['Sbin', 'mu', 'time', 'af']]
                 .groupby('Sbin'))

    # Avg over alleles from the same time point
    elif method == 'group':
        dataf = (data
                 .loc[:, ['Sbin', 'mut', 'mu', 'time', 'af']]
                 .groupby(['Sbin', 'mut', 'time'], as_index=False)
                 .mean()
                 .groupby('Sbin'))

    # Time bins
    elif method == 'binned':
        from hivwholeseq.utils.pandas import add_binned_column
        bins_t = np.logspace(np.log10(max(1, data['time'].min())),
                             np.log10(data['time'].max()),
                             8)
        binsc_t = np.sqrt(bins_t[1:] * bins_t[:-1])
        add_binned_column(data, 'tbin', 'time', bins_t, clip=True)

        dataf = (data
                 .loc[:, ['Sbin', 'mut', 'mu', 'tbin', 'af']]
                 .groupby(['Sbin', 'mut', 'tbin'], as_index=False)
                 .mean())
        dataf['time'] = binsc_t[dataf['tbin']]
        dataf = dataf.groupby('Sbin')

    # Time bins AND avg mutation rate instead of single ones
    elif method == 'binnedavg':
        from hivwholeseq.utils.pandas import add_binned_column
        bins_t = np.logspace(np.log10(max(1, data['time'].min())),
                             np.log10(data['time'].max()),
                             8)
        binsc_t = np.sqrt(bins_t[1:] * bins_t[:-1])
        add_binned_column(data, 'tbin', 'time', bins_t, clip=True)

        mutd = data.loc[:, ['mut', 'mu']].groupby('mut').mean()
        n_obs = (data
                 .loc[:, ['Sbin', 'mut', 'mu', 'tbin', 'af']]
                 .groupby(['Sbin', 'tbin', 'mut'], as_index=False)
                 .size())
        mu_obs = mutd.loc[n_obs.index.get_level_values('mut')].set_index(n_obs.index)
        mu_obs['n_obs'] = n_obs
        mu_avg = (mu_obs.prod(axis=1).mean(axis=0, level=['Sbin', 'tbin']) / 
                  mu_obs['n_obs'].mean(level=['Sbin', 'tbin']))
        dataf = (data
                 .loc[:, ['Sbin', 'tbin', 'af']]
                 .groupby(['Sbin', 'tbin'])
                 .mean())

        dataf['mu'] = mu_avg
        dataf['time'] = binsc_t[dataf.index.get_level_values('tbin')]
        dataf['Sbin'] = dataf.index.get_level_values('Sbin')
        dataf = dataf.groupby('Sbin')


    fits = []
    for iSbin, datum in dataf:
        x = np.array(datum['time'])
        y = np.array(datum['af'])
        mu = np.array(datum['mu'])

        ind = -(np.isnan(x) | np.isnan(y))
        x = x[ind]
        y = y[ind]
        mu = mu[ind]

        try:
            s = fit_fitness_cost(x, y, mu, VERBOSE=VERBOSE)
            if VERBOSE >= 3:
                regions = np.unique(data['region'])
                plot_fit_single(x, y, s, mu,
                                title=', '.join(regions)+', iSbin = '+str(iSbin))

        except RuntimeError:
            print 'Fit failed, Sbin = ', iSbin
            continue

        fits.append((iSbin, s))

    fits = pd.DataFrame(data=fits,
                        columns=['iSbin', 's'])
    fits['S'] = binsc_S[fits['iSbin']]
    fits['Smin'] = bins_S[fits['iSbin']]
    fits['Smax'] = bins_S[fits['iSbin'] + 1]
    
    return fits


def bootstrap_fit(data, fit_method='group', VERBOSE=0):
    '''Bootstrap fit for error bars'''
    pcodes = np.unique(data['pcode'])
    fits_bs = []
    for i in xrange(10):
        if VERBOSE >= 2:
            print 'Bootstrap n.'+str(i+1)
        pcodes_bs = pcodes[np.random.randint(len(pnames), size=len(pnames))]
        data_bs = pd.concat([data.loc[data['pcode'] == pc] for pc in pcodes_bs])
        fits_tmp = fit_saturation(data_bs, bins_S, binsc_S, method=fit_method, VERBOSE=VERBOSE)
        fits_tmp['bootstrap'] = i
        fits_bs.append(fits_tmp)
    fits_bs = pd.concat(fits_bs)
    return fits_bs[['iSbin', 's']].groupby('iSbin').std()['s']


def plot_function_minimization(x, y, params):
    '''Investigate inconsistencies in fits'''
    fun_min = lambda p: ((y - fun(x, p[0], p[1]))**2).sum()

    p1 = np.logspace(np.log10(params[0]) - 3, np.log10(params[0]) + 3, 10)
    p2 = np.logspace(np.log10(params[1]) - 3, np.log10(params[1]) + 3, 10)

    p1G = np.tile(p1, (len(p2), 1))
    p2G = np.tile(p2, (len(p1), 1)).T
    pG = np.dstack([p1G, p2G])
    z = np.log(np.array([[fun_min(ppp) for ppp in pp] for pp in pG]))

    fig, ax = plt.subplots()
    ax.imshow(z, interpolation='nearest')

    ax.set_xlabel('Log l')
    ax.set_ylabel('Log u')

    plt.ion()
    plt.show()

def plot_fit_single(x, y, s, mu, title=''):
    '''Investigate inconsistencies in fits'''
    fig, ax = plt.subplots()

    ax.scatter(x, y, s=50, color='k', label='data')

    yfit = fun(x, mu/s, mu)
    ax.scatter(x, yfit, s=30, marker='s', color='b', label='fits')

    if title:
        ax.set_title(title)
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Allele frequency')
    ax.set_ylim(1e-6, 1e-1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True)
    ax.legend(loc='upper left')

    plt.ion()
    plt.show()


def plot_fits(fits, title='', VERBOSE=0, data=None):
    '''Plot the fits for purifying selection'''
    if data is not None:
        mu = data.loc[:, 'mu'].max()
    else:
        mu = 5e-6

    ymin = 1e-5

    fig, axs = plt.subplots(1, 2, figsize=(13, 6))
    if title:
        fig.suptitle(title, fontsize=20)
    ax = axs[0]

    # Plot the time-averaged data for one 
    if data is not None:
        datap = (data.loc[data.loc[:, 'mu'] == mu]
                     .loc[:, ['Sbin', 'time', 'af']]
                     .groupby(['Sbin', 'time'])
                     .mean())
        datap['time'] = datap.index.get_level_values('time')
        datap['Sbin'] = datap.index.get_level_values('Sbin')
        datap = datap.groupby('Sbin')
        for iSbin, datum in datap:
            x = np.array(datum['time'])
            # Add pseudocounts to keep plot tidy
            y = np.array(datum['af']) + 1.1 * ymin
            color = cm.jet(1.0 * iSbin / len(fits))
            ax.scatter(x, y, s=40, color=color)

    # Plot the fits
    xfit = np.logspace(0, 3.5, 1000)
    for _, fit in fits.iterrows():
        iSbin = fit['iSbin']
        Smin = fit['Smin']
        Smax = fit['Smax']
        s = fit['s']
        yfit = fun(xfit, mu/s, mu)
        label = ('S e ['+'{:.1G}'.format(Smin)+', '+'{:.1G}'.format(Smax)+']'+
                 ', s = '+'{:.1G}'.format(s))
        color = cm.jet(1.0 * iSbin / len(fits))
        ax.plot(xfit, yfit, color=color, label=label, lw=2)

    ax.set_xlabel('Time [days from infection]')
    ax.set_ylabel('Allele frequency')
    ax.legend(loc='upper left', title='Entropy:', fontsize=14, ncol=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(ymin, 1)
    ax.set_xlim(50, 5000)
    ax.grid(True)

    # Plot the estimated fitness value
    ax3 = axs[1]
    
    if 'ds' in fits.columns:
        ax3.errorbar(fits['S'], fits['s'], yerr=fits['ds'], lw=2, c='k')
    else:
        ax3.plot(fits['S'], fits['s'], lw=2, c='k')

    ax3.set_xlabel('Entropy in subtype [bits]')
    ax3.set_ylabel('Fitness cost')
    ax3.set_ylim(1e-5, 5e-1)
    ax3.set_xlim(1e-3, 2)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.grid(True, which='both')

    plt.tight_layout(rect=(0, 0, 1, 0.96))


def plot_fits_regions(fits, title='', VERBOSE=0, ):
    '''Plot the fits for purifying selection'''
    fig, ax = plt.subplots()
    if title:
        ax.set_title(title, fontsize=20)
    
    regions = np.unique(fits['region']).tolist()

    for i, (region, fitsreg) in enumerate(fits.groupby('region')):
        color = cm.jet(1.0 * i / len(regions))

        if 'ds' in fitsreg.columns:
            ax.errorbar(fitsreg['S'], fitsreg['s'], yerr=fitsreg['ds'], lw=2,
                        color=color, label=region)
        else:
            ax.plot(fitsreg['S'], fitsreg['s'], lw=2, color=color, label=region)

    ax.set_xlabel('Entropy in subtype [bits]')
    ax.set_ylabel('Fitness cost')
    ax.set_ylim(1e-4, 5e-1)
    ax.set_xlim(1e-3, 2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, which='both')
    ax.legend(loc='upper right')

    plt.tight_layout(rect=(0, 0, 1, 0.96))


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Study accumulation of minor alleles for different kinds of mutations',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=pnames,
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+', default=regions,
                        help='Regions to analyze (e.g. p17 p24)')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot results')
    parser.add_argument('--method', default='binned',
                        choices=['single', 'group', 'binned', 'binnedavg'],
                        help='Fit method [group|single]')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot
    fit_method = args.method

    data = collect_data_fitness_cost(regions, pnames, VERBOSE=VERBOSE)


    # Bin by subtype entropy using quantiles
    from hivwholeseq.analysis.purifying_selection.joint_model import add_Sbins
    bins_S, binsc_S = add_Sbins(data, bins=5, VERBOSE=VERBOSE)

    # Get mutation rates from Abram 2010 PER DAY
    from hivwholeseq.analysis.mutation_rate.comparison_Abram import get_mu_Abram2010
    mu = get_mu_Abram2010() / 2
    data['mu'] = np.array(mu.loc[data['mut']])

    # GENOMEWIDE
    fitsgw = fit_saturation(data, bins_S, binsc_S, method=fit_method, VERBOSE=VERBOSE)
    fitsgw['ds'] = bootstrap_fit(data, fit_method=fit_method, VERBOSE=VERBOSE)

    # Store fitness cost to file
    if VERBOSE >= 1:
        print 'Save to file'
    fn_out = get_fitness_cost_entropy_filename('all_muAbram2010_'+fit_method)
    fitsgw.to_pickle(fn_out)

    if plot:
        plot_fits(fitsgw, title=', '.join(regions), VERBOSE=VERBOSE, data=data)

    # REGION BY REGION
    fits = []
    for region, datareg in data.groupby('region'):
        if VERBOSE >= 1:
            print region

        fitsreg = fit_saturation(datareg, bins_S, binsc_S, method=fit_method, VERBOSE=VERBOSE)
        fitsreg['ds'] = bootstrap_fit(datareg, fit_method=fit_method, VERBOSE=VERBOSE)
        fitsreg['region'] = region
        fits.append(fitsreg)

        # Store fitness cost to file
        if VERBOSE >= 1:
            print 'Save to file'
        fn_out = get_fitness_cost_entropy_filename(region+'_muAbram2010_'+fit_method)
        fitsreg.to_pickle(fn_out)

        if plot:
            plot_fits(fitsreg, title=region, VERBOSE=VERBOSE, data=datareg)

    fits = pd.concat(fits)
    fits.set_index(['region', 'iSbin'], drop=False, inplace=True)

    if plot:
        plot_fits_regions(fits, title=', '.join(regions), VERBOSE=VERBOSE)

        plt.ion()
        plt.show()


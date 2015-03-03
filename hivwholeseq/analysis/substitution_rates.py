# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/09/14
content:    Get substitution rate from divergence.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient




# Globals
pnames = ['20097', '15363', '15823', '15376', '20529', '9669', '15241', '15319']
regions = ['p17', 'PR', 'RT', 'vif', 'V3', 'RRE', 'nef']



# Functions
def plot_divergence_fit(data, VERBOSE=0):
    from operator import itemgetter

    fig, axs = plt.subplots(1, 2, figsize=(10, 5),
                            gridspec_kw={'width_ratios':[4, 1]})
    ax = axs[0]
    ax.set_xlabel('Time from transmission [days]')
    ax.set_ylabel('Divergence')
    #ax.set_yscale('log')

    for i, d in enumerate(data):
        pcode = d['pcode']
        region = d['region']
        dg = d['dg']
        times = d['t']
        r = d['r']

        ax.plot(times, dg, lw=2, label=', '.join([pcode, region]),
                color=cm.jet(1.0 * i / len(data)))

        ax.plot(times, times * r, lw=2, ls='--',
                color=cm.jet(1.0 * i / len(data)))

        axs[1].scatter(0, r * 365.24, color=cm.jet(1.0 * i / len(data)),
                       label=', '.join([pcode, region]))

    ax.legend(loc=4, fontsize=12)
    ax.grid(True)

    ax = axs[1]
    ax.xaxis.set_ticklabels('')
    ax.set_ylabel('Substitution rate [changes/site/yr]')
    ax.set_ylim(0.96 * 365.24 * min(data, key=itemgetter('r'))['r'],
                1.04 * 365.24 * max(data, key=itemgetter('r'))['r'])
    ax.set_xlim(-0.5, 0.5)
    ax.grid(True)
    plt.tight_layout()


def plot_substitution_rate(data, title='', VERBOSE=0):
    '''Plot the substitution rates'''

    if VERBOSE:
        print 'Plot substitution rates'
    
    # Sort patient codes
    pcodes = sorted(set(data['pcode']), key=lambda x: int(x[1:]))

    # Resort regions based on average substitution rate
    regions = (data[['region', 'rate']]
               .groupby('region')
               .mean()
               .sort('rate')
               .index
               .tolist())

    xfun = lambda region, pcode: regions.index(region) + 0.05 * pcodes.index(pcode)
    cfun = lambda pcode: cm.jet(1.0 * pcodes.index(pcode) / len(pcodes))

    fig, ax = plt.subplots(figsize=(2 * len(regions), 4))

    for pcode in pcodes:
        datum = (data
                 .loc[data['pcode'] == pcode]
                 .set_index('region', drop=False)
                 .loc[regions])

        y = datum['rate']
        x = map(lambda x: xfun(*x[1]), datum[['region', 'pcode']].iterrows())
        color = cfun(pcode)

        ax.scatter(x, y, color=color, s=90, label=pcode)
        if len(regions) > 1:
            ax.plot(x, y, color=color, lw=1, alpha=0.4)

    ylim = (0.5 * data['rate'].min(), 1.5 * data['rate'].max())
    xticksminor = [xfun(region, pcodes[len(pcodes) //2]) for region in regions]
    xticksmajor = [xt - 0.5 for xt in xticksminor] + [xticksminor[-1] + 0.5]
    xticklabels = regions
    xlim = (xticksminor[0] - 0.04 * (xticksminor[-1] - xticksminor[0]),
            xticksminor[0] + 1.04 * (xticksminor[-1] - xticksminor[0]))


    ax.set_ylim(*ylim)
    ax.set_xlim(*xlim)
    ax.set_xticks(xticksmajor)
    ax.set_xticklabels([])
    ax.set_xticks(xticksminor, minor=True)
    ax.set_xticklabels(xticklabels, fontsize=14, minor=True)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=14)
    ax.grid(True)
    ax.set_yscale('log')

    ax.set_xlabel('Genomic region', fontsize=14)
    ax.set_ylabel('Substitution rate\n[changes / year / site]', labelpad=10, fontsize=14)
    ax.legend(loc='lower right', title='Patients', ncol=len(pcodes) // 3,
              fontsize=14)

    if title:
        ax.set_title(title)

    plt.tight_layout(rect=(0, 0, 0.98, 1))




# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Measure substitution rates in different regions and patients',
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

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    pcodes = patients.code.tolist()

    data = []
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        for ifr, region in enumerate(regions):
            if VERBOSE >= 1:
                print pname, patient.code, region

            try:
                dg, ind = patient.get_divergence(region, cov_min=10)
            except ValueError:
                continue
            times = patient.times[ind]

            data.append({'pcode': patient.code, 'region': region, 'dg': dg, 't': times})

    if VERBOSE >= 1:
        print 'Fit slopes'

    for d in data:
        pcode = d['pcode']
        region = d['region']
        dg = d['dg']
        times = d['t']

        r = np.linalg.lstsq(times[:, np.newaxis], dg)[0][0]
        d['r'] = r


    if VERBOSE >= 1:
        print 'Prepare data for plot'

    datap = pd.DataFrame([{'pcode': d['pcode'],
                           'region': d['region'],
                           'rate': d['r'] * 365.24}
                          for d in data])

    if plot:
        if VERBOSE >= 1:
            print 'Plot'

        plot_divergence_fit(data, VERBOSE=VERBOSE)

        plot_substitution_rate(datap, VERBOSE=VERBOSE)

        plt.ion()
        plt.show()

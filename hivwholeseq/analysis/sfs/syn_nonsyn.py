# vim: fdm=indent
'''
author:     Fabio Zanini
date:       14/04/15
content:    Make the site frequency spectrum of syn/nonsyn derived alleles to
            plot in a linear scale for the paper.
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
from Bio.Seq import translate

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.patients.one_site_statistics import get_codons_n_polymorphic



# Globals
pnames = ['20097', '15363', '15823', '9669', '20529', '15241', '15376', '15319']
regions = ['p17', 'p24', 'PR', 'RT', 'IN', 'vif', 'vpu', 'nef', 'gp41']



# Functions
def collect_data(pnames, regions, VERBOSE=0):
    '''Collect data for the SFS'''
    data = []

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    # FIXME
    pbads = ('15313', '15107')
    patients = patients.loc[-patients.index.isin(pbads)]

    for region in regions:
        if VERBOSE >= 1:
            print region

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
                print 'Get initial consensus'
            icons = patient.get_initial_consensus_noinsertions(aft, VERBOSE=VERBOSE,
                                                               return_ind=True)
            consm = alpha[icons]

            if VERBOSE >= 2:
                print 'Get alleles'
            for posdna in xrange(aft.shape[-1]):
                if VERBOSE >= 3:
                    print posdna

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

                    # Find out whether it's syn or nonsyn
                    codanc = consm[posdna - posdna % 3: posdna - posdna % 3 + 3]
                    codder = codanc.copy()
                    codder[posdna % 3] = der
                    codanc = ''.join(codanc)
                    codder = ''.join(codder)
                    if ('-' in codanc) or ('-' in codder):
                        continue

                    if translate(codanc) == translate(codder):
                        mutclass = 'syn'
                    else:
                        mutclass = 'nonsyn'

                    for it, time in enumerate(timespos):
                        af = aft_der[it]
                        data.append((region, pcode,
                                     anc, der, mut,
                                     mutclass,
                                     time, af))

    data = pd.DataFrame(data=data, columns=('region', 'pcode',
                                            'anc', 'der', 'mut',
                                            'mutclass',
                                            'time', 'af'))

    return data


def get_sfs(data, bins_af='linear', attrnames=['tbin', 'mutclass'], VERBOSE=0,
            normalize='fraction'):
    '''Get SFS from data'''
    # Bin
    if bins_af == 'linear':
        bins_af = np.linspace(1e-3, 1, 10)
        binsc_af = 0.5 * (bins_af[1:] + bins_af[:-1])
    elif bins_af == 'logit':
        logistic_fun = lambda x: 1.0 / (1.0 + 10**(-x))
        bins_af_logit = np.linspace(-2.5, 1.5, 8)
        binsc_af_logit = 0.5 * (bins_af_logit[1:] + bins_af_logit[:-1])
        bins_af = logistic_fun(bins_af_logit)
        binsc_af = logistic_fun(binsc_af_logit)

    binsw_af = bins_af[1:] - bins_af[:-1]
    data['afbin'] = -1
    for b in bins_af:
        data.loc[data.loc[:, 'af'] >= b, 'afbin'] += 1

    # Classify
    datah = data.loc[(data.loc[:, 'afbin'] != -1) &
                     (data.loc[:, 'afbin'] != len(bins_af) - 1)]

    datah = (datah
             .loc[:, attrnames + ['afbin']]
             .groupby(attrnames + ['afbin'])
             .size()
             .unstack('afbin'))
    datah[np.isnan(datah)] = 0
    datah = datah.T

    # Normalize
    if normalize in ['fraction', 'density']:
        datah /= datah.sum(axis=0) # The sum of counts

    # (the bin widths)
    if normalize in ['density']:
        binsw_af = bins_af[1:] - bins_af[:-1]
        datah /= binsw_af

    # Add pseudocounts
    vmin = 0.1 * datah[datah > 0].min().min()
    datah[datah < vmin] = vmin

    # Add bin edges
    datah['afbin_left'] = bins_af[datah.index]
    datah['afbin_right'] = bins_af[datah.index + 1]
    datah['afbin_center'] = binsc_af[datah.index]

    return datah


def plot_sfs_synnonsyn(datah):
    plot_props = [{'mutclass': 'syn',
                   'color': 'steelblue',
                   'label': 'synonymous',
                  },
                  {'mutclass': 'nonsyn',
                   'color': 'darkred',                      
                   'label': 'nonsynonymous',
                  }]

    fig, ax = plt.subplots() 
    for ip, props in enumerate(plot_props):
        mutclass = props['mutclass']
        color = props['color']
        label = props['label']

        arr = datah.loc[:, mutclass]
        y = np.array(arr)
        x = datah.loc[:, 'afbin_left']
        w = datah.loc[:, 'afbin_right'] - x

        n_props = len(plot_props)
        ws = 1.0 / (n_props + 1)
        x += ws * w * ip

        bottom = 1e-3
        height = y - bottom
        ax.bar(x, height, width=ws * w, bottom=bottom,
               color=color,
               label=label)


    ax.set_xlabel('Allele frequency')
    ax.set_xlim(0, 1)
    ax.set_yscale('log')
    ax.grid(True)
    ax.legend(loc='upper right', ncol=1, fontsize=14)
    ax.set_ylabel('Spectrum')

    plt.ion()
    plt.show()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Accumulation of minor alleles stratified by abundance difference in subtype',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', action=PatientsAction, default=pnames,
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

    ## Bin by time
    #bins_t = 30.5 * np.array([-1, 12, 24, 48, 72, 120])
    #binsc_t = 0.5 * (bins_t[1:] + bins_t[:-1])
    #data['tbin'] = 0
    #for b in bins_t[1:]:
    #    data.loc[data.loc[:, 'time'] >= b, 'tbin'] += 1

    # SFS
    datah = get_sfs(data,
                    bins_af='linear',
                    normalize='fraction',
                    attrnames=['mutclass'], VERBOSE=VERBOSE)

    if plot:

        plot_sfs_synnonsyn(datah)

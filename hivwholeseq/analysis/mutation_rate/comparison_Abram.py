# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/01/15
content:    Study divergence at conserved/non- synonymous sites in different
            patients.
'''
# Modules
import os
import argparse
from itertools import izip
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.analysis.mutation_rate.explore_divergence_synonymous import (
    collect_data_mutation_rate, fit_mutation_rate)


# Globals



# Functions
def get_mu_Abram2010():
    '''Get the mutation rate matrix from Abram 2010'''
    muts = [a+'->'+b for a in alpha[:4] for b in alpha[:4] if a != b]

    muAbram = pd.Series(np.zeros(len(muts)), index=muts, dtype=float)
    muAbram.name = 'mutation rate Abram 2010'

    # Forward strand
    muAbram['C->A'] += 14
    muAbram['G->A'] += 146
    muAbram['T->A'] += 20
    muAbram['A->C'] += 1
    muAbram['G->C'] += 2
    muAbram['T->C'] += 18
    muAbram['A->G'] += 29
    muAbram['C->G'] += 0
    muAbram['T->G'] += 6
    muAbram['A->T'] += 3
    muAbram['C->T'] += 81
    muAbram['G->T'] += 4

    # Reverse strand
    muAbram['C->A'] += 24
    muAbram['G->A'] += 113
    muAbram['T->A'] += 32
    muAbram['A->C'] += 1
    muAbram['G->C'] += 2
    muAbram['T->C'] += 25
    muAbram['A->G'] += 13
    muAbram['C->G'] += 1
    muAbram['T->G'] += 8
    muAbram['A->T'] += 0
    muAbram['C->T'] += 61
    muAbram['G->T'] += 0

    # Normalize with the max (how do they calculate the average??)
    muAbramAv = 1e-5
    muAbram *= muAbramAv / muAbram.max()

    return muAbram


def comparison_Abram2010(mu, VERBOSE=2, title=''):
    '''Print a comparison with Abram 2010, J. Virol.'''
    muAbram = get_mu_Abram2010()
    data = pd.concat([muAbram, mu], axis=1)
    if VERBOSE >= 2:
        if title:
            print title
        print data

    data.rename(columns={0: 'Abram2010', 1: 'new'}, inplace=True)
    return data


def plot_mu_matrix(mu, time_unit='generation'):
    '''Plot the mutation rate matrix'''
    muM = np.ma.masked_all((4, 4))
    for ia1, a1 in enumerate(alpha[:4]):
        for ia2, a2 in enumerate(alpha[:4]):
            if a1+'->'+a2 in mu.index:
                muM[ia1, ia2] = mu.loc[a1+'->'+a2]

    if time_unit == 'generation':
        # Assume a generation time of 2 days
        muM *= 2

    fig, ax = plt.subplots()
    h = ax.imshow(np.log10(muM.T + 1e-10),
                  interpolation='nearest',
                  vmin=-9, vmax=-4)

    ax.set_xticks(np.arange(4))
    ax.set_yticks(np.arange(4))
    ax.set_xticklabels(alpha[:4])
    ax.set_yticklabels(alpha[:4])

    ax.set_xlabel('From:')
    ax.set_ylabel('To:')

    cb = fig.colorbar(h)
    cb.set_ticks(np.arange(-9, -3))
    cb.set_ticklabels(['$10^{'+str(x)+'}$' for x in xrange(-9, -3)])
    cb.set_label('changes / '+time_unit, rotation=270, labelpad=30)

    plt.tight_layout(rect=(-0.1, 0, 1, 1))

    plt.ion()
    plt.show()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Get mutation rate from the accumulation of synonymous minor alleles',
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

    data = defaultdict(list)

    patients = load_patients()
    if pnames is None:
        pnames = patients.index.tolist()

    data = collect_data_mutation_rate(regions, pnames, VERBOSE=VERBOSE)
    mu = fit_mutation_rate(data, VERBOSE=VERBOSE, Smin=0.1)

    print comparison_Abram2010(mu)
    
    if plot:
        plot_mu_matrix(mu)

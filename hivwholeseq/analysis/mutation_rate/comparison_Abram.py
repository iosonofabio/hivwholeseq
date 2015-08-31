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

from hivwholeseq.utils.miseq import alpha, alphal
from hivwholeseq.analysis.mutation_rate.explore_divergence_synonymous import (
    collect_data_mutation_rate, fit_mutation_rate)

from hivwholeseq.analysis.mutation_rate.mutation_rate import plot_mu_matrix

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

    # Normalize
    muAbramAv = 1.3e-5
    muAbram *= muAbramAv / (muAbram.sum() / 4.0)

    return muAbram


def add_Abram2010(mu, generation_time=2.0, VERBOSE=2, title=''):
    '''Print a comparison with Abram 2010, J. Virol.'''
    muAbram = get_mu_Abram2010() / generation_time
    data = pd.concat([muAbram, mu], axis=1)
    if VERBOSE >= 3:
        if title:
            print title
        print data

    if not title:
        title = 'new'

    data.rename(columns={0: 'Abram2010', 1: title}, inplace=True)
    return data


def plot_comparison(mu, method='joint'):
    '''Plot comparison between our estimate and Abram2010'''
    mus = add_Abram2010(mu)

    fig, ax = plt.subplots()
    x = np.array(mus['Abram2010'])
    y = np.array(mus['new'])

    r = np.corrcoef(np.log10(x), np.log10(y))[0, 1]

    ax.scatter(x, y, s=40, color='k', label='{:2.0%}'.format(r))
    xl = np.logspace(-8, -4, 100)
    ax.plot(xl, xl, lw=2, c='grey')
    ax.set_xlabel('Abram et al 2010 [changes / site / day]')
    ax.set_ylabel(method.capitalize()+' estimate [changes / site / day]')
    ax.set_xlim(1e-8, 1e-4)
    ax.set_ylim(1e-8, 1e-4)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper left', title='Correlation\ncoefficient:', fontsize=14)
    ax.grid(True)

    plt.tight_layout()


def plot_three_comparison():
    '''Plot three way comparison'''
    # Load mutation rate from joint model
    from hivwholeseq.analysis.filenames import analysis_data_folder
    mu_neu = pd.read_pickle(analysis_data_folder+'mu_neutralclass.pickle')
    mu_joi = pd.read_pickle(analysis_data_folder+'mu_joint.pickle')
    mu_Abr = get_mu_Abram2010() / 2.0

    data = pd.concat([mu_Abr, mu_neu, mu_joi], axis=1)
    data.rename(columns={0: 'Abram2010', 1: 'neutralclass', 2: 'joint'}, inplace=True)
    print data

    fig, axs = plt.subplots(2, 1, figsize=(6, 7))
    axfun = lambda x: axs[x >= len(data.index) // 2]
    colord = {'neutralclass': 'darkred', 'joint': 'steelblue'}


    for i, ax in enumerate(axs):
        width = 0.4
        xbase = np.arange(len(data.index) // 2) - width

        for j, key in enumerate(['neutralclass', 'joint']):
            xs = xbase + width * j
            ys = ((data[key] / data['Abram2010'] - 1)
                  .iloc[i * (len(data.index) // 2): (i + 1) * (len(data.index) // 2)])

            ax.bar(xs, ys, width=width, bottom=1,
                   color=colord[key],
                   label=key,
                  )

        ax.plot([-1, len(data.index) // 2], [1, 1], color='k', lw=1)
        muts = data.index[i * (len(data.index) // 2): (i + 1) * (len(data.index) // 2)].tolist()
        ax.set_xticks(np.arange(len(data.index) // 2))
        ax.set_xticklabels(muts)
        ax.set_yscale('log')

        ax.set_ylim(1e-1, 1e1)
        ax.grid(True, which='both', axis='y')

    ax.set_ylabel('Ratio to Abram et al 2010', position=(-0.02, 1))
    axs[0].legend(loc='upper right')

    plt.tight_layout(rect=(0.02, 0, 1, 1), h_pad=-4.5)





# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Get mutation rate from the accumulation of synonymous minor alleles',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--method', choices=['three', 'joint', 'neutralclass'], default='joint',
                        help='Which method to compare')

    args = parser.parse_args()
    pnames = args.patients
    VERBOSE = args.verbose
    method = args.method

    if method != 'three':
        # Load mutation rate from joint model
        from hivwholeseq.analysis.filenames import analysis_data_folder
        fn_mu = analysis_data_folder+'mu_'+method+'.pickle'
        mu = pd.read_pickle(fn_mu)
        plot_mu_matrix(mu, time_unit='days')
        plot_comparison(mu, method=method)
    else:
        plot_three_comparison()

    plt.ion()
    plt.show()

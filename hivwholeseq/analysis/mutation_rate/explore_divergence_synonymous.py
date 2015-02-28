# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/01/15
content:    Study divergence at conserved/non- synonymous sites in different
            patients.
'''
# Modules
import os
import sys
import argparse
from itertools import izip
from collections import defaultdict
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.sequence import translate_with_gaps
from hivwholeseq.utils.sequence import get_degeneracy_dict
import hivwholeseq.utils.plot
from hivwholeseq.cross_sectional.get_subtype_entropy_synonymous import \
    get_subtype_reference_alignment_entropy_syn
from hivwholeseq.cross_sectional.get_subtype_consensus import \
    get_subtype_reference_alignment_consensus
from hivwholeseq.patients.one_site_statistics import get_codons_n_polymorphic



# Functions
def translate_masked(cons):
    '''Translate a DNA sequence with Ns and weird gaps'''
    if len(cons) % 3:
        raise ValueError('The length of the DNA sequence is not a multiple of 3')

    from Bio.Seq import translate

    codons = cons.reshape((len(cons)//3, 3))

    # Mask codons that have 1 or 2 gaps, or any N
    mask_gaps = np.array((codons == '-').sum(axis=1) % 3, bool)
    mask_N = (codons == 'N').sum(axis=1) != 0

    prot = np.ma.empty(len(mask_gaps), 'S1')
    prot.mask = mask_gaps | mask_N

    for i, (m, codon) in enumerate(izip(prot.mask, codons)):
        if not m:
            if '-' in codon:
                prot[i] = '-'
            else:
                prot[i] = translate(''.join(codon))

    return prot


def collect_data_mutation_rate(regions, pnames, VERBOSE=0):
    '''Collect data for the mutation rate estimate'''

    degeneracy = get_degeneracy_dict()

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
        Ssub = get_subtype_reference_alignment_entropy_syn(region, VERBOSE=VERBOSE)
        # NOTE: Ssub is indexed by AMINO ACID, conssub by NUCLEOTIDE

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

            # Get only sites that are conserved or with only one site per codon changing
            ind_poly, pos_poly = get_codons_n_polymorphic(aft, icons, n=[0, 1], VERBOSE=VERBOSE)

            # NOTE: all 4-fold degenerate codons have all and only 3rd positions,
            # so they are easy to score, we ignore the rest for now
            # NOTE: the lambda is there for gaps
            prot_deg = np.array(map(lambda x: x if x else 0,
                                    map(degeneracy.get, protm)), int)
            ind_4fold = (prot_deg == 4).nonzero()[0]


            # Get only polymorphic codons that are also 4fold degenerate,
            # in which the 1st and 2nd positions are conserved
            inds = [ip for ip, pp in izip(ind_poly, pos_poly)
                    if (ip in ind_4fold) and (0 not in pp) and (1 not in pp)]

            # FIXME: deal better with depth (should be there already, I guess)
            aft[aft < 2e-3] = 0

            for pos in inds:
                if VERBOSE >= 3:
                    print pos

                # Only look at third position, so 
                posdna = pos * 3 + 2

                # Get the entropy
                if posdna not in coomapd['pat_to_subtype']:
                    continue
                pos_sub = coomapd['pat_to_subtype'][posdna]
                if (pos_sub not in Ssub) or (protm[pos] not in Ssub[pos_sub]):
                    continue
                Ssubpos = Ssub[pos_sub][protm[pos]]

                aft_pos = aft[:, :, posdna]

                # Get only non-masked time points
                indpost = -aft_pos[:, 0].mask
                if indpost.sum() == 0:
                    continue
                timespos = times[indpost]
                aft_pos = aft_pos[indpost]

                # Discard if the initial time point is already polymorphic
                aft_anc = aft_pos[:, icons[posdna]]
                if aft_anc[0] < 0.9:
                    continue

                # Iterate over the three derived nucleotides
                for inuc, aft_der in enumerate(aft_pos.T[:4]):
                    if inuc == icons[posdna]:
                        continue

                    anc = consm[posdna]
                    der = alpha[inuc]
                    mut = anc+'->'+der

                    for it, time in enumerate(timespos):
                        af = aft_der[it]
                        data.append((region, pcode,
                                     anc, der, mut,
                                     Ssubpos,
                                     time, af))


    data = pd.DataFrame(data=data, columns=('region', 'pcode',
                                            'anc', 'der', 'mut',
                                            'S',
                                            'time', 'af'))

    return data


def fit_mutation_rate_avgsample(data, VERBOSE=0, Smin=0.05):
    '''Fit the mutation rate from the data, each sample weighs 1'''
    muts = list(np.unique(data.loc[:, 'mut']))

    # Average within time points first
    datap = (data
             .loc[data.loc[:, 'S'] >= Smin]
             .loc[:, ['mut', 'time', 'af']]
             .groupby(['mut', 'time'])
             .mean()
             .loc[:, 'af']
             .unstack('time'))

    fits = {}
    for mut, datum in datap.iterrows():
        x = np.array(datum.index)
        y = np.array(datum)

        # Get rid of masked stuff (this happens if we miss data in the means)
        ind = -(np.isnan(x) | np.isnan(y))
        x = x[ind]
        y = y[ind]

        # Linear fit
        m = np.dot(y, x) / np.dot(x, x)
        fits[mut] = m

    return pd.Series(fits)


def fit_mutation_rate(data, VERBOSE=0, Smin=0.05):
    '''Fit the mutation rate from the data, each observation weighs 1'''
    muts = list(np.unique(data.loc[:, 'mut']))
    
    fits = {}
    for mut in muts:
        datum = data.loc[data.loc[:, 'mut'] == mut]

        x = np.array(datum['time'])
        y = np.array(datum['af'])

        # Get rid of masked stuff (this happens if we miss data in the means)
        ind = -(np.isnan(x) | np.isnan(y))
        x = x[ind]
        y = y[ind]

        # Linear fit
        m = np.dot(y, x) / np.dot(x, x)
        fits[mut] = m

    return pd.Series(fits)


def plot_mutation_rate(data, fits, VERBOSE=0):
    '''Plot mutation rate estimate'''
    muts = list(np.unique(data.loc[:, 'mut']))

    Sbins = [(0, 0.1), (0.1, 2)]
    Sbin = Sbins[-1]
    datap = (data
             .loc[(data.loc[:, 'S'] > Sbin[0]) & (data.loc[:, 'S'] <= Sbin[1])]
             .loc[:, ['mut', 'time', 'af']]
             .groupby(['mut', 'time'])
             .mean()
             .loc[:, 'af']
             .unstack('time'))


    fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 5),
                                       gridspec_kw={'width_ratios': [10, 12, 0.8]})

    for mut, datum in datap.iterrows():
        x = np.array(datum.index) / 30.5
        y = np.array(datum)

        # Get rid of masked stuff (this happens if we miss data in the means)
        ind = -(np.isnan(x) | np.isnan(y))
        x = x[ind]
        y = y[ind]

        color = cm.jet(1.0 * muts.index(mut) / len(muts))
        ax.scatter(x, y, color=color, label=mut)

        m = fits.loc[mut]
        xfit = np.linspace(0, 150)
        ax.plot(xfit, m * xfit, color=color, lw=1.5, alpha=0.5)

    ax.set_xlabel('Time from infection [months]')
    ax.set_ylabel('Allele frequency')
    ax.set_ylim(1e-4, 1)
    ax.set_yscale('log')
    ax.grid(True)
    ax.legend(loc='upper left', fontsize=10, ncol=2)

    mu = np.ma.masked_all((4, 4))
    for mut, fit in fits.iteritems():
        n1 = mut[0]
        n2 = mut[-1]
        # The mutation rate we want per generation (2 days)
        mu[alphal.index(n1), alphal.index(n2)] = fit * 2

    # Plot the log10 for better dynamic range
    vmin = np.floor(np.log10(mu.min()))
    vmax = np.ceil(np.log10(mu.max()))
    h = ax2.imshow(np.log10(mu), interpolation='nearest',
                   vmin=vmin, vmax=vmax)
    
    ax2.set_xticks(np.arange(4) + 0.0)
    ax2.set_yticks(np.arange(4) + 0.0)
    ax2.set_xticklabels(alphal[:4])
    ax2.set_yticklabels(alphal[:4])
    ax2.set_xlabel('to')
    ax2.set_ylabel('from')

    cticks = np.arange(vmin, vmax+1)

    from matplotlib.ticker import Formatter
    class LogTickFormatter(Formatter):
        def __call__(self, x, pos=None):
            '''Transform the logs into their 10**'''
            return '$10^{'+'{:1.0f}'.format(x)+'}$'

    cb = plt.colorbar(h, cax=ax3, ticks=cticks, format=LogTickFormatter())
    cb.set_label('Mutation rate\n[changes / generation]', rotation=90, labelpad=-100)

    plt.tight_layout(w_pad=0.05)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Explore conservation levels across patients and subtype',
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

    data = collect_data_mutation_rate(regions, pnames, VERBOSE=VERBOSE)
    fits = fit_mutation_rate(data, VERBOSE=VERBOSE)

    if plot:
        plot_mutation_rate(data, fits, VERBOSE=VERBOSE)

        plt.ion()
        plt.show()


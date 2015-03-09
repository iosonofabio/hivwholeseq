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
from hivwholeseq.cross_sectional.get_subtype_entropy_synonymous import \
    get_subtype_reference_alignment_entropy_syn
from hivwholeseq.cross_sectional.get_subtype_consensus import \
    get_subtype_reference_alignment_consensus
from hivwholeseq.patients.one_site_statistics import get_codons_n_polymorphic



# Globals
pnames = ['20097', '15376', '15823', '15241', '9669', '15319']
regions = ['p17', 'p24', 'nef', 'PR', 'RT', 'IN', 'vif']



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


def fit_mutation_rate(data, VERBOSE=0, Smin=0.05, method='group'):
    '''Fit the mutation rate from the data, each observation weighs 1'''
    muts = list(np.unique(data.loc[:, 'mut']))
    
    if method == 'single':
        dataf = (data
                 .loc[data.loc[:, 'S'] >= Smin]
                 .loc[:, ['mut', 'time', 'af']]
                 .groupby(['mut']))
    else:
        dataf = (data
                 .loc[data.loc[:, 'S'] >= Smin]
                 .loc[:, ['mut', 'time', 'af']]
                 .groupby(['mut', 'time'])
                 .mean())
        dataf['time'] = dataf.index.get_level_values('time')
        dataf['mut'] = dataf.index.get_level_values('mut')
        dataf = dataf.groupby(['mut'])

    fits = {}
    for mut, datum in dataf:
        x = np.array(datum['time'])
        y = np.array(datum['af'])

        # Get rid of masked stuff
        ind = -(np.isnan(x) | np.isnan(y))
        x = x[ind]
        y = y[ind]

        # Linear fit
        m = np.dot(y, x) / np.dot(x, x)
        fits[mut] = m

    fits =  pd.Series(fits)
    return fits


def plot_mu_matrix(mu, time_unit='days'):
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
    cb.set_label('changes / position / '+time_unit, rotation=270, labelpad=30)

    plt.tight_layout(rect=(-0.1, 0, 1, 1))

    plt.ion()
    plt.show()


def plot_fits(fits, title='', VERBOSE=0, data=None):
    '''Plot the fits for purifying selection'''
    ymin = 1e-5

    muts = fits.index.tolist()

    fig, ax = plt.subplots()
    ax.set_title(title, fontsize=20)

    # Plot the time-averaged data
    datap = (data.loc[:, ['mut', 'time', 'af']]
                 .groupby(['mut', 'time'])
                 .mean())
    datap['time'] = datap.index.get_level_values('time')
    datap['mut'] = datap.index.get_level_values('mut')
    datap = datap.groupby('mut')
    for mut, datum in datap:
        x = np.array(datum['time'])
        # Add pseudocounts to keep plot tidy
        y = np.array(datum['af']) + 1.1 * ymin

        color = cm.jet(1.0 * muts.index(mut) / len(fits))

        ax.scatter(x, y, s=40, color=color)

    # Plot the fits
    xfit = np.logspace(0, 3.5, 1000)
    for mut, mu in fits.iteritems():
        yfit = mu * xfit

        label = mut
        color = cm.jet(1.0 * muts.index(mut) / len(fits))

        ax.plot(xfit, yfit, color=color, label=label, lw=2)

    ax.set_xlabel('Time [days from infection]')
    ax.set_ylabel('Allele frequency')
    ax.legend(loc='upper left', fontsize=14, ncol=2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(ymin, 1)
    ax.set_xlim(100, 5000)
    ax.grid(True)

    plt.tight_layout()




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Explore conservation levels across patients and subtype',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=pnames,
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

    data = collect_data_mutation_rate(regions, pnames, VERBOSE=VERBOSE)

    mu_neu = fit_mutation_rate(data, VERBOSE=VERBOSE)
    # Add pseudocounts (we have little data)
    mu_neu = np.maximum(mu_neu, 4e-8)

    if plot:
        plot_fits(mu_neu, VERBOSE=VERBOSE, data=data)
        plot_mu_matrix(mu_neu)

        plt.ion()
        plt.show()

    #if VERBOSE >= 1:
    #    print 'Save to file'
    #from hivwholeseq.analysis.filenames import analysis_data_folder
    #fn_out_mu = analysis_data_folder+'mu_neutralclass.pickle'
    #mu_neu.to_pickle(fn_out_mu)


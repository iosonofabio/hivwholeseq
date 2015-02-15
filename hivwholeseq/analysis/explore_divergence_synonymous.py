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
            prot[i] = translate(''.join(codon))

    return prot



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

    degeneracy = get_degeneracy_dict()

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
                    print 'Skip'
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
                # Trim the stop codon if still there (why? why?)
                if protm[-1] == '*':
                    if VERBOSE >= 2:
                        print 'Ends with a stop, trim it'
                    icons = icons[:-3]
                    consm = consm[:-3]
                    protm = protm[:-1]
                    aft = aft[:, :, :-3]
                    coomap = coomap[:, coomap[:, 1] < len(consm)]

                else:
                    continue

            # Get the map as a dictionary from patient to subtype
            coomapd = dict(coomap[:, ::-1])

            # Get only sites that are conserved or with only one site per codon changing
            ind_poly, pos_poly = get_codons_n_polymorphic(aft, icons, n=[0, 1], VERBOSE=VERBOSE)

            # NOTE: all 4-fold degenerate codons have all and only 3rd positions,
            # so they are easy to score, we ignore the rest for now
            prot_deg = np.array(map(degeneracy.get, protm), int)
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
                pos_dna = pos * 3 + 2

                pos_sub = coomapd[pos]
                Ssubpos = Ssub[pos][protm[pos]]

                aft_pos = aft[:, :, pos_dna]

                # Get only non-masked time points
                indpost = -aft_pos[:, 0].mask
                if indpost.sum() == 0:
                    continue
                timespos = times[indpost]
                aft_pos = aft_pos[indpost]

                # Discard of the initial time point is already polymorphic
                aft_anc = aft_pos[:, icons[pos_dna]]
                if aft_anc[0] < 0.9:
                    continue

                # Iterate over the three derived nucleotides
                for inuc, aft_der in enumerate(aft_pos.T[:4]):
                    if inuc == icons[pos_dna]:
                        continue

                    anc = consm[pos_dna]
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


    if plot:

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

        fits = {}

        fig, ax = plt.subplots()
        for mut, datum in datap.iterrows():
            x = np.array(datum.index)
            y = np.array(datum)

            # Get rid of masked stuff (this happens if we miss data in the means)
            ind = -(np.isnan(x) | np.isnan(y))
            x = x[ind]
            y = y[ind]

            color = cm.jet(1.0 * muts.index(mut) / len(muts))
            ax.scatter(x, y, color=color, label=mut)

            # Linear fit
            m = np.dot(y, x) / np.dot(x, x)
            fits[mut] = m
            xfit = np.linspace(0, 3000)
            ax.plot(xfit, m * xfit, color=color, lw=1.5, alpha=0.5)

        ax.set_xlabel('Time')
        ax.set_ylabel('Allele frequency')
        ax.set_ylim(1e-4, 1)
        ax.set_yscale('log')
        ax.grid(True)
        ax.legend(loc='upper left', fontsize=10, ncol=2)

        plt.ion()
        plt.show()

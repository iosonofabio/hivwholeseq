# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/01/15
content:    Quantify positive selection by finding sweeps.
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
from Bio.Seq import translate

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


# Globals



# Functions



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Study accumulation of minor alleles for different kinds of mutations',
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

    regdata = []
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

        regdata.append({'name': region, 'consensus': conssub, 'S': Ssub, 'L': len(conssub)})

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

            # FIXME: deal better with depth (this should be already there?)
            aft[aft < 2e-3] = 0

            # Condition on fixation
            ind_sweep = zip(*((aft[0] < 0.05) & (aft[-1] > 0.95)).T.nonzero())

            for posdna, inuc in ind_sweep:
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
                nuc = alpha[inuc]
                mut = anc+'->'+nuc

                # Ignore indels
                if (inuc >= 4) or (ianc >= 4):
                    continue

                # Skip if the site is already polymorphic at the start
                if aftpos[ianc, 0] < 0.95:
                    continue

                # Define transition/transversion
                if frozenset(nuc+anc) in (frozenset('CT'), frozenset('AG')):
                    trclass = 'ts'
                else:
                    trclass = 'tv'

                # Find out whether it's syn or nonsyn
                codanc = consm[posdna - posdna % 3: posdna - posdna % 3 + 3]
                codder = codanc.copy()
                codder[posdna % 3] = nuc
                if translate(''.join(codanc)) == translate(''.join(codder)):
                    mutclass = 'syn'
                else:
                    mutclass = 'nonsyn'


                # Get the whole trajectory for plots against time
                for af, time in izip(aftpos[inuc], timespos):
                    data.append((region, pcode,
                                 posdna, pos_sub,
                                 anc, nuc, mut,
                                 mutclass, trclass,
                                 Ssubpos,
                                 time, af))

    data = pd.DataFrame(data=data,
                        columns=['region', 'pcode',
                                 'posdna', 'possub',
                                 'anc', 'der', 'mut',
                                 'class', 'tr',
                                 'Ssub',
                                 'time', 'af'])

    regdata = pd.DataFrame(regdata).set_index('name', drop=False)


    # Exclude synonymous ones
    # FIXME: better criterion
    data = data.loc[data['class'] == 'nonsyn']

    if plot:
        fig, ax = plt.subplots()

        # Fit sweep
        from scipy.optimize import curve_fit
        fun = lambda x, s, t0: 1.0 / (1.0 + np.exp(-s * (x - t0)))
        fits = []

        datap = data.groupby(['region', 'pcode', 'possub', 'mut'])
        for (region, pcode, pos_sub, mut), datum in datap:
            datum.sort('time', inplace=True)

            x = np.array(datum['time'])
            y = np.array(datum['af'])

            ind = -(np.isnan(x) | np.isnan(y))
            x = x[ind]
            y = y[ind]

            color = cm.jet(1.0 * pos_sub / regdata.loc[region].loc['L'])

            try:
                tmid = x[(y > 0.5)][0]
                s, t0 = curve_fit(fun, x, y, p0=(1e-3, tmid))[0]

                fits.append({'region': region,
                             'pcode': pcode,
                             'pos_sub': pos_sub,
                             'mut': mut,
                             's': s,
                             't0': t0,
                             'Ssub': datum.iloc[0]['Ssub']
                            })

                label=(', '.join(map(str, (region, pcode, pos_sub, mut)))+
                       ', s = '+'{:2.1e}'.format(s))
                xfit = np.linspace(0, x.max(), 1000)
                yfit = fun(xfit, s, t0)
                ax.plot(xfit, yfit, lw=1.5, color=color, alpha=0.5)
            except RuntimeError:
                label=', '.join(map(str, (region, pcode, pos_sub, mut)))

            ax.plot(x, y, color=color, lw=2, label=label)

        fits = pd.DataFrame(fits)
        
        ax.set_xlabel('Time [days from infection]')
        ax.set_ylabel('Allele frequency')
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True)
        ax.legend(loc='lower right', fontsize=10)


        # Plot the subtype entropy distribution of the fits
        fig, ax = plt.subplots()
        x = np.sort(np.array(fits['Ssub']))
        y = 1.0 - np.linspace(0, 1, len(x))
        Ssubmed = x[np.argmin(np.abs(y - 0.5))]
        ax.plot(x, y, color='k', lw=2)
        ax.grid(True)
        ax.set_xscale('log')
        ax.set_xlim(0.9 * x.min(), 1.1 * x.max())
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel('Entropy in subtype [bits]')
        ax.set_ylabel('Fraction of sweeps with Ssub > x')

        # Plot the fits for the selection coefficients
        fig, ax = plt.subplots()
        sh = np.sort(np.array(fits['s']))
        ax.plot(sh, 1.0 - np.linspace(0, 1, len(sh)), color='k',
                lw=2, label='all sweeps')

        # Plot separately the distributions for low/high entropy sweeps
        fitslow = fits.loc[fits['Ssub'] < 0.6 * Ssubmed]
        shlow = np.sort(np.array(fitslow['s']))
        ax.plot(shlow, 1.0 - np.linspace(0, 1, len(shlow)), color='darkblue',
                lw=2, label='Ssub < '+'{:2.2f}'.format(0.6 * Ssubmed))

        fitshigh = fits.loc[fits['Ssub'] > 1.4 * Ssubmed]
        shhigh = np.sort(np.array(fitshigh['s']))
        ax.plot(shhigh, 1.0 - np.linspace(0, 1, len(shhigh)), color='darkred',
                lw=2, label='Ssub > '+'{:2.2f}'.format(1.4 * Ssubmed))

        ax.grid(True)
        ax.set_xscale('log')
        ax.set_xlim(0.9 * sh.min(), 1.1 * sh.max())
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel('Fitness coefficient')
        ax.set_ylabel('Fraction of sweeps with s > x')
        ax.legend(loc='lower left')

        plt.ion()
        plt.show()

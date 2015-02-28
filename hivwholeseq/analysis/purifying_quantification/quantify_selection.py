# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/01/15
content:    Quantify purifying selection on different subtype entropy classes.
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

            # FIXME: deal better with depth (this should be already there?)
            aft[aft < 2e-3] = 0

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

    # Bin by subtype entropy
    bins_S = np.array([0, 0.03, 0.06, 0.1, 0.25, 0.7, 3])
    binsc_S = 0.5 * (bins_S[1:] + bins_S[:-1])
    data['Sbin'] = 0
    for b in bins_S[1:]:
        data.loc[data.loc[:, 'Ssub'] >= b, 'Sbin'] += 1


    # Fit exponential saturation
    from scipy.optimize import curve_fit
    fun = lambda x, l, u: l * (1 - np.exp(- u/l * x))
    fun_neu = lambda x, u: u * x
    fits = []
    dataf = (data
             .loc[data.loc[:, 'Ssub'] < bins_S[-2]]
             .loc[:, ['region', 'Sbin', 'time', 'af']]
             .groupby(['region', 'Sbin']))
    for (region, iSbin), datum in dataf:
        x = np.array(datum['time'])
        y = np.array(datum['af'])

        ind = -(np.isnan(x) | np.isnan(y))
        x = x[ind]
        y = y[ind]

        try:
            l, u = curve_fit(fun, x, y, p0=(y[-1], 1e-5))[0]
        except RuntimeError:
            continue

        ## If saturation is very high, fit unsaturated
        #if (l > 0.1) or (u <= 0):
        #    try:
        #        l = np.nan
        #        u = np.dot(x, y) / np.dot(x, x)
        #    except RuntimeError:
        #        continue

        fits.append((region, iSbin, l, u))

    fits = pd.DataFrame(data=fits,
                        columns=['region', 'iSbin', 'l', 'u'])
    fits['S'] = binsc_S[fits['iSbin']]
    
    # Estimate fitness cost
    mu = 0.5e-5
    fits['s'] = mu / fits['l']

    if plot:
        for (region, fitsreg) in fits.groupby('region'):

            fig, axs = plt.subplots(1, 2, figsize=(13, 6))
            fig.suptitle(region, fontsize=20)

            # Plot the fits
            ax = axs[0]
            xfit = np.logspace(0, 3.5, 1000)

            for _, fit in fitsreg.iterrows():
                iSbin = fit['iSbin']
                l = fit['l']
                u = fit['u']
                if np.isnan(l):
                    yfit = fun_neu(xfit, u)
                    label = ('S e ['+'{:2.2f}'.format(bins_S[iSbin])+', '+'{:2.2f}'.format(bins_S[iSbin+1])+']'+
                             ', u = '+'{:2.2e}'.format(u))
                else:
                    yfit = fun(xfit, l, u)
                    label = ('S e ['+'{:2.2f}'.format(bins_S[iSbin])+', '+'{:2.2f}'.format(bins_S[iSbin+1])+']'+
                             ', l = '+'{:2.2e}'.format(l)+
                             ', u = '+'{:2.2e}'.format(u))
                
                color = cm.jet(1.0 * iSbin / len(fitsreg))

                ax.plot(xfit, yfit, color=color, label=label, lw=2)

            ax.set_xlabel('Time [days from infection]')
            ax.set_ylabel('Allele frequency')
            ax.legend(loc='lower right', title='Entropy class', fontsize=10)
            ax.text(0.05, 0.9,
                    ('f(x) = l * [1 - e^(-ux/l)]'),
                    horizontalalignment='left',
                    verticalalignment='center',
                    transform=ax.transAxes)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.grid(True)


            ## Plot the value at some reasonable time
            #x0 = 2000
            #fitsreg['y0'] = fun(x0, fitsreg['l'], fitsreg['u'])
            #ax2 = axs[1]
            #ax2.plot(fitsreg['S'], fitsreg['y0'], lw=2, c='k')

            #ax2.set_xlabel('Entropy in subtype [bits]')
            #ax2.set_ylabel('Fit value at x0 = '+str(x0))
            #ax2.set_xscale('log')
            #ax2.set_yscale('log')
            #ax2.grid(True)


            # Plot the estimated fitness value
            ax3 = axs[1]
            ax3.plot(fitsreg['S'], fitsreg['s'], lw=2, c='k')
            ax3.set_xlabel('Entropy in subtype [bits]')
            ax3.set_ylabel('Fitness cost')
            ax3.set_ylim(5e-5, 1)
            ax3.set_xscale('log')
            ax3.set_yscale('log')
            ax3.grid(True)
            ax3.text(0.05, 0.12,
                     ('s := mu / l, with mu = 5 x 10^-6'),
                     horizontalalignment='left',
                     verticalalignment='center',
                     transform=ax3.transAxes)


            plt.tight_layout(rect=(0, 0, 1, 0.96))

        plt.ion()
        plt.show()


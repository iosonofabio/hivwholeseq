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
from collections import defaultdict
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.sequence import translate_with_gaps
from hivwholeseq.multipatient.explore_codon_usage_patsubtype import get_degenerate_dict
import hivwholeseq.utils.plot
from hivwholeseq.multipatient.explore_entropy_patsubtype import (
    get_subtype_reference_alignment, get_ali_entropy)


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

    degenerate_dict = get_degenerate_dict()

    data = defaultdict(list)

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for region in regions:
        if VERBOSE >= 1:
            print region

        # FIXME
#        if VERBOSE >= 2:
#            print 'Get subtype B reference sequence alignment to HXB2 of', region
#        ali_sub = get_subtype_reference_alignment(region, VERBOSE=VERBOSE)
#
#        if VERBOSE >= 2:
#            print 'Get entropy in subtype alignment'
#
        #Ssub = get_ali_entropy(ali_sub, VERBOSE=VERBOSE)
        Ssub = np.zeros(10000)

        for ipat, (pname, patient) in enumerate(patients.iterrows()):

            patient = Patient(patient)
            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=1000,
                                                                 depth_min=300,
                                                                 VERBOSE=VERBOSE)
            if len(ind) == 0:
                continue

            coomap = patient.get_map_coordinates_reference(region)
            # Shift HXB2 coordinates to alignment (region)
            # FIXME: this assumes we have the first position, should be fine for most cases
            coomap[:, 0] -= coomap[0, 0]
            coomap = dict(coomap[:, ::-1])

            times = patient.times[ind]

            icons = patient.get_initial_consensus_noinsertions(aft, VERBOSE=VERBOSE)
            cons = alpha[icons]
            
            if ('N' in cons) or ('-' in cons):
                continue

            prot = translate_with_gaps(cons)

            # Trim the stop codon if still there (why? why?)
            if prot[-1] == '*':
                icons = icons[:-3]
                cons = cons[:-3]
                prot = prot[:-1]
                aft = aft[:, :, :-3]

            # Premature stops in the initial consensus???
            if '*' in prot:
                continue

            prot_deg = np.array(map(degenerate_dict.get, prot), int)
            # NOTE: all 4-fold degenerate codons have all and only 3rd positions
            ind_4fold = (prot_deg == 4).nonzero()[0]
            ind_4fold_3rd = ind_4fold * 3 + 2

            # Exclude sites in which the other two sites change (sweeps or so)
            ind_good = ((aft[:, icons[ind_4fold_3rd - 2], ind_4fold_3rd - 2].min(axis=0) > 0.95) &
                        (aft[:, icons[ind_4fold_3rd - 1], ind_4fold_3rd - 1].min(axis=0) > 0.95))
            ind_4fold = ind_4fold[ind_good]
            ind_4fold_3rd = ind_4fold * 3 + 2

            aft4 = aft[:, :, ind_4fold_3rd]
            # FIXME: deal better with depth (should be there already, I guess)
            aft4[aft4 < 3e-3] = 0

            # Stratify by transitions/transversions
            for pos, aft4pos in enumerate(aft4.swapaxes(0, 2)):

                # Get only non-masked time points
                indpost = -aft4pos[0].mask
                if indpost.sum() == 0:
                    continue
                timespos = times[indpost]
                aft4pos = aft4pos[:, indpost]

                posdna = ind_4fold_3rd[pos]
                anc = cons[posdna]
                ianc = icons[posdna]

                # Skip if the site is already polymorphic at the start
                if aft4pos[ianc, 0] < 0.95:
                    continue

                Ssubpos = Ssub[coomap[posdna]]

                afkey = {'ts': [], 'tv': []}
                for inuc, af in enumerate(aft4pos[:4]):
                    nuc = alpha[inuc]
                    if nuc == anc:
                        continue

                    if frozenset(nuc+anc) in (frozenset('CT'), frozenset('AG')):
                        key = 'ts'
                    else:
                        key = 'tv'
                    afkey[key].append(af)

                # Get time points one by one
                for it, time in enumerate(timespos):
                    afts = afkey['ts'][0][it]
                    aftv = sorted([afkey['tv'][0][it], afkey['tv'][1][it]], reverse=True)

                    # Only take data points for which we see at least one mutations
                    if any(map(lambda x: x > 0, (afts, aftv[0]))):
                        data['af'].append((region, pname, time, anc, Ssubpos,
                                           afts, aftv[0], aftv[1]))

                # Get the whole trajectory for plots against time
                data['aft'].extend([(region, pname, time, anc, Ssubpos,
                                     afkey['ts'][0][it],
                                     afkey['tv'][0][it],
                                     afkey['tv'][1][it],
                                     )
                                    for it, time in enumerate(timespos)])


    data['af'] = pd.DataFrame(data=data['af'],
                              columns=['Region', 'Patient', 'Time', 'Anc', 'S subtype',
                                       'ts', 'tv1', 'tv2'])
    data['aft'] = pd.DataFrame(data=data['aft'],
                              columns=['Region', 'Patient', 'Time', 'Anc', 'S subtype',
                                       'ts', 'tv1', 'tv2'])

    # Print some basic output
    ts = np.array(data['af']['ts'])
    tv = np.array(data['af']['tv1'])
    print 'Only transition:', '{:1.0%}'.format(((ts > 0) & (tv == 0)).mean())
    print 'Only transversion:', '{:1.0%}'.format(((tv > 0) & (ts == 0)).mean())
    print 'Both:', '{:1.0%}'.format(((tv > 0) & (ts > 0)).mean())

    if plot:
        if pnames is not None:
            title = ' + '.join(pnames)+', '+' + '.join(regions)
        else:
            title = ' + '.join(regions)

        #fig, ax = plt.subplots()
        #ax.set_title(title)
        #ax.scatter(tv, ts)
        #ax.plot([1e-3, 1], [1e-3, 1], lw=1, c='k', ls='--')

        #indfit = (tv < 1e-2) & (ts < 1e-2) & (tv > 2e-3) & (ts > 2e-3)
        #tvfit = tv[indfit]
        #tsfit = ts[indfit]
        #m = np.inner(tvfit, tsfit) / np.inner(tvfit, tvfit)
        #ax.plot([1e-3, 1], [m * 1e-3, m], lw=1.5, c='k',
        #        label='Fit: m = '+'{:1.1f}'.format(m))


        #ax.grid(True)
        #ax.set_xlabel('Frequency highest transversion')
        #ax.set_ylabel('Frequency transition')
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        #ax.set_xlim(1e-3, 1)
        #ax.set_ylim(1e-3, 1)
        #plt.tight_layout()


        ts = np.array(data['aft']['ts'])
        tv = np.concatenate([np.array(data['aft']['tv1']),
                             np.array(data['aft']['tv2'])])
        times_ts = np.array(data['aft']['Time'])
        times_tv = np.tile(times_ts, 2)
        xs = {'ts': times_ts, 'tv': times_tv}
        ys = {'ts': ts, 'tv': tv}
        colors = {'ts': 'b', 'tv': 'g'}

        ## Plot by scatter
        #fig, ax = plt.subplots()
        #ax.set_title(title)
        #for key in ('ts', 'tv'):
        #    x = xs[key]
        #    y = ys[key]
        #    ax.scatter(x, y, c=colors[key])

        #    # Linear LS fit
        #    m = np.inner(x, y) / np.inner(x, x)
        #    ax.plot([x.min(), x.max()], [m * x.min(), m * x.max()],
        #            color=colors[key], label=key+' '+'{:1.1e}'.format(m)+' changes/day')

        #ax.legend(loc=2, fontsize=10)
        #ax.set_xlabel('Time [days from infection]')
        #ax.set_ylabel('Allele frequency')
        #ax.set_yscale('log')
        #ax.set_ylim(1e-3, 1)
        #plt.tight_layout()

        # Plot the mean only
        fig, ax = plt.subplots()
        ax.set_title(title)
        for key in ('ts', 'tv'):
            x = xs[key]
            y = ys[key]
            xunique = sorted(set(x))
            ymean = [y[x == xuniquei].mean() for xuniquei in xunique]
            ax.plot(xunique, ymean, lw=2, c=colors[key])

            # Linear LS fit
            m = np.inner(x, y) / np.inner(x, x)
            xfit = np.linspace(x.min(), x.max(), 1000)
            ax.plot(xfit, m * xfit,
                    ls='--',
                    color=colors[key], label=key+' '+'{:1.1e}'.format(m)+' changes/day')

        ax.grid(True)
        ax.legend(loc=2, fontsize=10)
        ax.set_xlabel('Time [days from infection]')
        ax.set_ylabel('Allele frequency')
        ax.set_yscale('log')
        ax.set_ylim(1e-4, 1)
        plt.tight_layout()


        ## Plot as a function of time, divided by subtype entropy
        #Ssub_ts = np.array(data['aft']['S subtype'])
        #Ssub_tv = np.tile(Ssub_ts, 2)
        #Ssubs = {'ts': Ssub_ts, 'tv': Ssub_tv}

        #fig, ax = plt.subplots()
        #ax.set_title(title)
        #Ssubbins = [[0, 0.05], [0.2, 2]]
        #markers = ['o', 'd']
        #for key in ('ts', 'tv'):
        #    xk = xs[key]
        #    yk = ys[key]
        #    for ibin, Ssubbin in enumerate(Ssubbins):
        #        indSsub = (Ssubs[key] > Ssubbin[0]) & (Ssubs[key] <= Ssubbin[1])
        #        x = xk[indSsub]
        #        y = yk[indSsub]
        #        ax.scatter(x, y, c=colors[key], marker=markers[ibin])

        #        # Linear LS fit
        #        m = np.inner(x, y) / np.inner(x, x)
        #        ax.plot([x.min(), x.max()], [m * x.min(), m * x.max()],
        #                color=colors[key],
        #                label=key+' '+'{:1.1e}'.format(m)+' changes/day, Ssub e '+str(Ssubbin))

        #ax.legend(loc=2, fontsize=10)
        #ax.set_xlabel('Time [days from infection]')
        #ax.set_ylabel('Allele frequency')
        #ax.set_yscale('log')
        #ax.set_ylim(1e-3, 1)
        #plt.tight_layout()



        plt.ion()
        plt.show()

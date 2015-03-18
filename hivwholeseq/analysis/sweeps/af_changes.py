# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/01/15
content:    Study allele frequency changes and classes of mutations.
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
from hivwholeseq.utils.sequence import get_coordinates_genomic_region



# Globals
pnames = ['20097', '15319', '9669', '15363', '15823', '15241', '15376']
regions = ['p17', 'p24', 'vif', 'nef', 'PR', 'RT', 'IN']



# Functions
def add_distance_from_sweeps(data):
    pos_sweeps = np.unique(data.loc[data['sweep'] == True, 'possubgw'])
    ds = 1e4 * np.ones(len(data))
    for pos in pos_sweeps:
        tmp = np.abs(data['possubgw'] - pos)
        ds = np.minimum(ds, tmp)
    data['distance'] = ds


def collect_data(pnames, regions, VERBOSE=0, plot=False):
    '''Collect data to study allele freqs around sweeps'''
    from hivwholeseq.reference import load_custom_reference

    regdata = []
    data = []

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    refseq = load_custom_reference('HXB2', 'gb')

    for region in regions:
        if VERBOSE >= 1:
            print region

        if VERBOSE >= 2:
            print 'Get subtype consensus (for checks only)'
        conssub = get_subtype_reference_alignment_consensus(region, VERBOSE=VERBOSE)

        if VERBOSE >= 2:
            print 'Get subtype entropy'
        Ssub = get_subtype_reference_alignment_entropy(region, VERBOSE=VERBOSE)

        if VERBOSE >= 2:
            print 'Get reference sequence'
        location = get_coordinates_genomic_region(refseq, region)
        start = location.nofuzzy_start
        end = location.nofuzzy_end

        regdata.append({'name': region, 'consensus': conssub, 'S': Ssub,
                        'location': (start, end), 'L': end - start})

        for ipat, (pname, patient) in enumerate(patients.iterrows()):
            pcode = patient.code
            if VERBOSE >= 2:
                print pname, pcode

            patient = Patient(patient)
            # We can call a sweep with little depth, but we need some depth for
            # the detailed dynamics of the neighbouring sites
            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=1000,
                                                                 depth_min=1,
                                                                 VERBOSE=VERBOSE)
            if len(ind) == 0:
                if VERBOSE >= 2:
                    print 'No time points: skip'
                continue

            times = patient.times[ind]
            ntemp = patient.get_n_templates_roi(region)[ind]

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

            # Condition on fixation
            pos_sweeps, inuc_sweeps = ((aft[0] < 0.05) & (aft[-1] > 0.95)).T.nonzero()
            pos_sweepsl = pos_sweeps.tolist()

            # Get all trajectories and then we'll filter by distance from sweeps
            for posdna, pos_sub in coomapd['pat_to_subtype'].iteritems():
                # Get the entropy
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
                ntpos = ntemp[indpost]
                aftpos = aftpos[:, indpost]

                anc = consm[posdna]
                ianc = icons[posdna]

                for inuc, nuc in enumerate(alpha[:4]):
                    if inuc == ianc:
                        continue

                    mut = anc+'->'+nuc

                    # Ignore indels
                    if (inuc >= 4) or (ianc >= 4):
                        continue

                    # Skip if the site is already polymorphic at the start
                    if aftpos[ianc, 0] < 0.95:
                        continue

                    # Is it a sweep?
                    if (posdna in pos_sweeps) and (inuc == inuc_sweeps[pos_sweepsl.index(posdna)]):
                        is_sweep = True
                    else:
                        is_sweep = False

                    # Define transition/transversion
                    if frozenset(nuc+anc) in (frozenset('CT'), frozenset('AG')):
                        trclass = 'ts'
                    else:
                        trclass = 'tv'

                    # Find out whether it's syn or nonsyn
                    codanc = consm[posdna - posdna % 3: posdna - posdna % 3 + 3]
                    codder = codanc.copy()
                    codder[posdna % 3] = nuc
                    codanc = ''.join(codanc)
                    codder = ''.join(codder)
                    if ('-' in codanc) or ('-' in codder):
                        continue

                    if translate(codanc) == translate(codder):
                        mutclass = 'syn'
                    else:
                        mutclass = 'nonsyn'

                    # Get allele frequency changes
                    daftpos = np.diff(aftpos)
                    dts = np.diff(timespos)
                    for it, (daf, dt, nt) in enumerate(izip(daftpos[inuc], dts, ntpos)):
                        data.append((region, pcode,
                                     posdna, pos_sub, pos_sub + start,
                                     anc, nuc, mut,
                                     codanc, codder,
                                     mutclass, trclass,
                                     is_sweep,
                                     Ssubpos,
                                     nt,
                                     daf, dt,
                                     timespos[it], timespos[it + 1],
                                     aftpos[inuc, it], aftpos[inuc, it + 1]))

    if len(data):
        data = pd.DataFrame(data=data,
                            columns=['region', 'pcode',
                                     'posdna', 'possub', 'possubgw',
                                     'anc', 'der', 'mut',
                                     'codanc', 'codder',
                                     'class', 'tr',
                                     'sweep',
                                     'Ssub',
                                     'n templates',
                                     'daf', 'dt',
                                     'time1', 'time2',
                                     'af1', 'af2'])

    regdata = pd.DataFrame(regdata).set_index('name', drop=False)

    return data, regdata



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Study accumulation of minor alleles for different kinds of mutations',
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

    data, regdata = collect_data(pnames, regions, VERBOSE=VERBOSE, plot=plot)

    # Bin by subtype entropy using quantiles
    from hivwholeseq.analysis.purifying_selection.joint_model import add_Sbins
    bins_S, binsc_S = add_Sbins(data, bins=5, VERBOSE=VERBOSE)
    data['|daf|'] = np.abs(data['daf'])
    data['dafdt'] = data['daf'] / data['dt']
    data['|dafdt|'] = data['|daf|'] / data['dt']


    if plot:

        # Plot distributions of daf for different entropy categories
        bins = np.linspace(-1, 1, 50)
        binsc = 0.5 * (bins[1:] + bins[:-1])
        fig, ax = plt.subplots()
        for iSbin, datum in data.groupby('Sbin'):
            S = binsc_S[iSbin]

            h = np.histogram(datum['daf'], bins=bins, density=True)[0]

            ax.plot(binsc, h, lw=2,
                    color=cm.jet(1.0 * iSbin / len(binsc_S)),
                    label='S = {:.1G}'.format(S),
                   )

        ax.set_xlabel('daf')
        ax.set_ylabel('density')
        ax.set_xlim(bins[0] -0.05, bins[-1] + 0.05)
        ax.set_yscale('log')
        ax.grid(True)
        ax.legend()

        plt.tight_layout()

        # Plot distributions of |daf| for different entropy categories
        bins = np.linspace(0, 1, 30)
        binsc = 0.5 * (bins[1:] + bins[:-1])
        fig, ax = plt.subplots()
        for iSbin, datum in data.groupby('Sbin'):
            S = binsc_S[iSbin]

            h = np.histogram(datum['|daf|'], bins=bins, density=True)[0]

            ax.plot(binsc, h, lw=2,
                    color=cm.jet(1.0 * iSbin / len(binsc_S)),
                    label='S = {:.1G}'.format(S),
                   )

        ax.set_xlabel('|daf|')
        ax.set_ylabel('density')
        ax.set_xlim(bins[0] -0.05, bins[-1] + 0.05)
        ax.set_yscale('log')
        ax.grid(True)
        ax.legend()

        plt.tight_layout()

        # Plot distributions of daf / dt for different entropy categories
        bins = np.linspace(-1, 1, 50) / 200.
        binsc = 0.5 * (bins[1:] + bins[:-1])
        fig, ax = plt.subplots()
        for iSbin, datum in data.groupby('Sbin'):
            Smin = bins_S[iSbin]
            Smax = bins_S[iSbin + 1]

            h = np.histogram(datum['dafdt'], bins=bins, density=True)[0]

            ax.plot(binsc, h, lw=2,
                    color=cm.jet(1.0 * iSbin / len(binsc_S)),
                    label='S = [{:.1G}, {:.1G}]'.format(Smin, Smax),
                   )

        ax.set_xlabel('daf / dt [1 / days]')
        ax.set_ylabel('density')
        ax.set_xlim(1.04 * bins[0], 1.04 * bins[-1])
        ax.set_yscale('log')
        ax.grid(True)
        ax.legend()

        plt.tight_layout()


        # Plot distributions of daf / dt for syn/nonsyn
        bins = np.linspace(-1, 1, 50) / 200.
        binsc = 0.5 * (bins[1:] + bins[:-1])
        fig, ax = plt.subplots()
        for synclass, datum in data.groupby('class'):
            h = np.histogram(datum['dafdt'], bins=bins, density=True)[0]

            ax.plot(binsc, h, lw=2,
                    label=synclass,
                   )

        ax.set_xlabel('daf / dt [1 / days]')
        ax.set_ylabel('density')
        ax.set_xlim(1.04 * bins[0], 1.04 * bins[-1])
        ax.set_yscale('log')
        ax.grid(True)
        ax.legend()

        plt.tight_layout()


        # Plot distributions of daf / dt for sweeps/nonsweeps
        bins = np.linspace(-1, 1, 50) / 200.
        binsc = 0.5 * (bins[1:] + bins[:-1])
        fig, ax = plt.subplots()
        for sweep, datum in data.groupby('sweep'):
            if sweep:
                label = 'sweep'
            else:
                label = 'nonsweep'

            h = np.histogram(datum['dafdt'], bins=bins, density=True)[0]

            ax.plot(binsc, h, lw=2,
                    label=label,
                   )

        ax.set_xlabel('daf / dt [1 / days]')
        ax.set_ylabel('density')
        ax.set_xlim(1.04 * bins[0], 1.04 * bins[-1])
        ax.set_yscale('log')
        ax.grid(True)
        ax.legend()

        plt.tight_layout()


        # Plot distributions of daf for entropy classes and distance to sweeps
        add_distance_from_sweeps(data)
        datanf = data.loc[data['distance'] != 0]
        fig, ax = plt.subplots()
        for iSbin, datum in datanf[datanf['af1'] < 0.05].groupby('Sbin'):
            S = binsc_S[iSbin]
            color = cm.jet(1.0 * iSbin / len(binsc_S))

            x = np.array(datum['distance'])
            y = np.array(datum['dafdt'])

            ax.scatter(x, y, lw=2,
                       color=color,
                       alpha=0.3,
                       label='S = {:.1G}'.format(S),
                       zorder=(10.-iSbin) / 10
                      )

        ax.set_ylabel('daf / dt [1 / days]')
        ax.set_xlabel('distance')
        ax.set_xlim(-1, 80)
        ax.grid(True)
        ax.legend()

        # And add exponential fit
        from scipy.optimize import curve_fit
        fun = lambda x, A, gamma: A * np.exp(-gamma * x)
        fits = []
        for iSbin, datum in datanf[datanf['af1'] < 0.05].groupby('Sbin'):
            x = np.array(datum['distance'])
            y = np.array(datum['dafdt'])
            A, gamma = curve_fit(fun, x, y, p0=(0.00002, 0.1))[0]
            fits.append({'Sbin': iSbin, 'A': A, 'gamma': gamma})

        fits = pd.DataFrame(fits)

        for _, fit in fits.iterrows():
            iSbin = fit['Sbin']
            S = binsc_S[iSbin]
            color = cm.jet(1.0 * iSbin / len(binsc_S))

            xfit = np.linspace(0, 80, 1000)
            yfit = fun(xfit, 0.001 * fit['A'] / fits['A'].max(), fit['gamma'])
            ax.plot(xfit, yfit, lw=2, color=color)


        plt.tight_layout()
 

        # Make histograms
        from hivwholeseq.utils.pandas import add_binned_column
        bins_d, binsc_d = add_binned_column(datanf, 'dbin', 'distance', bins=np.linspace(0, 30, 5))
        fig, ax = plt.subplots()
        for (iSbin, idbin), datum in datanf[datanf['af1'] < 0.05].groupby(('Sbin', 'dbin')):
            color = cm.jet(1.0 * iSbin / len(binsc_S))

            bins = np.linspace(0.01, 1, 50)
            h, bins = np.histogram(datum['daf'], bins=bins, density=True)
            binsc = 0.5 * (bins[:-1] + bins[1:])
            binsw = np.diff(bins)

            h = np.log10(h + 1e-3)
            h = 5 * (h - h.min()) / (h.max() - h.min())
            
            ax.barh(bins[:-1], h, height=binsw, left=binsc_d[idbin],
                    color=color, alpha=0.4)


        # Ask the inverse question, i.e. bin by daf and ask for entropy
        # (This gets rid of the noise)
        bins_daf, binsc_daf = add_binned_column(datanf, 'dafbin', 'dafdt', bins=np.linspace(0, 0.01, 20))
        fig, ax = plt.subplots()
        for (iSbin, idbin), datum in datanf[datanf['af1'] < 0.05].groupby(('Sbin', 'dbin')):
            color = cm.jet(1.0 * iSbin / len(binsc_S))

            bins = np.linspace(0.01, 1, 50)
            h, bins = np.histogram(datum['dafdt'], bins=bins, density=True)
            binsc = 0.5 * (bins[:-1] + bins[1:])
            binsw = np.diff(bins)

            h = np.log10(h + 1e-3)
            h = 5 * (h - h.min()) / (h.max() - h.min())
            
            ax.barh(bins[:-1], h, height=binsw, left=binsc_d[idbin],
                    color=color, alpha=0.4)


        # Hierarchical indexing by dbin, daf, and Sbin
        bins_d, binsc_d = add_binned_column(datanf, 'dbin', 'distance', bins=np.linspace(1, 11, 11))
        bins_daf, binsc_daf = add_binned_column(datanf, 'dafbin', 'daf', bins=np.logspace(-3, -1, 4))
        t = datanf.groupby(['dbin', 'dafbin', 'Sbin']).size().unstack('Sbin').fillna(0).T
        t = (1.0 * t / t.sum())
        t.set_index(binsc_S[t.index], inplace=True)
        t.index.name = 'S'
        t = t.T
       
        fig, ax = plt.subplots()
        ax.imshow(t, interpolation='nearest', aspect='auto')
        ax.set_yticks(np.arange(0, len(binsc_d)) * len(binsc_daf))
        ax.set_yticklabels(map(lambda x: str(int(x)), binsc_d))
        ax.set_xticks(np.arange(len(binsc_S)))
        ax.set_xticklabels(map('{:.1G}'.format, binsc_S))
        ax.set_xlabel('S')
        ax.set_ylabel('distance [bp] -> daf')

        plt.tight_layout()


        # Hierarchical indexing by daf, dbin, and Sbin
        bins_d, binsc_d = add_binned_column(datanf, 'dbin', 'distance', bins=np.arange(1, 21, 2))
        bins_daf, binsc_daf = add_binned_column(datanf, 'dafbin', 'daf', bins=np.logspace(-3, -1, 4))
        t = datanf.groupby(['dafbin', 'dbin', 'Sbin']).size().unstack('Sbin').fillna(0).T
        t = (1.0 * t / t.sum())
        t.set_index(binsc_S[t.index], inplace=True)
        t.index.name = 'S'
        t = t.T
       
        fig, ax = plt.subplots()
        ax.imshow(t, interpolation='nearest', aspect='auto')
        ax.set_yticks(np.arange(0, len(binsc_daf)) * len(binsc_d))
        ax.set_yticklabels(map('{:.1G}'.format, binsc_daf))
        ax.set_xticks(np.arange(len(binsc_S)))
        ax.set_xticklabels(map('{:.1G}'.format, binsc_S))
        ax.set_xlabel('S')
        ax.set_ylabel('daf -> distance [bp]')

        plt.tight_layout()



        # Plot histogram of distances
        fig, ax = plt.subplots()
        h =  np.bincount(data
                         .loc[:, ['posdna', 'pcode', 'region', 'distance']]
                         .groupby(('posdna', 'pcode', 'region'))
                         .mean()
                         .loc[:, 'distance']
                         .astype(int))
        ax.bar(np.arange(len(h)) - 0.5, h, width=1, facecolor='steelblue')

        ax.set_xlabel('Distance from next sweep [bp]')
        ax.set_ylabel('# alleles')
        ax.grid(True)

        plt.tight_layout()



        plt.ion()
        plt.show()

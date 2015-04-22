# vim: fdm=indent
'''
author:     Fabio Zanini
date:       08/04/15
content:    Compare the divergence of the consensus and the population.
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

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.cross_sectional.get_subtype_consensus import (
    get_subtype_reference_alignment_consensus)
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment_entropy)
from hivwholeseq.cross_sectional.get_subtype_allele_frequencies import (
    get_subtype_reference_alignment_allele_frequencies)
from hivwholeseq.patients.one_site_statistics import get_codons_n_polymorphic



# Globals
pnames = ['20097', '15363', '15823', '9669', '20529', '15241', '15376', '15319']
regions = ['p17', 'p24', 'p6', 'p7',
           'PR', 'RT', 'p15', 'IN',
           'vif', 'vpu', 'vpr', 'nef',
           'gp41', 'gp1201']



# Functions
def classify_syn(consm, posdna, der):
    '''Classify a mutation syn/nonsyn'''
    from hivwholeseq.utils.sequence import translate_with_gaps

    pos = posdna // 3
    pospoly = posdna % 3
    anccod = consm[pos * 3: (pos + 1) * 3]
    mutcod = anccod.copy()
    mutcod[pospoly] = der

    anccod = ''.join(anccod)
    mutcod = ''.join(mutcod)

    if ('-' in anccod) or ('-' in mutcod):
        raise ValueError('Gap found in codon')

    ancaa = translate_with_gaps(anccod)
    mutaa = translate_with_gaps(mutcod)

    # Define the type of mutations (syn, simliar aa, div aa)
    if ancaa == mutaa:
        mutclass = 'syn'
    else:
        mutclass = 'nonsyn'

    return mutclass


def collect_data(pnames, regions, VERBOSE=0):
    '''Collect data for the SFS'''
    data = []

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

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
                                                                 depth_min=30,
                                                                 VERBOSE=VERBOSE)
            if len(ind) == 0:
                if VERBOSE >= 2:
                    print 'No time points found: skip.'
                continue

            times = patient.times[ind]

            if VERBOSE >= 2:
                print 'Get coordinate map'
            coomap = patient.get_map_coordinates_reference(region, refname=('HXB2', region))

            icons = patient.get_initial_consensus_noinsertions(aft, VERBOSE=VERBOSE,
                                                               return_ind=True)
            consm = alpha[icons]

            if VERBOSE >= 2:
                print 'Collecting alleles'

            # Iterate over times
            for ii, it in enumerate(ind):
                time = times[ii]
                afs = aft[ii]
                
                # Get mask for allele frequencies
                mask = afs.mask.any(axis=0)

                # Get consensus from that time point
                consmt = alpha[afs.argmax(axis=0)]

                # Look site by site
                div = defaultdict(list)

                # Iterate over positions to distinguish syn/nonsyn
                for posdna in xrange(aft.shape[-1]):
                    pos = posdna // 3
                    pospoly = posdna % 3

                    if VERBOSE >= 3:
                        print posdna, pos

                    # Is this allele masked at this time?
                    if mask[posdna]:
                        continue

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

                        # Classify syn/nonsyn
                        try:
                            mutclass = classify_syn(consm, posdna, der)
                        except ValueError:
                            continue

                        div[('divergence', 'population', 'all')].append(afs[ider, posdna])
                        div[('divergence', 'consensus', 'all')].append(der == consmt[posdna])
                        div[('divergence', 'population', mutclass)].append(afs[ider, posdna])
                        div[('divergence', 'consensus', mutclass)].append(der == consmt[posdna])

                        # Diversity: <sum(i e ACGT) x_i(1-x_i)>
                        # averaged over positions
                        tmp = afs[ider, posdna] * (1.0 - afs[ider, posdna])
                        div[('diversity', 'population', 'all')].append(tmp)
                        div[('diversity', 'population', mutclass)].append(tmp)


                for (obs, ctype, cl), tmp in div.iteritems():
                    # NOTE: We count every derived allele, so each site contributes
                    # with three alleles. Divergence and diversity, however, are
                    # defined per site, so we sum over the three alleles
                    divtmp = np.mean(tmp) * 3
                    datum = {'patient': patient.code,
                             'region': region,
                             'time': time,
                             't0': patient.times[0],
                             'ctype': ctype,
                             'div': divtmp,
                             'class': cl,
                             'obs': obs,
                            }
                    data.append(datum)

    data = pd.DataFrame(data)

    return data


def plot_data(data, VERBOSE=0):
    '''Plot divergence'''
    # NOTE: gp41 has overlaps with tat/rev (hard to tell syn/nonsyn)
    # NOTE: gp120 has the V loops that are hard to call syn/nonsyn
    reg_groups = [['PR', 'IN', 'p15', 'RT'],
                  ['p17', 'p24', 'p6', 'p7'],
                  ['vif', 'vpu', 'vpr', 'nef'],
                  ['gp1201'],
                 ]

    pcodes = np.unique(data['patient']).tolist()
    pcodes.sort(key=lambda x: int(x[1:]))

    cls = np.unique(data['class'])

    fig, axg = plt.subplots(len(reg_groups), len(cls),
                            figsize=(2 + 4 * len(cls), 4 * len(reg_groups)),
                            sharey=True, sharex=True)

    fig.text(0.41, 0.035, 'Time [days from infection]', fontsize=16)
    fig.text(0.035, 0.6, 'Divergence [changes per site]', rotation=90, ha='center', fontsize=16)

    if len(cls) == 2:
        fig.text(0.32, 0.967, 'nonsyn', ha='center', fontsize=16)
        fig.text(0.75, 0.967, 'syn', ha='center', fontsize=16)
    else:
        fig.text(0.22, 0.967, 'all', ha='center', fontsize=16)
        fig.text(0.53, 0.967, 'nonsyn', ha='center', fontsize=16)
        fig.text(0.82, 0.967, 'syn', ha='center', fontsize=16)


    # Bin data in time
    from hivwholeseq.utils.pandas import add_binned_column
    _, times = add_binned_column(data, 'tbin', 'time',
                                 bins=np.linspace(0, 3300, 8),
                                 clip=True)
    data.loc[:, 'tbinned'] = times[data['tbin']]

    for irow, (axs, reg_group) in enumerate(izip(axg, reg_groups)):
        axs[-1].set_ylabel(', '.join(reg_group), rotation=270, labelpad=27,
                           fontsize=16)
        axs[-1].yaxis.set_label_position("right")
        for ax in axs:
            ax.grid(True)

        def plot_single(obs, ctype, cl, pcode, datum, groupby='time'):
            '''Plot single curve'''
            if obs == 'diversity':
                ls = '-'
                dashes = [8, 4, 2, 4, 2, 4]
            elif ctype == 'consensus':
                ls = '-'
                dashes = []
            else:
                ls = '--'
                dashes = [8, 6]

            if cl == 'all':
                ax = axs[-3]
            elif cl == 'nonsyn':
                ax = axs[-2]
            else:
                ax = axs[-1]

            if (cl == 'nonsyn') and (ctype == 'consensus'):
                label = pcode
            else:
                label = ''

            if pcode in pcodes:
                color = cm.jet(1.0 * pcodes.index(pcode) / len(pcodes))
                alpha = 0.2
            else:
                color = 'k'
                alpha = 1.0

            datump = datum.groupby(groupby, as_index=False).mean()
            x = datump[groupby]
            y = datump['div']

            ax.plot(x, y,
                    ls=ls,
                    dashes=dashes,
                    lw=2,
                    color=color,
                    alpha=alpha,
                    label=label,
                  )
        
        # Plot single patients
        datap = (data
                 .loc[data['region'].isin(reg_group)]
                 .groupby(['obs', 'ctype', 'class', 'patient']))
        for (obs, ctype, cl, pcode), datum in datap:
            plot_single(obs, ctype, cl, pcode, datum)

        # Plot average
        datap = (data
                 .loc[data['region'].isin(reg_group)]
                 .groupby(['obs', 'ctype', 'class']))
        for (obs, ctype, cl), datum in datap:
            plot_single(obs, ctype, cl, 'avg', datum, groupby='tbinned')

        if irow == 0:
            axs[-2].legend(loc=2, ncol=2, fontsize=14, title='Patients:')

        
    plt.tight_layout(rect=(0.05, 0.05, 0.98, 0.97))



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Divergence of consensi/population',
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


    data = collect_data(pnames, regions, VERBOSE=VERBOSE)


    if plot:
        tmp = plot_data(data, VERBOSE=VERBOSE)

        plt.ion()
        plt.show()

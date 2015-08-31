# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/01/15
content:    Plot bubbles for draft potential.
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
pnames = ['p1']
regions = ['p17']


# Functions
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

        #FIXME: use mocks for now
        #if VERBOSE >= 2:
        #    print 'Get subtype consensus (for checks only)'
        #conssub = get_subtype_reference_alignment_consensus(region, VERBOSE=VERBOSE)

        #if VERBOSE >= 2:
        #    print 'Get subtype entropy'
        #Ssub = get_subtype_reference_alignment_entropy(region, VERBOSE=VERBOSE)
        conssub = ''
        Ssub = np.zeros(1e4)

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

                    # Get the whole trajectory for plots against time
                    for af, time, nt in izip(aftpos[inuc], timespos, ntpos):
                        data.append((region, pcode,
                                     posdna, pos_sub,
                                     anc, nuc, mut,
                                     codanc, codder,
                                     mutclass, trclass,
                                     is_sweep,
                                     Ssubpos,
                                     nt,
                                     time, af))

    if len(data):
        data = pd.DataFrame(data=data,
                            columns=['region', 'pcode',
                                     'posdna', 'possub',
                                     'anc', 'der', 'mut',
                                     'codanc', 'codder',
                                     'class', 'tr',
                                     'sweep',
                                     'Ssub',
                                     'n templates',
                                     'time', 'af'])

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

    if plot:

        for pcode, datap in data.groupby('pcode'):
            fig, axs = plt.subplots(2, 1,
                                    sharex=True,
                                    figsize=(10, 7),
                                    gridspec_kw={'height_ratios':[4, 1]})
            ax = axs[0]

            tmax = datap['time'].max()
            posmin = regdata['location'].min()[0]
            posmax = regdata['location'].max()[1]
            L = posmax - posmin

            has_label = set()

            for region, datar in datap.groupby('region'):
                dmax = 80
                dex = 1.0 / 300
                pos0 = regdata.loc[region, 'location'][0]

                for (pos, der), datum in datar.groupby(['possub', 'der']):
                    x = np.array(datum['time'])
                    y = np.array(datum['af'])

                    if (y.max() < 0.1).all():
                        continue

                    if datum['sweep'].any():
                        color = 'steelblue'
                        alpha = 0.8
                        zorder = 5
                    else:
                        color = 'darkred'
                        alpha = 0.3
                        zorder = 3

                    xm = 0.5 * (x[1:] + x[:-1])
                    yd = np.diff(y) / np.diff(x)

                    w = dmax * np.abs(yd) / dex

                    ax.fill_betweenx(xm, pos0 + pos - w, pos0 + pos + w,
                                     edgecolor='none',
                                     facecolor=color,
                                     alpha=alpha,
                                     zorder=zorder,
                                    )

            # Proxy artists for legend
            from matplotlib.patches import Rectangle
            p1 = Rectangle((0, 0), 1, 1, fc='steelblue', alpha=0.8)
            p2 = Rectangle((0, 0), 1, 1, fc='darkred', alpha=0.5)
            ax.legend([p1, p2], ('fixed', 'nonfixed'),
                      loc='upper center', ncol=2)
            ax.set_ylabel('Time from infection [days]')
            ax.set_ylim(1.04 * tmax, -10)
            ax.set_title(pcode)
            ax.grid(True)

            ax = axs[1]
            ax.set_ylim(-4, 26)
            ax.set_xlabel('Position in HXB2 [bp]')
            ax.set_xlim(posmin - 0.04 * L, posmax + 0.04 * L)
            ax.set_yticks([])
            ax.grid(True, axis='x')

            from hivwholeseq.analysis.substitution_rate.rate_sliding_window import plot_region_boxes
            plot_region_boxes(regions, ax, VERBOSE=VERBOSE)
            plt.tight_layout(rect=(0, 0, 0.98, 1), h_pad=0.001)


        plt.ion()
        plt.show()

# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/01/15
content:    Quantify purifying selection on different subtype entropy classes.
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

from hivwholeseq.utils.miseq import alpha, alphal
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
from hivwholeseq.analysis.purifying_selection.filenames import get_fitness_cost_entropy_filename


# Globals



# Functions
def get_map_coordinates_reference(region, alim, refname='HXB2', VERBOSE=0):
    '''Get coordinate map between reference and a haplotype alignment'''
    from hivwholeseq.reference import load_custom_reference
    refseq = load_custom_reference(refname, 'gb')
    for feature in refseq.features:
        if feature.id == region:
            break

    regseq = feature.extract(refseq)

    from hivwholeseq.utils.sequence import get_consensus_from_MSA, align_pairwise
    consm = get_consensus_from_MSA(alim, alpha=alpha)

    # Align, masking/unmasking gaps as N
    consm[consm == '-'] = 'N'
    alipw = np.array(align_pairwise(consm, regseq, method='global', score_gapopen=-20))
    alipw[0, alipw[0] == 'N'] = '-'

    # Make map
    coomap = ((alipw != '-').cumsum(axis=1) - 1)[:, (alipw != '-').all(axis=0)]

    return coomap.T



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Study load by deleterious mutations',
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

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    data = []
    for region in regions:
        if VERBOSE >= 1:
            print region

        if VERBOSE >= 2:
            print 'Get subtype consensus'
        conssub = get_subtype_reference_alignment_consensus(region, VERBOSE=VERBOSE)

        if VERBOSE >= 2:
            print 'Get subtype entropy'
        Ssub = get_subtype_reference_alignment_entropy(region, VERBOSE=VERBOSE)

        if VERBOSE >= 2:
            print 'Get fitness costs'
        fn = get_fitness_cost_entropy_filename(region)
        fitness_cost = pd.read_pickle(fn)

        for ipat, (pname, patient) in enumerate(patients.iterrows()):
            pcode = patient.code
            if VERBOSE >= 2:
                print pname, pcode

            patient = Patient(patient)

            hct, ind, alim = patient.get_haplotype_count_trajectory(region)
            if len(ind) == 0:
                if VERBOSE >= 2:
                    print 'No time points: skip'
                continue

            times = patient.times[ind]

            if VERBOSE >= 2:
                print 'Get coordinate map'
            coomap = get_map_coordinates_reference(region, alim, VERBOSE=VERBOSE)

            consm = alim[hct[0].argmax()]
            iconsm = np.array(map(alphal.index, consm), int)
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


            if VERBOSE >= 2:
                print 'Calculating fitness'
            fitnesses = np.zeros(alim.shape[0])
            for iseq, seq in enumerate(alim):
                for pos, nuc in enumerate(seq):
                    # NOTE: Missing positions?
                    if pos not in coomapd['pat_to_subtype']:
                        continue

                    pos_sub = coomapd['pat_to_subtype'][pos]
                    nuc_sub = conssub[pos_sub]

                    if nuc != nuc_sub:
                        Spos = Ssub[pos_sub]
                        Sclass = (fitness_cost
                                  .loc[((fitness_cost['Smin'] < Spos) &
                                        (fitness_cost['Smax'] >= Spos))])

                        # If S is out of bounds, take the boundary value
                        if not len(Sclass):
                            if Spos < fitness_cost['Smin'].min():
                                s = fitness_cost['s'][fitness_cost['S'].argmin()]
                            else:
                                s = fitness_cost['s'][fitness_cost['S'].argmax()]

                        else:
                            s = Sclass.iloc[0].loc['s']

                        fitnesses[iseq] -= s


            data.append({'region': region, 'pcode': pcode, 'hct': hct,
                         'times': times, 'seqs': alim,
                         'fitness': fitnesses,
                        })

    data = pd.DataFrame(data)

    if plot: 
        if VERBOSE >= 1:
            print 'Plot'
        for _, datum in data.iterrows():
            fig, ax = plt.subplots()

            region = datum['region']
            pcode = datum['pcode']
            times = datum['times']
            hct = datum['hct']
            fitnesses = datum['fitness']

            for it, time in enumerate(times):
                hc = hct[it]
                # Get the absolute value for plotting with log scale
                fit_t = -np.concatenate([np.repeat(ftn_seq, hc_seq)
                                         for ftn_seq, hc_seq in izip(fitnesses, hc)])
                fit_t.sort()

                x = fit_t
                y = 1.0 - np.linspace(0, 1, len(x))
                color = cm.jet(1.0 * it / len(times))
                ax.plot(x, y, lw=2, color=color, label=str(int(time))+' days')

            ax.set_xlabel('Fitness cost')
            ax.set_ylabel('Fraction of haplotypes with fitness > x')
            ax.set_title(', '.join([region, pcode]))
            ax.grid(True)
            ax.set_ylim(-0.05, 1.05)
            ax.set_xlim(0, 0.2)

            plt.tight_layout()

        plt.ion()
        plt.show()


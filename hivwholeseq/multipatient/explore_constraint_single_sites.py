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
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.sequence_utils import translate_with_gaps
from hivwholeseq.multipatient.explore_codon_usage_patsubtype import get_degenerate_dict
import hivwholeseq.plot_utils
from hivwholeseq.multipatient.explore_entropy_patsubtype import (
    get_subtype_reference_alignment, get_ali_entropy)
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment_entropy)
from hivwholeseq.patients.one_site_statistics import get_codons_n_polymorphic



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Explore conservation in patients and subtype, site by site',
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

    data = defaultdict(list)

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for region in regions:
        if VERBOSE >= 1:
            print region

        try:
            Ssub = get_subtype_reference_alignment_entropy(region)
        except IOError:
            # FIXME: this is a hack
            Ssub = defaultdict(int)

        for ipat, (pname, patient) in enumerate(patients.iterrows()):

            patient = Patient(patient)
            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=1000,
                                                                 depth_min=300,
                                                                 VERBOSE=VERBOSE)
            if len(ind) == 0:
                continue

            # Get the coordinate map to HXB2 (and ref alignment)
            coomap = patient.get_map_coordinates_reference(region)
            coomap[:, 0] -= coomap[coomap[:, 1] == 0, 0][0]

            # NOTE: Pavel has sometimes cut the second half of the alignment (e.g. RT)
            if len(Ssub) < (coomap[:, 0]).max() + 1:
                coomap = coomap[coomap[:, 0] < len(Ssub)]
                aft = aft[:, :, :coomap[:, 1].max() + 1]

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

            ind_poly, pos_poly = get_codons_n_polymorphic(aft, icons, n=[0, 1], VERBOSE=VERBOSE)

            # FIXME: deal better with depth (this should be already there?)
            aft[aft < 2e-3] = 0

            for pos, pospoly in izip(ind_poly, pos_poly):

                # Forget the polymorphic positions, try all 3
                for pospoly  in xrange(3):
                
                    posdna = pos * 3 + pospoly
                    aftpos = aft[:, :, posdna].T


                    # Only keep positions that are also in HXB2 (so we get the subtype entropy)
                    if posdna not in coomap:
                        continue

                    posdnaref = coomap[posdna]
                    Ssubpos = Ssub[posdnaref]

                    # Get only non-masked time points
                    indpost = -aftpos[0].mask
                    if indpost.sum() == 0:
                        continue
                    timespos = times[indpost]
                    aftpos = aftpos[:, indpost]

                    anc = cons[posdna]
                    ianc = icons[posdna]

                    # Skip if the site is already polymorphic at the start
                    if aftpos[ianc, 0] < 0.95:
                        continue

                    for inuc, af in enumerate(aftpos[:4]):
                        nuc = alpha[inuc]
                        if nuc == anc:
                            continue

                        anccod = cons[pos * 3: (pos + 1) * 3]
                        mutcod = anccod.copy()
                        mutcod[pospoly] = nuc

                        anccod = ''.join(anccod)
                        mutcod = ''.join(mutcod)

                        ancaa = translate_with_gaps(anccod)
                        mutaa = translate_with_gaps(mutcod)

                        # Define the type of mutations (syn, simliar aa, div aa)
                        if ancaa == mutaa:
                            mutclass = 'syn'
                        else:
                            mutclass = 'nonsyn'

                        # Define transition/transversion
                        if frozenset(nuc+anc) in (frozenset('CT'), frozenset('AG')):
                            trclass = 'ts'
                        else:
                            trclass = 'tv'

                        mut = anc+'->'+nuc

                        # Get the whole trajectory for plots against time
                        data['aft'].extend([(region, pname, time,
                                             anccod, mutcod,
                                             posdna, posdnaref,
                                             anc, nuc, mut,
                                             Ssubpos,
                                             mutclass, trclass, af)
                                            for af, time in izip(aftpos[inuc], timespos)])

    data['aft'] = pd.DataFrame(data=data['aft'],
                              columns=['Region', 'Patient', 'Time',
                                       'Anccod', 'Mutcod',
                                       'Pos', 'PosRef',
                                       'Anc', 'Mut', 'Substitution',
                                       'Ssub',
                                       'Class', 'Trclass', 'af'])

    # Print some basic output
    if pnames is not None:
        title = ' + '.join(pnames)+',\n'+' + '.join(regions)
    else:
        title = ' + '.join(regions)

    g = data['aft'].groupby(['Region', 'PosRef'])

    # Correlate max frequency in any patient with subtype entropy
    # NOTE: this is the null of any more fancy model
    from scipy.stats import spearmanr
    data_maxf = g[['Ssub', 'af']].max()
    print 'Correlation between max frequency in patients and subtype entropy:',
    print spearmanr(*(np.array(data_maxf).T))


    # Ask about overcoming a fixed threshold at least once per patient
    threshold = 0.01
    g = data['aft'].groupby(['Region', 'PosRef', 'Patient'])
    datum_th = (g['af'].apply(lambda x: (x > threshold).any())
               .unstack('Patient')
               .mean(axis=1))
    tmp = [Ssub[i[1]] for i in datum_th.index]
    data_th = pd.DataFrame.from_dict({'frac_threshold': datum_th, 'Sub': tmp})
    print ('Correlation between fraction above frequency threshold '+
           '('+str(threshold)+')'+
           ' in patients and subtype entropy:'),
    print spearmanr(*(np.array(data_th).T))

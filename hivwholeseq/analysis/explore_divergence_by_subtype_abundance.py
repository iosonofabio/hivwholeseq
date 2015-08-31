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
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.cross_sectional.get_subtype_allele_frequencies import (
    get_subtype_reference_alignment_allele_frequencies)
from hivwholeseq.patients.one_site_statistics import get_codons_n_polymorphic



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Accumulation of minor alleles stratified by abundance difference in subtype',
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

    # FIXME
    pbads = ('p4', 'p7')
    patients = patients.loc[-patients.index.isin(pbads)]

    for region in regions:
        if VERBOSE >= 1:
            print region

        af_sub = get_subtype_reference_alignment_allele_frequencies(region)

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

            # Get the map as a dictionary from patient to subtype
            coomapd = dict(coomap[:, ::-1])

            # Get only sites that are conserved or with only one site per codon changing
            # NOTE: this is a more general concern about epistasis, but there's little to
            # do without moving to haplotypes (which might actually be fine)
            ind_poly, pos_poly = get_codons_n_polymorphic(aft, icons, n=[0, 1], VERBOSE=VERBOSE)

            # FIXME: deal better with depth (this should be already there?)
            aft[aft < 2e-3] = 0

            for pos, pospolys in izip(ind_poly, pos_poly):
                for poscod in xrange(3):
                    posdna = pos * 3 + poscod

                    if VERBOSE >= 3:
                        print posdna

                    # Look for this position in the subtype alignment
                    if posdna not in coomapd:
                        continue
                    pos_sub = coomapd[posdna]
                    if (pos_sub >= af_sub.shape[1]):
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

                        # Get only non-masked time points
                        aft_der = aft[:, ider, posdna]
                        indpost = -aft_der.mask
                        if indpost.sum() == 0:
                            continue
                        timespos = times[indpost]
                        aft_der = aft_der[indpost]
                        
                        mut = anc + '->' + der

                        # Get the difference in subtype abundances
                        afanc_sub = af_sub[ianc, pos_sub]
                        afder_sub = af_sub[ider, pos_sub]
                        daf_sub = afder_sub - afanc_sub

                        for it, time in enumerate(timespos):
                            af = aft_der[it]
                            data.append((region, pcode,
                                         anc, der, mut,
                                         afanc_sub, afder_sub, daf_sub,
                                         time, af))

    data = pd.DataFrame(data=data, columns=('region', 'pcode',
                                            'anc', 'der', 'mut',
                                            'afancsub', 'afdersub', 'dafsub',
                                            'time', 'af'))


    if plot:

        datap = data.loc[:, ['dafsub', 'time', 'af']]

        # Manually 
        bins_daf = np.array([-1, -0.5, -0.1, -0.02, 0.02, 0.1, 0.5, 1])
        bins_t = 30.5 * np.array([-1, 12, 24, 48, 72, 120])

        binsc_daf = 0.5 * (bins_daf[1:] + bins_daf[:-1])
        binsc_t = 0.5 * (bins_t[1:] + bins_t[:-1])

        afmean = np.ma.masked_all((len(binsc_daf), len(binsc_t)))
        for ib in xrange(len(binsc_daf)):
            for it in xrange(len(binsc_t)):
                datum = (datap.loc[(datap.loc[:, 'dafsub'] >= bins_daf[ib]) &
                                   (datap.loc[:, 'dafsub'] < bins_daf[ib+1]) &
                                   (datap.loc[:, 'time'] >= bins_t[it]) &
                                   (datap.loc[:, 'time'] < bins_t[it+1])]
                         .loc[:, 'af']
                         .mean())

                if not np.isnan(datum):
                    afmean[ib, it] = datum


        vmin = -4
        vmax = -1

        fig, ax = plt.subplots(figsize=(12, 8))
        h = ax.imshow(np.log10(afmean + 1e-5),
                      aspect='auto',
                      cmap=cm.jet,
                      interpolation='nearest',
                      vmin=vmin, vmax=vmax)

        ax.set_xticks(np.arange(len(binsc_t)))
        ax.set_xticklabels([str(int(time / 30.5)) for time in binsc_t], rotation=90)
        ax.set_xlabel('Time from infection [months]')
        
        ax.set_yticks(np.arange(len(bins_daf)) - 0.5)
        ax.set_yticklabels(['{:1.2f}'.format(be) for be in bins_daf])
        ax.set_ylabel('Abundance difference in subtype')

        cb = plt.colorbar(h)
        cb.set_label('log10(mean allele frequency)', rotation=270, labelpad=30)

        plt.tight_layout()

        plt.ion()
        plt.show()

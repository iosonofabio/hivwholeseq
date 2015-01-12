# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/07/14
content:    Explore sites that are very conserved within our patient set, across
            the whole infection, compared to their behaviour within the subtype.
'''
# Modules
import os
import argparse
from itertools import izip, combinations, chain
from collections import defaultdict
from operator import itemgetter
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO
from scipy.stats import pearsonr, spearmanr

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.filenames import root_patient_folder
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.sequence_utils import translate_alignment
from hivwholeseq.one_site_statistics import get_entropy
from hivwholeseq.multipatient.get_shared_alleles_trajectories import (
    get_shared_allele_frequencies, get_patient_indices)
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment,
    get_ali_entropy,
    get_subtype_reference_alignment_entropy)



# Globals
colors = {'p17': 'r',
          'p24': 'g',
          'PR': 'r',
          'RT': 'g',
          'p15': 'purple',
          'IN': 'orange',
          'gp41': 'r'}



# Functions
def get_allele_freqs_alignment(alim, positions=None, VERBOSE=0):
    '''Get allele frequencies from alignment'''
    from hivwholeseq.miseq import alpha, alphal

    # Make sure it's a numpy matrix, efficiently
    alim = np.asarray(alim, 'S1')

    if positions is None:
        positions = np.arange(alim.shape[-1])

    # Get allele freqs ignoring ambiguous: N, Y, R, etc. and renormalizing
    num = np.array([(alim == a).mean(axis=0) for a in alpha[:5]])
    num /= num.sum(axis=0)

    return num[:, positions]


def get_entropy_pats(afts, VERBOSE=0):
    '''Calculate entropies'''
    n_patients = afts.shape[0]
    lseq = afts.shape[1]

    if VERBOSE >= 1:
        print 'Calculate entropy for each patient'
    S = -np.ones((n_patients, lseq), float)
    S_time = np.ma.zeros((n_patients, lseq), object)
    times_covs = np.zeros(n_patients, object)
    for k, aft_pat in enumerate(afts):
        times_cov = aft_pat[0, 0] >= -0.1
        times_covs[k] = times_cov
        if times_cov.sum() < 2:
            S_time[k, :] = np.ma.masked
            continue
        
        for pos, aft_pos in enumerate(aft_pat):
            # FIXME: deal better with missing data (low coverage)
            Stmp = 0
            Stmp_time = np.zeros(len(times_cov))
            for j, aft_nuc in enumerate(aft_pos):
                aft_nuc_cov = aft_nuc[times_cov]
                Stmp -= ((aft_nuc_cov + 1e-8) * np.log2(aft_nuc_cov + 1e-8)).mean()
                Stmp_time -= ((aft_nuc + 1e-8) * np.log2(aft_nuc + 1e-8))
            S[k, pos] = Stmp
            S_time[k, pos] = Stmp_time

    return {'mean': S,
            'time': S_time}



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Explore conservation levels across patients and subtype',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot results')

    args = parser.parse_args()
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot

    plt.ioff()
    data_all = {}

    for ir, region in enumerate(regions):
        if VERBOSE >= 1:
            print region

        if VERBOSE >= 2:
            print 'Get subtype B reference sequence alignment to HXB2 of', region
        ali_sub = get_subtype_reference_alignment(region, VERBOSE=VERBOSE)

        if VERBOSE >= 2:
            print 'Get shared allele trajectories'
        # Low-coverage regions are bytecoded by -1
        data = get_shared_allele_frequencies(region, VERBOSE=VERBOSE, save=False)
        afts = data['afts']
        times = data['times']
        mapref = np.ma.masked_equal(data['mapref'], -1)
        refname = data['pnames'][-1]
        pnames = data['pnames'][:-1]
        samplenames = data['samplenames']
        n_patients = afts.shape[0]
        lseq = afts.shape[1]
        lalpha = afts.shape[2]
        patinds = get_patient_indices(samplenames, VERBOSE=VERBOSE)

        # Shrink coordinates (Pavel trimmed some protein e.g. RT)
        if mapref[-1] >= ali_sub.get_alignment_length():
            indpos = (mapref < ali_sub.get_alignment_length())
            mapref = mapref[indpos]
            afts = afts[:, :, indpos]
            lseq = afts.shape[1]


        ##############################################
        # SITE ENTROPY
        ##############################################
        if VERBOSE >= 2:
            print 'Get entropy of patients'
        Spats = get_entropy(afts, VERBOSE=VERBOSE)


        if VERBOSE >= 2:
            print 'Get entropy in subtype alignment'
        Ssub = get_ali_entropy(ali_sub, positions=mapref, VERBOSE=VERBOSE)


        if VERBOSE >= 2:
            print 'Comparing entropy in each patient and in the subtype'
        data = defaultdict(list)
        Ssub_shuffled = Ssub.copy(); np.random.shuffle(Ssub_shuffled)


        # Time resolved
        for ipat, pname in enumerate(pnames):
            Spat = Spats[patinds[ipat]]
            data['spearmanr-time'].append([spearmanr(Spat_t, Ssub) for Spat_t in Spat])
            data['spearmanr-time-shuffled'].append([spearmanr(Spat_t, Ssub_shuffled)
                                                    for Spat_t in Spat])
        
        if plot:
            if VERBOSE >= 2:
                print 'Plot entropy correlations with time'

            fig, ax = plt.subplots()
            for ipat, rps in enumerate(data['spearmanr-time']):
                ts = times[patinds[ipat]]
                ax.plot(ts, zip(*rps)[0], lw=2, label=pnames[ipat],
                        color=cm.jet(1.0 * ipat / len(pnames)))

            ax.set_xlabel('Time from infection [days]')
            ax.set_ylabel('Spearman r on entropy (pat VS subtype B)')

            ax.legend(loc=2, fontsize=10, ncol=2)
            ax.set_ylim(0, 0.5)
            ax.grid(True)
            ax.set_title(region)

            plt.tight_layout()
            plt.ion()
            plt.show()


        ## Time and entropy-class resolved
        #Sbins = [1e-3, 1e-1, 2]
        #ls = ['--', '-']
        #for Spat in Spats['time']:
        #    Spat = np.vstack(Spat).T

        #    Stmp = []
        #    for i in xrange(len(Sbins) - 1):
        #        ind = (Ssub >= Sbins[i]) & (Ssub < Sbins[i+1])
        #        Stmp.append([spearmanr(Spat_t[ind], Ssub[ind]) for Spat_t in Spat])
        #    data['spearmanr-time-entropy'].append(Stmp)


        ## Plot the various entropy classes
        #if plot:
        #    fig, ax = plt.subplots()
        #    for ipat, (ts, rpsS) in enumerate(izip(times, data['spearmanr-time-entropy'])):
        #        # FIXME: select only some patients for clarity
        #        if patients.iloc[ipat].name not in ['20097', '15319', '15107']:
        #            continue
        #        for iSbin, rps in enumerate(rpsS):
        #            ax.plot(ts, zip(*rps)[0], lw=2,
        #                    label=patients.iloc[ipat].name+', S_B e ['+str(Sbins[iSbin])+', '+str(Sbins[iSbin+1])+']',
        #                    ls=ls[iSbin], 
        #                    color=cm.jet(1.0 * ipat / len(patients)))

        #    ax.set_xlabel('Time from infection [days]')
        #    ax.set_ylabel('Spearman r on entropy (pat VS subtype B)')

        #    ax.legend(loc=2, fontsize=10, ncol=2)
        #    ax.set_ylim(0, 0.5)
        #    ax.grid(True)
        #    ax.set_title(region)

        #    plt.tight_layout()


        #    plt.ion()
        #    plt.show()


        ###############################################
        ## DEGENERATE CODONS
        ###############################################
        #if VERBOSE >= 2:
        #    print 'Get translated alignment'
        #ali_prot = get_subtype_reference_alignment(region, VERBOSE=VERBOSE, type='aa')
        #if len(ali_prot) != len(ali_sub):
        #    # Sometimes the first sequence is HXB2 itself, for some obscure reason
        #    if ((len(ali_prot) == len(ali_sub) + 1) and 
        #        ('HXB2' in ali_prot[0].name) and \
        #        ([seq.name for seq in ali_sub] == [seq.name for seq in ali_prot[1:]])):
        #        ali_prot = ali_prot[1:]

        #    # If things go bad, just remake the translated alignment
        #    else:
        #        ali_prot = translate_alignment(ali_sub, VERBOSE=VERBOSE)

        ## Get 4-fold degenerate codon positions
        #pos_4fold = get_degenerate_pos()


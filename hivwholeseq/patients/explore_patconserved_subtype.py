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
from hivwholeseq.sequencing.filenames import reference_folder
from hivwholeseq.patients.patients import load_patients, Patient



# Globals
tree_ali_foldername = reference_folder+'alignments_trees/typeM/'
colors = {'p17': 'r',
          'p24': 'g',
          'PR': 'r',
          'RT': 'g',
          'p15': 'purple',
          'IN': 'orange',
          'gp41': 'r'}



# Functions
def get_subtype_alignment(region, subtype='B', VERBOSE=0):
    '''Get the observables from subtype B alignments'''
    aliB_fn = tree_ali_foldername+region+'.'+subtype+'.nuc.aligned.fasta'
    aliB = AlignIO.read(aliB_fn, 'fasta')
    return aliB


def get_allele_freqs_alignment(alim, VERBOSE=0):
    '''Get allele frequencies from alignment'''
    from hivwholeseq.miseq import alpha, alphal

    # Make sure it's a numpy matrix, efficiently
    alim = np.asarray(alim, 'S1')

    # Get allele freqs ignoring ambiguous: N, Y, R, etc. and renormalizing
    num = np.array([(alim == a).mean(axis=0) for a in alpha[:5]])
    num /= num.sum(axis=0)

    return num


def get_entropy(num, VERBOSE=0):
    '''Get entropy from allele freqs'''
    S = np.zeros(num.shape[1])
    for nu in num:
        S -= (nu) * np.log2(nu + 1e-8)
    
    # Transform -0 in +0, for visual convenience
    S = np.abs(S)
    return S


def bridge_coordinate_maps(refcoos, ref, VERBOSE=0):
    '''Bridge patient alignment maps to subtype alignment'''
    refm = np.array(ref, 'S1')
    alics = (refm != '-').cumsum() - 1
    alicoos = np.array([(alics == refcoo).nonzero()[0][0] for refcoo in refcoos], int)
    return alicoos


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
    parser = argparse.ArgumentParser(description='Explore conservation levels across patients and subtype',
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
            print 'Get patient alignment coordinates, stripping gaps'
        fn_ali = root_patient_folder+'all/aft_shared_ali_'+region+'.fasta'
        ali = AlignIO.read(fn_ali, 'fasta')
        refseq = ali[-1]

        if VERBOSE >= 2:
            print 'Get map of aligned positions to reference (HXB2)'
        fn_map = root_patient_folder+'all/aft_shared_maps_'+region+'.npz'
        npdata = np.load(fn_map)
        mapref = np.ma.masked_equal(npdata['mapref'], -1)
        refname = npdata['pnames'][-1]
        pnames = npdata['pnames'][:-1]

        if ir == 0:
            patients = load_patients()
            patients = patients.loc[patients.index.isin(pnames)]


        if VERBOSE >= 2:
            print 'Get subtype B multiple sequence alignment of', region
        ali_sub = get_subtype_alignment(region, VERBOSE=VERBOSE)


        if VERBOSE >= 2:
            print 'Find reference in subtype alignment'
        for iref_sub, refseq_sub in enumerate(ali_sub):
            if refname in refseq_sub.name:
                break
        else:
            raise ValueError('Reference '+refname+' not found in subtype alignment')


        if VERBOSE >= 2:
            print 'Check integrity of reference from patient ali and subtype ali'
        if str(refseq.seq.ungap('-')) != str(refseq_sub.seq.ungap('-')):
            raise ValueError('Reference '+refname+' differs from the patient ali and subtype ali')


        if VERBOSE >= 2:
            print 'Build bridge map between aligned patients and subtype alignment via reference'
        map_ali = bridge_coordinate_maps(mapref, refseq_sub, VERBOSE=VERBOSE)


        if VERBOSE >= 2:
            print 'Get shared allele trajectories'
        
        # FIXME: we should use the function in the other script, that's what it's there for
        # Low-coverage regions are bytecoded by -1
        from hivwholeseq.patients.filenames import root_patient_folder
        npdata = np.load(root_patient_folder+'all/aft_shared_'+region+'.npz')
        # NOTE: this is not a rectangular tensor, so we keep the bytecode and will deal
        # with it as a mask all along
        afts = npdata['afts']
        times = npdata['times']
        n_patients = afts.shape[0]
        lseq = afts.shape[1]
        lalpha = afts.shape[2]


        if VERBOSE >= 2:
            print 'Get entropy of patients'
        Spats = get_entropy_pats(afts, VERBOSE=VERBOSE)


        if VERBOSE >= 2:
            print 'Get entropy in subtype alignment'
        alim_sub = np.array(ali_sub)
        af_sub = get_allele_freqs_alignment(alim_sub, VERBOSE=VERBOSE)[:, map_ali]
        Ssub = get_entropy(af_sub, VERBOSE=VERBOSE)

    
        if VERBOSE >= 2:
            print 'Comparing entropy in each patient and in the subtype'
        
        data = {}
        from scipy.stats import spearmanr, pearsonr
        data['spearmanr'] = [spearmanr(Spat, Ssub) for Spat in Spats['mean']]
        data['pearsonr'] = [pearsonr(Spat, Ssub) for Spat in Spats['mean']]

        data['spearmanr-time'] = [[spearmanr(Spat_t, Ssub) for Spat_t in np.vstack(Spat).T]
                                  for Spat in Spats['time']]
        

        if plot:
            if VERBOSE >= 2:
                print 'Plot correlations'

            fig, ax = plt.subplots()
            for ipat, (ts, rps) in enumerate(izip(times, data['spearmanr-time'])):
                ax.plot(ts, zip(*rps)[0], lw=2, label=patients.iloc[ipat].name,
                        color=cm.jet(1.0 * ipat / len(patients)))

            ax.set_xlabel('Time from infection [days]')
            ax.set_ylabel('Spearman r on entropy (pat VS subtype B)')

            ax.legend(loc=2, fontsize=10, ncol=2)
            ax.set_ylim(0, 0.5)
            ax.grid(True)
            ax.set_title(region)

            plt.tight_layout()
            plt.ion()
            plt.show()


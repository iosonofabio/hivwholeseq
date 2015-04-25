# vim: fdm=marker
'''
author:     Richard Neher, Fabio Zanini
date:       24/04/15
content:    Measure linkage disequilibrium.
'''
# Modules
import os
import sys
import glob
import argparse
import numpy as np
from itertools import izip
import pandas as pd
import hivwholeseq.utils.plot
from matplotlib import pyplot as plt

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.cross_sectional.get_subtype_consensus import (
    get_subtype_reference_alignment_consensus)
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment_entropy)
from hivwholeseq.cross_sectional.get_subtype_allele_frequencies import (
    get_subtype_reference_alignment_allele_frequencies)
from hivwholeseq.patients.one_site_statistics import get_codons_n_polymorphic



# Functions
def plot_LD(data, VERBOSE=0):
    '''Plot linkage disequilibrium'''
    import seaborn as sns
    sns.set_style('whitegrid')
    fs=16

    binc = data['binc']

    if 'LD' in data:
        LD_vs_distance = data['LD']

        fig = plt.figure(1, figsize=(5, 4))
        ax = plt.subplot(111)
        for frag, y in LD_vs_distance.iteritems():
            ax.plot(binc, y, label=frag, lw=3)
        ax.set_xticks(range(0, 401, 100))
        ax.set_yticks(np.arange(0., 1.01, 0.2))
        for item in ax.get_xticklabels() + ax.get_yticklabels():
            item.set_fontsize(fs)
        ax.set_ylabel('linkage disequilibrium r', fontsize=fs)
        ax.set_xlabel("distance [bp]", fontsize=fs)
        ax.legend(loc=1, fontsize=fs, title='PCR fragment:')

        plt.tight_layout(rect=(0, 0, 0.98, 1))
        #for ext in ['pdf','svg']:
        #    plt.savefig(figpath+'LD_rsq.'+ext)

    if 'Dp' in data:
        Dp_vs_distance = data['Dp']

        fig = plt.figure(2, figsize=(5, 4))
        ax = plt.subplot(111)
        for frag, y in Dp_vs_distance.iteritems():
            ax.plot(binc, y, label=frag, lw=3)

        ax.set_xticks(range(0, 401, 100))
        ax.set_yticks(np.arange(0., 1.01, 0.2))
        for item in ax.get_xticklabels() + ax.get_yticklabels():
            item.set_fontsize(fs)
        ax.set_ylabel("linkage disequilibrium D'", fontsize=fs)
        ax.set_xlabel("distance [bp]", fontsize=fs)
        ax.legend(loc=1, fontsize=fs, title='PCR fragment:')

        plt.tight_layout(rect=(0, 0, 0.98, 1))
        #for ext in ['pdf','svg']:
        #    plt.savefig(figpath+'LD_Dp.'+ext)


def collect_data(pnames, sample=None,
                 fragments=['F'+str(i) for i in xrange(1, 7)],
                 dmin=100,
                 varmin=0.2,
                 VERBOSE=0):
    '''Collect data on LD'''
    patients = load_patients(pnames=pnames)

    fn_ccs = {}
    cov_min=100
    LD_vs_distance = {}
    Dp_vs_distance = {}
    for frag in fragments:
        if VERBOSE >= 1:
            print frag

        from hivwholeseq.patients.filenames import root_patient_folder
        if sample is not None:
            ac_count_dir = root_patient_folder+pnames[0]+'/samples/'+sample+'/PCR1/'
            flist = sorted(glob.glob(ac_count_dir+'allele_cocounts*'))
        else:
            flist = []
            for pat, pdata in patients.iterrows():
                pdata = Patient(pdata)
                if VERBOSE >= 2:
                    print frag, pdata.code

                for samplename, sampl in pdata.samples.iterrows():
                    if VERBOSE >= 2:
                        print samplename, frag, sampl[frag], sampl[frag+'q']
                    if sampl[frag] == 'ok' and sampl[frag+'q'] > dmin:
                        flist.extend(glob.glob((root_patient_folder+pat+
                                                "/samples/"+samplename+
                                                '/PCR1/allele_cocounts*'+frag+
                                                '*.npz')))
                    else:
                        if VERBOSE >= 2:
                            print "EXCLUDED"
            flist.sort()

        fn_ccs[frag] = flist

    for frag, flist in fn_ccs.iteritems():

        # Scan the files with cocounts
        dists = []
        weights_LD = []
        weights_Dp = []
        bins = np.arange(0, 401, 40)
        binc = (bins[:-1] + bins[1:]) * 0.5
        for fname in flist:
            if VERBOSE >= 2:
                print fname
            
            cc = np.load(fname)['cocounts']

            L = cc.shape[-1]
            ac = cc[:, :,
                    np.arange(L), np.arange(L)].sum(axis=1)
            cov = ac[:5].sum(axis=0)

            af = 1.0 * ac / cov
            diverse_sites = ((1.0 - np.sum(af**2, axis=0)) > varmin)
            if VERBOSE >= 2:
                print np.where(diverse_sites)

            if np.sum(diverse_sites) <= 1:
                continue

            # Make reduced cocount matrix
            reduced_af = af[:, diverse_sites]
            positions = np.where(diverse_sites)[0]
            ndiv = reduced_af.shape[-1]
            reduced_cc = np.zeros((cc.shape[0], cc.shape[1], ndiv, ndiv))
            reduced_2paf = np.zeros((ndiv, ndiv), dtype=float)
            p = reduced_af.max(axis=0)
            pi = reduced_af.argmax(axis=0)
            q = 1 - p
            if VERBOSE >= 2:
                print p
            for di,dsite in enumerate(positions):
                reduced_cc[:,:,di,:] = cc[:,:,dsite, diverse_sites]
            # extract product count of major alleles
            for posi, (nuci, dsite) in enumerate(zip(pi,positions)):
                reduced_2paf[posi,:] = cc[nuci, pi, dsite,diverse_sites]

            reduced_cc_cov = reduced_cc[:5,:5].sum(axis=1).sum(axis=0)
            reduced_2paf /= (1e-10+reduced_cc_cov)
            reduced_2paf_LE = np.outer(p,p)
            reduced_2pqsq = np.outer(p*q,p*q)
            reduced_2pp = np.outer(p,p) 
            reduced_2pq = np.outer(p,q)

            # Calculate LD in various ways
            # r^2 = (2pfreq - p1*p2)**2/(p1*q1*p2*q2)
            LD = np.sqrt((reduced_2paf - reduced_2paf_LE)**2 / 
                         (1e-10 + reduced_2pqsq))

            # D' = (2pfreq - p1*p2)/min(p1(1-p2), p2(1-p1)) if D>0
            # D' = (2pfreq - p1*p2)/min(p1p2, p2p1) if D<0
            Dp = reduced_2paf - reduced_2paf_LE
            Dp[Dp>0] /= np.minimum(reduced_2pq, reduced_2pq.T)[Dp>0]
            Dp[Dp<0] /= np.minimum(reduced_2pp, reduced_2pp.T)[Dp<0]
            Dp = np.abs(Dp)
            np.fill_diagonal(LD,0)
            np.fill_diagonal(Dp,0)
            np.fill_diagonal(reduced_cc_cov,0)

            X,Y = np.meshgrid(positions, positions)
            dists.extend(np.abs(X-Y)[reduced_cc_cov>=cov_min])
            weights_LD.extend(LD[reduced_cc_cov>=cov_min])
            weights_Dp.extend(Dp[reduced_cc_cov>=cov_min])

        yn,xn = np.histogram(dists, bins=bins)
        
        y,x = np.histogram(dists, weights=weights_LD, bins=bins)
        LD_vs_distance[frag] = y / (1e-10 + yn)

        y,x = np.histogram(dists, weights=weights_Dp, bins=bins)
        Dp_vs_distance[frag] = y / (1e-10 + yn)

    return {'LD': LD_vs_distance,
            'Dp': Dp_vs_distance,
            'binc': binc,
           }



# Script
if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Calculate LD")
    parser.add_argument('--sample',
                        help='sample to look at')
    parser.add_argument('--patients', action=PatientsAction,
                        help='Patients to analyze')
    parser.add_argument('--fragments', default=['F1'], nargs='+',
                        help='Fragments to analyze')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-3]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot matrices')
    dmin = 100
    varmin = 0.2

    args = parser.parse_args()
    pnames = args.patients
    sample = args.sample
    fragments = args.fragments
    VERBOSE = args.verbose
    show_plot = args.plot

    data = collect_data(pnames, sample=sample, fragments=fragments,
                        dmin=dmin,
                        varmin=varmin,
                        VERBOSE=VERBOSE)


    if show_plot:
        plot_LD(data,
                VERBOSE=VERBOSE)

        plt.ion()
        plt.show()



# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/07/14
content:    Explore parallel evolution in different patients by studying shared
            allele frequency trajectories.
'''
# Modules
import os
import argparse
from operator import itemgetter
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from Bio import SeqIO, AlignIO
import matplotlib.pyplot as plt

from hivwholeseq.multipatient.get_shared_alleles_trajectories import (
    get_shared_allele_frequencies, get_patient_indices)
from hivwholeseq.one_site_statistics import get_entropy



# Functions
def classify_sites(afts):
    '''Classify sites: 0 = conserved, 1 = polymorphic, 2 = sweeping'''
    site_class = np.ma.zeros((afts.shape[0], afts.shape[1]), int)
    for i, aft in enumerate(afts):
        aft = reshape_single_aft(aft)
        for j, aftpos in enumerate(aft):
            if (-aftpos[0].mask).sum() < 2:
                site_class[i, j] = np.ma.masked
                continue

            aftpos = aftpos[:, -aftpos[0].mask].data

            if ((aftpos[:, 0] < aftpos[:, 0].max()) & (aftpos > 0.5).any(axis=1)).any():
                site_class[i, j] = 2

            # FIXME: make this dependent on depth
            elif ((aftpos[:, 0] < aftpos[:, 0].max()) & (aftpos > 0.02).any(axis=1)).any():
                site_class[i, j] = 1

    return site_class


def maxentropy_sites(afts):
    '''Maximal entropy per site during each infection'''
    M = np.ma.zeros((afts.shape[0], afts.shape[1]), float)
    for i, aft in enumerate(afts):
        aft = reshape_single_aft(aft)
        M[i] = -(aft * np.log2(aft + 1e-8)).sum(axis=1).min(axis=1)
    return M


def entropy_sites(afts):
    '''Entropy sample by sample per site'''
    n_samples = sum(len(aft[0, 0]) for aft in afts)
    M = np.ma.zeros((n_samples, afts.shape[1]), float)
    isam = 0
    for aft in afts:
        aft = reshape_single_aft(aft).swapaxes(0, 2)
        for isampat, af in enumerate(aft):
            M[isam] = np.maximum(0, -(af * np.log2(np.maximum(af, 1e-8))).sum(axis=0))
            isam += 1
    return M


def get_derived_frequency(afts, patinds):
    '''Derived allele frequencies'''
    M = np.ma.zeros((afts.shape[0], afts.shape[2]), float)
    L = M.shape[-1]
    for ipat, patind in enumerate(patinds):
        aft = afts[patind]
        ic0 = aft[0].argmax(axis=0)
        M[patind] = 1.0 - aft[:, ic0, np.arange(L)]
        
    return M


def get_distance(M, criterion='subtract'):
    '''Get distance matrix given a site class matrix'''
    d = np.zeros((M.shape[0], M.shape[0]), float)
    for i in xrange(d.shape[0]):
        for j in xrange(d.shape[0]):
            if criterion == 'sweeps':
                d[i, j] = ((M[i] != 2) != (M[j] != 2)).mean()
            elif criterion == 'conserved':
                d[i, j] = ((M[i] != 0) != (M[j] != 0)).mean()
            elif criterion == 'discrete':
                d[i, j] = (M[i]  != M[j]).mean()
            else:
                d[i, j] = (np.abs(M[i] - M[j])).mean()

    return d


def reshape_single_aft(aft):
    '''Reshape a single aft (it's square now)'''
    aftn = np.ma.zeros((aft.shape[0], aft.shape[1], aft[0, 0].shape[0]), float)
    for i, af in enumerate(aft):
        aftn[i] = np.vstack(af)
    aftn.mask = aftn == -1
    return aftn


def permute_matrix(M, inds):
    '''Permute matrix cols and rows according to inds'''
    N = np.zeros_like(M)
    for i, ind in enumerate(inds):
        N[i] = M[ind, inds]
    return N


# PLOT FUNCTIONS
def plot_sample_by_sample(data_single, VERBOSE=0, cluster=True):
    '''Plot some observable as a clustered heatmap, sample by sample'''

    (d, M, labels, label) = data_single

    if cluster:
        Z = linkage(d)
        figsize=(18, 12)
        fig, axs = plt.subplots(1, 3, figsize=figsize, gridspec_kw={'width_ratios': [1, 4, 4]})
        # sort axes for convenience
        axs = [axs[1], axs[2], axs[0]]
    else:
        figsize=(15, 12)
        fig, axs = plt.subplots(1, 2, figsize=figsize, gridspec_kw={'width_ratios': [1, 1]})

    
    # Plot as tree and image
    fig.suptitle(region)

    if cluster:
        # Plot dendrogram
        dg = dendrogram(Z, orientation='right', ax=axs[2], leaf_label_func=lambda x: '')
        axs[2].set_xticklabels('')
        leaves = dg['leaves']
    else:
        leaves = np.arange(d.shape[0])

    # Plot profiles
    ax = axs[0]
    ax.imshow(M[leaves][::-1], interpolation='nearest', aspect='auto')
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels('')
    ax.set_xlim(0, M.shape[1])
    
    # Add lines at patient boundaries to the central plot
    leafnames = labels[leaves][::-1]
    for isam, sname1 in enumerate(leafnames[:-1]):
        sname2 = leafnames[isam + 1]
        pname1 = sname1.split('_')[0]
        pname2 = sname2.split('_')[0]
        if pname1 != pname2:
            ax.plot([0, M.shape[1] - 1], [0.5 + isam] * 2, lw=2, c='k')

    # Plot distance matrix
    dshu = permute_matrix(d, leaves)
    ax = axs[1]
    ax.yaxis.tick_right()
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticklabels(labels[leaves][::-1].tolist(), fontsize=10)
    ax.set_xticklabels('')
    ax.imshow(dshu[::-1], interpolation='nearest', aspect='auto')

    # Add lines at patient boundaries to the distance matrix plots
    leafnames = labels[leaves][::-1]
    for isam, sname1 in enumerate(leafnames[:-1]):
        sname2 = leafnames[isam + 1]
        pname1 = sname1.split('_')[0]
        pname2 = sname2.split('_')[0]
        if pname1 != pname2:
            ax.plot(ax.get_xlim(), [0.5 + isam] * 2, lw=2, c='k')
            ax.plot([len(leafnames) - 1.5 - isam] * 2, ax.get_xlim(), lw=2, c='k')


    if cluster:
        plt.tight_layout(rect=(0, 0, 1, 0.96), w_pad=0.01)
    else:
        plt.tight_layout(rect=(0, 0, 1, 0.96), w_pad=0.5)
    plt.ion()
    plt.show()




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get shared allele trajectories',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. F1 V3)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--mask-fraction-max', type=float, default=0.2,
                        help='Maximal fraction of masked sites to accept sample')

    args = parser.parse_args()
    regions = args.regions
    VERBOSE = args.verbose
    mask_fraction_max = args.mask_fraction_max

    for region in regions:
        if VERBOSE >= 1:
            print region

        # Low-coverage regions are bytecoded by -1
        data = get_shared_allele_frequencies(region, VERBOSE=VERBOSE, save=False)
        afts = np.ma.masked_less(data['afts'], -0.99)
        pnames = data['pnames'][:-1]
        depthmaxs = data['depthmaxs']
        n_patients = len(pnames)
        samplenames = data['samplenames']
        patinds = get_patient_indices(samplenames, VERBOSE=VERBOSE)

        # SAMPLE BY SAMPLE
        data_plot_samples = []

        # Entropy sample-by-sample
        site_S1 = get_entropy(afts, VERBOSE=VERBOSE)

        d_S1 = get_distance(site_S1, criterion='subtract')
        ind_good = list((-np.isnan(d_S1)).nonzero()) + \
                   list((site_S1.mask.mean(axis=1) < mask_fraction_max).nonzero())
        ind_good = reduce(np.intersect1d, ind_good)
        snames = samplenames[ind_good]
        site_S1 = site_S1[ind_good]
        d_S1 = get_distance(site_S1, criterion='subtract')
        data_plot_samples.append((d_S1, site_S1, snames, 'entropy'))

        # Derived allele freuency sample-by-sample
        aft_der = get_derived_frequency(afts, patinds)
        d = get_distance(aft_der, criterion='subtract')
        ind_good = list((-np.isnan(d)).nonzero()) + \
                   list((aft_der.mask.mean(axis=1) < mask_fraction_max).nonzero())
        ind_good = reduce(np.intersect1d, ind_good)
        snames = samplenames[ind_good]
        aft_der = aft_der[ind_good]
        d = get_distance(aft_der, criterion='subtract')
        data_plot_samples.append((d, aft_der, snames, 'derived allele freq'))

        for data_single in data_plot_samples:
            plot_sample_by_sample(data_single, VERBOSE=VERBOSE, cluster=False)

            
            ## Parallel sweeps in different patients        
            #is_sweeping = site_class == 2
            #h = np.bincount(is_sweeping.sum(axis=0))
            #fig, ax = plt.subplots()
            #plt.plot(np.arange(len(h)), h, lw=2, label='Data, '+region)
            #from scipy.stats import poisson
            #mu = curve_fit(poisson.pmf, np.arange(len(h)), 1.0 * h / h.sum(), p0=1)[0][0]
            #hth = poisson.pmf(np.arange(len(h)), mu)
            #plt.plot(np.arange(len(hth)), hth * h.sum(), lw=2, c='k',
            #         label='mu = '+'{:1.1e}'.format(mu)+', '+region)
            #
            #plt.grid()
            #plt.ylabel('# sweeping sites')
            #plt.xlabel('# patients')
            #plt.yscale('log')
            #plt.legend(loc=1)
            #
            #plt.ion()
            #plt.show()

# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Extract the frequency of minor alleles.
'''
# Modules
import os
import sys
import cPickle as pickle
import argparse
import re
from collections import defaultdict
from itertools import izip
import pysam
import numpy as np
from Bio import SeqIO
import matplotlib.cm as cm

# Matplotlib parameters
import matplotlib
params = {'axes.labelsize': 20, 
          'text.fontsize': 20,
          'legend.fontsize': 18,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'text.usetex': False}
matplotlib.rcParams.update(params)
import matplotlib.pyplot as plt


# Globals
VERBOSE = 1

# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']

from mapping.miseq import alpha, read_types
from mapping.filenames import get_allele_counts_filename, get_coverage_filename
from mapping.mapping_utils import get_fragment_list



# Functions
def get_minor_allele_counts(counts, n_minor=1):
    '''Extract the minor allele counts
    
    Parameters:
       - counts: the complete allele counts
       - n_minor: how many minor alleles?
    '''
    # Copy input structure
    counts = counts.copy()

    # Prepare output data structures
    # The last axis is: [base index (A := 0, C := 1, ...), counts]
    # Note: the first array is the major allele
    all_sorted = np.zeros((n_minor + 1, counts.shape[0], counts.shape[2], 2), int)
    for js, ctype in enumerate(counts):
        for pos, cpos in enumerate(ctype.swapaxes(0, 1)):
            for i, all_counts in enumerate(all_sorted):
                imax = cpos.argmax()
                all_counts[js, pos] = (imax, cpos[imax])
                cpos[imax] = -1

    return all_sorted



# Script
if __name__ == '__main__':


    # Input arguments
    parser = argparse.ArgumentParser(description='Study minor allele frequency')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--subsample', action='store_true',
                        help='Apply only to a subsample of the reads')

    args = parser.parse_args()
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    subsample = args.subsample

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over all requested samples
    for adaID in adaIDs:

        nus_minor = {}
        for fragment in fragments:

            counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
            coverage = np.load(get_coverage_filename(data_folder, adaID, fragment))

            # TODO: study insertions

            (counts_major,
             counts_minor,
             counts_minor2) = get_minor_allele_counts(counts, n_minor=2)

            # Get minor allele frequencies
            nu_minor = 1.0 * counts_minor[:, :, 1] / (coverage + 1e-6)
            nus_minor[fragment] = nu_minor

        # Plot them
        fig, axs = plt.subplots(2, 3, figsize=(21, 12))
        fig.suptitle('adapterID '+'{:02d}'.format(adaID), fontsize=20)
        axs = axs.ravel()
        for i, fragment in enumerate(fragments):
            ax = axs[i]
            ax.set_yscale('log')
            ax.set_title(fragment)
            if i in [0, 4]:
                ax.set_ylabel(r'$\nu$')
            if i > 2:
                ax.set_xlabel('Position')

            for js, nu_minorjs in enumerate(nus_minor[fragment]):
                color = cm.jet(int(255.0 * js / len(read_types)))
                ax.plot(nu_minorjs, lw=1.5, c=color)
                ax.scatter(np.arange(len(nu_minorjs)), nu_minorjs, lw=1.5,
                           color=color)
    
        plt.tight_layout(rect=(0, 0, 1, 0.95))

        # Correct errors and plot final allele frequencies


        plt.show()
    
    
            ## Filter out sequencing errors
            #errors_seq = np.zeros(coverage.shape[1], bool)
    
            ### Exclude SNP calls in regions of low coverage
            ##errors_seq |= (coverage < 1000).all(axis=0)
    
            ### Sequencing errors are characterized by not being present in all four read
            ### types. We check all the types with enough coverage and test for homogeneously
            ### distributed counts
            ##for pos in (-errors_seq).nonzero()[0]:
            ##    read_types_cov = coverage[:, pos] > 1000
            ##    count_pos = counts_minor[read_types_cov, pos]
            ##    if count_pos.any():
            ##        cov = coverage[read_types_cov, pos]
            #        
    
            ### Sequencing errors tend to have this pattern: they are seen in both forward
            ### reads, but not reverse reads (or vice versa)
            ##nu_mav = nu_minor.mean(axis=0)
            ##corr_sense = (nu_minor[0] - nu_mav) * (nu_minor[2] - nu_mav) + (nu_minor[1] - nu_mav) * (nu_minor[3] - nu_mav)
            ##corr_antis = (nu_minor[0] - nu_mav) * (nu_minor[1] - nu_mav) + (nu_minor[2] - nu_mav) * (nu_minor[3] - nu_mav)
            ##errors_seq |= (corr_sense > 0) & (corr_antis < 0) & (nu_mav > 0.001)
            ##
            ##errors_seq |= (nu_minor > 0.01).sum(axis=0) < 2 # Ignore if they are not seen in at least 2 read types
    
            #### Tag as errors if their frequencies in different read types are too different
            ###for pos in xrange(len(errors_seq)):
            ###    nu_seen = nu_minor[nu_minor[:, pos] > 0.01, pos]
            ###    if len(nu_seen) and (nu_seen.std() / nu_seen.mean() > 0.1):
            ###        errors_seq[pos] = True
    
            ### Tag as errors if the coverage is real small
            ##errors_seq |= (coverage < 3000).any(axis=0)
    
            ##plt.scatter(corr_sense, corr_antis, s=30)
            ##plt.plot([-0.15, 0.15], [0] * 2, lw=1.5, ls='--', c='k')
            ##plt.plot([0] * 2, [-0.15, 0.15], lw=1.5, ls='--', c='k')
    
            ## Plot?
            ## Set special stuff for special plots
            #if adaID == 16: density = True
            #else: density = False
    
            ## Plot data
            #import matplotlib.cm as cm
            #fig, axs = plt.subplots(3, 1, figsize=(14, 14))
            #for i in xrange(len(read_types)):
            #    axs[0].plot(sequenced_region[0] + np.arange(coverage.shape[1]), coverage[i],
            #                c=cm.jet(int(255.0 * i / len(read_types))))
            #    read_type = read_types[i]
            #    y = nu_minor[i].copy()
            #    #y[errors_seq] = 0
            #    axs[1].plot(sequenced_region[0] + np.arange(coverage.shape[1]), y + 0.0000001,
            #                lw=1.5, c=cm.jet(int(255.0 * i / len(read_types))))
            #    h = np.histogram(y, bins=np.logspace(-3, 0, 100), density=density)
            #    axs[2].plot(np.sqrt(h[1][1:] * h[1][:-1]), h[0], lw=2, label=read_type, c=cm.jet(int(255.0 * i / len(read_types))))
    
            #axs[0].set_ylabel('Coverage [# reads]')
            #axs[1].set_ylim(1e-5, 1e0)
            #axs[1].set_yscale('log')
            #axs[2].set_xscale('log')
            #axs[2].set_xlim(0.95e-3, 1.05)
            #axs[2].set_ylim(-0.15, axs[2].get_ylim()[1])
            #axs[1].set_xlabel('HIV genome [b.p.]')
            #axs[1].set_ylabel('Minor allele frequency')
            #axs[2].set_xlabel('Minor allele frequency')
            #axs[2].set_ylabel('# alleles')
            #axs[2].legend(loc=1)
            #axs[0].set_title('Sample: '+sample_name)
    
            ## Set MORE special stuff
            #if adaID == 16:
            #    x = np.logspace(-3, 0, 1000)
            #    axs[2].plot(x, 1e-3 / x**2, c='k', lw=2, ls='--')
            #    axs[2].annotate(r'$p(\nu) \sim 1 / \nu^2$', (9e-2, 0.09), (8e-3, 0.04),
            #                    arrowprops=dict(arrowstyle="->"),
            #                    fontsize=20)
            #    axs[2].set_ylim(1e-2, 1e3)
            #    axs[2].set_yscale('log')
            #if adaID == 18:
            #    axs[1].plot([0, 10000], [2e-1] * 2, lw=1.5, ls='--', c='k', alpha=0.3)
            #    axs[2].plot([2e-1] * 2, [0, 380], lw=1.5, ls='--', c='k', alpha=0.3)
            #if adaID == 19:
            #    axs[1].plot([0, 10000], [6e-2] * 2, lw=1.5, ls='--', c='k', alpha=0.3)
            #    axs[1].plot([0, 10000], [1e-2] * 2, lw=1.5, ls='--', c='k', alpha=0.3)
            #    axs[2].plot([6e-2] * 2, [0, 380], lw=1.5, ls='--', c='k', alpha=0.3)
            #    axs[2].plot([1e-2] * 2, [0, 380], lw=1.5, ls='--', c='k', alpha=0.3)
    
            #plt.tight_layout()
    
            #plt.ion()
            #plt.show()

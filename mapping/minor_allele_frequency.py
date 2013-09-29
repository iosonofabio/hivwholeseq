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
from operator import itemgetter
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

from mapping.miseq import alpha, read_types
from mapping.filenames import get_allele_counts_filename, get_coverage_filename
from mapping.mapping_utils import get_fragment_list
from mapping.adapter_info import load_adapter_table



# Globals
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']



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


def filter_nus(counts, coverage):
    '''Filter allele frequencies from the four read types'''
    pvals = np.zeros((len(alpha), counts.shape[-1]))

    # Divide binarily
    nocounts = (coverage - counts.swapaxes(0, 1)).swapaxes(0, 1)

    # Sum read1 and read2
    counts_f = counts[0] + counts[2]
    counts_b = counts[1] + counts[3]
    nocounts_f = nocounts[0] + nocounts[2]
    nocounts_b = nocounts[1] + nocounts[3]
    nus_f = 1.0 * counts_f / (coverage[0] + coverage[2] + 1)
    nus_b = 1.0 * counts_b / (coverage[1] + coverage[3] + 1)

    # Test chi^2 for consistency across read 1 and read2,
    # using pseudocounts
    from scipy.stats import chi2_contingency
    for j in xrange(len(alpha)):
        for i in xrange(counts.shape[-1]):
            cm = np.array([[counts_f[j, i], nocounts_f[j, i]],
                           [counts_b[j, i], nocounts_b[j, i]]], int)
            chi2, pval = chi2_contingency(cm + 1)[:2]
            pvals[j, i] = pval

    errs = zip(*(((pvals < 1e-6) & (np.abs(nus_f - nus_b) > 1e-4)).nonzero()))
    errs.sort(key=itemgetter(1))

    # Take the mean of the two fed and rev for non-errors
    nu_filtered = np.ma.masked_all((len(alpha), counts.shape[-1]))
    covtot = coverage.sum(axis=0)
    ind = covtot > 0
    # Geometric mean? Not much changes
    #nu_filtered[:, ind] = 1.0 * (counts_f + counts_b)[:, ind] / covtot[ind]
    nu_filtered[:, ind] = np.sqrt(nus_f * nus_b)[:, ind]

    # Insert corrections
    if VERBOSE:
        print '{:4s}'.format('pos'), '{:3s}'.format('nuc'), \
                '{:5s}'.format('cov fw'), '{:10s}'.format('nu fw'), \
                '{:10s}'.format('nu re'), '{:5s}'.format('cov rv')
    for ai, pos in errs:
        # Take the allele farthest away from 0.5
        # (the sum to 1 is not guaranteed!)
        nu_tmp = np.array([nus_f[ai, pos], nus_b[ai, pos]])
        ind_tmp = np.argmax(np.abs(nu_tmp - 0.5))
        nu_filtered[ai, pos] = nu_tmp[ind_tmp]

        if VERBOSE and ((nus_f[ai, pos] < 0.5) or (nus_b[ai, pos] < 0.5)):
            print '{:4d}'.format(pos), alpha[ai], \
                    '{:7d}'.format(coverage[[0, 2], pos].sum()), \
                    '{:1.1e}'.format(nus_f[ai, pos]), \
                    '{:1.1e}'.format(nus_b[ai, pos]), \
                    '{:7d}'.format(coverage[[1, 3], pos].sum())

    # Get rid of the mask if not needed
    if not nu_filtered.mask.any():
        nu_filtered = nu_filtered.data

    return nu_filtered


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
        alls_minor = {}
        nus_filtered = {}
        nus_minor_filtered = {}
        for fragment in fragments:

            try:
                counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
                coverage = np.load(get_coverage_filename(data_folder, adaID, fragment))
    
    
                # TODO: study insertions
    
                (counts_major,
                 counts_minor,
                 counts_minor2) = get_minor_allele_counts(counts, n_minor=2)
    
                # Get minor allele frequencies and identities
                nu_minor = 1.0 * counts_minor[:, :, 1] / (coverage + 1e-6)
                nus_minor[fragment] = nu_minor
                all_minor = counts_minor[:, :, 0]
                alls_minor[fragment] = all_minor
    
                # Filter the minor frequencies by comparing the read types
                nu_filtered = filter_nus(counts, coverage)
                nut = np.zeros(nu_filtered.shape[-1])
                for pos, nupos in enumerate(nu_filtered.T):
                    nut[pos] = np.sort(nupos)[-2]
                nus_filtered[fragment] = nu_filtered
                nus_minor_filtered[fragment] = nut

            except IOError:
               pass 

        # Plot them
        fig, axs = plt.subplots(2, 3, figsize=(13, 8))
        fig.suptitle('adapterID '+'{:02d}'.format(adaID), fontsize=20)
        axs = axs.ravel()
        for i, fragment in enumerate(fragments):
            try:
                ax = axs[i]
                ax.set_yscale('log')
                ax.set_title(fragment)
                if i in [0, 3]:
                    ax.set_ylabel(r'$\nu$')
                if i > 2:
                    ax.set_xlabel('Position')
    
                for js, nu_minorjs in enumerate(nus_minor[fragment]):
                    color = cm.jet(int(255.0 * js / len(read_types)))
                    ax.plot(nu_minorjs, lw=1.5, c=color, label=read_types[js])
                    ax.scatter(np.arange(len(nu_minorjs)), nu_minorjs, lw=1.5,
                               color=color)
                
                # Plot filtered
                ax.plot(nus_minor_filtered[fragment] + 1e-6, lw=1.5, c='k',
                        alpha=0.5, label='Filtered')

                ax.set_xlim(-100, len(nu_minorjs) + 100)
            except:
                pass
    
        plt.legend(loc='upper center')
        plt.tight_layout(rect=(0, 0, 1, 0.95))

        plt.ion()
        plt.show() 

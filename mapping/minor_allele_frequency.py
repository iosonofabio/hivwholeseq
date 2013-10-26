# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Extract the frequency of minor alleles.
'''
# Modules
import argparse
from operator import itemgetter
import numpy as np
import matplotlib.cm as cm

# Matplotlib parameters
import matplotlib
params = {'axes.labelsize': 20, 
          'text.fontsize': 20,
          'legend.fontsize': 8,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'text.usetex': False}
matplotlib.rcParams.update(params)
import matplotlib.pyplot as plt

from mapping.datasets import MiSeq_runs
from mapping.miseq import alpha, read_types
from mapping.filenames import get_allele_counts_filename, get_coverage_filename
from mapping.adapter_info import load_adapter_table



# Globals
plot_grid = [(1, 1), (1, 2), (1, 3), (2, 2), (1, 5), (2, 3)]



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


def filter_nus(counts, coverage, VERBOSE=0):
    '''Filter allele frequencies from the four read types'''
    from scipy.stats import chi2_contingency

    # Divide binarily
    nocounts = (coverage - counts.swapaxes(0, 1)).swapaxes(0, 1)

    # Set counts and similia: sum read1 and read2
    counts_f = counts[0] + counts[2]
    counts_b = counts[1] + counts[3]
    nocounts_f = nocounts[0] + nocounts[2]
    nocounts_b = nocounts[1] + nocounts[3]
    cov_f = coverage[0] + coverage[2]
    cov_b = coverage[1] + coverage[3]
    ind_low_cov_f = cov_f < 10
    ind_low_cov_b = cov_b < 10
    ind_high_cov_both = (-ind_low_cov_f) & (-ind_low_cov_b)

    # Set filtered allele frequencies
    nu_filtered = np.ma.masked_all((len(alpha), counts.shape[-1]))

    # Iterate over positions
    for i in xrange(counts.shape[-1]):
        
        # 1. if we cover neither fwd nor rev, keep masked
        if ind_low_cov_f[i] and ind_low_cov_b[i]:
            if VERBOSE >= 4:
                print 'pos', i, 'base', alpha[j], 'not covered'
            pass

        # 2. if we cover only one of them, well, just take the
        # arithmetic sum of counts
        elif ind_low_cov_f[i] != ind_low_cov_b[i]:
            nu_filtered[:, i] = 1.0 * counts[:, :, i].sum(axis=0) / coverage[:, i].sum()
            if VERBOSE >= 4:
                print 'pos', i, 'base', alpha[j], 'covered only once'
                
        # 3. If we cover both, check whether the counts are significantly different
        else:

            # Check all alleles
            for j in xrange(len(alpha)):
                # To make a table, you must have coverage for both
                cm = np.array([[counts_f[j, i], nocounts_f[j, i]],
                               [counts_b[j, i], nocounts_b[j, i]]], int)
                chi2, pval = chi2_contingency(cm + 1)[:2]
    
                # If they are not significantly different, sum the counts
                if (pval > 1e-6):
                    nu_filtered[j, i] = 1.0 * counts[:, j, i].sum(axis=0) / coverage[:, i].sum()
                # If they are different by a significant and reasonable amount, take
                # the value further away from 0.5
                else:
                    nu_f = 1.0 * counts_f[j, i] / cov_f[i]
                    nu_b = 1.0 * counts_b[j, i] / cov_b[i]
                    if np.abs(nu_f - 0.5) > np.abs(nu_b - 0.5):
                        nu_filtered[j, i] = nu_f
                    else:
                        nu_filtered[j, i] = nu_b

                    if VERBOSE >= 3:
                        print 'pos', i, 'base', alpha[j], 'nu_f', nu_f, 'nu_b', nu_b

    # Renormalize to 1
    nu_filtered /= nu_filtered.sum(axis=0)

    # Get rid of the mask if not needed
    if not nu_filtered.mask.any():
        nu_filtered = nu_filtered.data

    return nu_filtered



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Study minor allele frequency')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

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

        covs = {}
        nus_minor = {}
        alls_minor = {}
        nus_filtered = {}
        nus_minor_filtered = {}
        for fragment in fragments:
            counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
            coverage = np.load(get_coverage_filename(data_folder, adaID, fragment))

            # Store coverage
            covs[fragment] = coverage
    
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

        # Plot them
        (n_plots_y, n_plots_x) = plot_grid[len(fragments) - 1]
        fig, axs = plt.subplots(n_plots_y, n_plots_x, figsize=(13, 8))
        if len(fragments) > 1:
            axs = axs.ravel()
        else:
            axs = [axs]
        fig.suptitle('adapterID '+'{:02d}'.format(adaID), fontsize=20)
        for i, fragment in enumerate(fragments):
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
            ax.plot(nus_minor_filtered[fragment], lw=1.5, c='k',
                    alpha=0.5, label='Filtered')
            ax.scatter(np.arange(len(nus_minor_filtered[fragment])),
                       nus_minor_filtered[fragment], lw=1.5, c='k',
                       alpha=0.5)

            # Plot 1/max(coverage)
            coverage = covs[fragment]
            cov_tot = coverage.sum(axis=0)
            ax.plot(1.0 / cov_tot, lw=1.2, c='r', label='1/sum(coverage)')

            ax.set_xlim(-100, len(nu_minorjs) + 100)
    
        plt.legend(loc='upper center')
        plt.tight_layout(rect=(0, 0, 1, 0.95))

        plt.ion()
        plt.show() 

        # Plot distribution of minor allele frequencies (filtered)
        (n_plots_y, n_plots_x) = plot_grid[len(fragments) - 1]
        fig, axs = plt.subplots(n_plots_y, n_plots_x, figsize=(13, 8))
        if len(fragments) > 1:
            axs = axs.ravel()
        else:
            axs = [axs]
        fig.suptitle('adapterID '+'{:02d}'.format(adaID), fontsize=20)
        for i, fragment in enumerate(fragments):
            ax = axs[i]
            ax.set_xscale('log')
            ax.set_title(fragment)
            if i >= (n_plots_y - 1) * n_plots_x:
                ax.set_xlabel(r'$\nu$')
    
            h = np.histogram(nus_minor_filtered[fragment] + 1e-6,
                             bins = np.logspace(-6, 0, 50),
                             # This can be switched
                             density=False)
            x = np.sqrt(h[1][1:] * h[1][:-1])
            y = 1e-6 + h[0]

            ax.plot(x, y,
                    lw=2,
                    color='grey',
                    label='Filtered')

        plt.legend(loc='upper right')
        plt.tight_layout(rect=(0, 0, 1, 0.95))

        plt.ion()
        plt.show() 

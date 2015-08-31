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

from hivwholeseq.sequencing.samples import load_sequencing_run, SampleSeq
from hivwholeseq.miseq import alpha, read_types
from hivwholeseq.sequencing.filenames import get_allele_counts_filename, get_coverage_filename, \
        get_allele_frequencies_filename
from hivwholeseq.sequencing.adapter_info import load_adapter_table
from hivwholeseq.utils.one_site_statistics import filter_nus



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


def plot_minor_allele_frequency(data_folder, adaID, fragments, VERBOSE=0,
                                savefig=False):
    '''Plot minor allele frequency along the genome'''
    from hivwholeseq.sequencing.filenames import get_minor_allele_frequency_figure_filename as gff
    import matplotlib
    params = {'axes.labelsize': 20, 
              'text.fontsize': 20,
              'legend.fontsize': 8,
              'xtick.labelsize': 16,
              'ytick.labelsize': 16,
              'text.usetex': False}
    matplotlib.rcParams.update(params)
    from matplotlib import cm
    import matplotlib.pyplot as plt

    plot_grid = [(1, 1), (1, 2), (1, 3), (2, 2), (1, 5), (2, 3)]

    # Store in globals structures
    covs = {}
    nus_minor = {}
    alls_minor = {}
    nus_filtered = {}
    nus_minor_filtered = {}

    for fragment in fragments:
        coverage = np.load(get_coverage_filename(data_folder, adaID, fragment))
        covs[fragment] = coverage
    
        counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
        (counts_major,
         counts_minor,
         counts_minor2) = get_minor_allele_counts(counts, n_minor=2)
    
        # Get minor allele frequencies and identities
        nu_minor = 1.0 * counts_minor[:, :, 1] / (coverage + 1e-6)
        nus_minor[fragment] = nu_minor
        all_minor = counts_minor[:, :, 0]
        alls_minor[fragment] = all_minor
    
        # Filter the minor frequencies by comparing the read types
        try:
            nu_filtered = np.load(get_allele_frequencies_filename(data_folder, adaID, fragment))
        except IOError:
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
    fig.suptitle('adapterID '+adaID, fontsize=20)
    labss = {'read1 f': 'read1 fwd', 'read1 r': 'read1 rev',
             'read2 f': 'read2 fwd', 'read2 r': 'read2 rev'}
    for i, fragment in enumerate(fragments):
        ax = axs[i]
        ax.set_yscale('log')
        ax.set_title(fragment)
        if i in [0, 3]:
            ax.set_ylabel(r'$\nu$')
        if i > 2:
            ax.set_xlabel('Position')
    
        # Plot divided by readtype
        for js, nu_minorjs in enumerate(nus_minor[fragment]):
            color = cm.jet(int(255.0 * js / len(read_types)))
            ax.plot(nu_minorjs, lw=1.5, c=color, label=labss[read_types[js]])
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
        ax.plot(1.0 / cov_tot, lw=1.2, c='r', label='Detection limit')

        ax.set_xlim(-100, len(nu_minorjs) + 100)
    
    plt.grid()
    plt.legend(loc='upper right')
    plt.tight_layout(rect=(0, 0, 1, 0.95))

    if savefig:
        outputfile = gff(data_folder, adaID, fragment)
        fig.savefig(outputfile)
        plt.close(fig)
    else:
        plt.ion()
        plt.show()


def plot_minor_allele_frequency_filtered(data_folder, adaID, fragments, VERBOSE=0,
                                         savefig=False):
    '''Plot minor allele frequency along the genome'''
    from hivwholeseq.sequencing.filenames import get_minor_allele_frequency_figure_filename as gff
    import matplotlib
    params = {'axes.labelsize': 20, 
              'text.fontsize': 20,
              'legend.fontsize': 8,
              'xtick.labelsize': 16,
              'ytick.labelsize': 16,
              'text.usetex': False}
    matplotlib.rcParams.update(params)
    from matplotlib import cm
    import matplotlib.pyplot as plt

    # Store in globals structures
    covs = {}
    nus_minor_filtered = {}
    for fragment in fragments:
        coverage = np.load(get_coverage_filename(data_folder, adaID, fragment))
        covs[fragment] = coverage

        try:
            nu_filtered = np.load(get_allele_frequencies_filename(data_folder,
                                                                  adaID, fragment))
        except IOError:
            counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
            nu_filtered = filter_nus(counts)

        nut = np.zeros(nu_filtered.shape[-1])
        for pos, nupos in enumerate(nu_filtered.T):
            nut[pos] = np.sort(nupos)[-2]
        
        nus_minor_filtered[fragment] = nut
 
    # Plot them
    plot_grid = [(1, 1), (1, 2), (1, 3), (2, 2), (1, 5), (2, 3)]
    (n_plots_y, n_plots_x) = plot_grid[len(fragments) - 1]
    fig, axs = plt.subplots(n_plots_y, n_plots_x, figsize=(13, 8))
    if len(fragments) > 1:
        axs = axs.ravel()
    else:
        axs = [axs]
    fig.suptitle('adapterID '+adaID, fontsize=20)
    for i, fragment in enumerate(fragments):
        ax = axs[i]
        ax.set_yscale('log')
        ax.set_title(fragment)
        if i in [0, 3]:
            ax.set_ylabel(r'$\nu$')
        if i > 2:
            ax.set_xlabel('Position')
        
        # Plot filtered
        ax.plot(nus_minor_filtered[fragment], lw=1.5, c='k',
                alpha=0.5, label='Filtered')
        ax.scatter(np.arange(len(nus_minor_filtered[fragment])),
                   nus_minor_filtered[fragment], lw=1.5, c='k',
                   alpha=0.5)

        ax.set_xlim(-100, len(nus_minor_filtered[fragment]) + 100)
    
    #plt.legend(loc='upper right')
    plt.tight_layout(rect=(0, 0, 1, 0.95))

    if savefig:
        outputfile = gff(data_folder, adaID, fragment, only_filtered=True)
        fig.savefig(outputfile)
        plt.close(fig)
    else:
        plt.ion()
        plt.show()






# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Study minor allele frequency')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--no-savefig', action='store_false', dest='savefig',
                        help='Show figure instead of saving it')
    parser.add_argument('--only-filtered', action='store_true', dest='only_filt',
                        help='Plot only filtered frequencies')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    savefig = args.savefig
    only_filt = args.only_filt

    # Specify the dataset
    dataset = load_sequencing_run(seq_run)
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = dataset['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Iterate over all requested samples
    for adaID in adaIDs:

        # If the script is called with no fragment, iterate over all
        sample = SampleSeq(dataset.samples.loc[dataset.samples.adapter == adaID].iloc[0])
        samplename = sample.name
        if not fragments:
            fragments_sample = sample.regions_generic
        else:
            fragments_sample = sorted(set(fragments) & set(sample.regions_generic))

        if not only_filt:
            plot_minor_allele_frequency(data_folder, adaID, fragments_sample,
                                        VERBOSE=VERBOSE, savefig=savefig)
        else:
            plot_minor_allele_frequency_filtered(data_folder, adaID, fragments_sample,
                                                 VERBOSE=VERBOSE, savefig=savefig)

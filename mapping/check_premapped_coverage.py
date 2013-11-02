# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/13
content:    After the preliminary mapping to reference, plot coverage and allele
            frequencies to spot major issues.
'''
# Modules
import argparse
import numpy as np

from mapping.datasets import MiSeq_runs
from mapping.miseq import read_types
from mapping.reference import load_HXB2
from mapping.filenames import get_premapped_file
from mapping.one_site_statistics import get_allele_counts_insertions_from_file_unfiltered
from mapping.minor_allele_frequency import get_minor_allele_counts
from mapping.primer_info import primers_coordinates_HXB2_inner as pcis_HXB2



# Functions
def plot_coverage_minor_allele(counts, frags_pos=None, suptitle=None):
    '''Plot the coverage and the minor allele frequency'''
    cov = counts.sum(axis=1)
    cov_tot = cov.sum(axis=0)
    counts_minor = get_minor_allele_counts(counts)[1, :, :, 1]
    # Use pseudocounts so-so (it is only rough)
    nus_minor = 1.0 * counts_minor / (1 + cov)

    import matplotlib.pyplot as plt
    from matplotlib import cm
    fig, axs = plt.subplots(1, 2, figsize=(13, 8))
    axs[0].plot(cov_tot.T, lw=2, c='k', label=read_types)
    axs[0].set_xlabel('Position [bases]')
    axs[0].set_ylabel('Coverage')

    # If the fragments positions are marked, plot them
    if frags_pos is not None:
        for i, frag_pos in enumerate(frags_pos.T):
            axs[0].plot(frag_pos, 2 * [(0.98 - 0.02 * i) * axs[0].get_ylim()[1]],
                        c=cm.jet(int(255.0 * i / len(frags_pos.T))), lw=2)


    for i, nu_minor in enumerate(nus_minor):
        color = cm.jet(int(255.0 * i / len(read_types)))
        axs[1].plot(nu_minor, label=read_types, c=color)
        axs[1].scatter(np.arange(counts.shape[-1]), nu_minor,
                       s=30, c=color,
                       label=read_types)
    axs[1].set_xlabel('Position [bases]')
    axs[1].set_ylabel('Minor allele frequency')
    axs[1].set_yscale('log')
    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=18)

    plt.tight_layout(rect=(0, 0, 1, 0.95))

    plt.ion()
    plt.show()
    


def check_premap(miseq_run, adaID, qual_min=30, match_len_min=10,
                   reference='HXB2', maxreads=-1, VERBOSE=0):
    '''Check premap to reference: coverage, etc.'''
    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']
    F5_primer = dataset['primerF5'][dataset['adapters'].index(adaID)]

    # Get fragments
    fragments = ['F1', 'F2', 'F3', 'F4', F5_primer, 'F6']

    if reference == 'HXB2':
        refseq = load_HXB2(cropped=True)
        # This structure contains the inner primers coordinates in cropped ref
        frags_pos = np.zeros((2, len(fragments)), int)
        for i, fragment in enumerate(fragments):
            pci = pcis_HXB2[fragment]
            frags_pos[:, i] = (pci[0][0], pci[1][1])
        frags_pos -= frags_pos.min()
        print frags_pos

    else:
        raise ValueError('Only HXB2 is implemented as a reference')

    # Open BAM and scan reads
    input_filename = get_premapped_file(data_folder, adaID, type='bam')

    # Get counts
    counts, inserts = get_allele_counts_insertions_from_file_unfiltered(input_filename,
                                                             len(refseq),
                                                             qual_min=qual_min,
                                                             match_len_min=match_len_min,
                                                             maxreads=maxreads,
                                                             VERBOSE=VERBOSE)

    # Plot results
    title=', '.join(map(lambda x: ' '.join([x[0], str(x[1])]),
                        [['run', miseq_run],
                         ['adaID', adaID],
                         ['n_reads', maxreads],
                        ]))
    plot_coverage_minor_allele(counts,
                               frags_pos=frags_pos,
                               suptitle=title)


                



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Check consensus')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('-n', type=int, default=1000,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    n_reads = args.n
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = MiSeq_runs[miseq_run]['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Iterate over samples
    for adaID in adaIDs:

            check_premap(miseq_run, adaID,
                         maxreads=n_reads,
                         VERBOSE=VERBOSE)



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

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.miseq import read_types
from hivwholeseq.reference import load_HXB2, load_custom_reference
from hivwholeseq.filenames import get_premapped_file
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file_unfiltered
from hivwholeseq.minor_allele_frequency import get_minor_allele_counts
from hivwholeseq.primer_info import primers_coordinates_HXB2_inner as pcis_HXB2
from hivwholeseq.primer_info import primers_coordinates_HXB2_outer as pcos_HXB2
from hivwholeseq.mapping_utils import get_number_reads



# Functions
def plot_coverage_minor_allele(counts, frags_pos=None, frags_pos_out=None,
                               suptitle=None):
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
    # Inner primers
    if frags_pos is not None:
        for i, frag_pos in enumerate(frags_pos.T):
            axs[0].plot(frag_pos, 2 * [(0.97 - 0.03 * (i % 2)) * axs[0].get_ylim()[1]],
                        c=cm.jet(int(255.0 * i / len(frags_pos.T))), lw=2)

    # Outer primers
    if frags_pos_out is not None:
        for i, frag_pos in enumerate(frags_pos_out.T):
            axs[0].plot(frag_pos, 2 * [(0.96 - 0.03 * (i % 2)) * axs[0].get_ylim()[1]],
                        c=cm.jet(int(255.0 * i / len(frags_pos_out.T))), lw=2)

    axs[0].set_xlim(-500, 9500)

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
        frags_pos_out = np.zeros((2, len(fragments)), int)
        for i, fragment in enumerate(fragments):
            pci = pcis_HXB2[fragment]
            frags_pos[:, i] = (pci[0][0], pci[1][1])
            pco = pcos_HXB2[fragment]
            frags_pos_out[:, i] = (pco[0][0], pco[1][1])
        start = frags_pos.min()
        frags_pos -= start
        frags_pos_out -= start
        print frags_pos, frags_pos_out

    else:
        refseq = load_custom_reference(reference)
        frags_pos = None
        frags_pos_out = None

    # Open BAM and scan reads
    input_filename = get_premapped_file(data_folder, adaID, type='bam')

    # Count reads if requested
    if VERBOSE:
        print 'N. of reads:', get_number_reads(input_filename)

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
                               frags_pos_out=frags_pos_out,
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
    parser.add_argument('--reference', default='HXB2',
                        help='Use alternative reference, e.g. chimeras (the file must exist)')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    n_reads = args.n
    VERBOSE = args.verbose
    refname = args.reference

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
                         reference=refname,
                         VERBOSE=VERBOSE)


# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/13
content:    After the preliminary mapping to reference, plot coverage and allele
            frequencies to spot major issues.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.miseq import read_types
from hivwholeseq.filenames import get_premapped_file, get_reference_premap_filename, \
        get_fragment_positions_filename
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file_unfiltered
from hivwholeseq.minor_allele_frequency import get_minor_allele_counts
from hivwholeseq.primer_info import primers_coordinates_HXB2_inner as pcis_HXB2
from hivwholeseq.primer_info import primers_coordinates_HXB2_outer as pcos_HXB2
from hivwholeseq.mapping_utils import get_number_reads
from hivwholeseq.samples import samples



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
    


def check_premap(seq_run, adaID, qual_min=30, match_len_min=10,
                 maxreads=-1, VERBOSE=0):
    '''Check premap to reference: coverage, etc.'''

    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    samplename = dataset['samples'][dataset['adapters'].index(adaID)]
    fragments = samples[samplename]['fragments']

    refseq = SeqIO.read(get_reference_premap_filename(data_folder, adaID), 'fasta')

    fragpos_filename = get_fragment_positions_filename(data_folder, adaID)
    if os.path.isfile(fragpos_filename):
        fragtmp = list(np.loadtxt(fragpos_filename, usecols=[0], dtype='S10'))
        postmp = np.loadtxt(fragpos_filename, usecols=[1, 4], dtype=int)
        frags_pos = np.array([postmp[fragtmp.index(fr)] for fr in fragments], int).T

    else:
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
                        [['run', seq_run],
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
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--maxreads', type=int, default=1000,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    n_reads = args.maxreads
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = MiSeq_runs[seq_run]['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Iterate over samples
    for adaID in adaIDs:

            check_premap(seq_run, adaID,
                         maxreads=n_reads,
                         VERBOSE=VERBOSE)



# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/10/13
content:    Check overlapping regions for PCR amplification biases.

            The simplest way of doing this is via allele frequencies. More
            thorough analyses on the reads could be performed.
'''
# Modules
import os
import numpy as np
import argparse
from operator import itemgetter
from itertools import izip
import Bio.SeqIO as SeqIO
import Bio.AlignIO as AlignIO

from hivwholeseq.sequencing.filenames import get_consensus_filename, \
        get_allele_counts_filename, get_coverage_filename, \
        get_overlap_nu_figure_filename
from hivwholeseq.sequencing.minor_allele_frequency import filter_nus


# Functions
def get_overlapping_fragments(fragments):
    '''Get the overlapping pairs from a list of fragments'''
    possible_pairs = [('F'+str(i+1), 'F'+str(i+2)) for i in xrange(5)]
    pairs = []
    for pair in possible_pairs:
        if (pair[0] in fragments) and (pair[1] in fragments):
            pairs.append(pair)
    return pairs


def get_overlap(data_folder, adaID, frag1, frag2, VERBOSE=0):
    '''Find the overlap coordinates for the two fragments'''
    from hivwholeseq.mapping_utils import align_muscle

    seq1 = SeqIO.read(get_consensus_filename(data_folder, adaID, frag1), 'fasta')
    seq2 = SeqIO.read(get_consensus_filename(data_folder, adaID, frag2), 'fasta')
    sm1 = np.array(seq1)
    sm2 = np.array(seq2)

    # Find the beginning of s2 in s1
    seed_len = 20
    matches_min = 16
    seed = sm2[:seed_len]
    found = False
    trials = 0
    while (not found) and (trials < 3):
        for pos in xrange(len(seq1) - 700, len(seq1) - seed_len):
            if (sm1[pos: pos + seed_len] == seed).sum() >= matches_min - trials:
                found = True
                start_s2 = pos
                break
        if not found:
            trials += 1

    if not found:
        return None

    if VERBOSE >= 3:
        print 'Beginning of '+frag2+' found in '+frag1

    # In an ideal world, the overlap is a holy place in which no indels happen.
    # We cannot assume that, sadly. However, we can search from the other side
    # and align: find the end of s1 in s2
    found = False
    seed = sm1[-seed_len:]
    trials = 0
    while (not found) and (trials < 3):
        for pos in xrange(700):
            if (sm2[pos: pos + seed_len] == seed).sum() >= matches_min - trials:
                found = True
                end_s1 = pos + seed_len
                break
        if not found:
            trials += 1
    if not found:
        return None

    if VERBOSE >= 3:
        print 'End of '+frag1+' found in '+frag2

    # Align
    ali = align_muscle(seq1[start_s2:], seq2[:end_s1])
    return (start_s2, end_s1, ali)


def check_overlap_consensus(data_folder, adaID, frag1, frag2, overlap, VERBOSE=0):
    '''Check the overlap of the consensi'''
    # Check for an overlap at all
    if overlap is None:
        print 'adaID', adaID, frag1, frag2, 'no overlap found!'
        return True

    # Check for indels in the consensus
    (start_s2, end_s1, ali) = overlap
    alim = np.array(ali)
    if (alim == '-').any():
        print 'adaID', adaID, frag1, frag2, 'Indels in the overlap'
        return True

    # Check for mismatches
    if (alim[0] != alim[1]).any():
        print 'adaID', adaID, frag1, frag2, 'Mismatches in the overlap'
        for pos2 in (alim[0] != alim[1]).nonzero()[0]:
            pos1 = start_s2 + pos2
            print str(frag1)+':', pos1, str(frag2)+':', pos2, \
                    alim[0, pos2], alim[1, pos2]
        return True

    return False


def check_overlap_allele_frequencies(data_folder, adaID, frag1, frag2, overlap,
                                     VERBOSE=0, ax=None):
    '''Check biases in allele frequencies in the overlap'''
    (start_s2, end_s1, ali) = overlap

    # Get allele counts and coverage
    cou1 = np.load(get_allele_counts_filename(data_folder, adaID, frag1))
    cov1 = np.load(get_coverage_filename(data_folder, adaID, frag1))
    cou2 = np.load(get_allele_counts_filename(data_folder, adaID, frag2))
    cov2 = np.load(get_coverage_filename(data_folder, adaID, frag2))

    # Cut the counts and coverage to the overlap region
    cou1 = cou1[:, :, start_s2:]
    cov1 = cov1[:, start_s2:]
    cou2 = cou2[:, :, :end_s1]
    cov2 = cov2[:, :end_s1]

    # Reduce the allele counts (fwd and rev have different status on most overlaps,
    # because of the uneven coverage)
    nu1 = filter_nus(cou1, cov1)
    nu2 = filter_nus(cou2, cov2)

    # FIXME
    if nu1.shape != nu2.shape:
        return

    # Print table of called polymorphisms
    print 'adaID', adaID, frag1, frag2, 'polymorphism matrix (NO | YES)'
    print 3 * ' ', '|', '{:^10s}'.format(frag1)
    print 15 * '-'
    print 3 * ' ', '|', \
            '{:3d}'.format(((nu1 < 1e-6) & (nu2 < 1e-6)).sum()), '|', \
            '{:3d}'.format(((nu1 > 3e-3) & (nu2 < 1e-6)).sum())
    print '{:3s}'.format(frag2), '+'+(5*'-')+'+'+(4*'-')
    print 3 * ' ', '|', \
            '{:3d}'.format(((nu1 < 1e-6) & (nu2 > 3e-3)).sum()), '|', \
            '{:3d}'.format(((nu1 > 3e-6) & (nu2 > 3e-3)).sum())
         

    # Plot scatter
    import matplotlib.pyplot as plt
    from matplotlib import cm
    if ax is None:
        show_plot = True
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        ax.set_title('allele frequencies, adaID '+str(adaID)+', '+\
                     str(frag1)+' - '+str(frag2),
                     fontsize=18)
    else:
        show_plot = False
        ax.set_title(str(frag1)+' - '+str(frag2), fontsize=18)
    ax.scatter(np.abs(nu1 - 1e-5), np.abs(nu2 - 1e-5), s=30,
               c=cm.jet([int(255.0 * i / len(nu1)) for i in xrange(len(nu1))]))
    # Plot diagonal
    ax.plot([1e-7, 2], [1e-7, 2], lw=1, c='k', ls='--')
    ax.set_xlabel(r'$\nu_1$', fontsize=20)
    ax.set_ylabel(r'$\nu_2$', fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.7e-5, 1.2)
    ax.set_ylim(0.7e-5, 1.2)

    if show_plot:
        plt.tight_layout(w_pad=0.05)
        plt.ion()
        plt.show()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Check overlap regions')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--savefig', action='store_true',
                        help='Save the figures to file')
    
    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    savefig = args.savefig

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = MiSeq_runs[seq_run]['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    for adaID in adaIDs:
        # Find overlap pairs
        pairs = get_overlapping_fragments(fragments)

        # Iterate over overlaps
        overlaps = []
        for (frag1, frag2) in pairs:

            # Determine the overlap
            overlap = get_overlap(data_folder, adaID,
                                  frag1, frag2, VERBOSE=VERBOSE)

            # Check consensus
            is_diff = check_overlap_consensus(data_folder, adaID, frag1, frag2, overlap,
                                              VERBOSE=VERBOSE)

            if overlap is None:
                continue

            overlaps.append(((frag1, frag2), overlap))

        # Check allele frequencies if present
        if not os.path.isfile(get_allele_counts_filename(data_folder, adaID, 'F1')):
            continue

        # Make unified figures
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, len(overlaps), figsize=(2+3.6*len(fragments), 6))
        fig.suptitle('allele frequencies, adaID '+str(adaID), fontsize=18)
        if len(overlaps) == 1:
            axs = [axs]
        for (ax, ((frag1, frag2), overlap)) in izip(axs, overlaps):

            # Check allele frequencies
            check_overlap_allele_frequencies(data_folder, adaID, frag1, frag2,
                                             overlap, VERBOSE=VERBOSE, ax=ax)

        plt.tight_layout(rect=(0, 0, 1, 0.95))
        plt.ion()
        plt.show()

        # Save the figure
        if savefig:
            fig_filename = get_overlap_nu_figure_filename(data_folder, adaID,
                                                          '-'.join(fragments))
            # Create the folder if necessary
            dirname = os.path.dirname(fig_filename)
            if not os.path.isdir(dirname): os.mkdir(dirname)
            fig.savefig(fig_filename)

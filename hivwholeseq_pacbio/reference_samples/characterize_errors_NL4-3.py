# vim: fdm=marker
'''
author:     Fabio Zanini
date:       04/02/14
content:    Characterize errors in the NL4-3 test sample.
'''
# Modules
import os
import argparse
import numpy as np
from itertools import izip
from Bio import SeqIO
import pysam
import matplotlib.pyplot as plt

from hivwholeseq_pacbio.datasets import data_folder_dict
from hivwholeseq_pacbio.filenames import get_premapped_file, get_reference_premap_filename
from hivwholeseq_pacbio.pacbio_rs_II import alpha, alphal



# Functions
def classify_mismatch_deletions_single(data_folder, samplename,
                                       qual_window=[0, 80],
                                       maxreads=-1, VERBOSE=0):
    '''Relate errors, phred score, and position in the read'''

    # (x, y) indicates: "x expected, y read" 
    errd = {(a1, a2): 0 for a1 in alpha for a2 in alpha if (a1 not in ('-', 'N')) and (a1 != a2)}

    # Get reference
    refseq = SeqIO.read(get_reference_premap_filename(data_folder, samplename), 'fasta')
    refmat = np.array(refseq)

    input_filename = get_premapped_file(data_folder, samplename, type='bam')
    with pysam.Samfile(input_filename, 'rb') as bamfile:
        for irp, read in enumerate(bamfile):

            if irp == maxreads:
                if VERBOSE:
                    print 'Maximal number of reads reached:', maxreads
                break

            if VERBOSE >= 2:
                if not ((irp+1) % 100):
                    print (irp+1)

            # If unmapped or unpaired, discard
            if read.is_unmapped or (not len(read.cigar)):
                if VERBOSE >= 3:
                    print 'Read unmapped/no CIGAR:', read.qname,
                continue

            # Exclude obvious mismaps
            if sum(bl for (bt, bl) in read.cigar if bt in (1, 2)) > 100:
                if VERBOSE >= 3:
                    print 'Obvious mismap:', read.qname
                continue

            pos_ref = read.pos
            pos_read = 0
            for (bt, bl) in read.cigar:
                # Skip inserts
                if bt == 1:
                    pos_read += bl

                # Deletions
                elif bt == 2:
                    ref = refseq[pos_ref: pos_ref + bl]
                    for c in ref:
                        errd[(c, '-')] += 1

                    pos_ref += bl

                # Only mismatches are kept
                else:
                    seq = np.fromstring(read.seq[pos_read: pos_read + bl], 'S1')
                    qual = np.fromstring(read.qual[pos_read: pos_read + bl], np.int8) - 33
                    ref = refmat[pos_ref: pos_ref + bl]
                    errs = seq != ref
                    qual_win = (qual >= qual_window[0]) & (qual < qual_window[1])

                    for trans in izip(ref[errs & qual_win], seq[errs & qual_win]):
                        errd[trans] += 1
                                       
                    pos_read += bl
                    pos_ref += bl

    return errd


def plot_mismatch_deletions_single(errd, title='', only_mismatches=False, VERBOSE=0):
    '''Plot the error transition matrix'''
    import matplotlib.pyplot as plt

    mat = np.zeros((len(alpha) - 2, len(alpha)))
    for (c1, c2), val in errd.iteritems():
        mat[alphal.index(c1), alphal.index(c2)] = val
    mat = np.log10(mat + 0.1)
    if only_mismatches:
        mat = mat[:, :-2]

    fig, ax = plt.subplots()
    ims = ax.imshow(mat, interpolation='nearest')
    ax.set_ylabel('From')
    ax.set_xlabel('To')
    ax.set_title(title)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(alphal[:mat.shape[1]])
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(alphal[:mat.shape[0]])

    plt.colorbar(mappable=ims, ax=ax)

    plt.tight_layout()

    plt.ion()
    plt.show()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Trim and divide reads into fragments')
    parser.add_argument('--verbose', default=0, type=int,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--qual_window', type=int, nargs=2, default=(0, 80),
                        help='Restrict to a phred quality window')
    parser.add_argument('--only_mismatches', action='store_true',
                        help='Ignore deletions')

    args = parser.parse_args()
    VERBOSE = args.verbose
    maxreads = args.maxreads
    qual_window = args.qual_window
    only_mismatches = args.only_mismatches

    # Specify the dataset
    seq_run = 'Upp23'
    samplename = 'S1'
    data_folder = data_folder_dict[seq_run]

    # Classify error
    errd = classify_mismatch_deletions_single(data_folder, samplename,
                                              qual_window=qual_window,
                                              maxreads=maxreads, VERBOSE=VERBOSE)

    plot_mismatch_deletions_single(errd, title='PacBio error matrix, NL4-3\nquality window '+\
                                   str(qual_window),
                                   only_mismatches=only_mismatches)


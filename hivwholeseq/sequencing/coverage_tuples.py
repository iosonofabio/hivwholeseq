# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/08/13
content:    Calculate the number of reads that cover a certain set of positions.
            This is used to normalize linkage analyses. For instance, suppose we
            find the pair of mutations A456T and G642T 10 times, how many reads
            actually covered both such that they could have given rise to such a
            signal?
'''
# Modules
import os
import argparse
from itertools import izip
import pysam
import numpy as np
from Bio import SeqIO

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.adapter_info import load_adapter_table
from hivwholeseq.sequencing.filenames import get_mapped_filename
from hivwholeseq.utils.mapping import get_ind_good_cigars, pair_generator



# Globals



# Functions
def format_tuples(tuples):
    '''Format input arguments for tules'''
    try:
        # Note: the first element is the fragment
        mtuples = map(lambda x: np.array(x.split(' '), int), tuples)
    except ValueError:
        raise ValueError('Tuple format not recognized. '\
                         "Example: "+os.path.basename(__file__)+\
                         " '2 56 78' '5 10 20 30'")
    return mtuples


from numba import autojit, int_, uint
@autojit(locals={"l": uint, "lm": uint, "mt": uint})
def add_read(ref_start, ref_end, mtuples, covs_pair):

    l = len(mtuples)
    for i in np.arange(l):
        mtuple = mtuples[i]
        cov_pair = covs_pair[i]
        lm = len(mtuple)
        for j in np.arange(lm):
            mt = mtuple[j]
            # Implement comparisons as div between integers
            if (mt >= ref_start) and (mt < ref_end):
                cov_pair[j] = True


def get_coverage_tuples(data_folder, adaID, fragment, mtuples,
                       maxreads=-1, VERBOSE=0):
    '''Get the joint coverage of a list of positions'''
    # Prepare data structures
    mtuples = [np.asarray(tup, int) for tup in mtuples]
    coverage = np.zeros(len(mtuples), int)

    # TODO: what to do if it is covered multiple times? or only some sites?
    covs_pair = [np.zeros(len(tup), bool) for tup in mtuples]

    # Open BAM
    bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                      filtered=True)
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)
    with pysam.Samfile(bamfilename, 'rb') as bamfile:
    
        # Iterate over all pairs
        for irp, reads in enumerate(pair_generator(bamfile)):

            # Limit to the first reads
            if irp == maxreads:
                if VERBOSE:
                    print 'Max reads reached:', maxreads
                break

            if VERBOSE >= 3:
                if not ((irp + 1) % 10000):
                    print irp + 1

            # Reinitialize temporary structure
            for cov_pair in covs_pair: cov_pair[:] = False

            # Look in both reads
            for read in reads:

                # NOTE: deletions count as covered, because in principle
                # we see that part of the reference
                cigar = read.cigar
                ref_start = read.pos
                ref_end = ref_start + sum(bl for (bt, bl) in cigar if bt in (0, 2))

                # Use numba to accelerate? better not
                if False:
                    add_read(ref_start, ref_end, mtuples, covs_pair)
                else:
                    for cov_pair, mtuple in izip(covs_pair, mtuples):
                        cov_pair[(mtuple >= ref_start) & (mtuple < ref_end)] = True

            # Check which tuples are fully covered
            for i, cov_pair in enumerate(covs_pair):
                if cov_pair.all():
                    coverage[i] += 1

    return coverage 



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Get coverage of multiple sites')
    parser.add_argument('--run', type=int, required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaID', type=int, required=True,
                        help='Adapter ID to analyze (e.g. 2)')
    parser.add_argument('--fragment', required=True,
                        help='Fragment to analyze (e.g. F1)')
    parser.add_argument('-n', type=int, default=-1,
                        help='Number of read pairs to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    # Example: coverage_tuples.py --adaID TS2 '5 30 56' '58 43 89'
    parser.add_argument('--tuples', nargs='+', required=True,
                        metavar=("'pos1 pos2 ...'", "'pos3 pos4 ...'"),
                        help="Tuples to analyze (e.g. '3 45 98' '54 62'")

    args = parser.parse_args()
    seq_run = args.run
    adaID = args.adaID
    fragment = args.fragment
    n_pairs = args.n
    VERBOSE = args.verbose
    mtuples = format_tuples(args.tuples)

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # Get the coverage
    coverage = get_coverage_tuples(data_folder, adaID, fragment, mtuples,
                                   maxreads=n_pairs, VERBOSE=VERBOSE)

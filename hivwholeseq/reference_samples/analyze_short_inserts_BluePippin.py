# vim: fdm=marker
'''
author:     Fabio Zanini
date:       10/02/14
content:    Despite the BluePippin size selection there seem to be a few short
            inserts left in some samples in Tue48, e.g. N2-S4 (p3). What is in there?
'''
# Modules
import os
import argparse
from collections import Counter, defaultdict
from itertools import izip
import numpy as np
import matplotlib.pyplot as plt
import pysam

from hivwholeseq.mapping_utils import pair_generator
from hivwholeseq.filenames import get_divided_filenames
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.samples import samples



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    args = parser.parse_args()
    VERBOSE = args.verbose
    maxreads = args.maxreads
    
    # Settings
    seq_run = 'Tue48'

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # p3
    adaID = 'N2-S4'

    fns = get_divided_filenames(data_folder, adaID, ['F1'], type='bam')
    unmapped_fn = filter(lambda x: 'unmapped' in x, fns)[0]

    # Scan the BAM file
    with pysam.Samfile(unmapped_fn, 'rb') as bamfile:
        for irp, reads in enumerate(pair_generator(bamfile)):
            if irp == maxreads:
                break

            if reads[0].is_unmapped or (not reads[0].is_proper_pair):
                if VERBOSE >=3:
                    print 'Unmapped', reads[0].qname
                continue

            print reads[0].is_reverse, reads[1].is_reverse

            i_fwd = reads[0].is_reverse
            i_rev = -i_fwd
            print reads[i_fwd].isize, reads[0].cigar
            print reads[i_rev].isize, reads[1].cigar
            print reads[i_fwd].seq
            print reads[i_rev].seq
            print

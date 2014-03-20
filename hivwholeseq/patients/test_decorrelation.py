# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/03/14
content:    Test script to measure decorrelation (recombination): pick pairs of
            alleles that are strongly correlated at the first time point, and see
            how fast the correlation disappears - if at all.
'''
# Modules
import sys
import os
import argparse
import numpy as np
import pysam
from Bio import SeqIO

from hivwholeseq.patients.patients import get_patient
from hivwholeseq.filter_mapped_reads import plot_distance_histogram, \
        plot_distance_histogram_sliding_window, get_distance_from_consensus, \
        check_overhanging_reads, trim_bad_cigar
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_mapped_to_initial_filename, get_filter_mapped_init_summary_filename
from hivwholeseq.mapping_utils import convert_sam_to_bam, pair_generator
from hivwholeseq.two_site_statistics import get_coallele_counts_from_file




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Filter mapped reads')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--tests', action='store_true',
                        help='Include consistency tests')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    maxreads = args.maxreads
    use_tests = args.tests

    patient = get_patient(pname)

    # If the script is called with no fragment, iterate over all
    # NOTE: This assumes that the fragment has been sequenced. If not, we should
    # take a later time point or similia, but let's not worry about it yet.
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Get the first sequenced sample
    sample_init = patient.initial_sample

    for fragment in fragments:

        reffilename = get_initial_consensus_filename(pname, fragment)
        refseq = SeqIO.read(reffilename, 'fasta')
        length = len(refseq)

        # FIXME: repeat for later time points
        sample = sample_init
        samplename = sample.name

        bamfilename = get_mapped_to_initial_filename(pname, samplename, fragment,
                                                     type='bam', filtered=True)

        cocounts = get_coallele_counts_from_file(bamfilename, length, 
                                                 maxreads=maxreads,
                                                 VERBOSE=VERBOSE,
                                                 use_tests=use_tests)

#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/03/14
content:    Get the joint counts at two sites for patient samples, after mapping.
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
        get_mapped_to_initial_filename, get_filter_mapped_init_summary_filename, \
        get_allele_cocounts_filename
from hivwholeseq.mapping_utils import convert_sam_to_bam, pair_generator
from hivwholeseq.two_site_statistics import get_coallele_counts_from_file
from hivwholeseq.fork_cluster import fork_get_cocounts_patient as fork_self



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Filter mapped reads')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--samples', nargs='*',
                        help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--tests', action='store_true',
                        help='Include consistency tests')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    samples = args.samples
    VERBOSE = args.verbose
    maxreads = args.maxreads
    use_tests = args.tests
    submit = args.submit

    patient = get_patient(pname)

    # If no samples are mentioned, use all sequenced ones
    if not samples:
        samples = patient.samples
    if VERBOSE >= 3:
        print 'samples', samples

    # If the script is called with no fragment, iterate over all
    # NOTE: This assumes that the fragment has been sequenced. If not, we should
    # take a later time point or similia, but let's not worry about it yet.
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    if submit:
        for fragment in fragments:
            for samplename in samples:
                fork_self(pname, samplename, fragment, VERBOSE=VERBOSE,
                          maxreads=maxreads, use_tests=use_tests)
        sys.exit()

    for fragment in fragments:
        reffilename = get_initial_consensus_filename(pname, fragment)
        refseq = SeqIO.read(reffilename, 'fasta')
        length = len(refseq)

        for samplename in samples:
            if VERBOSE >= 2:
                print pname, fragment, samplename

            bamfilename = get_mapped_to_initial_filename(pname, samplename, fragment,
                                                         type='bam', filtered=True)
            cocounts = get_coallele_counts_from_file(bamfilename, length, 
                                                     maxreads=maxreads,
                                                     VERBOSE=VERBOSE,
                                                     use_tests=use_tests)

            if VERBOSE >= 2:
                print 'Storing to file'
            cocounts.dump(get_allele_cocounts_filename(pname, samplename, fragment))

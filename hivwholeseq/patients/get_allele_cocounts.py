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

from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.sequencing.filter_mapped_reads import plot_distance_histogram, \
        plot_distance_histogram_sliding_window, get_distance_from_consensus, \
        check_overhanging_reads, trim_bad_cigar
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_filter_mapped_init_summary_filename, \
        get_allele_cocounts_filename
from hivwholeseq.mapping_utils import convert_sam_to_bam, pair_generator
from hivwholeseq.two_site_statistics import get_coallele_counts_from_file as gac
from hivwholeseq.cluster.fork_cluster import fork_get_cocounts_patient as fork_self



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele cocounts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele cocounts to file')
    parser.add_argument('--tests', action='store_true',
                        help='Include consistency tests')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--qualmin', type=int, default=30,
                        help='Minimal quality of base to call')
    parser.add_argument('--PCR', type=int, default=1,
                        help='Analyze only reads from this PCR (1 or 2)')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    VERBOSE = args.verbose
    maxreads = args.maxreads
    use_tests = args.tests
    submit = args.submit
    save_to_file = args.save
    qual_min = args.qualmin
    PCR = args.PCR

    samples = lssp()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    if submit:
        for fragment in fragments:
            for samplename, sample in samples.iterrows():
                fork_self(samplename, fragment, VERBOSE=VERBOSE,
                          qual_min=qual_min, PCR=PCR,
                          maxreads=maxreads, use_tests=use_tests)
        sys.exit()

    counts_all = []
    for fragment in fragments:
        counts = []
        for samplename, sample in samples.iterrows():
            sample = SamplePat(sample)
            pname = sample.patient

            if VERBOSE >= 2:
                print pname, fragment, samplename

            refseq = SeqIO.read(get_initial_reference_filename(pname, fragment), 'fasta')

            fn_out = sample.get_allele_cocounts_filename(fragment, PCR=PCR, qual_min=qual_min)
            fn = sample.get_mapped_filtered_filename(fragment, PCR=PCR, decontaminated=True) #FIXME
            if save_to_file:
                cocount = gac(fn, len(refseq), 
                              maxreads=maxreads,
                              VERBOSE=VERBOSE,
                              qual_min=qual_min,
                              use_tests=use_tests)

                cocount.dump(fn_out)

                if VERBOSE >= 2:
                    print 'Allele cocounts saved:', samplename, fragment

                counts.append(cocount)

            elif os.path.isfile(fn_out):
                cocount = np.load(fn_out)
                counts.append(cocount)

            elif os.path.isfile(fn):
                cocount = gac(fn, len(refseq), 
                              maxreads=maxreads,
                              VERBOSE=VERBOSE,
                              qual_min=qual_min,
                              use_tests=use_tests)
                counts.append(cocount)

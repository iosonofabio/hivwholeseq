#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       25/10/13
content:    Paste mapped reads in chunks together again.
'''
# Modules
import os
import time
import argparse
import subprocess as sp
import numpy as np
import pysam
import warnings

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.filenames import get_divided_filename
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.mapping_utils import stampy_bin, subsrate, \
        convert_sam_to_bam, convert_bam_to_sam, get_number_reads
from hivwholeseq.samples import samples as samples_seq
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.patients.filenames import get_initial_index_filename, \
        get_initial_hash_filename, get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_map_initial_summary_filename, \
        get_paste_mapped_chunks_initial_summary_filename
from hivwholeseq.fork_cluster import fork_paste_mapped_chunks_to_initial_consensus as fork_self
from hivwholeseq.clean_temp_files import remove_mapped_init_tempfiles
from hivwholeseq.patients.filter_mapped_reads import filter_mapped_reads



# Functions
def paste_chunks(patient, samplename, fragment, VERBOSE=0):
    '''Paste the mapped chunks together again'''
    output_filename = get_mapped_to_initial_filename(pname, samplename,
                                                     fragment,
                                                     type='bam')

    input_filename = get_mapped_to_initial_filename(pname, samplename,
                                                    fragment,
                                                    type='bam',
                                                    only_chunk=1)
    with pysam.Samfile(input_filename, 'rb') as bamfilein1:
        with pysam.Samfile(output_filename, 'wb', template=bamfilein1) as bamfileout:
            if VERBOSE >= 2:
                print 'Chunk 1'
            for read in bamfilein1:
                bamfileout.write(read)
            input_filename = get_mapped_to_initial_filename(pname, samplename,
                                                            fragment,
                                                            type='bam',
                                                            only_chunk=2)
            chunk_number = 2
            while os.path.isfile(input_filename):
                with pysam.Samfile(input_filename, 'rb') as bamfilein:
                    if VERBOSE >= 2:
                        print 'Chunk', chunk_number
                    for read in bamfilein:
                        bamfileout.write(read)

                chunk_number += 1
                input_filename = get_mapped_to_initial_filename(pname, samplename,
                                                                fragment,
                                                                type='bam',
                                                                only_chunk=chunk_number)

    if VERBOSE >= 1:
        print 'Remove temporary files: sample '+samplename, 'fragment', fragment
    remove_mapped_init_tempfiles(pname, samplename, fragment, VERBOSE=VERBOSE, remove_chunks=True)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Map to initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--samples', nargs='+',
                        help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--filter', action='store_true',
                        help='Filter reads immediately after mapping')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')

    args = parser.parse_args()
    pname = args.patient
    samples = args.samples
    fragments = args.fragments
    submit = args.submit
    VERBOSE = args.verbose
    filter_reads = args.filter
    summary = args.summary

    patient = get_patient(pname)

    # If no samples are mentioned, use all sequenced ones
    if not samples:
        samples = patient.samples
    if VERBOSE >= 3:
        print 'samples', samples

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    sample_init = patient.initial_sample

    for samplename in samples:
        for fragment in fragments:
            if VERBOSE:
                print pname, samplename, fragment
            if submit:
                fork_self(pname, samplename, fragment,
                          VERBOSE=VERBOSE, threads=threads,
                          n_pairs=n_pairs,
                          filter_reads=filter_reads,
                          summary=summary,
                          only_chunks=[only_chunk])
                continue

            if summary:
                sfn = get_paste_mapped_chunks_initial_summary_filename(pname, samplename, fragment)
                with open(sfn, 'w') as f:
                    f.write('Call: python paste_mapped_chunks.py'+\
                            ' --patient '+pname+\
                            ' --samples '+samplename+\
                            ' --fragments '+fragment+\
                            ' --verbose '+str(VERBOSE))
                    f.write('\n')

            paste_chunks(patient, samplename, fragment, VERBOSE=VERBOSE)

            if filter_reads:
                if summary:
                    sfn = get_filter_mapped_init_summary_filename(pname, samplename, fragment)
                    with open(sfn, 'w') as f:
                        f.write('Call: python filter_mapped_reads.py'+\
                                ' --patient '+pname+\
                                ' --samples '+samplename+\
                                ' --fragments '+fragment+\
                                ' --verbose '+str(VERBOSE))
                        f.write('\n')

                filter_mapped_reads(pname, samplename, fragment,
                                    VERBOSE=VERBOSE, summary=summary)




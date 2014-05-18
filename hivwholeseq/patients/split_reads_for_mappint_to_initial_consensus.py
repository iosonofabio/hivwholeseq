#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       18/05/14
content:    Split reads for mapping, from a patient by patient view.
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
        get_initial_hash_filename, get_initial_consensus_filename, \
        get_mapped_to_initial_filename, get_map_initial_summary_filename
from hivwholeseq.fork_cluster import fork_split_for_mapping as fork_split
from hivwholeseq.clean_temp_files import remove_mapped_init_tempfiles
from hivwholeseq.split_reads_for_mapping import split_reads


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
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--chunksize', type=int, default=10000,
                        help='Size of the chunks, has to fit 1 hr cluster time')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    samples = args.samples
    fragments = args.fragments
    submit = args.submit
    VERBOSE = args.verbose
    maxreads = args.maxreads
    chunksize = args.chunksize

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

    for samplename in samples:
        sample = patient.sample_table.loc[samplename]
        
        seq_run = sample['run']
        dataset = MiSeq_runs[seq_run]
        data_folder = dataset['folder']
        adaID = sample['adaID']
        fragments_sample = sample['fragments']

        for fragment in fragments:
            fragments_potential = filter(lambda x: fragment in x, fragments_sample)
            if len(fragments_potential) < 1:
                if VERBOSE:
                    print 'Fragment', fragment, 'not found in sample', samplename, 'of patient', pname
                continue

            if VERBOSE:
                print 'patient', pname, 'sample', samplename, 'fragment', fragment

            fragment_sample = fragments_potential[0]

            if submit:
                fork_split(seq_run, adaID, fragment_sample,
                           VERBOSE=VERBOSE,
                           maxreads=maxreads, chunk_size=chunksize)
                continue

            split_reads(data_folder, adaID, fragment_sample, chunk_size=chunksize,
                        maxreads=maxreads, VERBOSE=VERBOSE)


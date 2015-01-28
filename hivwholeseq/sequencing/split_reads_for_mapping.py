#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/05/14
content:    Split the reads for mapping into small aliquots, to exploit the
            short queue on the cluster.
'''
# Modules
import os
import sys
import time
import argparse
import pysam
import numpy as np
import subprocess as sp

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.adapter_info import load_adapter_table, foldername_adapter
from hivwholeseq.utils.mapping import stampy_bin, subsrate, bwa_bin, convert_sam_to_bam, \
        convert_bam_to_sam, get_number_reads, get_number_reads_open ,pair_generator
from hivwholeseq.sequencing.filenames import get_consensus_filename, get_mapped_filename,\
        get_read_filenames, get_divided_filename, get_map_summary_filename
from hivwholeseq.sequencing.filter_mapped_reads import match_len_min, trim_bad_cigars
from hivwholeseq.sequencing.filter_mapped_reads import filter_reads as filter_mapped_reads
from hivwholeseq.cluster.fork_cluster import fork_split_for_mapping as fork_self
from hivwholeseq.sequencing.samples import samples
from hivwholeseq.clean_temp_files import remove_mapped_tempfiles


# Functions
def split_reads(data_folder, adaID, fragment, chunk_size=10000, maxreads=-1, VERBOSE=0):
    '''Split reads into chunks for mapping'''

    input_filename = get_divided_filename(data_folder, adaID, fragment, type='bam')
    with pysam.Samfile(input_filename, 'rb') as bamfile:
        if VERBOSE:
            if maxreads == -1:
                n_reads = get_number_reads_open(bamfile) // 2
            else:
                n_reads = maxreads

            print 'Expected number of chunks:', 1 + (n_reads // chunk_size)

        chunk_number = 0
        chunkfile = None
        for irp, read_pair in enumerate(pair_generator(bamfile)):
            if irp == maxreads:
                break

            if VERBOSE >= 2:
                if not ((irp+1) % 10000):
                    print irp+1

            if not (irp % chunk_size):
                if chunkfile is not None:
                    chunkfile.close()
                chunk_number += 1
                chunk_filename = get_divided_filename(data_folder, adaID, fragment, type='bam', chunk=chunk_number)
                chunkfile = pysam.Samfile(chunk_filename, 'wb', template=bamfile)
                if VERBOSE >= 2:
                    print 'Chunk n', chunk_number, 'started'


            chunkfile.write(read_pair[0])
            chunkfile.write(read_pair[1])

        if chunkfile is not None:
            chunkfile.close()

    if VERBOSE:
        print 'Chunking finished'



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--chunksize', type=int, default=10000,
                        help='Size of the chunks, has to fit 1 hr cluster time')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    submit = args.submit
    maxreads = args.maxreads
    chunksize = args.chunksize

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = dataset['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Iterate over all requested samples
    for adaID in adaIDs:

        # If the script is called with no fragment, iterate over all
        samplename = dataset['samples'][dataset['adapters'].index(adaID)]
        if not fragments:
            fragments_sample = samples[samplename]['fragments']
        else:
            from re import findall
            fragments_all = samples[samplename]['fragments']
            fragments_sample = []
            for fragment in fragments:
                frs = filter(lambda x: fragment in x, fragments_all)
                if len(frs):
                    fragments_sample.append(frs[0])

        if VERBOSE >= 3:
            print 'adaID '+adaID+': fragments '+' '.join(fragments_sample)

        # Iterate over fragments
        for fragment in fragments_sample:

            if submit:
                fork_self(seq_run, adaID, fragment,
                          VERBOSE=VERBOSE,
                          maxreads=maxreads, chunk_size=chunksize)
                continue


            split_reads(data_folder, adaID, fragment, chunk_size=chunksize,
                        maxreads=maxreads, VERBOSE=VERBOSE)

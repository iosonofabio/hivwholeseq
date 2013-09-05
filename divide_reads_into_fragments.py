#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       28/08/13
content:    Classify the reads into the 6 fragments via a rough mapping before
            the final iterative mapping. This script reduces the chances of
            unpaired reads, unmapped, overhangs at fragments, and so on.

            This mapping is supposed to be fast and is a single run of stampy,
            so we run BWA before that to speed up.
'''
# Modules
import os
import sys
import subprocess as sp
import time
import argparse
import re
from operator import *
from itertools import izip
from collections import defaultdict
from collections import Counter
import numpy as np
import pysam
import Bio.SeqIO as SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import reverse_complement

# Horizontal import of modules from this folder
from mapping.adapter_info import load_adapter_table, foldername_adapter
from mapping.miseq import alpha, read_types, pair_generator
from mapping.mapping_utils import stampy_bin, subsrate, bwa_bin
from mapping.sequence_utils.annotate_HXB2 import load_HXB2
from mapping.filenames import get_HXB2_entire, get_HXB2_index_file, get_HXB2_hash_file
from mapping.filenames import get_read_filenames, get_premapped_file
from mapping.mapping_utils import convert_sam_to_bam, get_read_start_end
from mapping.mapping_utils import get_range_good_cigars



# Globals
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']

# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'divide_reads_into_fragments.py'
cluster_time = '0:59:59'
vmem = '8G'



# Functions
def make_index_and_hash(data_folder, cropped=True, VERBOSE=0):
    '''Make index and hash files for reference or consensus'''
    if VERBOSE:
        print 'Making index and hash files'

    # 1. Make genome index file for HXB2
    if not os.path.isfile(get_HXB2_index_file(data_folder, ext=True)):
        sp.call([stampy_bin,
                 '--species="HIV"',
                 '-G', get_HXB2_index_file(data_folder, ext=False),
                 get_HXB2_entire(data_folder, cropped=cropped),
                 ])
    
    # 2. Build a hash file for HXB2
    if not os.path.isfile(get_HXB2_hash_file(data_folder, ext=True)):
        sp.call([stampy_bin,
                 '-g', get_HXB2_index_file(data_folder, ext=False),
                 '-H', get_HXB2_hash_file(data_folder, ext=False),
                 ])


def make_bwa_hash(data_folder, cropped=True, VERBOSE=0):
    '''Make hash tables for BWA'''
    if VERBOSE:
        print 'Build BWA index'

    # Make folder if necessary
    dirname =  os.path.dirname(get_HXB2_hash_file(data_folder))
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # Call bwa index
    ref_filename = get_HXB2_entire(data_folder, cropped=cropped)
    sp.call([bwa_bin,
             'index',
             ref_filename,
             ])

    # Move the hashes into subfolder
    from glob import glob
    import shutil
    shutil.copy(ref_filename, dirname)
    hash_files = glob(ref_filename+'.*')
    for hash_file in hash_files:
        shutil.copy(hash_file, dirname)
        os.remove(hash_file)


def premap_bwa(data_folder, adaID, VERBOSE=0, subsample=False, cropped=True):
    '''Premap using BWA'''
    if VERBOSE:
        print 'Map via BWA: '+'{:02d}'.format(adaID)

    # Get input filenames
    index_prefix = os.path.dirname(get_HXB2_hash_file(data_folder, ext=False))+\
            '/'+os.path.basename(get_HXB2_entire(data_folder, cropped=cropped))
            
    readfiles = get_read_filenames(data_folder, adaID,
                                   subsample=subsample,
                                   filtered=True)

    # Get output filename
    output_filename = get_premapped_file(data_folder, adaID,
                                         type='sam', subsample=subsample,
                                         bwa=True)
    # Make folder if necessary
    dirname = os.path.dirname(output_filename)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # Map
    with open(output_filename, 'w') as f:
        call_list = [bwa_bin,
                     'mem',
                     index_prefix] + readfiles
        if VERBOSE >=2:
            print ' '.join(call_list)
        sp.call(call_list, stdout=f)


def premap_stampy(data_folder, adaID, VERBOSE=0, subsample=False, bwa=False):
    '''Call stampy for actual mapping'''
    if VERBOSE:
        print 'Map: adaID '+'{:02d}'.format(adaID)

    # Get input filenames
    if bwa:
        bwa_filename = get_premapped_file(data_folder, adaID,
                                        type='bam', subsample=subsample,
                                        bwa=True)
        if not os.path.isfile(bwa_filename):
            if VERBOSE >= 2:
                print 'Converting SAM to BAM'
            convert_sam_to_bam(bwa_filename)
        readfiles = [bwa_filename]

    else:
        readfiles = get_read_filenames(data_folder, adaID,
                                       subsample=subsample,
                                       filtered=True)
    # Get output filename
    output_filename =  get_premapped_file(data_folder, adaID, type='sam',
                                          subsample=subsample)
    # Make folder if necessary
    dirname =  os.path.dirname(output_filename)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    
    # Map
    # Note: no --solexa option as of 2013 (illumina v1.8)
    qsub_list = [stampy_bin,
                 '-g', get_HXB2_index_file(data_folder, ext=False),
                 '-h', get_HXB2_hash_file(data_folder, ext=False), 
                 '-o', output_filename,
                 '--substitutionrate='+subsrate]
    if bwa:
        qsub_list.append('--bamkeepgoodreads')
    qsub_list = qsub_list + ['-M'] + readfiles
    qsub_list = map(str,qsub_list)
    if VERBOSE >= 2:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def fork_self(data_folder, adaID, VERBOSE=0, subsample=False):
    '''Fork self for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+'{:02d}'.format(adaID)

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'premap_'+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--adaID', adaID,
                 '--verbose', VERBOSE,
                ]
    if subsample:
        qsub_list.append('--subsample')
    qsub_list = map(str, qsub_list)
    if VERBOSE >= 2:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def write_read_pair(reads, ranges, fileouts):
    '''Write read pair into FASTQ files'''

    (read1, read2) = reads

    # Find orientation
    if not reads[0].is_reverse:
        fq1 = read1.seq[ranges[0][0]: ranges[0][1]]
        qual1 = read1.qual[ranges[0][0]: ranges[0][1]]
        fq2 = reverse_complement(read2.seq[ranges[1][0]: ranges[1][1]])
        qual2 = read2.qual[ranges[1][0]: ranges[1][1]][::-1]
    else:
        fq1 = reverse_complement(read1.seq[ranges[0][0]: ranges[0][1]])
        qual1 = read1.qual[ranges[0][0]: ranges[0][1]][::-1]
        fq2 = read2.seq[ranges[1][0]: ranges[1][1]]
        qual2 = read2.qual[ranges[1][0]: ranges[1][1]]

    # Note: we lose the x:N:0 part of the qname after the
    # whitespace, but that's used for barcoding only
    fileouts[0].write("@%s\n%s\n+\n%s\n" % (reads[0].qname, fq1, qual1))
    fileouts[1].write("@%s\n%s\n+\n%s\n" % (reads[1].qname, fq2, qual2))


def divide_into_fragments(data_folder, adaID, VERBOSE=0, subsample=False):
    '''Divide mapped reads into fastq files for the fragments'''
    if VERBOSE:
        print 'Divide into fragments: adaID '+'{:02d}'.format(adaID)

    # Load annotated reference (via GIT submodule)
    seq = load_HXB2()

    # Extract fragments
    frags_pos = [(f.location.nofuzzy_start, f.location.nofuzzy_end)
                 for f in seq.features if f.type[:8] == 'fragment']
    frags_pos = np.array(frags_pos, int).T
    # Since the reference is cropped, subtract from the positions F1 start
    frags_pos -= frags_pos.min()
    # Note: now we are in the reference of the CROPPED HXB2!

    # Get the premapped reads
    bamfilename = get_premapped_file(data_folder, adaID, type='bam',
                                     subsample=subsample)
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)

    # Output fastq files
    fileouts = get_read_filenames(data_folder, adaID,
                                  subsample=subsample,
                                  premapped=True)
    
    # Iterate over the mapped reads and assign fragments
    fragments = []
    n_unmapped = 0
    with pysam.Samfile(bamfilename, 'rb') as bamfile,\
         open(fileouts[0][0], 'w') as fo1_F1,\
         open(fileouts[0][1], 'w') as fo1_F2,\
         open(fileouts[0][2], 'w') as fo1_F3,\
         open(fileouts[0][3], 'w') as fo1_F4,\
         open(fileouts[0][4], 'w') as fo1_F5,\
         open(fileouts[0][5], 'w') as fo1_F6,\
         open(fileouts[0][6], 'w') as fo1_um,\
         open(fileouts[1][0], 'w') as fo2_F1,\
         open(fileouts[1][1], 'w') as fo2_F2,\
         open(fileouts[1][2], 'w') as fo2_F3,\
         open(fileouts[1][3], 'w') as fo2_F4,\
         open(fileouts[1][4], 'w') as fo2_F5,\
         open(fileouts[1][5], 'w') as fo2_F6,\
         open(fileouts[1][6], 'w') as fo2_um:

        # Collect the file handles
        file_handles = ((fo1_F1, fo1_F2, fo1_F3, fo1_F4, fo1_F5, fo1_F6, fo1_um),
                        (fo2_F1, fo2_F2, fo2_F3, fo2_F4, fo2_F5, fo2_F6, fo2_um))

        for reads in pair_generator(bamfile):

            # Set read1 and read2
            (read1, read2) = reads

            # If it is unmapped or unpaired, keep in dumpyard
            if (read1.is_unmapped or read2.is_unmapped or
                (not read1.is_proper_pair) or (not read2.is_proper_pair)):
                n_unmapped += 1
                write_read_pair(reads,
                                ((0, read1.rlen), (0, read2.rlen)),
                                (file_handles[0][6],
                                 file_handles[1][6]))
                continue

            # Get the extremes of the good CIGARs (we do some trimming to
            # avoid reading back into the adapters and stuff)
            ranges_read = []
            ranges_ref = []
            for read in reads:
                (rread, rref) = get_range_good_cigars(read.cigar, read.pos)
                ranges_read.append(rread)
                ranges_ref.append(rref)

            # If no good CIGARs found in one of the reads, dump the pair
            if None in ranges_ref:
                n_unmapped += 1
                write_read_pair(reads,
                                ((0, read1.rlen), (0, read2.rlen)),
                                (file_handles[0][6],
                                 file_handles[1][6]))
                continue

            # Is the good CIGAR reagion is mapped within a single fragment?
            frags_reads = []
            for loc in ranges_ref:
                ind = ((loc[0] >= frags_pos[0]) &
                       (loc[0] < frags_pos[1]) &
                       (loc[1] > frags_pos[0]) &
                       (loc[1] <= frags_pos[1]))
                frags_reads.append(ind.nonzero()[0])
            
            # Find the fragment(s) good for the whole pair
            frags_pair = np.intersect1d(*frags_reads)
            
            # A few cases can happen
            # 1. If the intersection is a single fragment, good
            if len(frags_pair) == 1:
                n_frag = frags_pair[0]
                write_read_pair(reads,
                                ranges_read,
                                (file_handles[0][n_frag],
                                 file_handles[1][n_frag]))

            # 2. If 2+ fragments are possible (tie), pick one at random
            elif len(frags_pair) > 1:
                n_frag = np.random.choice(frags_pair)
                write_read_pair(reads,
                                ranges_read,
                                (file_handles[0][n_frag],
                                 file_handles[1][n_frag]))
            
            # 3. If no fragments are possible (e.g. one read crosses the
            # fragment boundary, they map to different fragments), dump it
            else:
                n_unmapped += 1
                write_read_pair(reads,
                                ((0, read1.rlen), (0, read2.rlen)),
                                (file_handles[0][6],
                                 file_handles[1][6]))


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--adaID', metavar='00', type=int, nargs='?',
                        default=0,
                        help='Adapter ID')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--subsample', action='store_true',
                        help='Apply only to a subsample of the reads')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--bwa', action='store_true',
                        help='Use BWA for premapping?')

    args = parser.parse_args()
    adaID = args.adaID
    VERBOSE = args.verbose
    subsample = args.subsample
    submit = args.submit
    bwa = args.bwa

    # If the script is called with no adaID, iterate over all
    if adaID == 0:
        adaIDs = load_adapter_table(data_folder)['ID']
    else:
        adaIDs = [adaID]

    # Make index and hash files for HXB2, using a cropped reference to
    # reduce the problems of LTRs
    if bwa:
        make_bwa_hash(data_folder, cropped=True, VERBOSE=VERBOSE)
    make_index_and_hash(data_folder, cropped=True, VERBOSE=VERBOSE)

    # Iterate over all adaIDs
    for adaID in adaIDs:

        # Submit to the cluster self if requested
        if submit:
            fork_self(data_folder, adaID, VERBOSE=VERBOSE, subsample=subsample)

        # or else, divide into fragments
        else:

            # Map roughly to HXB2
            if bwa:
                premap_bwa(data_folder, adaID,
                           VERBOSE=VERBOSE, subsample=subsample)
            premap_stampy(data_folder, adaID,
                          VERBOSE=VERBOSE, subsample=subsample, bwa=bwa)

            # Divide into fragments and unclassified
            divide_into_fragments(data_folder, adaID,
                                  VERBOSE=VERBOSE, subsample=subsample)



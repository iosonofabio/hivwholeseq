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
from mapping.mapping_utils import stampy_bin, subsrate
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
interval_check = 10
cluster_time = '0:59:59'
vmem = '8G'



# Functions
def make_index_and_hash(data_folder, cropped=False):
    '''Make index and hash files for reference or consensus'''
    # 1. Make genome index file for HXB2
    if not os.path.isfile(get_HXB2_index_file(data_folder, ext=True)):
        sp.call([stampy_bin,
                 '--species="HIV 6 fragments"',
                 '-G', get_HXB2_index_file(data_folder, ext=False),
                 get_HXB2_entire(data_folder, cropped=cropped),
                 ])
    
    # 2. Build a hash file for HXB2
    if not os.path.isfile(get_HXB2_hash_file(data_folder, ext=True)):
        sp.call([stampy_bin,
                 '-g', get_HXB2_index_file(data_folder, ext=False),
                 '-H', get_HXB2_hash_file(data_folder, ext=False),
                 ])


def call_stampy_for_mapping(data_folder, adaID, VERBOSE=3):
    '''Call stampy for actual mapping'''
    readfiles = get_read_filenames(data_folder, adaID,
                                   subsample=True, # FIXME
                                   filtered=True)

    output_filename =  get_premapped_file(data_folder, adaID, type='sam',
                                          subsample=True) #FIXME
    # If the output file exist, do not map
    if os.path.isfile(output_filename):
        return None
    
    # Stampy command line
    # Note: no --solexa option as of 2013 (illumina v1.8)
    qsub_list = ['qsub','-cwd',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'stampy_'+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 stampy_bin,
                 '-g', get_HXB2_index_file(data_folder, ext=False),
                 '-h', get_HXB2_hash_file(data_folder, ext=False), 
                 '-o', output_filename,
                 '--substitutionrate='+subsrate,
                 '-M'] + readfiles
    qsub_list = map(str,qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)

    # Call stampy and check output
    qsub_output = sp.check_output(qsub_list).rstrip('\n')
    if VERBOSE:
        print qsub_output
    jobID = qsub_output.split(' ')[2]
    if VERBOSE:
        print jobID

    return jobID


def fork_self(data_folder, adaID, VERBOSE=3):
    '''Fork self for each adapter ID'''
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
                 '--stage', 3,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def write_read_pair(reads, ranges, fileouts):
    '''Write read pair into FASTQ files'''
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




# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--adaID', metavar='00', type=int, nargs='?',
                        default=0,
                        help='Adapter ID')
    parser.add_argument('--verbose', type=int, default=0,
                        help=('Verbosity level [0-3]'))
    parser.add_argument('--stage', metavar='N', type=int, nargs='?',
                        default=1,
                        help=('1: initialize, 2: map, 3: divide'))

    args = parser.parse_args()
    stage = args.stage
    adaID = args.adaID
    VERBOSE = args.verbose

    # FIXME: extend to all adaIDs
    if adaID == 0:
        adaID = 2

    ###########################################################################
    # 1. PREPARE HXB2 FOR MAPPING
    ###########################################################################
    if stage == 1:

        # Make index and hash files for HXB2, using a cropped reference to
        # reduce the problems of LTRs
        if VERBOSE: print 'Making index and hash files...'
        make_index_and_hash(data_folder, cropped=True)
        if VERBOSE: print 'done.'

        # This was easy, move to the mapping part
        stage = 2

    ###########################################################################
    # 2. MAP AGAINST HXB2 OR CONSENSUS
    ###########################################################################
    if stage == 2:

        # If the script is called with no adaID, make it a master script for all
        # adapters: otherwise, keep it a master script for its own stampy only
        if adaID == 0:
            adaIDs = load_adapter_table(data_folder)['ID']
        else:
            adaIDs = [adaID]

        # Call stampy for mapping and get the jobIDs
        if VERBOSE: print 'Mapping...'
        jobIDs = [call_stampy_for_mapping(data_folder, adaID) for adaID in adaIDs]
        
        # Wait for all stampy children to finish, then either proceed or fork self
        mapping_is_done = np.zeros(len(jobIDs), bool)
        while not mapping_is_done.all():
            time.sleep(interval_check)
        
            if VERBOSE >= 3:
                print 'Mapping progress: ', list(map(int, mapping_is_done))
        
            # Ask qstat about our stampy jobs
            qstat_output = sp.check_output(['qstat'])
            if VERBOSE >= 3:
                print qstat_output

            # Check all unfinished mappings
            for i, adaID in enumerate(adaIDs):
                if mapping_is_done[i]:
                    continue

                # If the output file was there already, consider it done
                if jobIDs[i] is None:
                    mapping_is_done[i] = True
                else:
                    # If that mapping is done, proceed to consensus building
                    # FIXME: this check is rather lousy
                    is_in_qstat_output = jobIDs[i] in qstat_output
                    if not is_in_qstat_output:
                        mapping_is_done[i] = True

                # Fork self to a child unless there is only one process
                if (mapping_is_done[i]):
                    if len(adaIDs) > 1:
                        fork_self(data_folder, adaID)
                    else:
                        stage = 3

        if VERBOSE: print 'done.'

    ###########################################################################
    # 3. DIVIDE INTO FRAGMENTS
    ###########################################################################
    if stage == 3:

        # Divide the mapped reads into fragments
        seq = load_HXB2()

        # Extract fragments
        frags_pos = [(f.location.nofuzzy_start, f.location.nofuzzy_end)
                     for f in seq.features if f.type[:8] == 'fragment']
        frags_pos = np.array(frags_pos, int).T
        # Since the reference is cropped, subtract from the positions F1 start
        frags_pos -= frags_pos.min()
        # Note: now we are in the reference of the CROPPED HXB2!

        # Get the premapped reads
        samfilename = get_premapped_file(data_folder, adaID, type='sam',
                                         subsample=True) #FIXME
        bamfilename = get_premapped_file(data_folder, adaID, type='bam',
                                         subsample=True) #FIXME
        if not os.path.isfile(bamfilename):
            convert_sam_to_bam(bamfilename, samfilename)

        # Output fastq files
        fileouts = get_read_filenames(data_folder, adaID,
                                      subsample=True, # FIXME
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
                read1, read2 = reads

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

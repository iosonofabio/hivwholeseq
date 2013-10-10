#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       07/10/13
content:    Get the allele counts for phiX (control)
'''
# Modules
import os
import sys
import cPickle as pickle
from collections import defaultdict
import numpy as np
from Bio import SeqIO
import pysam
import argparse

from mapping.datasets import MiSeq_runs
from mapping.miseq import alpha, read_types
from mapping.filenames import get_phix_filename, get_mapped_phix_filename, \
        get_allele_counts_phix_filename, get_insert_counts_phix_filename
from mapping.mapping_utils import convert_sam_to_bam, get_ind_good_cigars



# Globals
maxreads = 5e10
match_len_min = 30
trim_bad_cigars = 3

# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'get_allele_counts_phiX.py'
cluster_time = '0:59:59'
vmem = '8G'



# Functions
def fork_self(miseq_run, VERBOSE=0):
    '''Fork self for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+'{:02d}'.format(adaID)

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'acn phiX',
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', miseq_run,
                 '--verbose', VERBOSE,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE >= 2:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def get_allele_counts(data_folder, VERBOSE=0, maxreads=1e10):
    '''Extract allele and insert counts from a bamfile'''

    # Read reference
    reffilename = get_phix_filename()
    refseq = SeqIO.read(reffilename, 'fasta')
    
    # Allele counts and inserts
    counts = np.zeros((len(read_types), len(alpha), len(refseq)), int)
    # Note: the data structure for inserts is a nested dict with:
    # position --> string --> read type --> count
    #  (dict)      (dict)       (list)      (int)
    inserts = defaultdict(lambda: defaultdict(lambda: np.zeros(len(read_types), int)))

    # Open BAM file
    # Note: the reads should already be filtered of unmapped stuff at this point
    bamfilename = get_mapped_phix_filename(data_folder, filtered=False)
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        # Iterate over single reads (no linkage info needed)
        for i, read in enumerate(bamfile):

            if i > maxreads:
                if VERBOSE >= 2:
                    print 'Max reads reached:', maxreads
                break
        
            # Print output
            if (VERBOSE >= 3) and (not ((i +1) % 10000)):
                print (i+1)
        
            # Exclude unmapped and unpaired
            if read.is_unmapped or (not read.is_proper_pair):
                continue

            # Divide by read 1/2 and forward/reverse
            js = 2 * read.is_read2 + read.is_reverse
        
            # Read CIGARs (they should be clean by now)
            cigar = read.cigar
            len_cig = len(cigar)
            (good_cigars, first_good_cigar, last_good_cigar) = \
                    get_ind_good_cigars(cigar, match_len_min=match_len_min,
                                        full_output=True)
            
            # Sequence and position
            # Note: stampy takes the reverse complement already
            seq = read.seq
            pos = read.pos

            # Iterate over CIGARs
            for ic, (block_type, block_len) in enumerate(cigar):

                # Check for pos: it should never exceed the length of the fragment
                if (block_type in [0, 1, 2]) and (pos > len(refseq)):
                    raise ValueError('Pos exceeded the length of the fragment')
            
                # Inline block
                if block_type == 0:
                    # Exclude bad CIGARs
                    if good_cigars[ic]: 
            
                        # The first and last good CIGARs are matches:
                        # trim them (unless they end the read)
                        if (ic == first_good_cigar) and (ic != 0):
                            trim_left = trim_bad_cigars
                        else:
                            trim_left = 0
                        if (ic == last_good_cigar) and (ic != len_cig - 1):
                            trim_right = trim_bad_cigars
                        else:
                            trim_right = 0
            
                        seqb = np.array(list(seq[trim_left:block_len - trim_right]), 'S1')
                        # Increment counts
                        for j, a in enumerate(alpha):
                            posa = (seqb == a).nonzero()[0]
                            if len(posa):
                                counts[js, j, pos + trim_left + posa] += 1
            
                    # Chop off this block
                    if ic != len_cig - 1:
                        seq = seq[block_len:]
                        pos += block_len
            
                # Deletion
                elif block_type == 2:
                    # Exclude bad CIGARs
                    if good_cigars[ic]: 
                        # Increment gap counts
                        counts[js, 4, pos:pos + block_len] += 1
            
                    # Chop off pos, but not sequence
                    pos += block_len
            
                # Insertion
                # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
                # THEN the insert, FINALLY comes seq[391:]
                elif block_type == 1:
                    # Exclude bad CIGARs
                    if good_cigars[ic]: 
                        seqb = seq[:block_len]
                        inserts[pos][seqb][js] += 1
            
                    # Chop off seq, but not pos
                    if ic != len_cig - 1:
                        seq = seq[block_len:]
            
                # Other types of cigar?
                else:
                    raise ValueError('CIGAR type '+str(block_type)+' not recognized')

    return counts, inserts



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    miseq_run = args.run
    VERBOSE = args.verbose
    submit = args.submit

    # Submit to the cluster self if requested
    if submit:
        fork_self(miseq_run, VERBOSE=VERBOSE)
        sys.exit()

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # Get counts
    counts, inserts = get_allele_counts(data_folder, VERBOSE=VERBOSE, maxreads=maxreads)
    
    # Save counts and inserts to file
    allele_count_filename = get_allele_counts_phix_filename(data_folder)
    insert_count_filename = get_insert_counts_phix_filename(data_folder)
    counts.dump(allele_count_filename)
    # Change defaultdict into plain dicts, for pickling
    inserts = {key: dict(value) for (key, value) in inserts.iteritems()}
    with open(insert_count_filename, 'w') as f:
        pickle.dump(inserts, f, protocol=-1)

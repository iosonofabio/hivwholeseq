#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       13/09/13
content:    Get the joint counts of two alleles (2-site statistics).
'''
# Modules
import os
import sys
import argparse
import subprocess as sp
import cPickle as pickle
from operator import itemgetter
from collections import defaultdict
from itertools import izip
import pysam
import numpy as np
from Bio import SeqIO

# Horizontal import of modules from this folder
from mapping.adapter_info import load_adapter_table
from mapping.miseq import alpha, alphal, read_pair_types
from mapping.filenames import get_mapped_filename, get_allele_counts_filename, \
        get_insert_counts_filename, get_coverage_filename, get_consensus_filename
from mapping.mapping_utils import get_ind_good_cigars, get_trims_from_good_cigars, \
        pair_generator


# Globals
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']

match_len_min = 30
trim_bad_cigars = 3

# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBSCRIPT = JOBDIR+'get_coallele_counts.py'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
# Different times based on subsample flag
cluster_time = ['23:59:59', '0:59:59']
vmem = '8G'



# Functions
def fork_self(data_folder, adaID, fragment, subsample=False, VERBOSE=3):
    '''Fork self for each adapter ID and fragment'''
    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'acn '+'{:02d}'.format(adaID)+' '+fragment,
                 '-l', 'h_rt='+cluster_time[subsample],
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                ]
    if subsample:
        qsub_list.append('--subsample')
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def get_coallele_counts(data_folder, adaID, fragment, subsample=False, VERBOSE=0):
    '''Extract allele and insert counts from a bamfile'''

    # Read reference
    reffilename = get_consensus_filename(data_folder, adaID, fragment,
                                         subsample=subsample)
    refseq = SeqIO.read(reffilename, 'fasta')
    ref = np.array(refseq)
    
    # Allele counts and inserts (TODO: compress this data?)
    # Note: the pair is of 2 types only, while the single reads usually are of 4
    counts = np.zeros((len(read_pair_types),
                       len(alpha), len(alpha),
                       len(refseq), len(refseq)), int)
    positions = np.zeros(501, int)
    ais = np.zeros_like(positions)
    # TODO: no inserts for now

    # Open BAM file
    # Note: the reads should already be filtered of unmapped stuff at this point
    bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                      filtered=True)
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        # Iterate over read pairs
        for i, reads in enumerate(pair_generator(bamfile)):
        
            # Print output
            if (VERBOSE >= 3) and (not ((i +1) % 10)):
                print (i+1) 

            # Divide by read 1/2 and forward/reverse
            js = reads[0].is_reverse
            count = counts[js]

            # List of mutations
            positions[:] = -1
            ais[:] = -1
            imut = 0

            # Collect from the pair of reads
            for read in reads:
        
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
 
                            # Get the mutations and add them
                            indb = map(alphal.index, seq[trim_left:block_len - trim_right])
                            positions[imut: imut + len(indb)] = \
                                    pos + trim_left + np.arange(len(indb))
                            ais[imut: imut + len(indb)] = indb
                            imut += len(indb)

                        # Chop off this block
                        if ic != len_cig - 1:
                            seq = seq[block_len:]
                            pos += block_len
 
                    # Deletion
                    elif block_type == 2:                
                        # Chop off pos, but not sequence
                        pos += block_len
                
                    # Insertion
                    # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
                    # THEN the insert, FINALLY comes seq[391:]
                    elif block_type == 1:
                        # Chop off seq, but not pos
                        if ic != len_cig - 1:
                            seq = seq[block_len:]
                
                    # Other types of cigar?
                    else:
                        raise ValueError('CIGAR type '+str(block_type)+' not recognized')

            # Put the mutations into the matrix
            # Transversal
            for ai1 in xrange(len(alpha)):
                for ai2 in xrange(len(alpha)):
                    coun = count[ai1, ai2]
                    pos1 = positions[ais == ai1]
                    pos2 = positions[ais == ai2]
                    coords = np.meshgrid(pos1, pos2)
                    ind = coords[0].ravel() * count.shape[0] + coords[1].ravel()
                    coun.ravel()[ind] += 1

            ## Longitudinal
            #for i1, (pos1, ia1) in enumerate(izip(positions, ais)):
            #    for (pos2, ia2) in izip(positions[:i1], ais[:i1]):
            #        count[ia1, pos1, ia2, pos2] += 1

        import ipdb; ipdb.set_trace()
                                        

    return counts

# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Get allele counts')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--subsample', action='store_true',
                        help='Apply only to a subsample of the reads')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    subsample = args.subsample
    submit = args.submit

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over all requested samples
    for adaID in adaIDs:
        for fragment in fragments:

            # Submit to the cluster self if requested
            if submit:
                fork_self(data_folder, adaID, fragment,
                          subsample=subsample, VERBOSE=VERBOSE)
                continue

            # Get cocounts
            counts = get_coallele_counts(data_folder, adaID, fragment,
                                         subsample=subsample, VERBOSE=VERBOSE)

            # Get coverage
            #coverage = counts.sum(axis=1)

            ## Save to file
            #write_output_files(data_folder, adaID, fragment,
            #                   counts, inserts, coverage,
            #                   subsample=subsample, VERBOSE=VERBOSE)


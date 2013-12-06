#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       13/09/13
content:    Get the joint counts of two alleles (2-site statistics).
'''
# Modules
import os
import argparse
import subprocess as sp
from itertools import izip
import pysam
import numpy as np
from Bio import SeqIO

# Horizontal import of modules from this folder
from hivwholeseq.adapter_info import load_adapter_table
from hivwholeseq.miseq import alpha, alphal, read_pair_types
from hivwholeseq.filenames import get_mapped_filename, get_allele_counts_filename, \
        get_coallele_counts_filename, get_consensus_filename
from hivwholeseq.mapping_utils import pair_generator
from hivwholeseq.get_allele_counts import get_allele_counts


# Globals
# FIXME
from hivwholeseq.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']

match_len_min = 30
trim_bad_cigars = 3
maxreads = 1e2 #FIXME

# Cluster submit
import hivwholeseq
JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
JOBSCRIPT = JOBDIR+'get_coallele_counts.py'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
cluster_time = '0:59:59'
vmem = '8G'



# Functions
def fork_self(data_folder, adaID, fragment, VERBOSE=3):
    '''Fork self for each adapter ID and fragment'''
    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'cac '+'{:02d}'.format(adaID)+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def get_coallele_counts(data_folder, adaID, fragment, VERBOSE=0):
    '''Extract allele and insert counts from a bamfile'''

    # Read reference
    reffilename = get_consensus_filename(data_folder, adaID, fragment,
                                         trim_primers=True)
    refseq = SeqIO.read(reffilename, 'fasta')
    
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

            # Limit to some reads for testing
            if i > maxreads:
                if VERBOSE:
                    print 'Max read number reached:', maxreads
                break
        
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
        
                # Sequence and position
                # Note: stampy takes the reverse complement already
                seq = read.seq
                pos = read.pos
    
                # Iterate over CIGARs
                len_cig = len(read.cigar)
                for ic, (block_type, block_len) in enumerate(read.cigar):
    
                    # Check for pos: it should never exceed the length of the fragment
                    if (block_type in [0, 1, 2]) and (pos > len(refseq)):
                        raise ValueError('Pos exceeded the length of the fragment')
                
                    # Inline block
                    if block_type == 0:
 
                        # Get the mutations and add them
                        indb = map(alphal.index, seq)
                        positions[imut: imut + len(indb)] = \
                                pos + np.arange(len(indb))
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

            if VERBOSE >= 4:
                for pos, ai in izip(positions, ais):
                    if pos == -1:
                        break
                    print pos, ai

            # Put the mutations into the matrix
            for ai1 in xrange(len(alpha)):
                for ai2 in xrange(len(alpha)):
                    coun = count[ai1, ai2]
                    pos1 = positions[ais == ai1]
                    if ai1 == ai2: pos2 = pos1
                    else: pos2 = positions[ais == ai2]
                    coords = np.meshgrid(pos1, pos2)
                    ind = coords[0].ravel() * coun.shape[0] + coords[1].ravel()
                    coun.ravel()[ind] += 1                                        

    return counts


def write_output_files(data_folder, adaID, fragment, counts, VERBOSE=0):
    '''Write coallele counts to file'''
    if VERBOSE >= 1:
        print 'Write to file: '+'{:02d}'.format(adaID)+' '+fragment

    # Save counts and coverage
    # TODO: use compressed files?
    counts.dump(get_coallele_counts_filename(data_folder, adaID, fragment))



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
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
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
                fork_self(data_folder, adaID, fragment, VERBOSE=VERBOSE)
                continue

            # Get cocounts
            cocounts = get_coallele_counts(data_folder, adaID, fragment,
                                           VERBOSE=VERBOSE)

            ## Check using the allele counts and the diagonal cocounts
            #counts, _ = get_allele_counts(data_folder, adaID, fragment,
            #                              VERBOSE=VERBOSE,
            #                              maxreads=2 * maxreads)

            #cocount = cocounts.sum(axis=0)
            #count = counts.sum(axis=0)

            ## Read reference
            #reffilename = get_consensus_filename(data_folder, adaID, fragment,
            #                                     trim_primers=True)
            #refseq = SeqIO.read(reffilename, 'fasta')
            #ref = np.array(refseq)
            #refi = np.array([(alpha == a).nonzero()[0][0] for a in ref], int)

            ## Check that counts and diagonal cocounts are the same thing
            #for i in xrange(100):
            #    j = refi[i]
            #    print cocount[j, j, i, i], count[j, i]

            # Save to file
            write_output_files(data_folder, adaID, fragment,
                               cocounts, VERBOSE=VERBOSE)


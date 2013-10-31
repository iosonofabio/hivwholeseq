#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Get the allele frequencies out of a BAM file and a reference.
'''
# Modules
import os
import argparse
import subprocess as sp
import cPickle as pickle
from collections import defaultdict
import pysam
import numpy as np
from Bio import SeqIO

from mapping.datasets import MiSeq_runs
from mapping.adapter_info import load_adapter_table
from mapping.miseq import alpha, read_types
from mapping.filenames import get_mapped_filename, get_allele_counts_filename, \
        get_insert_counts_filename, get_coverage_filename, get_consensus_filename
from mapping.mapping_utils import convert_sam_to_bam


# Globals
# Minimal quality required for a base to be considered trustful (i.e. added to 
# the allele counts), in phred score. Too high: lose coverage, too low: seq errors.
# Reasonable numbers are between 30 and 36.
qual_min = 35

# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBSCRIPT = JOBDIR+'get_allele_counts.py'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
cluster_time = '0:59:59'
vmem = '4G'



# Functions
def fork_self(miseq_run, adaID, fragment, VERBOSE=3):
    '''Fork self for each adapter ID and fragment'''
    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'acn '+'{:02d}'.format(adaID)+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', miseq_run,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def get_allele_counts_insertions_from_file(bamfilename, length,
                                           maxreads=-1, VERBOSE=0):
    '''Get the allele counts and insertions'''
    # Prepare output structures
    counts = np.zeros((len(read_types), len(alpha), length), int)
    # Note: the data structure for inserts is a nested dict with:
    # position --> string --> read type --> count
    #  (dict)      (dict)       (list)      (int)
    inserts = defaultdict(lambda: defaultdict(lambda: np.zeros(len(read_types), int)))

    # Open BAM file
    # Note: the reads should already be filtered of unmapped stuff at this point
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        # Iterate over single reads
        for i, read in enumerate(bamfile):

            # Max number of reads
            if i == maxreads:
                if VERBOSE >= 2:
                    print 'Max reads reached:', maxreads
                break
        
            # Print output
            if (VERBOSE >= 3) and (not ((i +1) % 10000)):
                print (i+1)
        
            # Divide by read 1/2 and forward/reverse
            js = 2 * read.is_read2 + read.is_reverse
        
            # Read CIGARs (they should be clean by now)
            seq = np.fromstring(read.seq, 'S1')
            qual = np.fromstring(read.qual, np.int8) - 33
            pos = read.pos

            # Iterate over CIGARs
            for ic, (block_type, block_len) in enumerate(read.cigar):

                # Check for pos: it should never exceed the length of the fragment
                if (block_type in [0, 1, 2]) and (pos >= length):
                    raise ValueError('Pos exceeded the length of the fragment')
            
                # Inline block
                if block_type == 0:
                    seqb = seq[:block_len]
                    qualb = qual[:block_len]
                    # Increment counts
                    for j, a in enumerate(alpha):
                        posa = ((seqb == a) & (qualb >= qual_min)).nonzero()[0]
                        if len(posa):
                            counts[js, j, pos + posa] += 1
            
                    # Chop off this block
                    if ic != len(read.cigar) - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
                        pos += block_len
            
                # Deletion
                elif block_type == 2:
                    # Increment gap counts
                    counts[js, 4, pos:pos + block_len] += 1
            
                    # Chop off pos, but not sequence
                    pos += block_len
            
                # Insertion
                # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
                # THEN the insert, FINALLY comes seq[391:]
                elif block_type == 1:
                    seqb = seq[:block_len]
                    qualb = qual[:block_len]
                    # Accept only high-quality inserts
                    if (qualb >= qual_min).all():
                        inserts[pos][seqb.tostring()][js] += 1
            
                    # Chop off seq, but not pos
                    if ic != len(read.cigar) - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
            
                # Other types of cigar?
                else:
                    raise ValueError('CIGAR type '+str(block_type)+' not recognized')

    return (counts, inserts)


def get_allele_counts(data_folder, adaID, fragment, VERBOSE=0,
                      maxreads=1e10):
    '''Extract allele and insert counts from a bamfile'''

    # Read reference
    reffilename = get_consensus_filename(data_folder, adaID, fragment,
                                         trim_primers=True)
    refseq = SeqIO.read(reffilename, 'fasta')

    # Open BAM file
    # Note: the reads should already be filtered of unmapped stuff at this point
    bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                      filtered=True)
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)

    # Call lower-level function
    return get_allele_counts_insertions_from_file(bamfilename, len(refseq),
                                                 maxreads=maxreads, VERBOSE=VERBOSE)


def write_output_files(data_folder, adaID, fragment,
                       counts, inserts, coverage, VERBOSE=0):
    '''Write allele counts, inserts, and coverage to file'''
    if VERBOSE >= 1:
        print 'Write to file: '+'{:02d}'.format(adaID)+' '+fragment

    # Save counts and coverage
    counts.dump(get_allele_counts_filename(data_folder, adaID, fragment))
    coverage.dump(get_coverage_filename(data_folder, adaID, fragment))

    # Convert inserts to normal nested dictionary for pickle
    inserts_dic = {k: dict(v) for (k, v) in inserts.iteritems()}
    with open(get_insert_counts_filename(data_folder, adaID, fragment), 'w') as f:
        pickle.dump(inserts_dic, f, protocol=-1)





# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Get allele counts')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    submit = args.submit

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

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
                fork_self(miseq_run, adaID, fragment, VERBOSE=VERBOSE)
                continue

            # Get counts
            counts, inserts = get_allele_counts(data_folder, adaID, fragment,
                                                VERBOSE=VERBOSE)

            # Get coverage
            coverage = counts.sum(axis=1)

            # Save to file
            write_output_files(data_folder, adaID, fragment,
                               counts, inserts, coverage,
                               VERBOSE=VERBOSE)

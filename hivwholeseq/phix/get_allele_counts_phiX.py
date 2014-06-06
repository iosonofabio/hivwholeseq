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
import subprocess as sp
import cPickle as pickle
from collections import defaultdict
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import ambiguous_dna
import pysam
import argparse

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.miseq import alpha, read_types
from hivwholeseq.filenames import get_phix_filename, get_mapped_phix_filename, \
        get_allele_counts_phix_filename, get_insert_counts_phix_filename, \
        get_consensus_phix_filename
from hivwholeseq.mapping_utils import convert_sam_to_bam, get_ind_good_cigars



# Globals
match_len_min = 30
trim_bad_cigars = 3

# Cluster submit
import hivwholeseq
JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'phix/get_allele_counts_phiX.py'
cluster_time = '0:59:59'
vmem = '8G'



# Functions
def fork_self(seq_run, maxreads=-1, VERBOSE=0, qual_min=None):
    '''Fork self for each adapter ID'''

    if VERBOSE:
        print 'Forking to the cluster'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'acpX'+seq_run,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--maxreads', maxreads,
                 '--verbose', VERBOSE,
                ]
    if qual_min is not None:
        qsub_list.extend(['--qual_min', qual_min])
    qsub_list = map(str, qsub_list)
    if VERBOSE >= 2:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def get_allele_counts(data_folder, qual_min=0, VERBOSE=0, maxreads=-1):
    '''Extract allele and insert counts from a bamfile'''

    reffilename = get_phix_filename()
    refseq = SeqIO.read(reffilename, 'fasta')

    from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file
    
    bamfilename = get_mapped_phix_filename(data_folder, filtered=True)
    (counts, inserts) = get_allele_counts_insertions_from_file(bamfilename,
                                                               len(refseq),
                                                               qual_min=qual_min,
                                                               maxreads=maxreads,
                                                               VERBOSE=VERBOSE)
    return (counts, inserts)


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--run', required=True,
                        help='MiSeq run to analyze (e.g. Tue28)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--qual_min', type=int, default=0,
                        help='Maximal quality of nucleotides')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose
    maxreads = args.maxreads
    qual_min = args.qual_min
    submit = args.submit

    # Submit to the cluster self if requested
    if submit:
        fork_self(seq_run, maxreads=maxreads, VERBOSE=VERBOSE, qual_min=qual_min)
        sys.exit()

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # Get counts
    counts, inserts = get_allele_counts(data_folder, VERBOSE=VERBOSE, maxreads=maxreads, qual_min=qual_min)
    consm = alpha[counts.sum(axis=0).argmax(axis=0)]
    cons_rec = SeqRecord(Seq(''.join(consm), ambiguous_dna),
                         id='PhiX_consensus_run_'+seq_run,
                         name='PhiX_consensus_run_'+seq_run,
                         description='PhiX consensus run '+seq_run)

    
    # Save counts and inserts to file
    allele_count_filename = get_allele_counts_phix_filename(data_folder, qual_min=qual_min)
    insert_count_filename = get_insert_counts_phix_filename(data_folder, qual_min=qual_min)
    counts.dump(allele_count_filename)
    SeqIO.write(cons_rec, get_consensus_phix_filename(data_folder), 'fasta')
    # Change defaultdict into plain dicts, for pickling
    inserts = {key: dict(value) for (key, value) in inserts.iteritems()}
    with open(insert_count_filename, 'w') as f:
        pickle.dump(inserts, f, protocol=-1)

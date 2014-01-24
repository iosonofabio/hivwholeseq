#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/08/13
content:    Build a subset of the mapped reads excluding mismappings.
'''
# Modules
import os
import argparse
from operator import itemgetter
import pysam
import numpy as np
from Bio import SeqIO


# Horizontal import of modules from this folder
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.adapter_info import load_adapter_table
from hivwholeseq.filenames import get_mapped_phix_filename
from hivwholeseq.mapping_utils import get_ind_good_cigars, convert_sam_to_bam,\
        pair_generator, get_range_good_cigars
from hivwholeseq.primer_info import primers_inner



# Globals
match_len_min = 30
trim_bad_cigars = 3

# Cluster submit
import hivwholeseq
JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr/'
JOBLOGOUT = JOBDIR+'logout/'
JOBSCRIPT = JOBDIR+'phix/filter_mapped_reads_phix.py'
cluster_time = '0:59:59'
vmem = '8G'



# Functions
def fork_self(seq_run, maxreads=-1, VERBOSE=0):
    '''Fork self'''
    import subprocess as sp

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'fmr phix',
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--verbose', VERBOSE,
                 '--maxreads', maxreads,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def trim_bad_cigar(read_pair, match_len_min=match_len_min,
                   trim_left=trim_bad_cigars, trim_right=trim_bad_cigars):
    '''Trim away bad CIGARs from the sides'''

    for read in read_pair:
        (good_cigars, first_good_cigar, last_good_cigar) = \
                get_ind_good_cigars(read.cigar, match_len_min=match_len_min,
                                    full_output=True)
    
        if not good_cigars.any():
            reads = None
            return
    
        # Get the good CIGARs coordinates
        ((start_read, end_read),
         (start_ref, end_ref)) = \
                get_range_good_cigars(read.cigar, read.pos,
                                      match_len_min=match_len_min,
                                      trim_left=trim_left,
                                      trim_right=trim_right)
    
        # Trim CIGAR because of bad CIGARs at the edges
        cigar = read.cigar[first_good_cigar: last_good_cigar + 1]
        # Trim cigar block lengths
        if first_good_cigar != 0:
            cigar[0] = (cigar[0][0],
                        cigar[0][1] - trim_left)
        if last_good_cigar != len(read.cigar) - 1:
            cigar[-1] = (cigar[-1][0],
                         cigar[-1][1] - trim_right)
    
        # Reset attributes
        seq = read.seq
        qual = read.qual
        read.seq = seq[start_read: end_read]
        read.qual = qual[start_read: end_read]
        read.pos = start_ref
        read.cigar = cigar    
        
        # Give up number of errors
        read.tags.pop(map(itemgetter(0), read.tags).index('NM'))

    # Fix mate pair
    i_fwd = read_pair[0].is_reverse
    i_rev = not i_fwd
    read_pair[i_fwd].mpos = read_pair[i_rev].pos
    read_pair[i_rev].mpos = read_pair[i_fwd].pos
    isize = read_pair[i_rev].pos + \
            sum(bl for bt, bl in read_pair[i_rev].cigar if bt in (0, 2)) -\
            read_pair[i_fwd].pos
    read_pair[i_fwd].isize = isize
    read_pair[i_rev].isize = -isize


def filter_reads(data_folder, maxreads=-1, VERBOSE=0):
    '''Filter the reads to good chunks'''

    bamfilename = get_mapped_phix_filename(data_folder, type='bam')
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)

    outfilename = get_mapped_phix_filename(data_folder, type='bam', filtered=True)
    trashfilename = outfilename[:-4]+'_trashed.bam'
 
    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        with pysam.Samfile(outfilename, 'wb', template=bamfile) as outfile,\
             pysam.Samfile(trashfilename, 'wb', template=bamfile) as trashfile:
 
            # Iterate over all pairs
            n_good = 0
            n_unmapped = 0
            n_unpaired = 0
            n_mutator = 0
            n_badcigar = 0
            for i_pairs, reads in enumerate(pair_generator(bamfile)):

                if i_pairs == maxreads:
                    if VERBOSE:
                        print 'Maximal number of read pairs reached:', maxreads
                    break

                if VERBOSE >= 2:
                    if not ((i_pairs+1) % 10000):
                        print i_pairs+1, maxreads
            
                (read1, read2) = reads

                # Basic checks (name, mapped, paired)
                if read1.qname != read2.qname:
                    map(trashfile.write, reads)
                    continue

                if read1.is_unmapped or read2.is_unmapped:
                    if VERBOSE >= 3:
                        print 'Read pair '+read1.qname+': unmapped'
                    n_unmapped += 1
                    map(trashfile.write, reads)
                    continue
            
                if (not read1.is_proper_pair) or (not read2.is_proper_pair):
                    if VERBOSE >= 3:
                        print 'Read pair '+read1.qname+': not properly paired'
                    n_unpaired += 1
                    map(trashfile.write, reads)
                    continue

                # Mismappings are often characterized by many mutations:
                # check the number of mismatches and skip reads with too many
                mm = (dict(read1.tags)['NM'], dict(read2.tags)['NM'])
                if (max(mm) > 50) or (sum(mm) > 50):
                    if VERBOSE >= 3:
                        print 'Read pair '+read1.qname+': too many mismatches '+\
                                '('+str(mm[0])+' + '+str(mm[1])+')'
                    n_mutator += 1
                    map(trashfile.write, reads)
                    continue

                # Trim the bad CIGARs from the sides (in place)
                trim_bad_cigar(reads, match_len_min=match_len_min,
                               trim_left=trim_bad_cigars,
                               trim_right=trim_bad_cigars)
                # If there are no good CIGARs, skip
                if reads is None:
                    map(trashfile.write, reads)
                    continue

                # Write the output: if we got here...
                n_good += 1
                map(outfile.write, reads)


    if VERBOSE >= 1:
        print 'Read pairs: '+str(n_good)+' good, '+str(n_unmapped)+' unmapped, '+\
                str(n_unpaired)+' unpaired, '+\
                str(n_mutator)+' many-mutations, '+str(n_badcigar)+' bad CIGAR.'



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Filter mapped reads')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--verbose', type=int, default=0,
                        help=('Verbosity level [0-3]'))
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--submit', action='store_true', default=False,
                        help='Submit the job to the cluster via qsub')

    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose
    maxreads = args.maxreads
    submit = args.submit

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    if submit:
        fork_self(seq_run, maxreads=maxreads, VERBOSE=VERBOSE)
    else:
        filter_reads(data_folder, maxreads=maxreads, VERBOSE=VERBOSE)

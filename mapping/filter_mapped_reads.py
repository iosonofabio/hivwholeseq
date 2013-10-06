#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/08/13
content:    Build a subset of the mapped reads excluding mismappings.
'''
# Modules
import os
import sys
import argparse
import cPickle as pickle
from operator import itemgetter
from collections import defaultdict
from itertools import izip
import pysam
import numpy as np
from Bio import SeqIO


# Horizontal import of modules from this folder
from mapping.adapter_info import load_adapter_table, foldername_adapter
from mapping.miseq import alpha, read_types
from mapping.filenames import get_consensus_filename, get_mapped_filename
from mapping.mapping_utils import get_ind_good_cigars, convert_sam_to_bam,\
        pair_generator, get_range_good_cigars
from mapping.primer_info import primers_inner



# Globals
# FIXME
from mapping.datasets import dataset_2 as dataset
data_folder = dataset['folder']

maxreads = 1e9
match_len_min = 30
trim_bad_cigars = 3

# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'filter_mapped_reads.py'
cluster_time = '0:59:59'
vmem = '8G'



# Functions
def fork_self(data_folder, adaID, fragment, VERBOSE=0):
    '''Fork self for each adapter ID'''
    import subprocess as sp

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'fmr '+'{:02d}'.format(adaID)+' '+fragment,
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


def trim_bad_cigar(read, match_len_min=match_len_min,
                   trim_left=trim_bad_cigars, trim_right=trim_bad_cigars):
    '''Trim away bad CIGARs from the sides'''

    # Get good CIGARs
    (good_cigars, first_good_cigar, last_good_cigar) = \
            get_ind_good_cigars(read.cigar, match_len_min=match_len_min,
                                full_output=True)

    if not good_cigars.any():
        read = None
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


def trim_primers(read, start_nonprimer, end_nonprimer):
    '''Trim the fwd and rev inner PCR primers'''

    # Is the read partially in the fwd primer?
    touches_fwd_primer = read.pos < start_nonprimer
    if not touches_fwd_primer:
        read.pos -= start_nonprimer
    else:
        # Trim the primer: the read starts with a long match (it has been CIGAR-
        # trimmed already), so we can just chop the primer off the first block
        start_read = start_nonprimer - read.pos
        if start_read > read.cigar[0][1]:
            raise ValueError('The fwd primer is not within the first CIGAR block!')
        cigar = read.cigar
        cigar[0] = (cigar[0][0], cigar[0][1] - start_read)
        
        # Reset attributes
        seq = read.seq
        qual = read.qual
        read.seq = seq[start_read:]
        read.qual = qual[start_read:]
        read.pos = start_nonprimer
        read.cigar = cigar
        
    # Is the read partially in the rev primer?
    end_ref = read.pos + sum(bl for (bt, bl) in read.cigar if bt in (0, 2))
    touches_rev_primer = end_ref > end_nonprimer - start_nonprimer
    if touches_rev_primer:

        # Trim the primer; the read ends with a long match (it has been CIGAR-
        # trimmed already), so we can just chop the primer off the last block
        end_read = read.rlen + end_nonprimer - start_nonprimer - end_ref
        if read.rlen - end_read > read.cigar[-1][1]:
            raise ValueError('The rev primer is not within the last CIGAR block!')
        cigar = read.cigar
        cigar[-1] = (cigar[-1][0], cigar[-1][1] - (read.rlen - end_read))

        # Reset attributes
        # Reset attributes
        seq = read.seq
        qual = read.qual
        read.seq = seq[:end_read]
        read.qual = qual[:end_read]
        read.cigar = cigar


def filter_reads(data_folder, adaID, fragment, VERBOSE=0):
    '''Filter the reads to good chunks'''
    from Bio import pairwise2

    # Resolve ambiguous fragment primers
    if fragment in ['F5a', 'F5b']:
        frag_gen = 'F5'
    else:
        frag_gen = fragment

    # Read reference (fragmented)
    reffilename = get_consensus_filename(data_folder, adaID, frag_gen)
    refseq = SeqIO.read(reffilename, 'fasta')
    ref = np.array(refseq)

    # Get the coordinate of the inner PCR primers in this consensus via local
    # alignment
    primer_fwd, primer_rev = primers_inner[fragment]
    len_localali = 30
    ref_fwd = str(refseq.seq)[:len_localali]
    ref_rev = str(refseq.seq)[-len_localali:]
    ali_fwd = pairwise2.align.localms(ref_fwd, primer_fwd, 2, -1, -0.5, -0.1)[0][:2]
    primer_fwd_end = len(ref_fwd) - ali_fwd[1][::-1].index(primer_fwd[-1])
    primer_fwd_end -= ali_fwd[0][:primer_fwd_end].count('-')

    # The reverse primer works, well, the other way around
    ali_rev = pairwise2.align.localms(ref_rev, primer_rev, 2, -1, -0.5, -0.1)[0][:2]
    primer_rev_start = len(ref) - len(ref_rev) + ali_rev[1].index(primer_rev[0])
    primer_rev_start += ref_rev[len(ref_rev) - len(ref) + primer_rev_start:].count('-')

    # Give decent names
    start_nonprimer = primer_fwd_end
    end_nonprimer = primer_rev_start

    # Get BAM files
    bamfilename = get_mapped_filename(data_folder, adaID, frag_gen, type='bam')
    # Try to convert to BAM if needed
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)

    outfilename = get_mapped_filename(data_folder, adaID, frag_gen, type='bam',
                                     filtered=True)
    trashfilename = outfilename[:-4]+'_trashed.bam'
 
    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        with pysam.Samfile(outfilename, 'wb', template=bamfile) as outfile,\
             pysam.Samfile(trashfilename, 'wb', template=bamfile) as trashfile:
 
            # Iterate over all pairs
            n_good = 0
            n_unmapped = 0
            n_unpaired = 0
            n_mutator = 0
            n_mismapped_edge = 0
            n_badcigar = 0
            for i_pairs, reads in enumerate(pair_generator(bamfile)):

                # Limit to the first reads
                if 2 * i_pairs >= maxreads: break
            
                # Assign names
                (read1, read2) = reads

                # Flag to decide on the read
                skip = False
            
                # Check a few things to make sure we are looking at paired reads
                if read1.qname != read2.qname:
                    raise ValueError('Read pair '+str(i_pairs)+': reads have different names!')
                # Ignore unmapped reads
                elif read1.is_unmapped or read2.is_unmapped:
                    if VERBOSE >= 2:
                        print 'Read pair '+read1.qname+': unmapped'
                    n_unmapped += 1
                    skip = True
            
                # Ignore not properly paired reads (this includes mates sitting on
                # different fragments)
                elif (not read1.is_proper_pair) or (not read2.is_proper_pair):
                    if VERBOSE >= 2:
                        print 'Read pair '+read1.qname+': not properly paired'
                    n_unpaired += 1
                    skip = True

                else:

                    # Mismappings are often characterized by many mutations:
                    # check the number of mismatches and skip reads with too many
                    mm = (dict(read1.tags)['NM'], dict(read2.tags)['NM'])
                    if (max(mm) > 50) or (sum(mm) > 50):
                        if VERBOSE >= 2:
                            print 'Read pair '+read1.qname+': too many mismatches '+\
                                    '('+str(mm[0])+' + '+str(mm[1])+')'
                        n_mutator += 1
                        skip = True

                    else:
    
                        # Mismappings are sometimes at fragment edges:
                        # Check for overhangs beyond the edge
                        for read in reads:
                            # Check overhangs
                            read_start = read.pos
                            read_end = read.pos + sum(x[1] for x in read.cigar if x[0] != 1)
                            if (((read_start == 0) and (read.cigar[0][0] == 1)) or
                                ((read_end == len(ref)) and (read.cigar[-1][0] == 1))):
                                n_mismapped_edge += 1
                                skip = True
                                break

                # If the read pair survived, check and trim good cigars
                if not skip:
                    for read in reads:

                        # Trim the bad CIGARs from the sides (in place)
                        trim_bad_cigar(read, match_len_min=match_len_min,
                                       trim_left=trim_bad_cigars,
                                       trim_right=trim_bad_cigars)

                        # If there are no good CIGARs, skip
                        if read is None:
                            n_badcigar += 1
                            skip = True
                            break

                        # Trim inner PCR primers (in place)
                        trim_primers(read, start_nonprimer, end_nonprimer)

                        # Give up the number of mismatches
                        read.tags.pop(map(itemgetter(0), read.tags).index('NM'))

                    # Mate pair stuff and insert size
                    if not skip:
                        read1.mpos = read2.pos
                        read2.mpos = read1.pos

                        # Insert size
                        readf = reads[read1.is_reverse]
                        readr = reads[read2.is_reverse]
                        isize = max([read1.pos + read1.rlen,
                                     read2.pos + read2.rlen]) -\
                                min([read1.pos, read2.pos])
                        readf.isize = isize
                        readr.isize = -isize

                # Write the output
                if skip:
                    map(trashfile.write, reads)
                else:
                    n_good += 1
                    map(outfile.write, reads)


    if VERBOSE >= 1:
        print 'Reads: '+str(n_good)+' good, '+str(n_unmapped)+' unmapped, '+\
                str(n_unpaired)+' unpaired, '+str(n_mismapped_edge)+' edge, '+\
                str(n_mutator)+' many-mutations, '+str(n_badcigar)+' bad CIGAR.'



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Extract linkage information')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help=('Verbosity level [0-3]'))
    parser.add_argument('--submit', action='store_true', default=False,
                        help='Submit the job to the cluster via qsub')
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
    for i, adaID in enumerate(adaIDs):
        for fragment in fragments:

            # Submit to the cluster self if requested
            if submit:
                fork_self(data_folder, adaID, fragment, VERBOSE=VERBOSE)
                continue

            # or else, perform the filtering
            if fragment == 'F5':
                fragment = dataset['primerF5'][i]
            filter_reads(data_folder, adaID, fragment, VERBOSE=VERBOSE)

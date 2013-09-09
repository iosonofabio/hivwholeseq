# vim: fdm=indent
'''
author:     Fabio Zanini
date:       06/08/13
content:    Get the histogram of insert sizes.

            This is particularly important because some of the maps contain
            contaminants (fosmid plasmid chunks), which give rise to spuriously
            short inserts.
'''
# Modules
import os
import sys
import cPickle as pickle
from collections import defaultdict
from itertools import izip
import pysam
import numpy as np
from Bio import SeqIO

from mapping.mapping_utils import pair_generator



# Globals
VERBOSE = 1
alpha = np.array(list('ACGT-N'), 'S1')
read_types = ['read1 f', 'read1 r', 'read2 f', 'read2 r']
maxreads = 5000
match_len_min = 30
trim_bad_cigars = 3
ref_filename = 'consensus_filtered_trimmed.fasta'
bam_filename = 'mapped_to_self_filtered_trimmed.bam'



# Script
if __name__ == '__main__':


    # Input arguments
    args = sys.argv
    if len(args) < 2:
        raise ValueError('This script takes the adapterID folder as input')
    data_folder = args[1].rstrip('/')+'/'

    # Read reference
    if os.path.isfile(data_folder+ref_filename): ref_file = data_folder+ref_filename
    else: ref_file = '/'.join(data_folder.split('/')[:-2]+['subsample/']+data_folder.split('/')[-2:])+ref_filename
    refseq = SeqIO.read(ref_file, 'fasta')
    ref = np.array(refseq)

    # Open BAM
    bamfilename = data_folder+bam_filename
    # Try to convert to BAM if needed
    if not os.path.isfile(bamfilename):
        samfile = pysam.Samfile(bamfilename[:-3]+'sam', 'r')
        bamfile = pysam.Samfile(bamfilename, 'wb', template=samfile)
        for s in samfile:
            bamfile.write(s)
    bamfile = pysam.Samfile(bamfilename, 'rb')


    # Iterate through reads
    for i, read in enumerate(bamfile):
    
        # Limit to the first reads
        if i >= maxreads: break

        # Print output
        if VERBOSE and not ((i +1) % 10000):
            print (i+1)
    
        # Ignore unmapped reads
        if read.is_unmapped or not read.is_proper_pair:
            continue
    
        # Divide by read 1/2 and forward/reverse
        if read.is_read1: js = 0
        else: js = 2
        if read.is_reverse: js += 1
    
        # Sequence and position
        # Note: stampy takes the reverse complement already
        seq = read.seq
        pos = read.pos
    
        # Read CIGAR code for indels, and anayze each block separately
        # Note: some reads are weird ligations of HIV and contaminants
        # (e.g. fosmid plasmids). Those always have crazy CIGARs, because
        # only the HIV part maps. We hence trim away short indels from the
        # end of reads (this is still unbiased).
        cigar = read.cigar
        len_cig = len(cigar)
        good_cigars = np.array(map(lambda x: (x[0] == 0) and (x[1] >= match_len_min), cigar), bool, ndmin=1)
        # If no long match, skip read
        # FIXME: we must skip the mate pair also? But what if it's gone already?
        # Note: they come in pairs: read1 first, read2 next, so we can just read two at a time
        if not (good_cigars).any():
            continue
        # If only one long match, no problem
        if (good_cigars).sum() == 1:
            first_good_cigar = last_good_cigar = good_cigars.nonzero()[0][0]
        # else include stuff between the extremes
        else:
            tmp = good_cigars.nonzero()[0]
            first_good_cigar = tmp[0]
            last_good_cigar = tmp[-1]
            good_cigars[first_good_cigar: last_good_cigar + 1] = True

        # Iterate over CIGARs
        for ic, block in enumerate(cigar):

            # Inline block
            if block[0] == 0:
                # Exclude bad CIGARs
                if good_cigars[ic]: 

                    # The first and last good CIGARs are matches: trim them (unless they end the read)
                    if (ic == first_good_cigar) and (ic != 0): trim_left = trim_bad_cigars
                    else: trim_left = 0
                    if (ic == last_good_cigar) and (ic != len_cig - 1): trim_right = trim_bad_cigars
                    else: trim_right = 0
    
                    seqb = np.array(list(seq[trim_left:block[1] - trim_right]), 'S1')
                    # Increment counts
                    for j, a in enumerate(alpha):
                        posa = (seqb == a).nonzero()[0]
                        if len(posa):
                            counts[js, j, pos + trim_left + posa] += 1
    
                # Chop off this block
                if ic != len_cig - 1:
                    seq = seq[block[1]:]
                    pos += block[1]
    
            # Deletion
            elif block[0] == 2:
                # Exclude bad CIGARs
                if good_cigars[ic]: 
                    # Increment gap counts
                    counts[js, 4, pos:pos + block[1]] += 1
    
                # Chop off pos, but not sequence
                pos += block[1]
    
            # Insertion
            elif block[0] == 1:
                # Exclude bad CIGARs
                if good_cigars[ic]: 
                    seqb = seq[:block[1]]
                    inserts[js][(pos, seqb)] += 1
    
                # Chop off seq, but not pos
                if ic != len_cig - 1:
                    seq = seq[block[1]:]
    
            # Other types of cigar?
            else:
                raise ValueError('CIGAR type '+str(block[0])+' not recognized')



# vim: fdm=indent
'''
author:     Fabio Zanini
date:       21/08/13
content:    Test script for extracting linkage information from the reads.
'''
# Modules
import os
import sys
import argparse
import cPickle as pickle
from collections import defaultdict
from itertools import izip
import pysam
import numpy as np
from Bio import SeqIO


# Horizontal import of modules from this folder
from mapping.adapter_info import load_adapter_table, foldername_adapter
from mapping.miseq import alpha, read_types, pair_generator


# Globals
VERBOSE = 1

# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']
HXB2_fragmented_file = 'HXB2_fragmented.fasta'

maxreads = 500000000000000000    #FIXME
match_len_min = 30
trim_bad_cigars = 3
ref_filename = 'consensus_filtered_trimmed.fasta'
bam_filename = 'mapped_to_self_filtered_trimmed.bam'



# Functions
def get_ind_good_cigars(cigar, match_len_min=30):
    '''Keep only CIGAR blocks between two long matches'''
    good_cigars = np.array(map(lambda x: (x[0] == 0) and (x[1] >= match_len_min), cigar), bool, ndmin=1)
    # If there are no or one good CIGAR, keep that
    if (good_cigars).sum() < 2:
        return good_cigars
    # If there are 2+, keep stuff in the middle
    else:
        tmp = good_cigars.nonzero()[0]
        first_good_cigar = tmp[0]
        last_good_cigar = tmp[-1]
        good_cigars[first_good_cigar: last_good_cigar + 1] = True
    return good_cigars





# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Map HIV reads recursively')
    parser.add_argument('--adaID', metavar='00', type=int, required=True,
                        help='Adapter ID sample to analyze')

    args = parser.parse_args()
    adaID = args.adaID


    # Read reference
    ref_file = data_folder+'subsample/'+foldername_adapter(adaID)+ref_filename
    refseq = SeqIO.read(ref_file, 'fasta')
    ref = np.array(refseq)

    # Open BAM
    # FIXME (work on full data)
    bamfilename = data_folder+'subsample/'+foldername_adapter(adaID)+bam_filename
    # Try to convert to BAM if needed
    if not os.path.isfile(bamfilename):
        samfile = pysam.Samfile(bamfilename[:-3]+'sam', 'r')
        bamfile = pysam.Samfile(bamfilename, 'wb', template=samfile)
        for s in samfile:
            bamfile.write(s)
    bamfile = pysam.Samfile(bamfilename, 'rb')

    # TODO: there will be chromosomes

    # Prepare data structures
    count_pairs = [{} for read_type in read_types]

    # Iterate over all pairs
    for i_pairs, (read1, read2) in pair_generator(bamfile):

        # Check a few things to make sure we are looking at paired reads
        if read1.qname != read2.qname:
            raise ValueError('Read pair '+str(i_pairs)+': reads have different names!')

        # Ignore unmapped reads
        if (read1.is_unmapped or (not read1.is_proper_pair)
            or read2.is_unmapped or (not read2.is_proper_pair)):
            continue)

        # Limit to the first reads
        if 2 * i_pairs >= maxreads: break

        # Take both reads at the same time
        if read2.is_reverse:
            read_pair = {'fwd': read1, 'rev': read2}
        else:
            read_pair = {'fwd': read2, 'rev': read1}

        # Store positions of all cigars in the reference of the read
        poss = {}
        for key, read in read_pair.iteritems():
            pos = [0]
            for (block_type, block_len) in enumerate(read.cigar):
                pos.append(pos[-1] + block_len)
            poss[key] = pos

        # Keep only good match CIGARs (TODO: what about indels?)
        good_cigars = {}
        for key, read in read_pair.iteritems():
            is_good = get_ind_good_cigars(read.cigar, match_len_min=match_len_min)
            is_match = [True if block_type == 0 else False for (block_type, block_len) in read.cigar]
            good_cigars[key] = is_good & is_match

        # OK, now it's about storing all pairs of mutations within and between the mates
        for (dir1, dir2) in [('fwd', 'fwd'),
                             ('rev', 'rev'),
                             ('fwd', 'rev')]:
            good_cigar1 = good_cigars[dir1]
            good_cigar2 = good_cigars[dir2]
            
            read1 = read_pair[dir1]
            read2 = read_pair[dir2]

            # Iterate only over good cigars
            for ic1, (block_type1, block_len1) in enumerate(read1.cigar):
                if not good_cigar1[ic1]:
                    continue

                seq1 = np.array(list(read1.seq[poss[dir1][ic1]: poss[dir1][ic1] + block_len1]), 'S1')
                ref1 = ref[read1.pos: read1.pos + block_len1]

                # Look for mutations
                mutpos1 = (seq1 != ref1).nonzero()[0]
                if not len(mutpos1):
                    continue

                for ic1, (block_type2, block_len2) in enumerate(read2.cigar):
                    if not good_cigar2[ic2]:
                        continue
                
                    seq2 = np.array(list(read2.seq[poss[dir2][ic2]: poss[dir2][ic2] + block_len2]), 'S1')
                    ref2 = ref[read2.pos: read2.pos + block_len2]

                # Look for mutations
                mutpos1 = (seq1 != ref1).nonzero()[0]
                if not len(mutpos1):
                    

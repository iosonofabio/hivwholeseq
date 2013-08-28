# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/08/13
content:    Calculate the number of reads that cover a certain set of positions.
            This is used to normalize linkage analyses. For instance, suppose we
            find the pair of mutations A456T and G642T 10 times, how many reads
            actually covered both such that they could have given rise to such a
            signal?
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
from mapping.filenames import get_last_reference, get_last_mapped
from mapping.mapping_utils import get_ind_good_cigars



# Globals
VERBOSE = 1

# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']

maxreads = 5000000    #FIXME
match_len_min = 30
trim_bad_cigars = 3



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Extract linkage information')
    parser.add_argument('--adaID', metavar='00', type=int, required=True,
                        help='Adapter ID sample to analyze')
    # Example: coverage_tuples.py --adaID 02 '5 30 56' '58 43 89'
    parser.add_argument('tuples', nargs='+',
                        help="'fragment1 tuple1' ['fragment2 tuple2' ...]")
    args = parser.parse_args()
    adaID = args.adaID
    try:
        # Note: the first element is the fragment
        mtuples = map(lambda x: np.array(x.split(' '), int), args.tuples)
    except ValueError:
        raise ValueError('Tuple format not recognized. '\
                         "Example: "+os.path.basename(__file__)+\
                         " '2 56 78' '5 10 20 30'")

    # Open BAM
    bamfilename = get_last_mapped(data_folder, adaID, type='bam')
    # Try to convert to BAM if needed
    if not os.path.isfile(bamfilename):
        samfile = pysam.Samfile(bamfilename[:-3]+'sam', 'r')
        bamfile = pysam.Samfile(bamfilename, 'wb', template=samfile)
        for s in samfile:
            bamfile.write(s)
    bamfile = pysam.Samfile(bamfilename, 'rb')

    # Chromosome list
    chromosomes = bamfile.references
    
    # Read reference (fragmented)
    refseqs_raw = list(SeqIO.parse(get_last_reference(data_folder, adaID, ext=True),
                                   'fasta'))
    # Sort according to the chromosomal ordering
    refseqs = []
    for chromosome in chromosomes:
        for seq in refseqs_raw:
            if chromosome == seq.id:
                refseqs.append(seq)
                break
    refs = [np.array(refseq) for refseq in refseqs]

    # Prepare data structures
    coverage = np.zeros(len(mtuples), int)

    # Iterate over all pairs
    for i_pairs, reads in enumerate(pair_generator(bamfile)):

        # Limit to the first reads
        if 2 * i_pairs >= maxreads: break

        # Assign names
        read1 = reads[0]
        read2 = reads[1]

        # Check a few things to make sure we are looking at paired reads
        if read1.qname != read2.qname:
            raise ValueError('Read pair '+str(i_pairs)+': reads have different names!')

        # Ignore unmapped reads
        if (read1.is_unmapped or (not read1.is_proper_pair)
            or read2.is_unmapped or (not read2.is_proper_pair)):
            continue

        # Make a list of covered sites for the read_pair (some might be covered
        # more than once)
        covered_read = [np.zeros(len(x) - 1, int) for x in mtuples]
        for read in reads:

            # Find out on what chromosome the read has been mapped
            fragment = read.tid
            ref = refs[fragment]

            seq = read.seq
            good_cigar = get_ind_good_cigars(read.cigar,
                                             match_len_min=match_len_min)

            # The following two indices indicate the block position in the read
            # and in the reference sequence. Because of indels, they are updated
            # separately
            pos_read = 0
            pos_ref = read.pos

            # TODO: include indels as 'mutations'
            # TODO: include CIGAR trimming (we should really filter them out!)
            for (block_type, block_len), is_good in izip(read.cigar, good_cigar):
                # Match
                if block_type == 0:
                    if is_good:
                        for it, (mtuple, covt) in enumerate(izip(mtuples, covered_read)):
                            covt += (fragment == mtuple[0]) &\
                                    (pos_ref <= mtuple[1:]) &\
                                    (pos_ref + block_len > mtuple[1:])

                    pos_read += block_len
                    pos_ref += block_len

                # Deletion
                elif block_type == 2:
                    pos_ref += block_len

                # Insert
                elif block_type == 1:
                    pos_read += block_len

                # Other types of cigar?
                else:
                    raise ValueError('CIGAR type '+str(block_type)+' not recognized')

            # TODO: the problem is mismapped reads!


        # Check which tuples are fully covered
        for it, covt in enumerate(covered_read):
            coverage[it] += covt.min()

    print coverage 

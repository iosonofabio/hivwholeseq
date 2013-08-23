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
from mapping.filenames import get_last_reference, get_last_mapped
from mapping.mapping_utils import get_ind_good_cigars



# Globals
VERBOSE = 1

# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']

maxreads = 50000    #FIXME
match_len_min = 30
trim_bad_cigars = 3



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Extract linkage information')
    parser.add_argument('--adaID', metavar='00', type=int, required=True,
                        help='Adapter ID sample to analyze')
    args = parser.parse_args()
    adaID = args.adaID

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
    muts_all = []

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

        # If the reads are mapped to different fragments, that's mismapping
        if read1.tid != read2.tid:
            print read1.qname
            continue

        # Make a list of mutations for the read_pair
        muts = []
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
                        reftmp = ref[pos_ref: pos_ref + block_len]
                        seqtmp = seq[pos_read: pos_read + block_len]
                        seqtmp = np.array(list(seqtmp), 'S1')
                        mut_pos = (reftmp != seqtmp).nonzero()[0]

                        ## FIXME: this is mismapping at the beginning of the reference
                        ## (the insert length is wrong by 2 bases!)
                        #if read.qname == 'HWI-M01346:28:000000000-A53RP:1:1101:11993:2529':
                        #    import pdb; pdb.set_trace()

                        if len(mut_pos):
                            mut_der_all = seqtmp[mut_pos]
                            muts.extend(zip(mut_pos + pos_ref, mut_der_all))
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


        if len(muts):
            muts_all.append((read1.qname, fragment, muts))

    # Get rid of mismapped stuff (no read has 50 or more SNPs, not even in the)
    mismapped = [x[0] for x in muts_all if len(x[2]) > 50]
    #muts_all = [x for x in muts_all if len(x[1]) < 50]

    # Write results to file
    import cPickle as pickle
    mut_file = 'mutations.pickle'
    with open(data_folder+'subsample/'+foldername_adapter(adaID)+mut_file, 'w') as f:
        pickle.dump(muts_all, f, protocol=-1)
        
    # Now we can have a look at the statistics
    # 1. See the SFS
    from operator import *
    muts_flat = []
    for mut in muts_all:
        muts_flat.extend(((mut[1], m[0], m[1]) for m in mut[2]))
    from collections import Counter
    allele_counts = Counter(muts_flat)
    allele_counts_list = list(allele_counts.iteritems())
    allele_counts_list.sort(key=itemgetter(1), reverse=True)

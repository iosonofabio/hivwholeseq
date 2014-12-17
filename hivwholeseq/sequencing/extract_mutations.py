#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       21/08/13
content:    Call SNPs from reads.
'''
# TODO: finish up this mess!
# Modules
import os
import sys
import argparse
from itertools import izip
import pysam
import numpy as np
from Bio import SeqIO

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.adapter_info import load_adapter_table
from hivwholeseq.sequencing.filenames import get_last_reference, get_last_mapped
from hivwholeseq.mapping_utils import get_ind_good_cigars, pair_generator
from hivwholeseq.cluster.fork_cluster import fork_extract_mutations as fork_self



# Globals
maxreads = 50000
match_len_min = 30
trim_bad_cigars = 3



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Extract linkage information')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaID', metavar='00', required=True,
                        help='Adapter ID sample to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true', default=False,
                        help='Submit the job to the cluster via qsub')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')

    args = parser.parse_args()
    seq_run = args.run
    adaID = args.adaID
    VERBOSE = args.verbose
    submit = args.submit
    summary = args.summary

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # Branch to the cluster if required
    if submit:
        # If no adaID is specified, use all
        if adaID == 0:
            adaIDs = load_adapter_table(data_folder)['ID']
        else:
            adaIDs = [adaID]
            for adaID in adaIDs:
                fork_self(seq_run, adaID, VERBOSE=VERBOSE, summary=summary) 
        sys.exit()

    ###########################################################################
    # The actual script starts here
    ###########################################################################
    # Open BAM
    bamfilename = get_last_mapped(data_folder, adaID, type='bam', filtered=True)
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

        # Log
        if VERBOSE and (not (i_pairs % 10000)): print i_pairs,
    
        # Assign names
        read1 = reads[0]
        read2 = reads[1]
    
        # Check a few things to make sure we are looking at paired reads
        if read1.qname != read2.qname:
            raise ValueError('Read pair '+str(i_pairs)+': reads have different names!')
    
        # Ignore unmapped reads
        if (read1.is_unmapped or (not read1.is_proper_pair)
            or read2.is_unmapped or (not read2.is_proper_pair)):
            raise ValueError('Read pair '+str(i_pairs)+': not filtered properly')
    
        # If the reads are mapped to different fragments, that's mismapping
        if read1.tid != read2.tid:
            raise ValueError('Read pair '+str(i_pairs)+': not filtered properly')
    
        # Find out on what chromosome the read has been mapped
        fragment = read1.tid
        ref = refs[fragment]

        # Make a list of mutations for the read_pair
        muts = []
        for read in reads:
    
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
                #if read.is_read2:
                #    print pos_read, pos_ref
                # Match
                if block_type == 0:
                    if is_good:
                        reftmp = ref[pos_ref: pos_ref + block_len]
                        seqtmp = seq[pos_read: pos_read + block_len]
                        seqtmp = np.array(list(seqtmp), 'S1')
                        mut_pos = (reftmp != seqtmp).nonzero()[0]

                        #if (read.qname == 'HWI-M01346:28:000000000-A53RP:1:1101:5648:7688') and (read.is_read2):
                        #    raise ValueError
    

    
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
    
    
        # Guard against mismapped reads, which contain a lot of mutations
        if 0 < len(muts) <= 50:
            muts_all.append((read1.qname, fragment, muts))
        elif len(muts) > 50:
            raise ValueError('Read pair '+str(i_pairs)+': not filtered properly')

        if VERBOSE and (not (i_pairs % 10000)): print len(muts_all)
    
    
    
    ## Write results to file (~1 Gb per 1.5x10^6 reads in the patient sample)
    #import cPickle as pickle
    #mut_file = 'mutations.pickle'
    #with open(data_folder+foldername_adapter(adaID)+mut_file, 'w') as f:
    #    pickle.dump(muts_all, f, protocol=-1)

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
from collections import defaultdict
from itertools import izip
import pysam
import numpy as np
from Bio import SeqIO


# Horizontal import of modules from this folder
from mapping.adapter_info import load_adapter_table, foldername_adapter
from mapping.miseq import alpha, read_types, pair_generator
from mapping.filenames import get_last_reference, get_last_mapped
from mapping.mapping_utils import get_ind_good_cigars, convert_sam_to_bam



# Globals
VERBOSE = 1

# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']

maxreads = 500000000000    #FIXME
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
def fork_self(data_folder, adaID):
    '''Fork self for each adapter ID'''
    import subprocess as sp

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'exm_'+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--adaID', adaID,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Extract linkage information')
    parser.add_argument('--adaID', metavar='00', type=int, required=True,
                        help='Adapter ID sample to analyze')
    parser.add_argument('--submit', action='store_true', default=False,
                        help='Submit the job to the cluster via qsub')
    args = parser.parse_args()
    adaID = args.adaID
    submit = args.submit

    # Branch to the cluster if required
    if submit:
        # If no adaID is specified, use all
        if adaID == 0:
            adaIDs = load_adapter_table(data_folder)['ID']
        else:
            adaIDs = [adaID]
        for adaID in adaIDs:
            fork_self(data_folder, adaID) 
        sys.exit()

    ###########################################################################
    # The actual script starts here
    ###########################################################################
    # Open BAM
    bamfilename = get_last_mapped(data_folder, adaID, type='bam')
    # Try to convert to BAM if needed
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)
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
    
    # Open output BAM file
    bamoutfilename = get_last_mapped(data_folder, adaID, type='bam', filtered=True)
    bamoutfile = pysam.Samfile(bamoutfilename, 'wb', template=bamfile)

    # Iterate over all pairs
    n_unmapped = 0
    n_unpaired = 0
    n_mismapped_edge = 0
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
        if read1.is_unmapped or read2.is_unmapped:
            if VERBOSE >= 2:
                print 'Read pair '+read1.qname+': unmapped'
            n_unmapped += 1
            continue
    
        # Ignore not properly paired reads (this includes mates sitting on
        # different fragments)
        if (not read1.is_proper_pair) or (not read2.is_proper_pair):
            if VERBOSE >= 2:
                print 'Read pair '+read1.qname+': not properly paired'
            n_unpaired += 1
            continue
    
        # Find out on what chromosome the read has been mapped
        fragment = read1.tid
        ref = refs[fragment]

        # Mismappings are often characterized by a large number of mutations,
        # and/or an overhang before the beginning of the fragment or beyond 
        # its end
        muts = []
        skip = False
        for read in reads:
    

            # Check overhangs
            read_start = read.pos
            read_end = read.pos + sum(x[1] for x in read.cigar if x[0] != 1)
            if (((read_start == 0) and (read.cigar[0][0] == 1)) or
                ((read_end == len(ref)) and (read.cigar[-1][0] == 1))):
                skip = True
                n_mismapped_edge += 1
                break

        # Write the output
        if not skip:
            bamoutfile.write(read1)
            bamoutfile.write(read2)



#            # Check mutations
#            seq = read.seq
#            good_cigar = get_ind_good_cigars(read.cigar,
#                                             match_len_min=match_len_min)
#    
#            # The following two indices indicate the block position in the read
#            # and in the reference sequence. Because of indels, they are updated
#            # separately
#            pos_read = 0
#            pos_ref = read.pos
#    
#            # TODO: include indels as 'mutations'
#            # TODO: include CIGAR trimming (we should really filter them out!)
#            for (block_type, block_len), is_good in izip(read.cigar, good_cigar):
#                # Match
#                if block_type == 0:
#                    if is_good:
#                        reftmp = ref[pos_ref: pos_ref + block_len]
#                        seqtmp = seq[pos_read: pos_read + block_len]
#                        seqtmp = np.array(list(seqtmp), 'S1')
#                        mut_pos = (reftmp != seqtmp).nonzero()[0]
#    
#                        ## FIXME: this is mismapping at the beginning of the reference
#                        ## (the insert length is wrong by 2 bases!)
#                        #if read.qname == 'HWI-M01346:28:000000000-A53RP:1:1101:11993:2529':
#                        #    import pdb; pdb.set_trace()
#    
#                        if len(mut_pos):
#                            mut_der_all = seqtmp[mut_pos]
#                            muts.extend(zip(mut_pos + pos_ref, mut_der_all))
#                    pos_read += block_len
#                    pos_ref += block_len
#    
#                # Deletion
#                elif block_type == 2:
#                    pos_ref += block_len
#    
#                # Insert
#                elif block_type == 1:
#                    pos_read += block_len
#    
#                # Other types of cigar?
#                else:
#                    raise ValueError('CIGAR type '+str(block_type)+' not recognized')
#    
#    
#        if len(muts):
#            # Guard against mismapped reads, which appear as containing a lot of
#            # mutations
#            if len(muts) <= 50:
#                muts_all.append((read1.qname, fragment, muts))
#    
#        # Log
#        if VERBOSE and (not (i_pairs % 10000)): print i_pairs, len(muts_all)
#    
#    
#    # Get rid of mismapped stuff (no read has 50 or more SNPs, not even in the)
#    mismapped = [x[0] for x in muts_all if len(x[2]) > 50]
#    muts_all_red = [x for x in muts_all if len(x[2]) < 50]
#    
#    ## Write results to file (~1 Gb per 1.5x10^6 reads in the patient sample)
#    #import cPickle as pickle
#    #mut_file = 'mutations.pickle'
#    #with open(data_folder+foldername_adapter(adaID)+mut_file, 'w') as f:
#    #    pickle.dump(muts_all_red, f, protocol=-1)
#

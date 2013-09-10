#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Get the allele frequencies out of a BAM file and a reference.
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
from mapping.adapter_info import load_adapter_table
from mapping.miseq import alpha, read_types
from mapping.filenames import get_last_reference, get_last_mapped
from mapping.filenames import get_allele_counts_filename, get_insert_counts_filename, get_coverage_filename
from mapping.mapping_utils import get_ind_good_cigars, get_trims_from_good_cigars


# Globals
VERBOSE = 1

# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']

maxreads = 500000000000000    #FIXME
match_len_min = 30
trim_bad_cigars = 3



def get_allele_counts(refs, bamfile):
    '''Extract allele and insert counts from a bamfile'''
    
    # Prepare allele counts
    counts = [np.zeros((len(read_types), len(alpha), len(ref)), int)
              for ref in refs]
    
    # Insertions must be dealt with separately
    inserts = []
    for ref in refs:
        inserts.append([defaultdict(int) for k in read_types])

    # Count alleles and inserts
    # Note: the reads should already be filtered of unmapped stuff at this point
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

        # Fragment
        fragment = read.tid
        ref = refs[fragment]
        count = counts[fragment]
        insert = inserts[fragment]
    
        # Read CIGAR code for indels, and anayze each block separately
        # Note: some reads are weird ligations of HIV and contaminants
        # (e.g. fosmid plasmids). Those always have crazy CIGARs, because
        # only the HIV part maps. We hence trim away short indels from the
        # end of reads (this is still unbiased).
        cigar = read.cigar
        len_cigar = len(cigar)
        good_cigars = get_ind_good_cigars(cigar, match_len_min=match_len_min)
        trims = get_trims_from_good_cigars(good_cigars,
                                           trim_left=trim_bad_cigars,
                                           trim_right=trim_bad_cigars)

        # Iterate over CIGARs
        for ic, ((block_type, block_len),
                 good_cigar,
                 (trim_left, trim_right)) in enumerate(izip(cigar, good_cigars, trims)):

            # Inline block
            if block_type == 0:
                # Exclude bad CIGARs
                if good_cigar: 
                    # Get subsequence
                    seqb = np.array(list(seq[trim_left:block_len - trim_right]), 'S1')
                    # Increment counts
                    for j, a in enumerate(alpha):
                        posa = (seqb == a).nonzero()[0]
                        if len(posa):
                            count[js, j, pos + trim_left + posa] += 1
    
                # Chop off this block
                if ic != len_cigar - 1:
                    seq = seq[block_len:]
                    pos += block_len
    
            # Deletion
            elif block_type == 2:
                # Exclude bad CIGARs
                if good_cigar: 
                    # Increment gap counts
                    count[js, 4, pos:pos + block_len] += 1
    
                # Chop off pos, but not sequence
                pos += block_len
    
            # Insertion
            elif block_type == 1:
                # Exclude bad CIGARs
                if good_cigar: 
                    seqb = seq[:block_len]
                    insert[js][(pos, seqb)] += 1
    
                # Chop off seq, but not pos
                if ic != len_cigar - 1:
                    seq = seq[block_len:]
    
            # Other types of cigar?
            else:
                raise ValueError('CIGAR type '+str(block[0])+' not recognized')

    return counts, inserts



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
    
    # Get counts
    counts, inserts = get_allele_counts(refs, bamfile)

    # Get coverage
    coverage = [c.sum(axis=1) for c in counts]

    # Save counts and inserts to file
    with open(get_allele_counts_filename(data_folder, adaID), 'w') as f:
        pickle.dump(counts, f, protocol=-1)
    with open(get_coverage_filename(data_folder, adaID), 'w') as f:
        pickle.dump(coverage, f, protocol=-1)
    with open(get_insert_counts_filename(data_folder, adaID), 'w') as f:
        pickle.dump(inserts, f, protocol=-1)

    ## Raw estimate of diversity
    #coverage = counts.sum(axis=0).sum(axis=0)
    #nu_major = counts.sum(axis=0).max(axis=0) / (coverage + 0.000001)
    #nu = counts.sum(axis=0) / (coverage + 0.000001)

    ### Plot info about the major allele
    #import matplotlib.pyplot as plt
    ##fig, axs = plt.subplots(2, 1, figsize=(14.3, 11.9))
    ##plt.sca(axs[0])
    ##plt.plot(nu_major, lw=2, c='k')
    ##plt.xlabel('position (b.p.)')
    ##plt.ylabel(r'$\nu_{maj}$')
    ##plt.ylim(-0.05, 1.05)
    ### Maximal entropy line
    ##plt.plot([0, len(ref)], [1.0 / len(alpha)] * 2, lw=1.5, ls='--')
    ##plt.title('Major allele frequency along the HIV genome')
    ### Plot histogram of major allele frequency
    ##h = np.histogram(nu_major, bins=np.linspace(1.0 / len(alpha), 1, 20),
    ##                 density=True)
    ##x = 0.5 * (h[1][1:] + h[1][:-1])
    ##y = h[0]
    ##plt.sca(axs[1])
    ##plt.plot(x, y + 1e-6, lw=1.5, c='k')
    ##plt.xlabel(r'$\nu_{maj}$')
    ##plt.title('Distribution of major allele frequency')
    ##plt.yscale('log')
    ##plt.tight_layout()

    ## Per read-type analysis
    #fig, axs = plt.subplots(2, 1, figsize=(14.3, 11.9))
    #import matplotlib.cm as cm
    #for i, (read_type, count, insert) in enumerate(izip(read_types, counts, inserts)):
    #    coverage = count.sum(axis=0)
    #    nu_major = count.max(axis=0) / (coverage + 0.000001)
    #    nu = count / (coverage + 0.000001)

    #    # Plot info about the major allele
    #    plt.sca(axs[0])
    #    plt.plot(nu_major + len(read_types) - i - 1, lw=2, c=cm.jet(int(255.0 * i / len(read_types))))
    #    plt.plot([0, len(ref)], [1.0 / len(alpha) + len(read_types) - i - 1] * 2, lw=1.5, ls='--', c='k')
    #    plt.xlabel('position (b.p.)')
    #    plt.ylabel(r'$\nu_{maj}$')
    #    plt.ylim(-0.05, 4.05)
    #    # Maximal entropy line
    #    plt.title('Major allele frequency along the HIV genome')
    #    # Plot histogram of major allele frequency
    #    h = np.histogram(nu_major, bins=np.linspace(1.0 / len(alpha), 1, 20),
    #                     density=True)
    #    x = 0.5 * (h[1][1:] + h[1][:-1])
    #    y = h[0]
    #    plt.sca(axs[1])
    #    plt.plot(x, y + 1e-4, lw=1.5,
    #             c=cm.jet(int(255.0 * i / len(read_types))),
    #             label=read_type)
    #    plt.xlabel(r'$\nu_{maj}$')
    #    plt.title('Distribution of major allele frequency')
    #    plt.yscale('log')

    #plt.legend(loc=2)
    #plt.tight_layout()


    #plt.ion()
    #plt.show()

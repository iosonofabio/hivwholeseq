#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Get the allele frequencies out of a SAM/BAM file and a reference.
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



# Globals
VERBOSE = 1
alpha = np.array(list('ACGT-N'), 'S1')
read_types = ['read1 f', 'read1 r', 'read2 f', 'read2 r']
maxreads = 500000000000000000    #FIXME
match_len_min = 30
trim_bad_cigars = 3
ref_filename = 'consensus_filtered_trimmed.fasta'
bam_filename = 'mapped_to_self_filtered_trimmed.bam'
allele_count_filename = 'allele_counts.npy'
insert_count_filename = 'insert_counts.pickle'



def get_allele_counts(ref, bamfile):
    '''Extract allele and insert counts from a bamfile'''
    
    # Prepare allele counts
    counts = np.zeros((len(read_types), len(alpha), len(ref)), int)
    
    # Insertions must be dealt with separately
    inserts = [defaultdict(int) for k in read_types]

    # Count alleles
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

    return counts, inserts



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
    
    # Get counts
    counts, inserts = get_allele_counts(ref, bamfile)

    # Save counts and inserts to file
    np.save(data_folder+allele_count_filename, counts)
    with open(data_folder+insert_count_filename, 'w') as f:
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

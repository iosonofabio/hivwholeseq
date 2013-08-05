# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Get the allele frequencies out of a SAM/BAM file and a reference.
'''
# Modules
import os
import sys
from collections import defaultdict
import pysam
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO



# Globals
VERBOSE = 1
alpha = np.array(list('ACGT-N'), 'S1')
read_types = ['read1 f', 'read1 r', 'read2 f', 'read2 r']
maxreads = 20000
ref_filename = 'consensus_filtered_trimmed.fasta'
bam_filename = 'mapped_to_self_filtered_trimmed.bam'



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
    
        # Ignore unmapped reads
        if read.is_unmapped:
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
        cigar = read.cigar
        len_cig = len(cigar)
        for ic, block in enumerate(cigar):
            # Inline block
            if block[0] == 0:
                seqb = np.array(list(seq[:block[1]]), 'S1')
                # Increment counts
                for j, a in enumerate(alpha):
                    posa = (seqb == a).nonzero()[0]
                    if not len(posa):
                        continue
                    counts[js, j, pos + posa] += 1
    
                # Chop off this block
                if ic != len_cig - 1:
                    seq = seq[block[1]:]
                    pos += block[1]
    
            # Deletion
            elif block[0] == 2:
                # Increment gap counts
                counts[js, 4, pos:pos + block[1]] += 1
    
                # Chop off pos, but not sequence
                pos += block[1]
    
            # Insertion
            elif block[0] == 1:
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
    ref_file = data_folder+ref_filename
    bamfilename = data_folder+bam_filename

    # Read reference
    refseq = SeqIO.read(ref_file, 'fasta')
    ref = np.array(refseq)

    # Open BAM
    bamfile = pysam.Samfile(bamfilename, 'rb')
    
    # Get counts
    counts, inserts = get_allele_counts(ref, bamfile)

    # Raw estimate of diversity
    coverage = counts.sum(axis=0).sum(axis=0)
    nu_major = counts.sum(axis=0).max(axis=0) / (coverage + 0.000001)
    nu = counts.sum(axis=0) / (coverage + 0.000001)

    # Plot info about the major allele
    fig, axs = plt.subplots(2, 1, figsize=(14.3, 11.9))
    plt.sca(axs[0])
    plt.plot(nu_major, lw=2, c='k')
    plt.xlabel('position (b.p.)')
    plt.ylabel(r'$\nu_{maj}$')
    plt.ylim(-0.05, 1.05)
    # Maximal entropy line
    plt.plot([0, len(ref)], [1.0 / len(alpha)] * 2, lw=1.5, ls='--')
    plt.title('Major allele frequency along the HIV genome')

    # Plot histogram of major allele frequency
    h = np.histogram(nu_major, bins=np.linspace(1.0 / len(alpha), 1, 20),
                     density=True)
    x = 0.5 * (h[1][1:] + h[1][:-1])
    y = h[0]
    plt.sca(axs[1])
    plt.plot(x, y + 1e-6, lw=1.5, c='k')
    plt.xlabel(r'$\nu_{maj}$')
    plt.title('Distribution of major allele frequency')
    plt.yscale('log')

    plt.tight_layout()
    plt.ion()
    plt.show()

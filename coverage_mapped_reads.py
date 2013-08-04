# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/08/13
content:    Check the coverage of mapped reads on the HIV genome.
'''
# Modules
import os
import sys
from collections import defaultdict
from operator import itemgetter
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pysam
from map_HIV_HXB2 import load_adapter_table


# Globals
VERBOSE = 1
maxreads = 20000
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/subsample/'
adapters_table_file = 'adapters_table.dat'
ref_file = 'HXB2.fa.gz'
mapped_file_sam = 'mapped_to_HXB2.sam'
mapped_file_bam = 'mapped_to_HXB2.bam'
consensus_file = 'consensus.fasta'



# Script
if __name__ == '__main__':

    # Directory to read
    adapter_table = load_adapter_table(data_folder)
    adaID = 2
    dirname = 'adapterID_'+'{:02d}'.format(adaID)+'/'

    # Convert SAM (output of stampy) to BAM (more efficient)
    if not os.path.isfile(data_folder+dirname+mapped_file_bam):
        samfile = pysam.Samfile(data_folder+dirname+mapped_file_sam, 'r')
        bamfile = pysam.Samfile(data_folder+dirname+mapped_file_bam, 'wb', template=samfile)
        for s in samfile:
            bamfile.write(s)

    # Open BAM
    bamfile = pysam.Samfile(data_folder+dirname+mapped_file_bam, 'rb')

    # Read reference
    with gzip.open(ref_file, "r") as ref_fileh:
        refseq = SeqIO.read(ref_fileh, 'fasta')
    ref = np.array(refseq)

    # Prepare allele counts
    alpha = np.array(list('ACGT-N'), 'S1')
    read_types = ['read1 f', 'read1 r', 'read2 f', 'read2 r']
    counts = [np.zeros((len(alpha), len(refseq)), int) for k in read_types]

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
        if read.is_read1:
            js = 0
        else:
            js = 2
        if read.is_reverse:
            js += 1

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
                    counts[js][j, pos + posa] += 1

                # Chop off this block
                if ic != len_cig - 1:
                    seq = seq[block[1]:]
                    pos += block[1]

            # Deletion
            elif block[0] == 2:
                # Increment gap counts
                counts[js][4, pos:pos + block[1]] += 1

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

    # Make consensi for each of the four categories
    consensi = []
    for count in counts:
        # Positions without reads are considered N
        # (this should happen only at the ends)
        count.T[(count).sum(axis=0) == 0] = np.array([0, 0, 0, 0, 0, 1])
        consensi.append(alpha[count.argmax(axis=0)])

    # Add inserts that are present in more than half the reads
    inserts_consensi = [[k for (k, i) in insert.iteritems()
                         if ((i > 10) and
                             (i > 0.5 * counts[js][:, max(0, k[0]-5):min(len(ref), k[0]+5)].sum(axis=0).mean()))]
                        for js, insert in enumerate(inserts)]

    # Make final consensus
    count_all = np.sum(counts, axis=0)
    consensus = alpha[count_all.argmax(axis=0)]
    keys_all = set([k for insert in inserts for k in insert])
    inserts_consensus = [k for k in keys_all
                         if sum([insert[k] for insert in inserts]) > 0.5 * count_all[:, k[0]].sum()]
    consensus_final = []
    pos = 0
    for k in inserts_consensus:
        consensus_final.append(''.join(consensus[pos:k[0]]))
        consensus_final.append(k[1])
        pos = k[0]
    consensus_final.append(''.join(consensus[pos:]))
    consensus_final = ''.join(consensus_final)

    # Print divergence to reference excluding inserts
    print 'Divergence to ref:', (ref != consensus)[consensus != 'N'].mean()

    # Save consensus
    consensusseq = SeqRecord(Seq(consensus_final),
                             id='{:02d}'.format(adaID)+'_consensus',
                             name='{:02d}'.format(adaID)+'_consensus')
    SeqIO.write(consensusseq, data_folder+dirname+consensus_file, 'fasta')


# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/08/13
content:    Check the coverage of mapped reads on the HIV genome.
'''
# Modules
import os
import sys
import re
from collections import defaultdict
from collections import Counter
from itertools import izip
from operator import itemgetter
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pysam
from map_HIV_HXB2 import load_adapter_table
import matplotlib.pyplot as plt



# Globals
VERBOSE = 1
maxreads = 20000
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/subsample/'
adapters_table_file = 'adapters_table.dat'
ref_file = 'HXB2.fa.gz'
mapped_file_sam = 'mapped_to_HXB2.sam'
mapped_file_bam = 'mapped_to_HXB2.bam'
consensus_file = 'consensus.fasta'
coverage_file = 'coverage.pdf'



# Script
if __name__ == '__main__':

    # Directory to read
    adapter_table = load_adapter_table(data_folder)
    
    # Iterate for all adapterIDs
    for adaID in adapter_table['ID']:
        dirname = 'adapterID_'+'{:02d}'.format(adaID)+'/'
    
        # Convert SAM (output of stampy) to BAM (more efficient)
        if not os.path.isfile(data_folder+dirname+mapped_file_bam):
            print 'Converting SAM to BAM...',
            samfile = pysam.Samfile(data_folder+dirname+mapped_file_sam, 'r')
            bamfile = pysam.Samfile(data_folder+dirname+mapped_file_bam, 'wb', template=samfile)
            for s in samfile:
                bamfile.write(s)
            print 'DONE.'
    
        # Open BAM
        bamfile = pysam.Samfile(data_folder+dirname+mapped_file_bam, 'rb')
    
        # Read reference
        with gzip.open(ref_file, "r") as ref_fileh:
            refseq = SeqIO.read(ref_fileh, 'fasta')
        ref = np.array(refseq)
    
        # Prepare allele counts
        alpha = np.array(list('ACGT-N'), 'S1')
        read_types = ['read1 f', 'read1 r', 'read2 f', 'read2 r']
        counts = np.zeros((len(read_types), len(alpha), len(refseq)), int)
    
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
    
        # Make consensi for each of the four categories
        consensi = np.zeros((len(read_types), len(refseq)), 'S1')
        for js, count in enumerate(counts):
            # Positions without reads are considered N
            # (this should happen only at the ends)
            count.T[(count).sum(axis=0) == 0] = np.array([0, 0, 0, 0, 0, 1])
            consensi[js] = alpha[count.argmax(axis=0)]
    
        # Add inserts that are present in more than half the reads
        inserts_consensi = [[k for (k, i) in insert.iteritems()
                             if ((i > 10) and
                                 (i > 0.5 * counts[js, :, max(0, k[0]-5):min(len(ref), k[0]+5)].sum(axis=0).mean()))]
                            for js, insert in enumerate(inserts)]
    
        # Make final consensus
        # This has two steps: 1. counts; 2. inserts
        # We should check that different read types agree (e.g. forward/reverse) on
        # both counts and inserts if they are covered
        # 1.0: Prepare everything as 'N': ambiguous
        consensus = np.repeat('N', len(refseq))
        # 1.1: Put safe stuff
        ind_agree = (consensi == consensi[0]).all(axis=0)
        consensus[ind_agree] = consensi[0, ind_agree]
        # 1.2: Ambiguous stuff requires more care
        polymorphic = []
        # 1.2.1: If is covered by some read types only, look among those
        amb_pos = (-ind_agree).nonzero()[0]
        for pos in amb_pos:
            cons_pos = consensi[:, pos]
            cons_pos = cons_pos[cons_pos != 'N']
            # If there is unanimity, ok
            if (cons_pos == cons_pos[0]).all():
                consensus[pos] = cons_pos[0]
            # else, assign only likely mismapped deletions
            else:
                # Restrict to nongapped things (caused by bad mapping)
                cons_pos = cons_pos[cons_pos != '-']
                if (cons_pos == cons_pos[0]).all():
                    consensus[pos] = cons_pos[0]
                # In case of polymorphisms, take any most abundant nucleotide
                else:
                    polymorphic.append(pos)
                    tmp = zip(*Counter(cons_pos).items())
                    consensus[pos] = tmp[0][np.argmax(tmp[1])]
    
        # 2.0 get all inserts
        insert_names = set()
        for insert in inserts:
            insert_names |= set(insert.keys())
        # 2.1 iterate over inserts and check all read types
        insert_consensus = []
        for insert_name in insert_names:
            # Get counts for the insert and local coverage
            ins_counts = np.array([insert[insert_name] for insert in inserts])
            cov_loc = counts[:, :, insert_name[0]: min(len(refseq), insert_name[0] + 2)].mean(axis=2).sum(axis=1)
            # Get read types in which the insert is called
            types_ins_called = ins_counts > 0.5 * cov_loc
            # Get read types which actually cover this region
            types_cov = cov_loc > 3
            # Check the insert counts compared to the local coverage
            # if any read type has coverage and all read types with coverage agree, ok
            if types_cov.sum() and (types_ins_called[types_cov]).all():
                insert_consensus.append(insert_name)
        insert_consensus.sort(key=itemgetter(0))
        # 3. put inserts in
        consensus_final = []
        pos = 0
        for insert_name in insert_consensus:
            # Indices should be fine...
            # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
            # THEN the insert, FINALLY comes seq[391:]
            consensus_final.append(''.join(consensus[pos:insert_name[0]]))
            consensus_final.append(insert_name[1])
            pos = insert_name[0]
        consensus_final.append(''.join(consensus[pos:]))
        # Trip initial and final Ns, and strip gaps
        consensus_final = ''.join(consensus_final).strip('N')
        consensus_final = re.sub('-', '', consensus_final)
    
        # Plot coverage along the genome
        plt.figure(figsize=(14, 10))
        for read_type, count in izip(read_types, counts):
            plt.plot(count.sum(axis=0), lw=1.5, label=read_type)
        plt.xlabel('Position (HXB2)')
        plt.ylabel('# reads (tot '+str(maxreads)+')')
        plt.legend(loc=1)
        plt.title('Coverage for adapterID '+'{:02d}'.format(adaID))
        plt.tight_layout()
        plt.savefig(data_folder+dirname+coverage_file)
    
        #plt.ion()
        #plt.show()
    
        # Print divergence to reference excluding inserts
        print 'Adapter ID'+'{:02d}'.format(adaID)+': divergence to ref:', (ref != consensus)[consensus != 'N'].mean()
    
        # Save consensus
        consensusseq = SeqRecord(Seq(consensus_final),
                                 id='{:02d}'.format(adaID)+'_consensus',
                                 name='{:02d}'.format(adaID)+'_consensus')
        SeqIO.write(consensusseq, data_folder+dirname+consensus_file, 'fasta')
    
    

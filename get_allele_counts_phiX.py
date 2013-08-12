import os
import sys
import cPickle as pickle
from collections import defaultdict
from itertools import izip
import numpy as np
from Bio import SeqIO
import pysam
import get_allele_counts


# Globals
VERBOSE = 1
alpha = np.array(list('ACGT-N'), 'S1')
read_types = ['read1 f', 'read1 r', 'read2 f', 'read2 r']
maxreads = 500000000000000000    #FIXME
match_len_min = 30
trim_bad_cigars = 3
ref_file = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/phiX_genome.fasta'
bam_filename = 'mapped_to_phiX_filtered_trimmed.bam'
allele_count_filename = 'allele_counts.npy'
insert_count_filename = 'insert_counts.pickle'

# Script
if __name__ == '__main__':

    # Input arguments
    data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/unclassified_reads/'

    # Read reference
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
    counts, inserts = get_allele_counts.get_allele_counts(ref, bamfile)

    # Save counts and inserts to file
    np.save(data_folder+allele_count_filename, counts)
    with open(data_folder+insert_count_filename, 'w') as f:
        pickle.dump(inserts, f, protocol=-1)

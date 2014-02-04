# vim: fdm=marker
'''
author:     Fabio Zanini
date:       30/01/14
content:    Check the mapped reads.
'''
# Modules
import os
import argparse
from itertools import izip
from Bio import SeqIO
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hivwholeseq_pacbio.samples import samples
from hivwholeseq_pacbio.map_reads import Read
from hivwholeseq_pacbio.coverage import load_mapped_reads


# Globals
data_folder = '/ebio/ag-neher/share/data/PacBio_HIV_Karolinska/run23/'



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get coverage of PacBio reads')
    parser.add_argument('--sample', required=True,
                        help='Sample to analyze (e.g. S1)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    samplename = args.sample
    VERBOSE = args.verbose

    # Parse consensus reads
    sample = samples.set_index('name').loc[samplename]

    reads = load_mapped_reads(data_folder, samplename)
    
    # Get NL4-3 reference
    from hivwholeseq.reference import load_NL43  
    refseq = load_NL43()
    refm = np.array(refseq)
    refs = refm.tostring()

    # Start checking random reads
    n_reads = len(reads)
    for i in xrange(30):
        ind = np.random.randint(n_reads)
        read = reads[ind]

        if read is None:
            continue

        pos_ref = read.pos
        pos_read = 0
        ali_ref = []
        ali_read = []
        for (bt, bl) in read.cigar:
            if bt == 0:
                ali_ref.append(refs[pos_ref: pos_ref + bl])
                ali_read.append(read.seq[pos_read: pos_read + bl])
                pos_ref += bl
                pos_read += bl

            elif bt == 1:
                ali_ref.append('-' * bl)
                ali_read.append(read.seq[pos_read: pos_read + bl])
                pos_read += bl
                
            else:
                ali_ref.append(refs[pos_ref: pos_ref + bl])
                ali_read.append('-' * bl)
                pos_ref += bl

        print ''.join(ali_ref)[:50]
        print ''.join(ali_read)[:50]
        print

            

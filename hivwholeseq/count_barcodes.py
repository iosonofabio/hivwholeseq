# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       02/08/13
content:    Count barcode abundance. We expect a high abundance for the true ones,
            plus sequencing errors around them.
'''
# Modules
import argparse
from collections import defaultdict
from Bio import SeqIO
import matplotlib.pylab as plt

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.filenames import get_raw_read_files


# Globals
maxreads = -1



# Functions
def count_barcodes(dataset, VERBOSE=0):
    '''Count the abundance of each barcode'''

    # Get the read filenames
    data_filenames = get_raw_read_files(dataset)
    datafile = data_filenames['adapter']

    # Count the abundance of each barcode
    bc_counts = defaultdict(int)
    rc = 0
    with open(datafile, 'r') as infile:
        for read in SeqIO.parse(infile, 'fastq'):
            bc_counts[read.seq.tostring()] += 1
            rc += 1
            if rc == maxreads:
                break
    
    print sorted(bc_counts.items(), key=lambda x:x[1], reverse=True)[:20]
    
    # Plot results
    plt.figure()
    ax=plt.subplot(111)
    plt.plot(range(1,len(bc_counts)+1), sorted(bc_counts.values(), reverse=True))
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.xlabel('barcode rank')
    plt.ylabel('abundance')

    plt.ion()
    plt.show()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    miseq_run = args.run
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]

    # Count the barcodes
    count_barcodes(dataset, VERBOSE=VERBOSE)

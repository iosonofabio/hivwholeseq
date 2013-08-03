# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       01/08/13
content:    Plot the distribution of read lengths as a first quality assessment
            for the Illumina sequencing on the MiSeq.
'''
# Modules
from Bio import SeqIO
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


# Globals
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'
datafile_read1 = data_folder+'copies/lane1_NoIndex_L001_R1_001.fastq'
datafile_read2 = data_folder+'copies/lane1_NoIndex_L001_R3_001.fastq'

# Quality threshold
phred_min = 20

# Max reads
max_reads = 30000



# Script
if __name__ == '__main__':
    
    bins = np.arange(255) - 0.5

    # Repeat both for read 1 and read 2
    datafiles = {'read 1': datafile_read1,
                 'read 2': datafile_read2}
    histograms = {}
    for readname, datafile in datafiles.iteritems():
    
        # Result data structures
        read_histogram = np.zeros(len(bins) - 1, int)
    
        # Read data
        with open(datafile, 'r') as f: 
            seq_iter = SeqIO.parse(f, 'fastq')
    
            for i, seq in enumerate(seq_iter):
    
                # Analyze only so many reads
                if i >= max_reads:
                    break
    
                phred = np.asarray(seq.letter_annotations['phred_quality'], int)
                ind = phred >= phred_min
    
                # If the read is good, take its length
                if ind.all():
                    read_histogram[len(seq)] += 1
                # else take the length until the first bad spot
                else:
                    read_histogram[(-ind).nonzero()[0][0]] += 1

        # Store data
        histograms[readname] = read_histogram
            
        # Plot
        center_bins = (bins[:-1]+bins[1:]) * 0.5
        n_reads = sum(read_histogram)
        mean_length = 1.0 * sum(center_bins * read_histogram) / n_reads
        var_length = 1.0 * sum(center_bins**2 * read_histogram) / n_reads - mean_length**2
        lb = '#reads: '+str(n_reads)
        plt.plot(center_bins, read_histogram, lw=1.5,
                 label=readname+': '+str(100 * read_histogram[150:].sum() / read_histogram.sum())+'% > 150 b.p.')

    plt.title(lb)
    plt.xlabel('readlength')
    plt.xlim(-1, 255)
    plt.ylim(1, plt.ylim()[1] * 1.1)
    plt.yscale('log')
    plt.legend(loc='upper center')
    plt.savefig('readlength.pdf')
    
    plt.ion()
    plt.show()

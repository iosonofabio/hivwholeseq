# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       01/08/13
content:    Plot the distribution of read lengths as a first quality assessment
            for the Illumina sequencing on the MiSeq.
'''
# Modules
import numpy as np
import scipy as sp
import pylab as pl


# Globals
data_folder = '/ebio/ag-neher/share/data/HIV_illumina_data/run28/'
datafile = data_folder+'s_8_2_sequence.txt'



# Script
if __name__ == '__main__':
    
    # Read data
    with open(datafile, 'r') as f:
        count = 0
        count_max = 10000
        outer_count = 0
        
        read_histogram = np.zeros(200)
        bins = np.arange(201) - 0.5
        
        line = f.readline()
        while line != '' and outer_count < 20000:
            read_length = []
            count = 0
            while count < count_max and line != '':
                if line[0] == '+':
                    line = f.readline().strip()
                    try:
                        read_length.append(line.index('B'))
                    except:
                        read_length.append(len(line))
                    count += 1
                else:
                    line = f.readline()
            
            y,x = sp.histogram(read_length, bins=bins)
            read_histogram += y
            line = f.readline()
            outer_count += 1
            print outer_count
        
        
    # Plot
    center_bins = (bins[:-1]+bins[1:]) * 0.5
    n_reads = sum(read_histogram)
    mean_length = sum(center_bins * read_histogram) / n_reads
    var_length = sum(center_bins**2 * read_histogram) / n_reads-mean_length**2
    lb = '#reads: '+str(n_reads)
    lb += '  mean: '+str(round(mean_length, 2))
    lb += '  std: '+str(round(np.sqrt(var_length), 2))
    pl.plot(center_bins, read_histogram)
    pl.title(lb)
    pl.xlabel('readlength')
    pl.savefig('read2_readlength.pdf')
    
    pl.ion()
    pl.show()

#!/usr/bin/env python
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
import matplotlib
from matplotlib import cm
import sys
import os

params = {'backend': 'pdf',  
          'axes.labelsize': 20, 
          'text.fontsize': 20,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 18,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': False}
matplotlib.rcParams.update(params)
import matplotlib.pyplot as plt


# Globals
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'

# Quality threshold
phred_min = 20

# Max reads
max_reads = 5000000
# length of reads and list of positions (needed for rapid updating of the phred distribution)
L=250
positions = np.arange(L)

# Script
if __name__ == '__main__':
    
    bins = np.arange(L+2) - 0.5
    if len(sys.argv)>1:
        barcode = sys.argv[1]
        if barcode[-1]!='/':
            barcode+='/'
    else:
        barcode = 'adapterID_02/'
        
    outdir= os.getcwd()+'/q_control_'+barcode[:-1]
    if not os.path.exists(outdir):
        try:
            os.mkdir(outdir)
        except:
            print "cannot make directory:",outdir
            sys.exit(1)

    # Repeat both for read 1 and read 2
    datafiles = {'read 1': data_folder+barcode+'read1.fastq',
                 'read 2': data_folder+barcode+'read2.fastq'}
    histograms = {}


    for readname, datafile in datafiles.iteritems():
    
        # Result data structures
        longest_good_block = np.zeros(L+1, int)
        first_bad_nucleotide = np.zeros(L+1, int)
        phred_score_dis = np.zeros((50,L))
    
        # Read data
        with open(datafile, 'r') as f: 
            seq_iter = SeqIO.parse(f, 'fastq')
    
            for i, seq in enumerate(seq_iter):
    
                # Analyze only so many reads
                if i >= max_reads:
                    break
    
                # read the phred score, make int array
                phred = np.asarray(seq.letter_annotations['phred_quality'], int)
                # update the phred score distribution
                phred_score_dis[(phred,positions)]+=1

                # get all positions above the cut-off
                ind = np.asarray(phred >= phred_min, int)
                #convert to str 00111111011111 and find longest stretch of 1s
                good_blocks=np.asarray([len(x) for x in "".join(map(str,ind)).split('0')],int)
                lgb= np.max(good_blocks)
                longest_good_block[lgb]+=1
                
                # If the read is good, take its length
                if lgb==len(seq):
                    first_bad_nucleotide[len(seq)] += 1
                # else take the length until the first bad spot
                else:
                    first_bad_nucleotide[(1-ind).nonzero()[0][0]] += 1

        # Store data
        histograms[readname] = longest_good_block, first_bad_nucleotide, phred_score_dis



    # Plot distribution of the first bad nucleotide for each data file
    fig_first_bad = plt.figure()
    ax1_first_bad = plt.subplot(111)
    center_bins = (bins[:-1]+bins[1:]) * 0.5
    for readname in datafiles:
        fbn = histograms[readname][1]
        n_reads = sum(fbn)
        mean_length = 1.0 * sum(center_bins * fbn) / n_reads
        var_length = 1.0 * sum(center_bins**2 * fbn) / n_reads - mean_length**2
        lb = '#reads: '+str(n_reads)
        plt.plot(center_bins, fbn, lw=1.5,
                 label=readname+': '+str(100 * fbn[150:].sum() / fbn.sum())+'% > 150 b.p.')
        
    plt.title(lb)
    plt.xlabel('position of first bad nucleotide')
    plt.xlim(-1, 255)
    plt.ylim(1, plt.ylim()[1] * 1.1)
    plt.yscale('log')
    plt.legend(loc='upper center')
    plt.savefig(outdir+'/readlength_fbn.pdf')


    # plot the distribution of the longest block with continuous quality above cutoff
    fig_blocks = plt.figure()
    ax1_blocks = plt.subplot(111)
    for readname in datafiles:
        lgb = histograms[readname][0]
        n_reads = sum(lgb)
        mean_length = 1.0 * sum(center_bins * lgb) / n_reads
        var_length = 1.0 * sum(center_bins**2 * lgb) / n_reads - mean_length**2
        lb = '#reads: '+str(n_reads)
        plt.plot(center_bins, lgb, lw=1.5,
                 label=readname+': '+str(100 * lgb[150:].sum() / lgb.sum())+'% > 150 b.p.')

         
    plt.title(lb)
    plt.xlabel('length of longest good block')
    plt.xlim(-1, 255)
    plt.ylim(1, plt.ylim()[1] * 1.1)
    plt.yscale('log')
    plt.legend(loc='upper center')
    plt.savefig(outdir+'/readlength_lgb.pdf')

    # make a figure on the cumulative phred score distribution 
    # at different positions in the read
    fig_cum_phred_score = plt.figure(figsize = (12,6))
    for ri,readname in enumerate(datafiles):
        cum_phred_score_dis = histograms[readname][2].cumsum(axis=0)
        plt.subplot(1,2,ri)
        plt.title(readname)
        for pos in positions[::5]:
            plt.plot(1.0-cum_phred_score_dis[:,pos]/cum_phred_score_dis[-1,pos], c=cm.jet(pos))

        sm = plt.cm.ScalarMappable(cmap=cm.jet, norm=plt.normalize(0, 255))
        sm._A = []
        cbar = plt.colorbar(sm)
        plt.xlabel('phred score q')

    plt.ylabel('fraction better than q')
    plt.savefig(outdir+'/quality_along_read.pdf')
   

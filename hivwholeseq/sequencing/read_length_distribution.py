#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       01/08/13
content:    Plot the distribution of read lengths as a first quality assessment
            for the Illumina sequencing on the MiSeq.
'''
# Modules
import argparse
from Bio import SeqIO
import numpy as np
import matplotlib
from matplotlib import cm
import sys
import os

params = {'axes.labelsize': 20, 
          'text.fontsize': 20,
          'font.sans-serif': 'Helvetica',
          'legend.fontsize': 18,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'text.usetex': False}
matplotlib.rcParams.update(params)
import matplotlib.pyplot as plt

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.filenames import get_read_filenames
from hivwholeseq.sequencing.adapter_info import load_adapter_table



# Globals
# Quality threshold
phred_min = 20

# Max reads
max_reads = 5000
# length of reads and list of positions (needed for rapid updating of the phred distribution)
L=250
positions = np.arange(L)
bins = np.arange(L+2) - 0.5



# Script
if __name__ == '__main__':

    # Parse input arguments
    parser = argparse.ArgumentParser(description='Filter & trim demultiplexed reads.')
    parser.add_argument('--run', type=int, required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--verbose', type=int, default=0,
                        help=('Verbosity level [0-3]'))
    parser.add_argument('--submit', action='store_true', default=False,
                        help='Submit the job to the cluster via qsub')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    VERBOSE = args.verbose
    submit = args.submit

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Iterate over adaIDs
    for adaID in adaIDs:

        outdir = os.getcwd()+'/q_control_'+adaID
        if not os.path.exists(outdir):
            try:
                os.mkdir(outdir)
            except:
                print "cannot make directory:",outdir
                sys.exit(1)
    
        # Repeat both for read 1 and read 2
        datafiles = get_read_filenames(data_folder, adaID, filtered=False)
        datafiles = {'read 1': datafiles[0],
                     'read 2': datafiles[1]}
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
    
        # Plot
        # Prepare four plots: cumulative read1, read2, length of longest block,
        # and first bad nucleotide
        fig, axs = plt.subplots(2, 2, figsize=(11, 11))
        axs = axs.ravel()
    
        # Plot distribution of the first bad nucleotide for each data file
        ax = axs[3]
        center_bins = (bins[:-1]+bins[1:]) * 0.5
        for readname in datafiles:
            fbn = histograms[readname][1]
            n_reads = sum(fbn)
            mean_length = 1.0 * sum(center_bins * fbn) / n_reads
            var_length = 1.0 * sum(center_bins**2 * fbn) / n_reads - mean_length**2
            lb = '#reads: '+str(n_reads)
            ax.plot(center_bins, fbn, lw=1.5,
                    label=readname+': '+str(100 * fbn[150:].sum() / fbn.sum())+'% > 150 b.p.')
        ax.set_title(lb)
        ax.set_xlabel('position of first bad nucleotide')
        ax.set_xlim(-1, 255)
        ax.set_ylim(1, plt.ylim()[1] * 1.1)
        ax.set_yscale('log')
        ax.legend(loc='upper center')
    
        # plot the distribution of the longest block with continuous quality above cutoff
        ax = axs[2]
        for readname in datafiles:
            lgb = histograms[readname][0]
            n_reads = sum(lgb)
            mean_length = 1.0 * sum(center_bins * lgb) / n_reads
            var_length = 1.0 * sum(center_bins**2 * lgb) / n_reads - mean_length**2
            lb = '#reads: '+str(n_reads)
            ax.plot(center_bins, lgb, lw=1.5,
                    label=readname+': '+str(100 * lgb[150:].sum() / lgb.sum())+'% > 150 b.p.')
        ax.set_title(lb)
        ax.set_xlabel('length of longest good block')
        ax.set_xlim(-1, 255)
        ax.set_ylim(1, plt.ylim()[1] * 1.1)
        ax.set_yscale('log')
        ax.legend(loc='upper center')
    
        # make a figure on the cumulative phred score distribution 
        # at different positions in the read
        for ri, readname in enumerate(datafiles):
            cum_phred_score_dis = histograms[readname][2].cumsum(axis=0)
            ax = axs[ri]
            plt.sca(ax)
            for pos in positions[::5]:
                plt.plot(1.0-cum_phred_score_dis[:,pos]/cum_phred_score_dis[-1,pos], c=cm.jet(pos))
    
            sm = plt.cm.ScalarMappable(cmap=cm.jet, norm=plt.normalize(0, 255))
            sm._A = []
            cbar = plt.colorbar(sm)
            ax.set_xlabel('phred score q')
            ax.set_title(readname)
    
        axs[0].set_ylabel('fraction better than q')

        fig.suptitle('adaID '+adaID, fontsize=18)
        plt.tight_layout(rect=(0, 0, 0.99, 0.95))

        plt.ion()
        plt.show()

        #fig.savefig(outdir+'/read_quality.pdf')
       


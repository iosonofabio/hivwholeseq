import os
import sys
from collections import defaultdict
from itertools import izip
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy import stats

# Globals
VERBOSE = 1
alpha = np.array(list('ACGT-N'), 'S1')
read_types = ['read1_f', 'read1_r', 'read2_f', 'read2_r']
fwd = np.array(np.arange(4)%2==0)
rev = np.array(np.arange(4)%2==1)
ref_filename = 'consensus_filtered_trimmed.fasta'
allele_count_filename = 'allele_counts.npy'

col=['g', 'b', 'r', 'm', 'c', 'y']
# Script
if __name__ == '__main__':
    if len(sys.argv)==2:
        sample_folder = 'adapterID_'+'{:02d}'.format(int(sys.argv[1]))
    else:
        sample_folder = 'unclassified_reads'


    # Get data
    out_filename = 'major_minor_sequences.fa'
    data_folder = ('/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'+
                   sample_folder+'/')
    ref_file = data_folder+ref_filename


    # Get counts
    # The first dimension is read type, the second alphabet, the third position
    counts = np.load(data_folder+allele_count_filename)
    coverage = np.sum(counts, axis=1)
    good_region_5p = np.where(coverage[1,:]>3000)[0]
    good_region_3p = np.where(coverage[0,:]>3000)[0]
    counts = counts[:,:,good_region_3p[0]+1000:good_region_5p[-1]]
    coverage = coverage[:,good_region_3p[0]+1000:good_region_5p[-1]]

    L=counts.shape[2]
    freqs = np.zeros_like(counts,float)
    for ri, rt in enumerate(read_types):
        freqs[ri,:,:] = 1.0*counts[ri,:,:]\
            /np.repeat(coverage[ri:(ri+1),:],alpha.shape[0],axis=0)
    minor_freqs = 1.0-np.max(freqs,axis=1)

    ##############################################################################
    # make coverage, minor variant, and cumulative frequency plot
    ##############################################################################
    plt.figure()
    for ri, rt in enumerate(read_types):
        plt.subplot(311)
        plt.plot(coverage[ri,:], label=rt)
        plt.subplot(312)
        plt.plot(minor_freqs[ri,:], label=rt)
        plt.subplot(313)
        plt.plot(sorted(minor_freqs[ri,:]), np.linspace(1,0,minor_freqs.shape[1]), label=rt)

    plt.subplot(311)
    plt.ylabel('coverage')
    plt.subplot(312)
    plt.yscale('log')
    plt.ylabel('minor allele freq')
    plt.subplot(313)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('minor allele freq')
    plt.xlabel('cumulative distribution')

    ##############################################################################
    # Scatter frequencies of every read type against every other
    ##############################################################################
    plt.figure(figsize=(15,10))
    panel_count=0
    region_to_analyse = np.min(coverage, axis=0)>1000
    for ri1, rt1 in enumerate(read_types):
        for ri2,rt2 in enumerate(read_types[:ri1]):
            panel_count+=1
            plt.subplot(2,3,panel_count)
            plt.yscale('log')
            plt.xscale('log')
            plt.ylim([1e-5,2])
            plt.xlim([1e-5,2])
            for ni,nuc in enumerate(alpha):
                plt.scatter(freqs[ri1,ni,region_to_analyse]+1e-10, freqs[ri2,ni,region_to_analyse]+1e-10, c=col[ni], label=nuc)
            plt.ylabel(rt2)
            plt.xlabel(rt1)
    plt.legend(loc=2)

    ##############################################################################
    # correct frequencies with reads different frequencies in forward and reverse direction
    ##############################################################################
    corrected_freqs = np.zeros(freqs.shape[1:])
    pval_distribution=[]
    for ni,nuc in enumerate(alpha[:5]):
        non_nuc_counts = [np.sum(coverage[fwd,:] - counts[fwd,ni,:], axis=0), 
                          np.sum(coverage[rev,:] - counts[rev,ni,:], axis=0)]
        nuc_counts =  [ np.sum(counts[fwd,ni,:], axis=0),  np.sum(counts[rev,ni,:], axis=0)]
        fwd_rev_coverage =[np.sum(coverage[fwd,:], axis=0), np.sum(coverage[rev,:], axis=0)]
        
        for pos in xrange(L):
            if (fwd_rev_coverage[0][pos]<100 or fwd_rev_coverage[1][pos]<100):
                corrected_freqs[ni,pos]=np.nan
            else:
                cm = np.array([[non_nuc_counts[0][pos], nuc_counts[0][pos]],
                               [non_nuc_counts[1][pos], nuc_counts[1][pos]]])

                chi2, pval = stats.chi2_contingency(cm+1)[:2]
                tmp_freqs = np.array([1.0*nuc_counts[0][pos]/fwd_rev_coverage[0][pos], 1.0*nuc_counts[1][pos]/fwd_rev_coverage[1][pos]])
                pval_distribution.append(pval)
                if pval<1e-6 and np.abs(tmp_freqs[1]-tmp_freqs[0])>0.001:
                    tmp_choice = np.argmax((2*tmp_freqs-1)**2)
                    corrected_freqs[ni,pos] = tmp_freqs[tmp_choice]
                    print pos, "corrections", corrected_freqs[ni,pos], tmp_freqs
                else:
                    corrected_freqs[ni,pos] = np.mean(tmp_freqs)

    corrected_minor_freq = 1.0-np.max(corrected_freqs, axis=0)
    plt.figure()
    bins = np.logspace(-4,np.log10(0.5),30)
    bc = 0.5*(bins[1:]+bins[:-1])
    for ri, rt in enumerate(read_types):    
        y,x = np.histogram(minor_freqs[ri,:], bins)
        plt.plot(sorted(minor_freqs[ri,region_to_analyse]), 
                 np.linspace(1,0,np.sum(region_to_analyse)), label=rt)

    plt.plot(sorted(corrected_minor_freq), np.linspace(1,0,corrected_minor_freq.shape[0]), label='corrected')
    plt.legend()
    plt.yscale('log')
    plt.xlim([0,1e-2])
    #plt.xscale('log')
    plt.ylabel('cumulative distribution')
    plt.xlabel('frequency')

    plt.figure()
    bins = np.logspace(-4,np.log10(1.0),50)
    bc = 0.5*(bins[1:]+bins[:-1])
    bw = bins[1:]-bins[:-1]
    y,x = np.histogram(corrected_minor_freq, bins)
    plt.plot(bc, y/bw)
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('error frequency distribution')
    plt.xlabel('frequency')

    plt.figure()
    plt.plot(sorted(pval_distribution), np.linspace(0,1,len(pval_distribution)))
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('cumulative pvalue distribution')
    plt.xlabel('p-value')



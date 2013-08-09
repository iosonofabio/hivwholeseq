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

sample_names = {2:'NL4-3', 4:'SF162', 7:'F10', 16:'patient', 18:'SF162_NL4-3', 19:'SF162_NL4-3_F10'}
col=['g', 'b', 'r', 'm', 'c', 'y']
# Script
if __name__ == '__main__':

    for adaID in [2,4,7,16,18,19]:
        # Get data
        out_filename = 'major_minor_sequences.fa'
        data_folder = ('/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'+
                       'adapterID_'+'{:02d}'.format(adaID)+'/')
        ref_file = data_folder+ref_filename

        # Read reference
        # Read reference
        if os.path.isfile(data_folder+ref_filename): ref_file = data_folder+ref_filename
        else: ref_file = '/'.join(data_folder.split('/')[:-2]+['subsample/']+data_folder.split('/')[-2:])+ref_filename
        refseq = SeqIO.read(ref_file, 'fasta')[500:9500]
        ref = np.array(refseq)

        # Get counts
        # The first dimension is read type, the second alphabet, the third position
        counts = np.load(data_folder+allele_count_filename)[:, :, 1000:8500]
        L=counts.shape[2]
        freqs = np.zeros_like(counts,float)
        coverage = np.sum(counts, axis=1)
        minor_freqs=np.zeros_like(coverage, float)
        corrected_freqs = np.zeros(freqs.shape[1:])

        for ri, rt in enumerate(read_types):
            freqs[ri,:,:] = 1.0*counts[ri,:,:]\
                /np.repeat(coverage[ri:(ri+1),:],alpha.shape[0],axis=0)

        minor_freqs = 1.0-np.max(freqs,axis=1)

        pval_distribution=[]
        if (adaID>10):
            plt.figure(figsize=(12,6))
        for ni,nuc in enumerate(alpha[:5]):
            non_nuc_counts = [np.sum(coverage[fwd,:] - counts[fwd,ni,:], axis=0), 
                              np.sum(coverage[rev,:] - counts[rev,ni,:], axis=0)]
            nuc_counts =  [ np.sum(counts[fwd,ni,:], axis=0),  np.sum(counts[rev,ni,:], axis=0)]
            fwd_rev_coverage =[np.sum(coverage[fwd,:], axis=0), np.sum(coverage[rev,:], axis=0)]

            overlap = (fwd_rev_coverage[0]>50000)*(fwd_rev_coverage[1]>50000)
            within_fragment = (fwd_rev_coverage[0]<50000)*(fwd_rev_coverage[1]<50000)
            if adaID>10 and ni<4:
                plt.subplot(121)
                plt.scatter(1.0*nuc_counts[0][overlap]/fwd_rev_coverage[0][overlap], 1.0*nuc_counts[1][overlap]/fwd_rev_coverage[1][overlap], c=col[ni])
                plt.subplot(122)
                plt.scatter(1.0*nuc_counts[0][within_fragment]/fwd_rev_coverage[0][within_fragment], 1.0*nuc_counts[1][within_fragment]/fwd_rev_coverage[1][within_fragment], c=col[ni], label=nuc)

            for pos in xrange(L):
                if (fwd_rev_coverage[0][pos]<100 or fwd_rev_coverage[1][pos]<100):
                    corrected_freqs[ni,pos]=np.nan
                else:
                    cm = np.array([[non_nuc_counts[0][pos], nuc_counts[0][pos]],
                                   [non_nuc_counts[1][pos], nuc_counts[1][pos]]])
                #print cm
                    chi2, pval = stats.chi2_contingency(cm+1)[:2]
                    tmp_freqs = np.array([1.0*nuc_counts[0][pos]/fwd_rev_coverage[0][pos], 1.0*nuc_counts[1][pos]/fwd_rev_coverage[1][pos]])
                    pval_distribution.append(pval)
                    if pval<1e-6 and np.abs(tmp_freqs[1]-tmp_freqs[0])>0.001:
                        tmp_choice = np.argmax((2*tmp_freqs-1)**2)
                        corrected_freqs[ni,pos] = tmp_freqs[tmp_choice]
                        print pos, "corrections", corrected_freqs[ni,pos], tmp_freqs
                    else:
                        corrected_freqs[ni,pos] = np.mean(tmp_freqs)
        if (adaID>10):
            plt.suptitle(sample_names[adaID])
            plt.subplot(121)
            plt.title('on different fragments')
            plt.yscale('log')
            plt.xscale('log')
            plt.xlim([1e-4,2])
            plt.ylim([1e-4,2])
            plt.xlabel('forward frequency')
            plt.ylabel('reverse frequency')
            plt.subplot(122)
            plt.title('on the same fragments')
            plt.yscale('log')
            plt.xscale('log')
            plt.xlim([1e-4,2])
            plt.ylim([1e-4,2])
            plt.xlabel('forward frequency')
            plt.ylabel('reverse frequency')
            plt.legend(loc=2)
            plt.savefig('fwd_rev_scatter_'+sample_names[adaID]+'.png')
        corrected_minor_freq = 1.0-np.max(corrected_freqs, axis=0)
        #plt.figure()
        #plt.title('adapterID_'+'{:02d}'.format(adaID))
        #bins = np.logspace(-4,np.log10(0.5),30)
        #bc = 0.5*(bins[1:]+bins[:-1])
        #for ri, rt in enumerate(read_types):    
        #    y,x = np.histogram(minor_freqs[ri,:], bins)
        #    plt.plot(sorted(minor_freqs[ri,:]), 
        #             np.linspace(1,0,minor_freqs.shape[1]), label=rt)

        #y,x = np.histogram(corrected_minor_freq, bins)
        #plt.plot(sorted(corrected_minor_freq), np.linspace(1,0,corrected_minor_freq.shape[0]), label='corrected')
        #plt.yscale('log')
        #plt.xscale('log')

        #plt.figure()
        #plt.title('adapterID_'+'{:02d}'.format(adaID))
        #plt.plot(coverage.T)

        #plt.figure()
        #plt.title('adapterID_'+'{:02d}'.format(adaID))
        #bins = np.logspace(-4,np.log10(1.0),50)
        #bc = 0.5*(bins[1:]+bins[:-1])
        #bw = bins[1:]-bins[:-1]
        #y,x = np.histogram(corrected_minor_freq, bins)
        #plt.plot(bc, y/bw)
        #plt.yscale('log')
        #plt.xscale('log')

        #plt.figure()
        #plt.title('adapterID_'+'{:02d}'.format(adaID))
        #plt.plot(sorted(pval_distribution), np.linspace(0,1,len(pval_distribution)))
        #plt.xscale('log')
        #plt.yscale('log')


        # pull out minor variants
        consensus = np.zeros_like(L, dtype='|S1')
        consensus = alpha[np.argmax(corrected_freqs, axis=0)]

        minor_variant = np.array(consensus)
        rare_variant = np.array(consensus)
        for ai,a in enumerate(alpha):
            ind = (corrected_freqs[ai,:]>=0.02)*(corrected_freqs[ai,:]<0.5)
            print np.sum(ind), np.where(ind)
            minor_variant[ind]=a
            ind = (corrected_freqs[ai,:]>=0.001)*(corrected_freqs[ai,:]<0.02)
            print np.sum(ind), np.where(ind)
            rare_variant[ind]=a

        with open(out_filename, 'a') as outfile:
            outfile.write('>'+sample_names[adaID]+'_adapterID_'+'{:02d}'.format(adaID)+'_major'+'\n')
            outfile.write("".join(consensus[:])+'\n')
            if adaID>10:
                outfile.write('>'+sample_names[adaID]+'_adapterID_'+'{:02d}'.format(adaID)+'_minor'+'\n')
                outfile.write("".join(minor_variant[:])+'\n')
                outfile.write('>'+sample_names[adaID]+'_adapterID_'+'{:02d}'.format(adaID)+'_rare'+'\n')
                outfile.write("".join(rare_variant[:])+'\n')


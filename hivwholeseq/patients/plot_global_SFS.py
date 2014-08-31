# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/05/14
content:    Plot site frequency spectra for derived alleles.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO

from hivwholeseq.miseq import alpha
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file, \
        get_allele_counts_insertions_from_file_unfiltered, \
        filter_nus
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories as plot_nus
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories_3d as plot_nus_3d
from hivwholeseq.patients.one_site_statistics import get_allele_frequency_trajectories


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele frequency trajectories')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele frequency trajectories to file')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--PCR1', action='store_true',
                        help='Show only PCR1 samples where possible (still computes all)')
    parser.add_argument('--BSC', action='store_true',
                        help='compare to BSC site frequency spectra')

    args = parser.parse_args()
    pnames = map(str, [9669, 15107, 15241, 15313,15376,15823,20097,20529])
    fragments = args.fragments
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot
    submit = args.submit
    use_PCR1 = args.PCR1
    add_bsc = args.BSC

    if add_bsc:
        import glob, pickle
        '''LOAD SFS OF DIRECT BSC SIMULATIONS'''                                        
        filelist = glob.glob('/ebio/ag-neher/share/users/rneher/WeakSelectionCoalescent/BS_SFS/SFS*.pickle')
        sfs_bsc = {}
        sfs_bc = {}
        sfs_bw={}
        sfs_bsc_filecount ={}
        for fname in filelist:                                                          
            entries=fname.split('_')                                                    
            N=int(entries[-2])                                                          
            file=open(fname, 'r')                                                       
            temp, temp_bc=pickle.load(file)                                           
            if (N in sfs_bsc):                                                            
                sfs_bsc[N]+=temp/sfs_bw[N]
                sfs_bsc_filecount[N]+=1                                                     
            else:                                                                       
                sfs_bc[N]=temp_bc                                                 
                sfs_bw[N] = sfs_bc[N][1:]-sfs_bc[N][:-1]
                sfs_bsc_filecount[N]=1                                                      
                sfs_bsc[N]=temp/sfs_bw[N]

        #Produce centered bins                                                          
        for N in sfs_bc:
            sfs_bc[N]=np.sqrt(sfs_bc[N][1:]*sfs_bc[N][:-1])
        #divide the spectra by the file count to normalize                              
        for N in sfs_bsc:                                                                 
            sfs_bsc[N]/=sfs_bsc_filecount[N]                                                  


        #calculate pi for the spectra and normalize to equal pi (== equal T_2)
        for N in sfs_bsc:
            print np.sum(sfs_bc[N]*(1-sfs_bc[N])*sfs_bw[N]*sfs_bsc[N])
            # normalize bsc SFS to unit pairwise difference
            sfs_bsc[N]/=np.sum(sfs_bc[N]*(1-sfs_bc[N])*sfs_bw[N]*sfs_bsc[N])

    #bins = np.logspace(-3, 0, 11)
    tbins = np.linspace(-7,7,21) 
    bins = np.exp(tbins)/(1+np.exp(tbins))
    binsc = np.sqrt(bins[1:] * bins[:-1])
    binw = np.diff(bins)
    hist = np.zeros_like(binsc)

    delta_af_bins = np.concatenate([-np.logspace(-2,0,11)[::-1], np.logspace(-2,0,11)])
    delta_af_binsc = 0.5*(delta_af_bins[1:] + delta_af_bins[:-1])
    delta_af_binw = np.diff(delta_af_bins)
    delta_af_hist = np.zeros_like(delta_af_binsc)

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    for pname in pnames:
        patient = get_patient(pname)
        times = patient.times()
        samplenames = patient.samples
        if use_PCR1:
            # Keep PCR2 only if PCR1 is absent
            ind = np.nonzero(map(lambda x: ('PCR1' in x[1]) or ((times == times[x[0]]).sum() == 1), enumerate(samplenames)))[0]
            times = times[ind]


        # Iterate over samples and fragments
        for fragment in fragments:

            if VERBOSE >= 1:
                print fragment

            act_filename = get_allele_count_trajectories_filename(pname, fragment)
            aft_filename = get_allele_frequency_trajectories_filename(pname, fragment)

            if os.path.isfile(aft_filename):
                aft = np.load(aft_filename)
                act = np.load(act_filename)

                aft[np.isnan(aft)] = 0
                aft[(aft < 1e-4) | (aft > 1)] = 0

                # Get rid of gaps and low-coverage regions
                is_gap = ((aft.argmax(axis=1) == 6) | (act.sum(axis=1) < 100)).any(axis=0)
                if VERBOSE >= 2:
                    print 'Fraction of gap sites (excluded):', is_gap.mean()

                if use_PCR1:
                    aft = aft[ind]
                    act = act[ind]

                # Get rid of ancestral alleles
                aft_der = aft.copy()
                aft_der[:, :, is_gap] = 0
                for i, ai in enumerate(aft[0].argmax(axis=0)):
                    aft_der[:, ai, i] = 0
                # take out everything at high frequency in first sample to improve polarization
                for i, ai in enumerate(aft[0].argmax(axis=0)):
                    aft_der[:, aft_der[0,:,i] > 0.1, i] = 0
                #remove first time sample
                aft_der[0]=0
                hist += np.histogram(aft_der, bins=bins, density=False)[0]/binw
                delta_af = np.diff(aft_der, axis=0)
                delta_af_hist += np.histogram(delta_af[np.where( (aft_der[:-1,:,:] > 0.1)*(aft_der[:-1,:,:] < 0.2))], 
                                              bins=delta_af_bins, density=False)[0]/delta_af_binw

    if plot:
        N=10000
        import matplotlib.pyplot as plt
        trfun = lambda x: np.log10(x / (1 - x))

        fig, ax = plt.subplots()
        ax.plot(trfun(binsc), hist, lw=2, c='k',marker='o', label = 'HIV polymorphisms')
        alpha = hist[0]
        if add_bsc:
            ax.plot(trfun(sfs_bc[N]), alpha*sfs_bsc[N]/sfs_bsc[N][np.argmax(sfs_bc[N]>1.3e-3)], 
                    lw=2, c='r', ls = '-', label = 'Strong selection')
        else:
            ax.plot(trfun(binsc[:8]), alpha*binsc[0]**2/binsc[:8]**2, lw=2, c='r')
        ax.plot(trfun(binsc), alpha*binsc[0]/binsc, lw=2, c='g', label = 'Neutral')

        ax.set_xlabel('derived allele frequency')
        ax.set_ylabel('SFS [density = counts / sum / binsize]')
        ax.set_xlim(-3.1, 3.1)
        ax.set_ylim([3e2,3e7])
        tickloc = np.array([0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999])
        ax.set_xticks(trfun(tickloc))
        ax.set_xticklabels(map(str, tickloc))
        from matplotlib.ticker import FixedLocator
        ticklocminor = np.concatenate([[10**po * x for x in xrange(2 , 10)] for po in xrange(-4, -1)] + \
                                      [[0.1 * x for x in xrange(2 , 9)]] + \
                                      [[1 - 10**po * (10 - x) for x in xrange(2, 10)] for po in xrange(-2, -5, -1)])
        ax.xaxis.set_minor_locator(FixedLocator(trfun(ticklocminor)))
        #ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title('All patients, fragments '+str(fragments))
        ax.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.ion()
        plt.show()



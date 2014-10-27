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
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories as plot_nus
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories_3d as plot_nus_3d
from hivwholeseq.patients.one_site_statistics import get_allele_frequency_trajectories



# Functions
def load_BSC_SFS(VERBOSE=0):
    '''Load sfs of direct bsc simulations'''
    import glob
    import cPickle as pickle

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
        if VERBOSE >= 3:
            print np.sum(sfs_bc[N]*(1-sfs_bc[N])*sfs_bw[N]*sfs_bsc[N])
        # normalize bsc SFS to unit pairwise difference
        sfs_bsc[N]/=np.sum(sfs_bc[N]*(1-sfs_bc[N])*sfs_bw[N]*sfs_bsc[N])

    return (sfs_bc, sfs_bsc)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get site frequency spectra',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--BSC', action='store_true',
                        help='compare to BSC site frequency spectra')
    parser.add_argument('--min-depth', type=int, default=500, dest='min_depth',
                        help='Minimal depth to consider the site')
    parser.add_argument('--saveplot', action='store_true',
                        help='Save the plot')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    use_plot = args.plot
    add_bsc = args.BSC
    depth_min = args.min_depth
    use_saveplot = args.saveplot

    # Prepare histogram data structures
    # Bin either logarithmically or logit
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

    if add_bsc:
        (sfs_bc, sfs_bsc) = load_BSC_SFS()

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        for fragment in fragments:
            if VERBOSE >= 1:
                print patient.name, fragment
    
            # TODO: include low depth due to low viral load
            aft, ind = patient.get_allele_frequency_trajectories(fragment,
                                                                 cov_min=depth_min)

            if VERBOSE >= 2:
                print 'Filter out masked positions'
            ind_nonmasked = -aft.mask.any(axis=0).any(axis=0)
            aft = aft[:, :, ind_nonmasked].data

            if VERBOSE >= 2:
                print 'Remove first time sample'
            aft_der = aft[1:].copy()

            if VERBOSE >= 2:
                print 'Filter out ancestral alleles'
            for i, ai in enumerate(aft[0].argmax(axis=0)):
                aft_der[:, ai, i] = 0
                # take out everything at high frequency in first sample to
                # improve polarization
                aft_der[:, aft[0, :, i] > 0.1, i] = 0

            # Add to the histogram
            hist += np.histogram(aft_der, bins=bins, density=False)[0]/binw
            delta_af = np.diff(aft_der, axis=0)
            h = np.histogram(\
                delta_af[np.where((aft_der[:-1,:,:] > 0.1)*(aft_der[:-1,:,:] < 0.2))], 
                bins=delta_af_bins, density=False)
            delta_af_hist += h[0] / delta_af_binw

    
    if use_plot:
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
        ax.set_ylim([3e1,3e7])
        tickloc = np.array([0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999])
        ax.set_xticks(trfun(tickloc))
        ax.set_xticklabels(map(str, tickloc))
        from matplotlib.ticker import FixedLocator
        ticklocminor = np.concatenate([[10**po * x for x in xrange(2 , 10)]
                                       for po in xrange(-4, -1)] + \
                                      [[0.1 * x for x in xrange(2 , 9)]] + \
                                      [[1 - 10**po * (10 - x) for x in xrange(2, 10)]
                                       for po in xrange(-2, -5, -1)])
        ax.xaxis.set_minor_locator(FixedLocator(trfun(ticklocminor)))
        #ax.set_xscale('log')
        ax.set_yscale('log')
        if len(fragments) == 6:
            ax.set_title('All patients, all fragments')
        else:
            ax.set_title('All patients, fragments '+str(fragments))
        ax.grid(True)
        plt.legend()

        plt.tight_layout()


        if use_saveplot:
            from hivwholeseq.patients.filenames import get_SFS_figure_filename
            if len(fragments) == 6:
                fn_out = get_SFS_figure_filename('all', 'all-fragments')
            else:
                fn_out = get_SFS_figure_filename('all', fragments)

        plt.ion()
        plt.show()



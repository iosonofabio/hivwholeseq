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
from hivwholeseq.utils.one_site_statistics import get_allele_counts_insertions_from_file, \
        get_allele_counts_insertions_from_file_unfiltered, \
        filter_nus
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories as plot_nus
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories_3d as plot_nus_3d
from hivwholeseq.patients.one_site_statistics import get_allele_frequency_trajectories



# Functions
def load_beta_SFS(VERBOSE=0, alpha=1, bins=20):
    '''Load sfs of direct bsc simulations'''
    sample_size = 3000

    import os
    from hivwholeseq.theory.filenames import get_sfs_betatree_filename
    fn = get_sfs_betatree_filename(sample_size, alpha)
    if os.path.isfile(fn):
        if VERBOSE >= 2:
            print 'Recycling beta coalescent SFS with N = '+str(sample_size)+\
                    ' and alpha = '+str(alpha)
        data = np.loadtxt(fn, unpack=True)
        sfs_bc = {(sample_size, alpha): data[0]}
        sfs_bsc = {(sample_size, alpha): data[1]}

    else:
        from hivwholeseq.theory.betatree.src.sfs import SFS
        if VERBOSE >= 2:
            print 'Generating beta coalescent SFS with N = '+str(sample_size)+\
                    ' and alpha = '+str(alpha)
        sfs_beta = SFS(sample_size=sample_size, alpha=alpha)
        sfs_beta.getSFS(ntrees=1000)
        sfs_beta.binSFS(mode='logit', bins=bins)
        sfs_beta.bin_center = np.sqrt(bins[1:] * bins[:-1])

        sfs_bc = {(sample_size, alpha): sfs_beta.bin_center}
        sfs_bsc = {(sample_size, alpha): sfs_beta.binned_sfs}

    return (sfs_bc, sfs_bsc)



# Script
if __name__ == '__main__':

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
    parser.add_argument('--save', action='store_true',
                        help='Save the SFS to file')
    parser.add_argument('--saveplot', action='store_true',
                        help='Save the plot')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    use_plot = args.plot
    add_bsc = args.BSC
    depth_min = args.min_depth
    use_save = args.save
    use_saveplot = args.saveplot

    # Prepare histogram data structures
    # Bin either logarithmically or logit
    #bins = np.logspace(-3, 0, 11)
    tbins = np.linspace(-7,7,21) 
    bins = np.exp(tbins)/(1+np.exp(tbins))
    binsc = np.sqrt(bins[1:] * bins[:-1])
    binw = np.diff(bins)
    hist = np.zeros_like(binsc)

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

            if VERBOSE >= 2:
                print 'Get initial allele frequencies'
            af0 = patient.get_initial_allele_frequencies(fragment, cov_min=depth_min)

            if VERBOSE >= 2:
                print 'Get allele frequencies'
            aft, ind = patient.get_allele_frequency_trajectories(fragment,
                                                                 depth_min=depth_min)

            if VERBOSE >= 2:
                print 'Filter out masked positions'
            ind_nonmasked = -aft.mask.any(axis=0).any(axis=0)
            af0 = af0[:, ind_nonmasked]
            aft = aft[:, :, ind_nonmasked].data

            if VERBOSE >= 2:
                print 'Remove first time sample (if still there)'
            aft_der = aft[int(0 in ind):].copy()

            if VERBOSE >= 2:
                print 'Filter out ancestral alleles'
            for i, ai in enumerate(af0.argmax(axis=0)):
                aft_der[:, ai, i] = 0
                # take out everything at high frequency in first sample to
                # improve polarization
                aft_der[:, af0[:, i] > 0.1, i] = 0

            hist += np.histogram(aft_der, bins=bins, density=False)[0] / binw
    
    # Add neutral spectrum
    sfs_neu = hist[0] * binsc[0]/binsc

    # Add selection spectra
    if add_bsc:
        (sfs_bc, sfs_bsc) = load_beta_SFS(bins=bins, VERBOSE=VERBOSE, alpha=1)
        tmp = load_beta_SFS(bins=bins, VERBOSE=VERBOSE, alpha=1.5)
        sfs_bc.update(tmp[0])
        sfs_bsc.update(tmp[1])
        for (N, alpha) in sfs_bsc:
            ind = (sfs_bsc[(N, alpha)] > 0).nonzero()[0]
            sfs_bc[(N, alpha)] = sfs_bc[(N, alpha)][ind]
            sfs_bsc[(N, alpha)] = sfs_bsc[(N, alpha)][ind] / sfs_bsc[(N, alpha)][ind[0]] * hist[ind[0]]
    else:
        indmax = (binsc < 0.1).nonzero()[0][-1] + 1
        sfs_sel = hist[0] * binsc[0]**2 / binsc[:indmax]**2

    if use_save:
        from hivwholeseq.patients.filenames import get_SFS_filename
        if pnames is None:
            fn_out = get_SFS_filename(['all'], fragments)
        else:
            fn_out = get_SFS_filename(pnames, fragments)
        # NOTE: do NOT make the call below recursive by default
        if not os.path.isdir(os.path.dirname(fn_out)):
            os.mkdir(os.path.dirname(fn_out))

        d_out = {'HIV_bin_centers': binsc, 'HIV_sfs': hist,
                 'neutral_bin_centers': binsc, 'neutral_sfs': sfs_neu}

        if add_bsc:
            for (N, alpha) in sfs_bsc:
                if alpha == 1:
                    key = 'bsc_'+str(N)
                else:
                    key= 'betatree_alpha_'+str(alpha)+'_'+str(N)
                d_out[key+'_bin_centers'] = sfs_bc[(N, alpha)]
                d_out[key+'_sfs'] = sfs_bsc[(N, alpha)]
        else:
            d_out['sel_bin_centers'] = binsc
            d_out['sel_sfs'] = sfs_sel

        np.savez(fn_out, **d_out)

    if use_plot:
        from hivwholeseq.utils.import plot
        from matplotlib import cm
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot(binsc, hist, lw=2, c='k',marker='o', label = 'HIV, depth >= '+str(depth_min))
        if add_bsc:
            for (N, alpha) in sfs_bc:
                ax.plot(sfs_bc[(N, alpha)], sfs_bsc[(N, alpha)],
                        lw=2, ls = '-',
                        color=cm.jet_r(1.0 * (alpha - 1)),
                        label = 'Beta coalescent, $\\alpha = '+str(alpha)+'$')
        else:
            ax.plot(binsc[:len(sfs_sel)], sfs_sel, lw=2, c='r')
        ax.plot(binsc, sfs_neu, lw=2, c='b', label = 'Neutral, $\\alpha = 2$')

        ax.set_xlabel('derived allele frequency')
        ax.set_ylabel('SFS [density = counts / sum / binsize]')
        ax.set_xlim(10**(-3.1), 1 - 10**(-3.1))
        ax.set_ylim([3e1,3e7])
        ax.set_xscale('logit')
        ax.set_yscale('log')
        if len(fragments) == 6:
            ax.set_title('All patients, all fragments')
        else:
            ax.set_title('All patients, fragments '+str(fragments))
        ax.grid(True)
        ax.set_ylim(1e2, 1e8)
        ax.legend(loc=(0.53, 0.68), fontsize=14)

        plt.tight_layout()

        if use_saveplot:
            from hivwholeseq.patients.filenames import get_SFS_figure_filename
            if len(fragments) == 6:
                fn_out = get_SFS_figure_filename('all', 'all-fragments')
            else:
                fn_out = get_SFS_figure_filename('all', fragments)

        plt.ion()
        plt.show()



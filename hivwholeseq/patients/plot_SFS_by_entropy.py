# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/05/14
content:    Plot site frequency spectra for derived alleles.
'''
# Modules
import os, sys
import argparse
import numpy as np
from Bio import SeqIO

from hivwholeseq.reference import load_custom_reference, load_custom_alignment
from hivwholeseq.miseq import alpha
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file, \
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


def get_coordinate_map(refname, ali):
    '''Get coordinate map of sequence in alignment'''
    for seq in ali:
        if refname in seq.name:
            break
    else:
        raise ValueError('Reference not found in alignment')

    seqm = np.array(seq, 'S1')
    ind = (seqm != '-').nonzero()[0]
    mapali = np.vstack([np.arange(len(ind)), ind])
    return mapali



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

    # Select type M entropy strata
    S_bins = np.array([0, 0.01, 0.05, 0.1, 0.5, 1, 2])

    # Prepare histogram data structures
    # Bin either logarithmically or logit
    #bins = np.logspace(-3, 0, 11)
    tbins = np.linspace(-7,7,21) 
    bins = np.exp(tbins)/(1+np.exp(tbins))
    binsc = np.sqrt(bins[1:] * bins[:-1])
    binw = np.diff(bins)
    hists = np.zeros((len(S_bins) - 1, len(binsc)))

    if VERBOSE >= 1:
        print 'Load alignment, reference, and coordinate map'
    ali = load_custom_alignment('HIV1_FLT_2013_genome_DNA')
    alim = np.array(ali, 'S1')
    S = np.zeros(alim.shape[1])
    for a in alpha[:5]:
        nu = (alim == a).mean(axis=0)
        S -= nu * np.log2(nu + 1e-8)

    refname = 'HXB2'
    refseq = load_custom_reference('HXB2', format='gb')
    mapali = get_coordinate_map(refname, ali)
    if len(refseq) != mapali[0, -1] + 1:
        raise ValueError('Reference '+refname+' in alignment is not complete')
    Sref = S[mapali[1]]

    Srefind = - np.ones(len(Sref), int)
    for i, Sbin in enumerate(S_bins[:-2]):
        ind = (Sref >= Sbin) & (Sref < S_bins[i + 1])
        Srefind[ind] = i
    Srefind[Srefind < 0] = len(S_bins) - 2
    
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    if VERBOSE >= 1:
        print 'Analyze patients'
    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        for fragment in fragments:
            if VERBOSE >= 1:
                print patient.name, fragment

            mapco = patient.get_map_coordinates_reference(fragment, refname=refname)
    
            # TODO: include low depth due to low viral load
            aft, ind = patient.get_allele_frequency_trajectories(fragment,
                                                                 cov_min=depth_min)

            if VERBOSE >= 2:
                print 'Filter out masked positions'
            ind_nonmasked = -aft.mask.any(axis=0).any(axis=0)

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

            # Add to the histograms
            for ih in xrange(len(S_bins) - 1):
                ind = mapco[Srefind[mapco[:, 0]] == ih, 1]
                # Discard masked sites (low coverage and similia)
                ind = ind[ind_nonmasked[ind]]
                hists[ih] += np.histogram(aft_der[:, :, ind], bins=bins, density=False)[0]/binw

    if add_bsc:
        (sfs_bc, sfs_bsc) = load_beta_SFS(bins=bins, VERBOSE=VERBOSE, alpha=1)
        tmp = load_beta_SFS(bins=bins, VERBOSE=VERBOSE, alpha=1.5)
        sfs_bc.update(tmp[0])
        sfs_bsc.update(tmp[1])
        for (N, alpha) in sfs_bsc:
            ind = (sfs_bsc[(N, alpha)] > 0).nonzero()[0]
            sfs_bc[(N, alpha)] = sfs_bc[(N, alpha)][ind]
            sfs_bsc[(N, alpha)] = sfs_bsc[(N, alpha)][ind] / sfs_bsc[(N, alpha)][ind[0]] * hists[1, ind[0]]
    
    if use_plot:
        from hivwholeseq import plot_utils
        from matplotlib import cm
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()

        # Plot HIV
        for ih in xrange(len(S_bins) - 1):
            ax.plot(binsc, hists[ih] / hists[ih, 0] * hists[-1, 0], lw=2, c='k', marker='o',
                    label='HIV, $S_{type M} \in ['+str(S_bins[ih])+', '+str(S_bins[ih + 1])+']$',
                    color=cm.jet(1.0 * ih / hists.shape[0]))

        # Plot theory
        al = hists[1, 0]
        if add_bsc:
            for (N, alpha) in sfs_bc:
                ax.plot(sfs_bc[(N, alpha)], sfs_bsc[(N, alpha)],
                        lw=2, ls = '-',
                        color=cm.jet_r(1.0 * (alpha - 1)),
                        label = 'Beta coalescent, $\\alpha = '+str(alpha)+'$')
        else:
            ax.plot(binsc[:8], al*binsc[0]**2/binsc[:8]**2, lw=2, c='r')
        ax.plot(binsc, al*binsc[0]/binsc, lw=2, c='b', label = 'Neutral, $\\alpha = 2$')

        ax.set_xlabel('derived allele frequency')
        ax.set_ylabel('SFS [density = counts / sum / binsize]')
        ax.set_xlim(10**(-3.1), 1.0 - 10**(-3.1))
        ax.set_ylim([3e1,3e7])
        ax.set_xscale('logit')
        ax.set_yscale('log')
        if len(fragments) == 6:
            ax.set_title('All patients, all fragments')
        else:
            ax.set_title('All patients, fragments '+str(fragments))
        ax.grid(True)
        ax.set_ylim(1e2, 1e8)
        ax.legend(loc=1, fontsize=10)

        plt.tight_layout()


        # Plot HIV divided by one another
        den = hists[-1]
        fig, ax = plt.subplots()
        for ih in xrange(len(S_bins) - 1):
            y = hists[ih] / den
            y /= y[0]
            ax.plot(binsc, y, lw=2, c='k', marker='o',
                    label='HIV, $S_{type M} \in ['+str(S_bins[ih])+', '+str(S_bins[ih + 1])+']$',
                    color=cm.jet(1.0 * ih / hists.shape[0]))

        ax.set_xlabel('derived allele frequency')
        ax.set_ylabel('$SFS / SFS_{S_{max}}$')
        ax.set_xlim(10**(-3.1), 1.0 - 10**(-3.1))
        ax.set_ylim([1e-3, 1.2])
        tickloc = np.array([0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999])
        ax.set_xscale('logit')
        ax.set_yscale('log')
        if len(fragments) == 6:
            ax.set_title('All patients, all fragments')
        else:
            ax.set_title('All patients, fragments '+str(fragments))
        ax.grid(True)
        ax.set_ylim(1e-3, 1e2)
        ax.legend(loc=1, fontsize=10)

        plt.tight_layout()

        # Plot correlation between entropy and "fixation" probability (high-nu tail)
        x = np.sqrt(S_bins[1:-1] * S_bins[2:])
        y = [(h[h > 0][-3:] / hists[-1, h > 0][-3:]).mean() for h in hists[1:]]
        xerr = np.vstack([x - S_bins[1:-1], S_bins[2:] - x])
        r = np.corrcoef(x, y)[0, 1]
        fig, ax = plt.subplots()
        if len(fragments) == 6:
            ax.set_title('All patients, all fragments')
        else:
            ax.set_title('All patients, fragments '+str(fragments))
        ax.set_ylabel('$<SFS / SFS_{S_{max}}>_{last bins}$')
        ax.set_xlabel('$<S>$ [bits]')
        ax.errorbar(x, y, xerr=xerr, fmt='o', ms=5, lw=2, label='HIV, r = '+str(np.round(r, 2)))

        Tc = 1.0
        xx = np.linspace(0, x[-1], 1000)
        yy = np.exp((xx - x[-1]) * Tc)
        ax.plot(xx, yy, c='k', lw=2, label='$\exp[ (S - S_{max}) \cdot T_c ]$, $T_c = '+str(Tc)+'$')

        ax.legend(loc=4)
        ax.grid(True)
        ax.set_xlim(-0.05, 1.1 * (x[-1] + xerr[1, -1]))
        ax.set_ylim(-0.05, 1.05)

        plt.tight_layout()


        if use_saveplot:
            from hivwholeseq.patients.filenames import get_SFS_figure_filename
            if len(fragments) == 6:
                fn_out = get_SFS_figure_filename('all', 'all-fragments')
            else:
                fn_out = get_SFS_figure_filename('all', fragments)

        plt.ion()
        plt.show()



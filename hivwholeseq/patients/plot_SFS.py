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
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories as plot_nus
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories_3d as plot_nus_3d
from hivwholeseq.patients.one_site_statistics import get_allele_frequency_trajectories



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele frequency trajectories')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
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

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot
    submit = args.submit
    use_PCR1 = args.PCR1

    patient = get_patient(pname)
    times = patient.times()
    samplenames = patient.samples
    if use_PCR1:
        # Keep PCR2 only if PCR1 is absent
        ind = np.nonzero(map(lambda x: ('PCR1' in x[1]) or ((times == times[x[0]]).sum() == 1), enumerate(samplenames)))[0]
        times = times[ind]

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    # Iterate over samples and fragments
    for fragment in fragments:

        if VERBOSE >= 1:
            print fragment

        act_filename = get_allele_count_trajectories_filename(pname, fragment)
        aft_filename = get_allele_frequency_trajectories_filename(pname, fragment)

        aft = np.load(aft_filename)
        act = np.load(act_filename)

        aft[np.isnan(aft)] = 0
        aft[(aft < 1e-5) | (aft > 1)] = 0

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

        bins = np.logspace(-2, -0.5, 11)
        binsc = np.sqrt(bins[1:] * bins[:-1])

        hist = np.histogram(aft_der, bins=bins, density=True)[0]

        if plot:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots()
            ax.plot(binsc, hist, lw=2, c='k')
            alpha = hist[-3]
            ax.plot([binsc[0], binsc[-1]], [alpha * ((binsc[-3] / binsc[0])**2),
                                            alpha * ((binsc[-3] / binsc[-1])**2)], lw=2, c='r')
            ax.plot([binsc[0], binsc[-1]], [alpha * ((binsc[-3] / binsc[0])),
                                            alpha * ((binsc[-3] / binsc[-1]))], lw=2, c='g')

            ax.set_xlabel('Freq')
            ax.set_ylabel('SFS [density = counts / sum / binsize]')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_title(pname+', '+fragment)
            ax.grid(True)

            plt.tight_layout()
            plt.ion()
            plt.show()



# vim: fdm=marker
'''
author:     Fabio Zanini
date:       22/05/14
content:    Plot coverage across time, for a fragment or genomewide.
'''
# Modules
import argparse
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.patients.filenames import get_allele_count_trajectories_filename, \
        get_coverage_to_initial_figure_filename, get_figure_folder



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
    parser.add_argument('--PCR1', action='store_true',
                        help='Show only PCR1 samples where possible (still computes all)')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele frequency trajectories to file')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    save_to_file = args.save
    VERBOSE = args.verbose
    use_PCR1 = args.PCR1

    patient = get_patient(pname)
    times = patient.times()
    samplenames = np.array(patient.samples)
    if use_PCR1:
        # Keep PCR2 only if PCR1 is absent
        ind = np.nonzero(map(lambda x: ('PCR1' in x[1]) or ((times == times[x[0]]).sum() == 1), enumerate(samplenames)))[0]
        times = times[ind]
        samplenames = samplenames[ind]

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    # Iterate over samples and fragments
    for fragment in fragments:
        if VERBOSE >= 1:
            print pname, fragment

        act_filename = get_allele_count_trajectories_filename(pname, fragment)
        act = np.load(act_filename)
        cov = act.sum(axis=1)

        if use_PCR1:
            cov = cov[ind]

        # Plot
        fig, ax = plt.subplots()
        for i in xrange(cov.shape[0]):
            PCRtype = samplenames[i].split('_')[-1]
            color = cm.jet(1.0 * i / cov.shape[0])
            ax.plot(cov[i] + 0.1, lw=2, c=color, label=str(times[i])+' ('+PCRtype+')')
        ax.set_xlabel('Position [bp]')
        ax.set_ylabel('Coverage')
        ax.set_title(' '.join([pname, fragment]))
        ax.set_yscale('log')
        ax.set_ylim(0.1, 1e6)
        ax.set_xlim(-30, cov.shape[-1] + 30)

        ax.legend(loc='lower center', title='Time [days]:', fontsize=12)

        if save_to_file:
            out_fn = get_coverage_to_initial_figure_filename(pname, fragment)
            mkdirs(get_figure_folder(pname))
            fig.savefig(out_fn)

    plt.ion()
    plt.show()


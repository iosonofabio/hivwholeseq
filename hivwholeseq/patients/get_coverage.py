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

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.patients.samples import itersample
from hivwholeseq.patients.filenames import get_allele_count_trajectories_filename, \
        get_coverage_to_initial_figure_filename, get_figure_folder



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get coverage')
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--regions', nargs='+', default=['genomewide'],
                        help='Regions to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the allele frequency trajectories')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    regions = args.regions
    VERBOSE = args.verbose
    use_plot = args.plot

    samples = lssp()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()

    if use_plot:
        fig, ax = plt.subplots()

    for ir, region in enumerate(regions):

        for isam, (samplename, sample) in enumerate(itersample(samples)):
            if VERBOSE >= 1:
                print samplename, region

            cov = sample.get_coverage(region)
        
            if use_plot:
                color = cm.jet(1.0 * (isam + ir * len(samples)) / (len(regions) * len(samples)))
                ax.plot(cov + 0.1, lw=2, c=color, label=samplename+', '+region)

        if use_plot:
            ax.set_xlabel('Position [bp]')
            ax.set_ylabel('Coverage')
            ax.set_yscale('log')
            ax.set_ylim(0.1, 1e6)

    if use_plot:
        ax.legend(loc='lower center', fontsize=12)
        ax.grid(True)
        ax.set_xlim(-30, 2100)

        plt.ion()
        plt.show()


#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/06/14
content:    Calculate the allele counts and frequencies for each patient sample
            (and both PCR1 and PCR2 if present).
'''
# Modules
import argparse
import numpy as np
import matplotlib.pyplot as plt


from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.patients.samples import SamplePat



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele counts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')
    parser.add_argument('--qualmin', type=int, default=30,
                        help='Minimal quality of base to call')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    regions = args.regions
    VERBOSE = args.verbose
    qual_min = args.qualmin
    use_plot = args.plot

    samples = lssp()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()

    for region in regions:
        for samplename, sample in samples.iterrows():
            sample = SamplePat(sample)

            if VERBOSE >= 1:
                print region, samplename

            count = sample.get_allele_counts(region, qual_min=qual_min)

            if use_plot:
                x = np.tile(np.arange(count.shape[1]), (count.shape[0], 1))
                color = np.tile(np.arange(count.shape[0]), (count.shape[1], 1)).T

                fig, ax = plt.subplots(figsize=(12, 6))
                
                ax.scatter(x, count + 0.1, lw=2, c=color)
                ax.set_xlabel('Position [bp]')
                ax.set_ylabel('Coverage')
                ax.set_xlim(-1, count.shape[-1])
                ax.set_ylim(ymin=0.09)
                ax.set_yscale('log')
                ax.grid(True)
                ax.set_title(samplename+', '+region)

    if use_plot:
        plt.ion()
        plt.show()

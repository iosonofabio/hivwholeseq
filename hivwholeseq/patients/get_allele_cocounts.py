#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/03/14
content:    Get the joint counts at two sites for patient samples, after mapping.
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
    parser = argparse.ArgumentParser(description='Get allele cocounts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map')
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
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

            cocount = np.load(fn_out)['cocounts']

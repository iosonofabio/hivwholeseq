#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/12/14
content:    Compress cocount matrices for faster IO.
'''
# Modules
import os
import sys
import argparse
import numpy as np

from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.cluster.fork_cluster import fork_compress_cocounts_patient as fork_self



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele cocounts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--qualmin', type=int, default=30,
                        help='Minimal quality of base to call')
    parser.add_argument('--PCR', type=int, default=1,
                        help='Analyze only reads from this PCR (1 or 2)')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    VERBOSE = args.verbose
    submit = args.submit
    qual_min = args.qualmin
    PCR = args.PCR

    samples = lssp()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    if submit:
        for fragment in fragments:
            for samplename, sample in samples.iterrows():
                fork_self(samplename, fragment, VERBOSE=VERBOSE,
                          qual_min=qual_min, PCR=PCR)
        sys.exit()


    for samplename, sample in samples.iterrows():
        sample = SamplePat(sample)
        pname = sample.patient

        for fragment in fragments:

            if VERBOSE >= 1:
                print pname, samplename, fragment

            fn = sample.get_allele_cocounts_filename(fragment, PCR=PCR,
                                                     qual_min=qual_min,
                                                     compressed=False)
            
            fn_out = sample.get_allele_cocounts_filename(fragment, PCR=PCR,
                                                         qual_min=qual_min,
                                                         compressed=True)

            if not os.path.isfile(fn):
                if VERBOSE >= 2:
                    print 'Input file not found, skipping'
                continue

            if VERBOSE >= 2:
                print 'Loading cocounts'
            cocount = np.load(fn)

            if VERBOSE >= 2:
                print 'Storing compressed cocounts'
            np.savez_compressed(fn_out, cocounts=cocount)

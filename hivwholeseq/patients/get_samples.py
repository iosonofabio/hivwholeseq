# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/02/14
content:    Get patient samples.
'''
# Modules
import argparse

from hivwholeseq.patients.samples import load_samples_sequenced as lssp



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Load patient samples')
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    VERBOSE = args.verbose

    samples = lssp()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()

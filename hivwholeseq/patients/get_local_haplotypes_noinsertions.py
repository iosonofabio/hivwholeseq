# vim: fdm=marker
'''
author:     Fabio Zanini
date:       04/09/14
content:    Get local haplotypes from the reads fully covering a small region
            of interest, cluster them somewhat and pack them into trajectories.

            NOTE: this script ignores insertions, which require a multiple
            sequence alignment, for the sake of simplicity. Large regions are
            going to fail, because read pairs cover them with a hole in the middle.
'''
# Modules
import os
import argparse
from operator import itemgetter
import numpy as np
from Bio import SeqIO



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get coverage trajectories',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--roi', required=True, nargs='+',
                        help='Region of interest (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot
    use_PCR1 = args.PCR1

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()
    samplenames = patient.samples.index

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    for fragment in fragments:
        if VERBOSE >= 1:
            print fragment


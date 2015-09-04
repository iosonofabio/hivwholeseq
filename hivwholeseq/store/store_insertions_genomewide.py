# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/03/14
content:    Merge allele frequency trajectories of all fragments.
'''
# Modules
import os
import argparse
import numpy as np
from collections import Counter
from Bio import SeqIO

from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.utils.miseq import alpha, read_types
from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.patients.filenames import get_initial_reference_filename
from hivwholeseq.utils.sequence import find_annotation
from hivwholeseq.store.store_insertions import save_insertions



# Functions
def merge_insertions(ics, VERBOSE=0):
    '''Merge the insertions of all fragments'''
    from collections import Counter
    ic = [Counter() for rt in read_types]

    for (fragment, start), icts in ics.iteritems():
        if VERBOSE >= 2:
            print fragment, start
        for irt, ict in enumerate(icts):
            for (position, insertion), value in ict.iteritems():
                ic[irt][(position + start, insertion)] += value

    return ic



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Store genomewide insertions',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', action=PatientsAction,
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele count trajectories to file')
    parser.add_argument('--PCR', type=int, default=1,
                        help='Analyze only reads from this PCR (e.g. 1)')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    VERBOSE = args.verbose
    save_to_file = args.save
    PCR = args.PCR

    samples = lssp()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()

    for samplename, sample in samples.iterrows():
        if VERBOSE >= 1:
            print samplename

        sample = SamplePat(sample)
        pname = sample.patient
        ref = sample.get_reference('genomewide', 'gb')

        # Collect the insertions (where possible)
        ics = {}
        for fragment in ['F'+str(i) for i in xrange(1, 7)]:
            try:
                ic = sample.get_insertions(fragment, merge_read_types=False)
            except IOError:
                continue
            start = find_annotation(ref, fragment).location.nofuzzy_start
            ics[(fragment, start)] = ic

        if not len(ics):
            if VERBOSE >= 1:
                print 'No data found: skipping'
            continue

        # Merge insertions
        ic = merge_insertions(ics, VERBOSE=VERBOSE)
        if save_to_file:
            fn_out = sample.get_insertions_filename('genomewide')
            save_insertions(fn_out, ic)
            if VERBOSE >= 1:
                print 'Genomewide insertions saved to:', fn_out

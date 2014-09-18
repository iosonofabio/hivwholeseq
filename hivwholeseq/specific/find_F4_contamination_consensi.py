# vim: fdm=marker
'''
author:     Fabio Zanini
date:       16/09/14
content:    Three patient samples have F4 consensus sequences that are far away
            from anybody else and similar to each other. They are very similar
            to an RNA culture standard, 38304, i.e. Tue59 N6-S3.
'''
# Modules
import os
import sys
import argparse
import numpy as np
from seqanpy import align_global, align_local

from hivwholeseq.reference import load_custom_reference
from hivwholeseq.sequencing.samples import load_samples_sequenced, SampleSeq
from hivwholeseq.sequence_utils import pretty_print_pairwise_ali
from hivwholeseq.generic_utils import getchar


# Scripts
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Check contamination from RNA sample',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--verbose', type=int, default=1,
                        help='Verbosity level [0-3]')
    
    args = parser.parse_args()
    VERBOSE = args.verbose
    pnames = None

    contseq = load_custom_reference('38304_F4')

    # Take mapped and filtered seqs and look for residual contamination
    from hivwholeseq.patients.patients import load_samples_sequenced as lssp
    from hivwholeseq.patients.patients import SamplePat
    if pnames is not None:
        samples = lssp(patients=pnames)
    else:
        samples = lssp()

    for samplename, sample in samples.iterrows():
        sample = SamplePat(sample)

        if VERBOSE >= 1:
            print samplename,
    
        try:
            consseq = sample.get_consensus('F4')
        except IOError:
            print ''
            continue

        ali = align_global(contseq, consseq, score_gapext=0)
        score = ali[0]
        if VERBOSE >= 1:
            print score

        # Less than 10 mismatches from the contamination is suspicious
        if score > 3 * len(ali[1]) - (6 * 10):
            is_contam = True
            print 'SUSPECT CONTAMINATION!'
        else:
            is_contam = False

        if (is_contam and (VERBOSE >= 2)) or (VERBOSE >= 3):
            pretty_print_pairwise_ali(ali[1:], '38304', samplename, width=90)

        if is_contam:
            import ipdb; ipdb.set_trace()

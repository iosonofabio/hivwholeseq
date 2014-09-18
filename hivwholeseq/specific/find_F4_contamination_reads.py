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



# Functions
def fish_contamination(bamfilename, contseq, VERBOSE=0, deltascore=60,
                       output_alignments=False,
                       **kwargs):
    '''Fish contaminated reads from mapped reads
    
    Args:
      **kwargs: passed down to the pairwise alignment function
    '''
    import pysam
    from seqanpy import align_overlap

    if 'score_match' in kwargs:
        score_match = kwargs['score_match']
    else:
        score_match = 3

    conts = ''.join(contseq)
    alignments = []

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        n_good = 0
        n_cont = 0

        # We could advance pair by pair, or read by read (should not change much)
        for ir, read in enumerate(bamfile):
            if VERBOSE >= 2:
                if not ((ir + 1) % 100):
                    if not ((ir + 1) == 100):
                        sys.stdout.write('\x1b[1A')
                    print ir + 1

            if ir >= 10000:
                break

            (score, ali1, ali2) = align_overlap(conts, read.seq, **kwargs)
            
            start = len(ali2) - len(ali2.lstrip('-'))
            end = len(ali2.rstrip('-'))
            ali1 = ali1[start: end]
            ali2 = ali2[start: end]

            scoremax = len(ali1) * score_match
            if score > scoremax - deltascore:
                n_cont += 1
                if output_alignments:
                    alignments.append([score, scoremax, ali1, ali2])

                if VERBOSE >= 2:
                    print 'Good:', n_good, 'cont:', n_cont

                if VERBOSE >= 3:
                    print scoremax, score
                    pretty_print_pairwise_ali([ali1, ali2], width=90, name1='ref', name2='read')

            else:
                n_good += 1
                if VERBOSE >= 4:
                    print scoremax, score
                    pretty_print_pairwise_ali([ali1, ali2], width=90, name1='ref', name2='read')

    if output_alignments:
        return (n_good, n_cont, alignments)
    else:
        return (n_good, n_cont)



# Scripts
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Check contamination from RNA sample 38304',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--fragment', default='F4',
                        help='Fragment to check')
    parser.add_argument('--verbose', type=int, default=1,
                        help='Verbosity level [0-3]')
    parser.add_argument('--allsamples', action='store_true',
                        help='Analyze all samples, not only the 6 suspicious ones')
    
    args = parser.parse_args()
    VERBOSE = args.verbose
    pnames = None
    fragment = args.fragment
    allsamples = args.allsamples

    contseq = load_custom_reference('38304_'+fragment)

    # Take mapped and filtered seqs and look for residual contamination
    from hivwholeseq.patients.patients import load_samples_sequenced as lssp
    from hivwholeseq.patients.patients import SamplePat
    if pnames is not None:
        samples = lssp(patients=pnames)
    else:
        samples = lssp()

    # Take only the samples that Lina is suspecting
    if not allsamples:
        samplenames = ['VK04-4187', '14908', '18113', '12879', '6154', '18798']
        samples = samples.loc[samples.index.isin(samplenames)]

    for samplename, sample in samples.iterrows():
        sample = SamplePat(sample)

        if VERBOSE >= 1:
            print samplename,
            if VERBOSE >= 2:
                print ''
    
        bamfilename = sample.get_mapped_filtered_filename(fragment)
        if not os.path.isfile(bamfilename):
            print ''
            continue

        (n_good, n_cont) = fish_contamination(bamfilename, contseq, VERBOSE=VERBOSE,
                                              deltascore=30)


        if VERBOSE >= 1:
            if VERBOSE >= 2:
                print samplename,
            print n_good, n_cont

        #ali = align_global(contseq, consseq, score_gapext=0)
        #score = ali[0]
        #if VERBOSE >= 1:
        #    print score

        ## Less than 10 mismatches from the contamination is suspicious
        #if score > 3 * len(ali[1]) - (6 * 10):
        #    is_contam = True
        #    print 'SUSPECT CONTAMINATION!'
        #else:
        #    is_contam = False

        #if (is_contam and (VERBOSE >= 2)) or (VERBOSE >= 3):
        #    pretty_print_pairwise_ali(ali[1:], '38304', samplename, width=90)

        #if is_contam:
        #    import ipdb; ipdb.set_trace()

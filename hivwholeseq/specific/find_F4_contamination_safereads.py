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
from hivwholeseq.utils.sequence import pretty_print_pairwise_ali
from hivwholeseq.utils.generic import getchar
from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.patients.patients import SamplePat
from hivwholeseq._secret import (
    samplenames_suspected_contamination_consensus as samples_cont)



# Functions
def fish_original(bamfilename, bamfilename_out,
                  contseq, consensi, VERBOSE=0, deltascore=60,
                  maxreads=10000,
                  output_alignments=False,
                  **kwargs):
    '''Fish non contaminated reads from mapped reads of broken consensi
    
    Args:
      **kwargs: passed down to the pairwise alignment function
    '''
    import pysam
    from seqanpy import align_overlap

    from hivwholeseq.utils.mapping import pair_generator

    if 'score_match' in kwargs:
        score_match = kwargs['score_match']
    else:
        score_match = 3

    conts = ''.join(contseq)
    conss = {sn: ''.join(consseq) for sn, consseq in consensi.iteritems()}
    alignments = []

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        with pysam.Samfile(bamfilename_out, 'wb', template=bamfile) as bamfile_out:
            n_good = 0
            n_cont = 0

            # We could advance pair by pair, or read by read (should not change much)
            for irp, reads in enumerate(pair_generator(bamfile)):
                if VERBOSE >= 2:
                    if not ((irp + 1) % 100):
                        if not ((irp + 1) == 100):
                            sys.stdout.write('\x1b[1A')
                        print irp + 1

                if irp == maxreads:
                    break

                skip = False
                for read in reads:
                    (score, alic1, alic2) = align_overlap(conts, read.seq, **kwargs)
                    start = len(alic2) - len(alic2.lstrip('-'))
                    end = len(alic2.rstrip('-'))
                    alic1 = alic1[start: end]
                    alic2 = alic2[start: end]
                    scoremax = len(alic1) * score_match
                    delta_read = scoremax - score
                    if delta_read <= deltascore:
                        skip = True
                        break

                    for consname, cons in conss.iteritems():
                        (score, ali1, ali2) = align_overlap(cons, read.seq, **kwargs)
                        start = len(ali2) - len(ali2.lstrip('-'))
                        end = len(ali2.rstrip('-'))
                        ali1 = ali1[start: end]
                        ali2 = ali2[start: end]
                        scoremax = len(ali1) * score_match
                        delta_read = scoremax - score

                        if delta_read <= deltascore:
                            if output_alignments:
                                alignments.append([score, scoremax, ali1, ali2])

                            if VERBOSE >= 2:
                                print 'Good:', n_good, 'cont:', n_cont

                            if VERBOSE >= 3:
                                print scoremax, score
                                pretty_print_pairwise_ali([alic1, alic2], width=90,
                                                          name1='38304', name2='read')
                                pretty_print_pairwise_ali([ali1, ali2], width=90,
                                                          name1=consname, name2='read')

                            break
                    else:
                        skip = True
                    if skip:
                        break

                if skip:
                    n_cont += 1
                else:
                    n_good += 1
                    bamfile_out.write(reads[0])
                    bamfile_out.write(reads[1])

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
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    
    args = parser.parse_args()
    VERBOSE = args.verbose
    fragment = args.fragment
    maxreads = args.maxreads

    contseq = load_custom_reference('38304_'+fragment)

    # Take only the samples with contaminated consensus
    samples = lssp()
    samples_focal = samples.loc[samples_cont]

    for samplename, sample in samples_focal.iterrows():
        sample = SamplePat(sample)
        samples_pat = samples.loc[samples.patient == sample.patient]
        consensi = {}
        for sn, s in samples_pat.iterrows():
            try:
                consensi[sn] = SamplePat(s).get_consensus(fragment)
            except IOError:
                continue
        del consensi[samplename]

        if VERBOSE >= 1:
            print samplename,
            if VERBOSE >= 2:
                print ''
    
        bamfilename = sample.get_mapped_filtered_filename(fragment)
        if not os.path.isfile(bamfilename):
            print ''
            continue

        bamfilename_out = sample.get_mapped_filtered_filename(fragment, decontaminated=True)

        (n_good, n_cont) = fish_original(bamfilename, bamfilename_out,
                                         contseq, consensi,
                                         VERBOSE=VERBOSE,
                                         maxreads=maxreads,
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

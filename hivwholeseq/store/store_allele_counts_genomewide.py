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
from Bio import SeqIO

from hivwholeseq.miseq import alpha, read_types
from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.patients.filenames import get_initial_reference_filename



# Functions
def merge_allele_counts(ref_genomewide, acs, VERBOSE=0):
    '''Merge the allele counts of all fragments
    
    Note: we do not require full coverage of all fragments, the missing
          ones will just have zero counts. Sometimes, cherry-picking the data
          fragment by fragment might be a better choice.
    '''
    from hivwholeseq.miseq import alpha, read_types
    from seqanpy import align_overlap

    ac = np.zeros((len(read_types), len(alpha), len(ref_genomewide)), int)

    pos_ref = 1000
    for (fr, ref, acsi) in acs:

        # Find the coordinates
        (score, ali1, ali2) = align_overlap(ref_genomewide[pos_ref - 1000:],
                                            ref,
                                            #score_gapopen=-20,
                                           )
        fr_start = len(ali2) - len(ali2.lstrip('-'))
        fr_end = len(ali2.rstrip('-'))

        if VERBOSE:
            print fr, pos_ref - 1000 + fr_start, pos_ref - 1000 + fr_end

        # Scan the alignment
        pos_ref = pos_ref - 1000 + fr_start
        fr_start_ref = pos_ref
        fr_end_ref = pos_ref + fr_end - fr_start
        pos_fr = 0
        for pos_ali in xrange(fr_start, fr_end):
            # Gap in genomewise, ignore position
            if ali1[pos_ali] == '-':
                pos_fr += 1
                continue

            # Gap in fragment, ignore FIXME: probably we should put deletions
            elif ali2[pos_ali] == '-':
                pos_ref += 1
                continue

            # Add the counts
            # NOTE: all fragments are treated the same, even in case of coverage
            # differences of orders of magnitude. This means, larger coverage
            # always wins. Maybe we want to implement this somewhat differently
            ac[:, :, pos_ref] += acsi[:, :, pos_fr]
            pos_fr += 1
            pos_ref += 1

        if VERBOSE >= 3:
            from hivwholeseq.utils.sequence import pretty_print_pairwise_ali
            cons = alpha[ac.sum(axis=0).argmax(axis=0)]
            pretty_print_pairwise_ali((ali1[fr_start: fr_end],
                                       cons[fr_start: fr_end]),
                                      name1='gw',
                                      name2=fr,
                                      width=100)

    return ac



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Store genomewide allele counts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze (e.g. VL98-1253 VK03-4298)')
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
        conss_genomewide = SeqIO.read(get_initial_reference_filename(pname, 'genomewide'), 'fasta')

        # Collect the allele counts (where possible)
        acs = []
        for fragment in ['F'+str(i) for i in xrange(1, 7)]:
            try:
                ref = ''.join(SeqIO.read(get_initial_reference_filename(pname, fragment), 'fasta'))
                ac = sample.get_allele_counts(fragment, merge_read_types=False)
                acs.append((fragment, ref, ac))
            except IOError:
                continue

        if not len(acs):
            if VERBOSE >= 1:
                print 'No data found: skipping'
            continue

        # Merge allele counts
        ac = merge_allele_counts(conss_genomewide, acs, VERBOSE=VERBOSE)
        if save_to_file:
            fn_out = sample.get_allele_counts_filename('genomewide')
            np.save(fn_out, ac)
            if VERBOSE >= 1:
                print 'Genomewide allele counts saved to:', fn_out

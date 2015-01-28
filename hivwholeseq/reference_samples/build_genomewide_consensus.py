# vim: fdm=marker
'''
author:     Fabio Zanini
date:       16/01/15
content:    Join fragments for a reference (consensus) sequence into a
            genomewide sequence.
'''
# Modules
import os
import argparse

from hivwholeseq.reference import (
    load_custom_reference, get_custom_reference_filename)



# Functions
def merge_sequences(seqs, skip_initial=30, accept_gaps=False, VERBOSE=0):
    '''Merge sequences with overlaps
    
    Parameters:
       seqs (list): sequences to merge
       skip_initial (int): trim from the beginning of overlaps because we do not
       really trust those bases
       accept_gaps (bool): accept gaps in the overlaps
    '''
    from itertools import izip
    from seqanpy import align_ladder
    import numpy as np

    seqs = map(''.join, seqs)

    left_trim = 0
    seqs_all = []
    for iov, (seq1, seq2) in enumerate(izip(seqs[:-1], seqs[1:])):
        if VERBOSE >= 1:
            print 'Overlap n', iov+1

        (score, ali1, ali2) = align_ladder(seq1[left_trim:], seq2, score_gapopen=-20)
        start2 = len(ali2) - len(ali2.lstrip('-'))
        end1 = len(ali1.rstrip('-'))

        # Append first sequence until overlap
        seqs_all.append(ali1[:start2 + skip_initial])

        # Check overlap
        ov1 = ali1[start2 + skip_initial: end1 - skip_initial]
        ov2 = ali2[start2 + skip_initial: end1 - skip_initial]

        if VERBOSE >= 2:
            from hivwholeseq.utils.sequence import pretty_print_pairwise_ali
            pretty_print_pairwise_ali((ov1, ov2), width=100,
                                      name1='seq1', name2='seq2')

        if (not accept_gaps) and (('-' in ov1) or ('-' in ov2)):
            raise ValueError('Gaps in the overlap n. '+str(iov+1))

        # Trust the first sequence until half, then the other one
        i_mid = len(ov1) // 2
        seqs_all.append(ov1[:i_mid])
        seqs_all.append(ov2[i_mid:])

        # Set the left trim for the trailing sequence
        left_trim = len(ali2[: end1 - skip_initial].replace('-', ''))

    if VERBOSE >= 1:
        print 'Add last sequence'
    seqs_all.append(seq2[left_trim:])

    return ''.join(seqs_all)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Merge fragments in reference')
    parser.add_argument('--reference', required=True,
                        help='Reference to analyze (e.g. LAI-III)')
    parser.add_argument('--save', action='store_true',
                        help='Save to file')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    refname = args.reference
    use_save = args.save
    VERBOSE = args.verbose


    fragments = ['F'+str(i) for i in xrange(1, 7)]
    consensi = [load_custom_reference(refname+'_'+fragment) for fragment in fragments]

    consensus = merge_sequences(consensi, VERBOSE=VERBOSE)
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    consrec = SeqRecord(Seq(consensus, alphabet=consensi[0].seq.alphabet),
                        id=refname+'_genomewide',
                        name=refname+'_genomewide',
                        description=refname+', genomewide reference (merged)',
                       )

    if use_save:
        if VERBOSE >= 1:
            print 'Save to file'
        from Bio import SeqIO
        fn = get_custom_reference_filename(refname, format='fasta')
        if os.path.isfile(fn):
            raise IOError('Destination file already exists')
        SeqIO.write(consrec, fn, 'fasta')


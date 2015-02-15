# vim: fdm=marker
'''
author:     Fabio Zanini
date:       15/01/15
content:    Build a reference alignment from LANL sequences.
'''
# Modules
import os
import argparse
from collections import defaultdict
from Bio import AlignIO

from hivwholeseq.cross_sectional.filenames import (
    get_raw_LANL_sequences_filename,
    get_subtype_reference_alignment_filename)
from hivwholeseq.utils.sequence import align_codon_pairwise


# Globals



# Functions
def align_to_reference(seq, refstr, VERBOSE=0, codon_align=False):
    '''Align sequence to refernce, stripping reference gaps'''
    import numpy as np
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from seqanpy import align_overlap
    from hivwholeseq.utils.sequence import pretty_print_pairwise_ali

    seqstr = ''.join(seq).upper()

    n_amb = len(seqstr) - sum(map(seqstr.count, ('A', 'C', 'G', 'T', '-')))
    if n_amb > 2:
        raise ValueError('Too many ambiguous sites')

    def align_dna(seqstr, refstr):
        (score, alis, alir) = align_overlap(seqstr, refstr)
        start = len(alir) - len(alir.lstrip('-'))
        end = len(alir.rstrip('-'))
        alist = alis[start: end]
        alirt = alir[start: end]
        return (alist, alirt)

    (alis, alir) = align_dna(seqstr, refstr)

    if codon_align:
        (alis, alir) = align_codon_pairwise(alis.replace('-', ''), alir.replace('-', ''))

    # If the sequence is shorter than HXB2, skip
    if '-' in (alis[0], alis[-1]):
        raise ValueError('The sequence does not fully cover the region')

    # If the sequence has too much gapping close to the edges, it's also short
    if (alis[:15].count('-') > 5) or (alis[-15:].count('-') > 5):
        raise ValueError('The sequence does not fully cover the region')

    if VERBOSE >= 2:
        pretty_print_pairwise_ali((alis, alir), width=100,
                                  name2=refname, name1=seq.name)


    # Strip gaps in HXB2
    alism = np.fromstring(alis, 'S1')
    alirm = np.fromstring(alir, 'S1')
    ind = (alirm != '-')
    seq_aliref = ''.join(alism[ind])

    rec = SeqRecord(Seq(seq_aliref, seq.seq.alphabet),
                    id=seq.id,
                    name=seq.name,
                    description=seq.description)

    return rec


def build_reference_alignments(region, refname,
                               VERBOSE=0,
                               subtypes=['B', 'C', 'A', 'AE', 'F1', 'D', 'O', 'H'],
                               codon_align=False,
                              ):
    '''Build reference alignment by subtype'''
    from hivwholeseq.reference import load_custom_reference
    from Bio import SeqIO
    from Bio.Align import MultipleSeqAlignment

    ref = load_custom_reference(refname, region=region)
    refstr = ''.join(ref)

    seq_by_subtype = defaultdict(list)

    fn_in = get_raw_LANL_sequences_filename(region)
    if VERBOSE >= 2:
        print fn_in

    seq_iter = SeqIO.parse(fn_in, 'fasta')

    for i, seq in enumerate(seq_iter):
        if VERBOSE >= 1:
            if not ((i+1) % 100):
                print i+1

        subtype = seq.id.split('.')[0]

        if subtype not in subtypes:
            continue

        if VERBOSE >= 3:
            print subtype

        try:
            rec = align_to_reference(seq, refstr, VERBOSE=VERBOSE,
                                     codon_align=codon_align)
        except ValueError:
            continue

        seq_by_subtype[subtype].append(rec)


    for subtype, seqs in seq_by_subtype.iteritems():
        seq_by_subtype[subtype] = MultipleSeqAlignment(seqs)

    return seq_by_subtype



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Align to reference',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--region', required=True,
                        help='Region to anign (e.g. V3)')
    parser.add_argument('--reference', default='HXB2',
                        help='Reference to use for alignment')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--subtypes', nargs='+', default=['B'],
                        help='Subtypes to keep')
    parser.add_argument('--codonalign', action='store_true',
                        help='Align codon by codon')

    args = parser.parse_args()
    region = args.region
    refname = args.reference
    VERBOSE = args.verbose
    subtypes = args.subtypes
    codalign = args.codonalign

    alis = build_reference_alignments(region, refname,
                                      subtypes=subtypes,
                                      codon_align=codalign,
                                      VERBOSE=VERBOSE)

    for subtype, ali in alis.iteritems():
        if VERBOSE >= 1:
            print subtype
        fn = get_subtype_reference_alignment_filename(region, subtype=subtype,
                                                      refname=refname,
                                                      VERBOSE=VERBOSE)
        AlignIO.write(ali, fn, 'fasta')


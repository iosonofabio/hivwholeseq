# vim: fdm=marker
'''
author:     Fabio Zanini
date:       15/01/15
content:    Clean a reference alignment from LANL sequences from problematic seqs.
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
def filter_sequence(seq, VERBOSE=0):
    '''Align sequence to refernce, stripping reference gaps'''
    seqstr = ''.join(seq).upper()
    n_amb = len(seqstr) - sum(map(seqstr.count, ('A', 'C', 'G', 'T', '-')))
    if n_amb > 2:
        return False

    return True


def clean_reference_alignment(region, refname,
                              VERBOSE=0,
                              subtype='B',
                             ):
    '''Clean reference alignment'''
    from hivwholeseq.reference import load_custom_reference
    from Bio import SeqIO
    from Bio.Align import MultipleSeqAlignment

    from hivwholeseq.cross_sectional.filenames import (
        get_subtype_reference_alignment_filename)


    fn = get_subtype_reference_alignment_filename(region, subtype=subtype,
                                                  refname=refname,
                                                  VERBOSE=VERBOSE)
    ali = AlignIO.read(fn, 'fasta')
    nseqs = len(ali)

    ali = MultipleSeqAlignment(filter(filter_sequence, ali))
    nseqsnew = len(ali)

    if VERBOSE >= 2:
        print refname, region, subtype+':', nseqsnew, 'of', nseqs, 'seqs kept'

    return ali



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
    parser.add_argument('--subtype', default='B',
                        help='Subtype to analyze')

    args = parser.parse_args()
    region = args.region
    refname = args.reference
    VERBOSE = args.verbose
    subtype = args.subtype

    ali = clean_reference_alignment(region, refname,
                                    subtype=subtype,
                                    VERBOSE=VERBOSE)

    fn = get_subtype_reference_alignment_filename(region, subtype=subtype,
                                                  refname=refname,
                                                  VERBOSE=VERBOSE)
    AlignIO.write(ali, fn, 'fasta')


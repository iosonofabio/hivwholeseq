# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/15
content:    Translate a reference alignment into amino acids.
'''
# Modules
import argparse
from Bio import AlignIO

from hivwholeseq.sequence_utils import translate_alignment
from hivwholeseq.cross_sectional.get_subtype_reference_alignment import (
    get_subtype_reference_alignment)
from hivwholeseq.cross_sectional.filenames import get_subtype_reference_alignment_filename



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

    args = parser.parse_args()
    region = args.region
    refname = args.reference
    VERBOSE = args.verbose
    subtypes = args.subtypes

    for subtype in subtypes:
        if VERBOSE >= 1:
            print subtype

        if VERBOSE >= 2:
            print 'Get nucleotide alignment'
        ali = get_subtype_reference_alignment(region, subtype=subtype,
                                              refname=refname,
                                              type='nuc',
                                              VERBOSE=VERBOSE)
 
        if VERBOSE >= 2:
            print 'Translating'
        aliaa = translate_alignment(ali, VERBOSE=VERBOSE)

        if VERBOSE >= 2:
            print 'Save to file'
        fn = get_subtype_reference_alignment_filename(region, subtype=subtype,
                                                      refname=refname,
                                                      type='aa',
                                                      VERBOSE=VERBOSE)
        AlignIO.write(aliaa, fn, 'fasta')


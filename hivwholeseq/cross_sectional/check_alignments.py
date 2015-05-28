# vim: fdm=indent
'''
author:     Fabio Zanini
date:       28/05/15
content:    Check that the pairwise alignments are fine.
'''
# Modules
import os
import argparse
import numpy as np
from collections import defaultdict
from Bio import AlignIO

from hivwholeseq.cross_sectional.filenames import (
    get_subtype_reference_alignment_filename)




# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Align to reference',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--regions', required=True, nargs='+',
                        help='Region to anign (e.g. V3)')
    parser.add_argument('--reference', default='HXB2',
                        help='Reference to use for alignment')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--subtypes', nargs='+', default=['B'],
                        help='Subtypes to keep')

    args = parser.parse_args()
    regions = args.regions
    refname = args.reference
    VERBOSE = args.verbose
    subtypes = args.subtypes


    from hivwholeseq.reference import load_custom_reference
    from hivwholeseq.utils.sequence import find_annotation
    ref = load_custom_reference('HXB2', 'gb')

    for region in regions:
        regm = np.array(find_annotation(ref, region).extract(ref), 'S1')
        for subtype in subtypes:
            fn = get_subtype_reference_alignment_filename(region,
                                                          subtype=subtype,
                                                          refname=refname,
                                                          VERBOSE=VERBOSE)
            alim = np.array(AlignIO.read(fn, 'fasta'), 'S1')
            weird = ((alim != regm).mean(axis=1) > 0.2)
            print region, subtype, weird.sum()
                


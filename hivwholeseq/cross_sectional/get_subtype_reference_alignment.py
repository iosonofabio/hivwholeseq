# vim: fdm=marker
'''
author:     Fabio Zanini
date:       12/01/15
content:    Get subtype entropy from alignments, since it's used so often.
'''
# Modules
import os
import argparse
import numpy as np

from hivwholeseq.miseq import alpha
from hivwholeseq.one_site_statistics import get_entropy
from hivwholeseq.cross_sectional.filenames import (
    get_subtype_reference_alignment_filename)



# Functions
def get_subtype_reference_alignment(region, subtype='B',
                                    refname='HXB2',
                                    type='nuc',
                                    VERBOSE=0):
    '''Get the observables from subtype B reference alignments'''
    from Bio import AlignIO
    ali_fn = get_subtype_reference_alignment_filename(region,
                                                      subtype=subtype,
                                                      refname=refname,
                                                      type=type,
                                                      VERBOSE=VERBOSE)
    ali = AlignIO.read(ali_fn, 'fasta')
    return ali



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Get subtype alignment',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--subtype', default='B',
                        help='Subtype of the alignment')
    parser.add_argument('--reference', default='HXB2',
                        help='Reference of the alignment')

    args = parser.parse_args()
    regions = args.regions
    VERBOSE = args.verbose
    subtype = args.subtype
    refname = args.reference


    alis = {}
    for region in regions:
        if VERBOSE >= 1:
            print region

        if VERBOSE >= 2:
            print 'Get alignment'
        ali = get_subtype_reference_alignment(region, subtype=subtype,
                                              refname=refname,
                                                  VERBOSE=VERBOSE)
        alis[region] = ali


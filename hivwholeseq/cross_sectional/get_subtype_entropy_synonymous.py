# vim: fdm=marker
'''
author:     Fabio Zanini
date:       12/01/15
content:    Get subtype entropy from alignments, since it's used so often.
'''
# Modules
import os
import argparse
import cPickle as pickle
import numpy as np

from hivwholeseq.miseq import alpha
from hivwholeseq.one_site_statistics import get_entropy
from hivwholeseq.cross_sectional.filenames import (
    get_subtype_reference_alignment_filename,
    get_subtype_reference_alignment_entropy_syn_filename)
from hivwholeseq.cross_sectional.get_subtype_reference_alignment import get_subtype_reference_alignment



# Functions
def get_ali_entropy_syn(alim, positions=None, alpha=alpha[:5], VERBOSE=0):
    '''Get entropy of alignment at some positions'''
    from collections import defaultdict
    from Bio.Seq import translate

    if len(ali[0]) % 3:
        raise ValueError('The alignment length is not a multiple of 3')

    if positions is None:
        positions = np.arange(len(ali[0]) // 3)

    # The data structure is a nested dict by position and amino acid
    S = {}
    # Iterate over codons
    for pos in positions:
        if VERBOSE >= 3:
            print pos

        asub = alim[:, pos * 3: (pos + 1) * 3]
        aacount = defaultdict(lambda: defaultdict(int))
        for cod in asub:
            cod = ''.join(cod)
            aacount[translate(cod)][cod] += 1

        Spos = {}
        for aa, codd in aacount.iteritems():
            af = np.array(codd.values(), float)
            af /= af.sum()

            Spos[aa] = get_entropy(af)
        S[pos] = Spos

    return S


def get_subtype_reference_alignment_entropy_syn(region, subtype='B', VERBOSE=0,
                                                refname='HXB2',
                                                type='nuc'):
    '''Get the entropy of a large subtype reference alignment'''
    import cPickle as pickle

    fn = get_subtype_reference_alignment_entropy_syn_filename(region,
                                                              subtype=subtype,
                                                              refname=refname,
                                                              type=type,
                                                              VERBOSE=VERBOSE)

    with open(fn, 'r') as f:
        data = pickle.load(f)

    return data



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Calculate entropy of subtype alignment',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--subtype', default='B',
                        help='Subtype of the alignment')
    parser.add_argument('--reference', default='HXB2',
                        help='Reference of the alignment')
    parser.add_argument('--save', action='store_true',
                        help='Save to file')

    args = parser.parse_args()
    regions = args.regions
    VERBOSE = args.verbose
    subtype = args.subtype
    refname = args.reference
    use_save = args.save


    Ss = {}
    for region in regions:
        if VERBOSE >= 1:
            print region

        if use_save:

            if VERBOSE >= 2:
                print 'Get alignment'
            ali = get_subtype_reference_alignment(region, subtype=subtype,
                                                  refname=refname,
                                                  VERBOSE=VERBOSE)
            alim = np.array(ali, 'S1')

            if VERBOSE >= 2:
                print 'Compute entropy'
            S = get_ali_entropy_syn(alim, VERBOSE=VERBOSE)

            if VERBOSE >= 2:
                print 'Store to file'
            fn_out = get_subtype_reference_alignment_entropy_syn_filename(region,
                                                                      subtype=subtype,
                                                                      refname=refname,
                                                                      VERBOSE=VERBOSE)

            with open(fn_out, 'wb') as f:
                pickle.dump(S, f)

        else:
            if VERBOSE >= 2:
                print 'Get entropy from file'
            S = get_subtype_reference_alignment_entropy_syn(region,
                                                        subtype=subtype,
                                                        refname=refname,
                                                        VERBOSE=VERBOSE)

        Ss[region] = S


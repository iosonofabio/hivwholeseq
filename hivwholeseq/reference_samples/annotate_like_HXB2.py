# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/08/15
content:    Annotate another sequence using features from HXB2.
'''
# Modules
import os
import argparse
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from hivwholeseq.reference import load_custom_reference, save_custom_reference
from hivwholeseq.utils.sequence import trim_to_refseq


# Functions
def annotate_like_HXB2(refname, VERBOSE=0):
    '''Annotate copying from HXB2'''
    hxb2 = load_custom_reference('HXB2', 'gb')
    ref = load_custom_reference(refname, 'fasta')
    refs = str(ref.seq)

    def get_sublocation(sublocation):
        hxb2_seq = sublocation.extract(hxb2)
        ref_seq = trim_to_refseq(refs, hxb2_seq).replace('-', '')
        start = refs.find(ref_seq)
        end = start + len(ref_seq)
        return FeatureLocation(start, end, strand=+1)

    for fea in hxb2.features:
        if VERBOSE >= 1:
            print fea.id
        loc = [get_sublocation(loc) for loc in fea.location.parts]
        if len(loc) == 1:
            loc = loc[0]
        else:
            loc = CompoundLocation(loc)

        feature = SeqFeature(loc, type=fea.type, id=fea.id)

        # Test length of old and new
        if fea.id not in ["LTR5'", "LTR3'", 'V4']:
            L1 = len(fea.extract(hxb2))
            L2 = len(feature.extract(ref))
            s = str(L2)+' vs '+str(L1)
            if 1.0 * L2 / L1 < 0.9:
                raise ValueError('Feature: '+fea.id+' is too short: '+s)
            elif 1.0 * L2 / L1 > 1.1:
                raise ValueError('Feature: '+fea.id+' is too long: '+s)

        ref.features.append(feature)

    return ref



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Annotate like HXB2',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--refname', required=True,
                        help='Name of reference to annotate')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save annotated sequence to file')

    args = parser.parse_args()
    refname = args.refname
    VERBOSE = args.verbose
    save = args.save

    ref = annotate_like_HXB2(refname, VERBOSE=VERBOSE)

    if save:
        save_custom_reference(ref, refname, format='gb')

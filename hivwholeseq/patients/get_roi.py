# vim: fdm=marker
'''
author:     Fabio Zanini
date:       03/12/14
content:    Get genomic region of interest from a number of possible inputs.
'''
# Modules
import argparse
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.argparse_utils import RoiAction



# Functions
def get_fragmented_roi(patient, roi, VERBOSE=0, include_genomewide=False):
    '''From a Region Of Interest, get fragment(s), start and end coordinates'''
    def find_fragment(fea_frags, start, end, include_genomewide=False):
        for fea in fea_frags:
            fr_start = fea.location.nofuzzy_start
            fr_end = fea.location.nofuzzy_end

            if (fr_start <= start) and (fr_end >= end):
                return (fea.id, start - fr_start, end - fr_start)
        
        if not include_genomewide:
            raise ValueError('No fragment found that fully covers this roi')
        else:
            return ('genomewide', start, end)

    refseq = patient.get_reference('genomewide', 'gb')
    fea_frags = [fea for fea in refseq.features if fea.type == 'fragment']
    feanames = [fea.id for fea in refseq.features]

    if roi[0] in ['F'+str(i) for i in xrange(1, 7)]:
        start = roi[1]
        if roi[2] != '+oo':
            end = roi[2]
        else:
            fea = refseq.features[feanames.index(roi[0])]
            end = fea.location.nofuzzy_end - fea.location.nofuzzy_start
        roi = (roi[0], start, end)

        if VERBOSE >= 3:
            print 'Fragment selected', roi
        return roi

    if roi[0] == 'genomewide':
        start = roi[1]
        if roi[2] != '+oo':
            end = roi[2]
        else:
            end = len(refseq)

        roi = find_fragment(fea_frags, start, end, include_genomewide=include_genomewide)
        if VERBOSE >= 3:
            print 'Genomewide selected', roi
        return roi

    elif roi[0] in feanames:
        fea = refseq.features[feanames.index(roi[0])]
        start = roi[1] + fea.location.nofuzzy_start
        if roi[2] != '+oo':
            end = roi[2] + fea.location.nofuzzy_start
        else:
            end = fea.location.nofuzzy_end

        roi = find_fragment(fea_frags, start, end, include_genomewide=include_genomewide)
        if VERBOSE >= 3:
            print 'Feature selected', roi
        return roi

    raise ValueError('Roi not understood')


def get_fragments_covered(patient, roi, VERBOSE=0):
    '''Get the list of fragments interested by this roi'''
    if roi[0] in ['F'+str(ifr) for ifr in xrange(1, 7)]:
        return [roi[0]]

    (fragment, start, end) = get_fragmented_roi(patient, roi, VERBOSE=VERBOSE,
                                                include_genomewide=True)
    if fragment != 'genomewide':
        return [fragment]

    refseq = patient.get_reference('genomewide', 'gb')
    frags_cov = []
    for fea in refseq.features:
        if fea.type == 'fragment':
            if (fea.location.nofuzzy_start < end) & (fea.location.nofuzzy_end > start):
                frags_cov.append(fea.id)

    return frags_cov



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get ROI',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--roi', required=True, action=RoiAction,
                        help='Region of interest (e.g. F1 300 350)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--include-genomewide', action='store_true',
                        help='Accept genomewide as a fragment')


    args = parser.parse_args()
    pname = args.patient
    roi = args.roi
    VERBOSE = args.verbose
    include_genomewide = args.include_genomewide

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()

    (fragment, start, end) = get_fragmented_roi(patient, roi, VERBOSE=VERBOSE,
                                                include_genomewide=include_genomewide)


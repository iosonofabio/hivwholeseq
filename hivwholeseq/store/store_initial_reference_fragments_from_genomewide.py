# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/01/15
content:    Rebuild initial reference for single fragments from the genomewide
            one, making sure we use F3b and F5a.
'''
# Modules
import os
import argparse
from Bio import SeqIO
from seqanpy import align_overlap

from hivwholeseq.patients.patients import load_patient
from hivwholeseq.data.primers import primers_outer



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Build initial references from genomewide',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='+',
                        default=['F'+str(i) for i in xrange(1, 7)],
                        help='Fragments to taylor')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    VERBOSE = args.verbose
    fragments = args.fragments

    primers_outer['F3'] = primers_outer['F3b']
    primers_outer['F5'] = primers_outer['F5a']

    patient = load_patient(pname)
    refseqgw = patient.get_reference('genomewide')

    for fragment in fragments:
        if VERBOSE >= 1:
            print pname, fragment

        if VERBOSE >= 2:
            print 'Cutting out fragment', fragment

        # Get start coordinate
        if fragment == 'F1':
            start = 0
        else:
            prfwd = primers_outer[fragment][0]
            (score, ali1, ali2) = align_overlap(refseqgw, prfwd, score_gapopen=-20)
            start = len(ali2) - len(ali2.lstrip('-')) + len(prfwd)

        # Get end coordinate
        if fragment == 'F6':
            end = len(refseqgw)
        else:
            prrev = primers_outer[fragment][1]
            (score, ali1, ali2) = align_overlap(refseqgw, prrev, score_gapopen=-20)
            end = len(ali2) - len(ali2.lstrip('-'))

        refseq = refseqgw[start: end]

        refseq.id = patient.code+'_ref_'+fragment
        refseq.name = refseq.id
        refseq.description = 'Patient '+patient.code+', initial reference '+fragment

        if VERBOSE >= 2:
            print 'Store to file'
        fn_out = patient.get_reference_filename(fragment, format='fasta')
        SeqIO.write(refseq, fn_out, 'fasta')
        print 'Reference changed, remap!'


# vim: fdm=marker
'''
author:     Fabio Zanini
date:       16/12/14
content:    Merge initial reference from the 6 fragments into one sequence.
'''
# Modules
import os
import argparse
from Bio import SeqIO

from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.build_genomewide_consensus import merge_fragments



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Build initial genomewide reference',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    VERBOSE = args.verbose

    patient = load_patient(pname)

    fragments = ['F'+str(i) for i in xrange(1, 7)]
    refseqs = {fr: patient.get_reference(fr) for fr in fragments}

    if VERBOSE >= 1:
        print 'Merge fragments'
    refseq = merge_fragments(refseqs, VERBOSE=VERBOSE)

    refseq.id = patient.name+'_ref'
    refseq.name = refseq.id
    refseq.description = 'Patient '+patient.name+', initial genomewide reference'

    if VERBOSE >= 1:
        print 'Store to file'
    fn_out = patient.get_reference_filename('genomewide', format='fasta')
    SeqIO.write(refseq, fn_out, 'fasta')


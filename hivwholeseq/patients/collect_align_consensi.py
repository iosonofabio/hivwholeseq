# vim: fdm=marker
'''
author:     Fabio Zanini
date:       25/10/13
content:    Collect and align all consensi from a patient.
'''
# Modules
import os
import argparse
from Bio import SeqIO
from Bio import AlignIO

from mapping.samples import samples as samples_seq
from mapping.datasets import MiSeq_runs
from mapping.filenames import get_merged_consensus_filename, get_consensus_filename
from mapping.patients.filenames import get_consensi_alignment_genomewide_filename, \
        get_consensi_alignment_filename
from mapping.patients.patients import get_patient, get_sequenced_samples
from mapping.mapping_utils import align_muscle



# Functions
def get_consensi_frag(patient, fragment, VERBOSE=0):
    '''Collect all consensi for a single fragment'''
    consensi = []
    for sample in patient['samples']:
        # Find run and adaID of the sample
        miseq_run = samples_seq[sample]['run']
        dataset = MiSeq_runs[miseq_run]
        data_folder = dataset['folder']
        adaID = dataset['adapters'][dataset['samples'].index(sample)]

        # Get input filename
        input_filename = get_consensus_filename(data_folder, adaID, fragment)
        # If present, add it
        if os.path.isfile(input_filename):
            seq = SeqIO.read(input_filename, 'fasta')

            # Change name
            seq.description = ', '.join([seq.id, seq.description])
            seq.id = seq.name = sample

            consensi.append(seq)

    return consensi


def get_consensi_genomewide(patient, VERBOSE=0):
    '''Collect the consensi'''
    # FIXME: assume all fragments are there for now (no patchy stuff yet)
    fragments = ['F'+str(i) for i in xrange(1, 7)]
    consensi = []
    for sample in patient['samples']:
        # Find run and adaID of the sample
        miseq_run = samples_seq[sample]['run']
        dataset = MiSeq_runs[miseq_run]
        data_folder = dataset['folder']
        adaID = dataset['adapters'][dataset['samples'].index(sample)]

        # Get input filename
        input_filename = get_merged_consensus_filename(data_folder, adaID, fragments)
        seq = SeqIO.read(input_filename, 'fasta')

        # Change name
        seq.description = ', '.join([seq.id, seq.description])
        seq.id = seq.name = sample

        consensi.append(seq)

    return consensi



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose

    # Get patient and the sequenced samples
    patient = get_patient(pname)
    patient['samples'] = get_sequenced_samples(patient)

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Collect and align fragment by fragment
    for fragment in fragments:
        consensi = get_consensi_frag(patient, fragment, VERBOSE=VERBOSE)
        ali = align_muscle(*consensi)
        output_filename = get_consensi_alignment_filename(pname, fragment)
        AlignIO.write(ali, output_filename, 'fasta')

    # Collect and align genome wide consensi
    consensi = get_consensi_genomewide(patient, VERBOSE=VERBOSE)

    # Align them
    ali = align_muscle(*consensi)

    # Write alignment
    output_filename = get_consensi_alignment_genomewide_filename(pname)
    AlignIO.write(ali, output_filename, 'fasta')

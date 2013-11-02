# vim: fdm=marker
'''
author:     Fabio Zanini
date:       25/10/13
content:    Collect and align all consensi from a patient.
'''
# Modules
import argparse
from Bio import SeqIO
from Bio import AlignIO

from mapping.datasets import MiSeq_runs
from mapping.filenames import get_merged_consensus_filename
from mapping.patients.filenames import get_consensi_alignment_genomewide_filename
from mapping.patients.patients import get_patient, get_sequenced_samples
from mapping.mapping_utils import align_muscle



# Functions
def get_consensi_genomewide(patient):
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

    # FIXME: skip the fragments for now, only do the genome wide consensus
    consensi_genomewide = get_consensi_genomewide(patient)

    # Align them
    ali_genomewide = align_muscle(*consensi_genomewide)

    # Write alignment
    output_filename_genomewide = get_consensi_alignment_genomewide_filename(pname)
    AlignIO.write(ali_genomewide, output_filename_genomewide, 'fasta')

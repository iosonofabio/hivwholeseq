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

from hivwholeseq.sequencing.samples import samples as samples_seq
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.filenames import get_merged_consensus_filename, get_consensus_filename
from hivwholeseq.patients.filenames import get_consensi_alignment_genomewide_filename, \
        get_consensi_alignment_filename
from hivwholeseq.patients.patients import patients, get_patient
from hivwholeseq.mapping_utils import align_muscle



# Functions
def make_output_folders(pname, VERBOSE=0):
    '''Make the patient folder and similia'''
    from hivwholeseq.generic_utils import mkdirs
    from hivwholeseq.patients.filenames import get_foldername
    mkdirs(get_foldername(pname))
    if VERBOSE >= 1:
        print 'Folder created:', get_foldername(pname)


def get_consensi_frag(patient, fragment, VERBOSE=0):
    '''Collect all consensi for a single fragment'''
    consensi = []
    for (index, sample) in patient.sample_table.iterrows():
        seq_run = sample['run']
        adaID = sample['adaID']
        date = sample['date']

        dataset = MiSeq_runs[seq_run]
        data_folder = dataset['folder']

        input_filename = get_consensus_filename(data_folder, adaID, fragment)
        if os.path.isfile(input_filename):
            seq = SeqIO.read(input_filename, 'fasta')
            seq.description = ', '.join([seq.id, seq.description])
            seq.id = seq.name = sample['name']+'_'+date

            consensi.append(seq)

    return consensi


def get_consensi_genomewide(patient, VERBOSE=0):
    '''Collect the consensi'''
    # FIXME: assume all fragments are there for now (no patchy stuff yet)
    fragments = ['F'+str(i) for i in xrange(1, 7)]
    consensi = []

    for (index, sample) in patient.sample_table.iterrows():
        seq_run = sample['run']
        adaID = sample['adaID']
        date = sample['date']

        dataset = MiSeq_runs[seq_run]
        data_folder = dataset['folder']

        input_filename = get_merged_consensus_filename(data_folder, adaID, fragments)
        if os.path.isfile(input_filename):
            seq = SeqIO.read(input_filename, 'fasta')
            seq.description = ', '.join([seq.id, seq.description])
            seq.id = seq.name = sample['name']+'_'+date

            consensi.append(seq)

    return consensi


def align_consensi(consensi, VERBOSE=0):
    '''Align with muscle and sort'''
    from Bio.Align import MultipleSeqAlignment as MSA

    ali = align_muscle(*consensi)

    ali_sorted = []
    for cons in consensi:
        for row in ali:
            if row.id == cons.id:
                ali_sorted.append(row)
                break

    ali_sorted = MSA(ali_sorted)
    return ali_sorted



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Align consensi from all time points')
    parser.add_argument('--patients', default=None, nargs='*',
                        help='Patients to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose

    if pnames is None:
        pnames = [p.id for p in patients]

    for pname in pnames:
        if VERBOSE >= 1:
            print pname,

        # Get patient and the sequenced samples
        patient = get_patient(pname)

        make_output_folders(pname)

        # If the script is called with no fragment, iterate over all
        if not fragments:
            fragments = ['F'+str(i) for i in xrange(1, 7)] + ['genomewide']
        if VERBOSE >= 3:
            print 'fragments', fragments

        # Collect and align fragment by fragment
        for fragment in fragments:
            if fragment == 'genomewide':
                continue
            if VERBOSE >= 2:
                print fragment,
            consensi = get_consensi_frag(patient, fragment, VERBOSE=VERBOSE)
            if len(consensi) < 2:
                continue

            ali = align_consensi(consensi)
            output_filename = get_consensi_alignment_filename(pname, fragment)
            AlignIO.write(ali, output_filename, 'fasta')
        if (VERBOSE >= 2) and ('genomewide' not in fragments):
            print

        # Collect and align genome wide consensi
        if 'genomewide' in fragments:
            if VERBOSE >= 2:
                print 'genomewide'
            consensi = get_consensi_genomewide(patient, VERBOSE=VERBOSE)
            if len(consensi) < 2:
                continue
        
            # Align them
            ali = align_consensi(consensi)
        
            # Write alignment
            output_filename = get_consensi_alignment_genomewide_filename(pname)
            AlignIO.write(ali, output_filename, 'fasta')

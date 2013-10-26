# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/13
content:    Update the initial consensus for the paitent (when we get more data).
'''
# Modules
import os
import argparse
import shutil
from Bio import SeqIO

from mapping.filenames import get_patient_initial_consensus_filename, \
        get_patient_foldername, get_consensus_filename
from mapping.datasets import MiSeq_runs
from mapping.patients.patients import get_patient, get_initial_sequenced_sample
from mapping.samples import samples, date_to_integer



# Functions



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

    # Get the patient
    patient = get_patient(pname)

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Make dir for the patient if absent
    pfolder = get_patient_foldername(pname)
    if not os.path.isdir(pfolder):
        os.mkdir(pfolder)
        if VERBOSE >= 1:
            print pname+': folder created.'
    
    # Get the first sequenced sample
    sample_init = get_initial_sequenced_sample(patient)
    miseq_run = samples[sample_init]['run']
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']
    adaID = dataset['adapters'][dataset['samples'].index(sample_init)]

    # Check for the existence of an initial consensus already
    for fragment in fragments:
        input_filename = get_consensus_filename(data_folder, adaID, fragment)
        output_filename = get_patient_initial_consensus_filename(pname, fragment)
        # If absent, just copy the thing over
        if not os.path.isfile(output_filename):
            shutil.copy(input_filename, output_filename)
            if VERBOSE >= 1:
                print pname+': initial consensus file created for sample', sample_init

        # if present, check whether the sequences are the same (if so, no remapping
        # is needed!). Overwrite the file anyway, because single samples carry
        # their consensus (mapping reference) with them in the folder (not much
        # overhead and MUUUCH cleaner than otherwise).
        else:
            seq_in = SeqIO.read(input_filename, 'fasta')
            seq_out = SeqIO.read(output_filename, 'fasta')

            if str(seq_in.seq) != str(seq_out.seq):
                print 'NOTE: initial consensus updated to '+sample_init+', remap!'
            
            SeqIO.write(seq_in, output_filename, 'fasta')

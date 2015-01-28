# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/11/14
content:    Store annotated initial references and other sequences for the
            website.
'''
# Modules
import os
import sys
from Bio import SeqIO

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.website.filenames import get_patient_reference_filename
from hivwholeseq.utils.sequence import correct_genbank_features_save
from hivwholeseq.reference import load_custom_reference



# Script
if __name__ == '__main__':

    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name

        ref = patient.get_reference('genomewide', format='gb')
        ref.id = 'ref_'+patient.code
        ref.name = 'ref_'+patient.code
        ref.description = 'reference for patient '+patient.code+', genomewide'

        # NOTE: the genbank Biopython writer is buggy, or the Genbank format is
        # buggy, anyway we need to trick feature types into notes.
        correct_genbank_features_save(ref)

        fn_out = get_patient_reference_filename(patient.code, format='gb')
        SeqIO.write(ref, fn_out, 'gb')

    # Store HXB2
    print 'HXB2'
    ref = load_custom_reference('HXB2', format='gb')
    correct_genbank_features_save(ref)
    fn_out = get_patient_reference_filename('HXB2', format='gb')
    SeqIO.write(ref, fn_out, 'gb')

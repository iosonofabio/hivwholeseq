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
from Bio import AlignIO

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.website.filenames import get_consensi_alignment_filename
from hivwholeseq.sequence_utils import correct_genbank_features_save



# Script
if __name__ == '__main__':

    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        for fragment in ['F'+str(i) for i in xrange(1, 7)] + ['genomewide']:
            print patient.code, patient.name, fragment

            try:
                ali = patient.get_consensi_alignment(fragment)
            except IOError:
                continue

            # Relabel the sequences
            for i in xrange(len(ali)):
                # The initial reference has the patient name
                if 'init' in ali[i].name:
                    ali[i].id = '_'.join([patient.code, fragment, 'initial_reference'])
                    ali[i].name = ali[i].id
                    ali[i].description = ', '.join([patient.code, fragment,
                                                    'initial reference'])
                
                # The other sequences have the sample name
                else:
                    tstr = (ali[i].id).split('_')[0]
                    ali[i].id = '_'.join([patient.code, fragment, tstr, 'days'])
                    ali[i].name = ali[i].id
                    ali[i].description = ', '.join([patient.code, fragment, tstr+' days'])

                print ali[i].id, ali[i].description


            # Write output
            fn_out = get_consensi_alignment_filename(patient.code, fragment,
                                                     format='fasta')
            AlignIO.write(ali, fn_out, 'fasta')


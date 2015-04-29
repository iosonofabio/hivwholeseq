# vim: fdm=indent
'''
author:     Fabio Zanini
date:       29/03/15
content:    Store sample tables for website.
'''
# Modules
import os
import sys
import numpy as np

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.website.filenames import get_sample_table_filename


# Script
if __name__ == '__main__':

    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name

        # Samples table
        tables = [{'kind': 'quantitative', 'suffix': 'q'},
                  {'kind': 'qualitative', 'suffix': ''}]

        for tdata in tables:

            table = (patient.samples[['days since infection'] +
                                     ['F'+str(i)+tdata['suffix'] for i in xrange(1, 7)]]
                     .set_index('days since infection'))

            table['RNA templates'] = np.array(patient.n_templates, int)

            # Write output
            fn_out = get_sample_table_filename(patient.code, kind=tdata['kind'])
            table.to_csv(fn_out, sep='\t')

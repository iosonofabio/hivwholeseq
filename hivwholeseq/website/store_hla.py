# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/11/14
content:    Store HLA of patients in a format suitable for the website.
'''
# Modules
import os
import sys
import numpy as np

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.website.filenames import get_HLA_filename



# Script
if __name__ == '__main__':

    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name


        

        # Write output
        fn_out = get_HLA_filename(patient.code)
        header = ['Locus', 'Type', 'Not excluded']
        hla = []
        for t in ('A', 'B', 'C', 'DRB1', 'DQB1'):
            d = [t, None, None]
            d[1] = '/'.join([patient['HLA-'+t],
                             patient['HLA-'+t+'-2']])
            d[2] = []            
            for suffix in ('', '-2'):
                datum = patient['HLA-'+t+suffix+' not excluded']
                if isinstance(datum, basestring):
                    d[2].append(str(datum))
            d[2] = ' '.join(d[2])

            hla.append(d)

        np.savetxt(fn_out,
                   hla,
                   fmt='%s',
                   delimiter='\t',
                   header='\t'.join(header),
                   comments='#')

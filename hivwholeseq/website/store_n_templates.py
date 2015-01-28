# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/11/14
content:    Store number of templates to PCR, in a format suitable for the website.
'''
# Modules
import os
import sys
import numpy as np

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.website.filenames import get_ntemplates_filename



# Script
if __name__ == '__main__':

    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name

        n_templates = patient.n_templates
        ind = -(np.ma.getmaskarray(n_templates))
        n_templates = np.array(n_templates[ind])
        times = patient.times[ind]

        # Write output
        fn_out = get_ntemplates_filename(patient.code)
        np.savetxt(fn_out,
                   np.vstack([times, n_templates]).T,
                   delimiter='\t',
                   header='t\tn_templates',
                   comments='')


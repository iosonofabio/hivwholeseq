# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/11/14
content:    Store physiological data, e.g. viral load and CD4+, in a format
            suitable for the website.
'''
# Modules
import os
import sys
import numpy as np

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.website.filenames import get_viral_load_filename, get_cell_count_filename



# Script
if __name__ == '__main__':

    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name

        # Viral load
        viral_load = np.array(patient.viral_load)
        ind = viral_load > 0
        times = patient.times[ind]
        viral_load = viral_load[ind]

        # Write output
        fn_out = get_viral_load_filename(patient.code)
        np.savetxt(fn_out,
                   np.vstack([times, viral_load]).T,
                   delimiter='\t',
                   header='t\tvl',
                   comments='')

        # CD4+ counts
        counts = np.array(patient.cell_count)
        ind = counts > 0
        times = patient.times[ind]
        counts = counts[ind]

        # Write output
        fn_out = get_cell_count_filename(patient.code)
        np.savetxt(fn_out,
                   np.vstack([times, counts]).T,
                   delimiter='\t',
                   header='t\tcc',
                   comments='')


# vim: fdm=marker
'''
author:     Fabio Zanini
date:       12/11/14
content:    Store reads in a way that's appropriate for the website (only
            symlinks).
'''
# Modules
import os
import shutil
import sys
import numpy as np


from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.website.filenames import get_reads_filename, get_timeline_filename
from hivwholeseq.patients.patients import load_patients, Patient



# Script
if __name__ == '__main__':

    fragments = ['F'+str(i) for i in xrange(1, 7)]
    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name
        
        # reads
        for fragment in fragments:
            print fragment
            for it, sample in enumerate(patient.itersamples()):
                print it, sample.name

                fn = sample.get_mapped_filtered_filename(fragment, decontaminated=True)
                if not os.path.isfile(fn):
                    continue

                fn_out = get_reads_filename(patient.code, fragment, it)
                if not os.path.isfile(fn_out):
                    mkdirs(os.path.dirname(fn_out))
                    os.symlink(fn, fn_out)

        # time points
        fn_out = get_timeline_filename(patient.code)
        np.savetxt(fn_out, patient.times, header='Time from infection [days]')

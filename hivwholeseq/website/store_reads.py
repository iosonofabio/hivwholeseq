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


from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.website.filenames import get_reads_filename, get_timeline_filename
from hivwholeseq.patients.patients import load_patients, Patient


# Functions
def copy_or_symlink_reads(fn_in, fn_out, maxsize=20e6):
    '''Symlink reads or, if too many, copy a subsample'''
    import os

    # Remove previous file/link
    if os.path.isfile(fn_out):
        os.remove(fn_out)
    elif os.path.islink(fn_out):
        os.unlink(fn_out)

    if os.stat(fn_in).st_size <= maxsize:
        os.symlink(fn_in, fn_out)
    else:
        from hivwholeseq.utils.mapping import extract_mapped_reads_subsample
        extract_mapped_reads_subsample(fn_in, fn_out, 50000, VERBOSE=2)



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
                print it, sample.name, fragment

                fn = sample.get_mapped_filtered_filename(fragment, decontaminated=True)
                if not os.path.isfile(fn):
                    continue

                fn_out = get_reads_filename(patient.code, fragment, it)

                mkdirs(os.path.dirname(fn_out))

                copy_or_symlink_reads(fn, fn_out)

        # time points
        fn_out = get_timeline_filename(patient.code)
        np.savetxt(fn_out, patient.times, header='Time from infection [days]')

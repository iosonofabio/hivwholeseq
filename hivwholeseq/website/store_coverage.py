# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/11/14
content:    Store coverage for the website.
'''
# Modules
import os
import sys
import numpy as np

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.website.filenames import get_coverage_filename
from hivwholeseq.patients.filenames import get_allele_count_trajectories_filename



# Script
if __name__ == '__main__':

    fragments = ['F'+str(i) for i in xrange(1, 7)]
    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name

        # Coverage
        fn_out = get_allele_count_trajectories_filename(pname, 'genomewide')
        npz = np.load(fn_out)
        ind = npz['ind']
        act = npz['act']

        times = patient.times[ind]
        cov = act.sum(axis=1)

        # Write output
        fn_out = get_coverage_filename(patient.code, 'genomewide')
        np.savez(fn_out, times=times, cov=cov)

# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/11/14
content:    Store divergence and diversity for the website.
'''
# Modules
import os
import sys
import numpy as np

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.website.filenames import get_divergence_filename, get_diversity_filename



# Script
if __name__ == '__main__':

    fragments = ['F'+str(i) for i in xrange(1, 7)]
    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name

        for fragment in fragments:
            # Divergence
            dg, ind = patient.get_divergence(fragment,
                                             cov_min=100)
            times = patient.times[ind]
            ind = -dg.mask
            dg = dg[ind].data
            times = times[ind]
            
            # Write output
            fn_out = get_divergence_filename(patient.code, fragment)
            np.savetxt(fn_out,
                       np.vstack([times, dg]).T,
                       delimiter='\t',
                       header='t\tdg',
                       comments='')

            # Diversity
            ds, ind = patient.get_diversity(fragment,
                                            cov_min=100)
            times = patient.times[ind]
            ind = -ds.mask
            ds = ds[ind].data
            times = times[ind]

            # Write output
            fn_out = get_diversity_filename(patient.code, fragment)
            np.savetxt(fn_out,
                       np.vstack([times, ds]).T,
                       delimiter='\t',
                       header='t\tds',
                       comments='')

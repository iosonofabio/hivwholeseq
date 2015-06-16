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

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, iterpatient
from hivwholeseq.website.filenames import get_divergence_filename, get_diversity_filename
from hivwholeseq.website import _regions as regions_all



# Script
if __name__ == '__main__':

    fragments = ['F'+str(i) for i in xrange(1, 7)]
    regions = fragments + regions_all
    patients = load_patients()
    for pname, patient in iterpatient(patients):
        for region in regions:
            print patient.code, patient.name, region

            # Divergence
            dg, ind = patient.get_divergence(region,
                                             cov_min=100)
            times = patient.times[ind]
            ind = -dg.mask
            dg = dg[ind].data
            times = times[ind]
            
            # Write output
            fn_out = get_divergence_filename(patient.code, region)
            np.savetxt(fn_out,
                       np.vstack([times, dg]).T,
                       delimiter='\t',
                       header='t [days]\tdg',
                       comments='')


            # Diversity
            ds, ind = patient.get_diversity(region,
                                            cov_min=100)
            times = patient.times[ind]
            ind = -ds.mask
            ds = ds[ind].data
            times = times[ind]

            # Write output
            fn_out = get_diversity_filename(patient.code, region)
            np.savetxt(fn_out,
                       np.vstack([times, ds]).T,
                       delimiter='\t',
                       header='t [days]\tds',
                       comments='')

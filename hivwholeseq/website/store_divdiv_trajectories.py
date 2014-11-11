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
from hivwholeseq.website.filenames import get_divergence_trajectories_filename, \
        get_diversity_trajectories_filename
from hivwholeseq.patients.filenames import get_divergence_trajectories_local_filename, \
        get_diversity_trajectories_local_filename



# Script
if __name__ == '__main__':

    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name

        # Divergence
        npz = np.load(get_divergence_trajectories_local_filename(patient.name,
                                                                 'genomewide'))
        ind = npz['ind']
        dg = npz['dg']
        block_length = npz['block_length']
        times = patient.times[ind]

        # Write output
        fn_out = get_divergence_trajectories_filename(patient.code, 'genomewide')
        np.savez(fn_out, times=times, dg=dg, block_length=block_length)

        # Diversity
        npz = np.load(get_diversity_trajectories_local_filename(patient.name,
                                                                'genomewide'))
        ind = npz['ind']
        ds = npz['ds']
        block_length = npz['block_length']
        times = patient.times[ind]

        # Write output
        fn_out = get_diversity_trajectories_filename(patient.code, 'genomewide')
        np.savez(fn_out, times=times, ds=ds, block_length=block_length)


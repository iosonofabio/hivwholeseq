# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/11/14
content:    Store allele counts for the website. We need two formats, one binary
            fast-access for the plots, and one annotated and standard for
            downloads.

            We store both allele count trajectories for whole patients and sample
            by sample.
'''
# Modules
import os
import sys
import shutil
import numpy as np
import pandas as pd

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.website.filenames import get_insertion_trajectories_filename as get_fn_out_traj
from hivwholeseq.website.filenames import get_insertions_filename as get_fn_out_sample


# Globals



# Functions



# Script
if __name__ == '__main__':

    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name

        # Allele count trajectories
        (inse, ind) = patient.get_insertion_trajectories('genomewide')
        if not ind:
            continue
        inse = pd.Series(inse, name='insertion counts')
        inse.index.names = ['DSI', 'position', 'insertion']

        # Write to file
        fn_out = get_fn_out_traj(patient.code, 'genomewide')
        mkdirs(os.path.dirname(fn_out))
        inse.to_pickle(fn_out)

        # Sample by sample
        print 'Sample by sample'
        for i, sample in enumerate(patient.itersamples()):
            samplename = patient.code+'_sample_'+str(i+1)
            print "Sample time:", sample['days since infection'], sample.name

            for region in ['F'+str(j) for j in xrange(1, 7)] + ['genomewide']:
                try:
                    inse = sample.get_insertions(region)
                except IOError:
                    continue
                if not inse:
                    continue
                inse = pd.Series(inse, name='insertion counts')
                inse.index.names = ['position', 'insertion']

                # Write to file
                fn_out = get_fn_out_sample(samplename, region)
                mkdirs(os.path.dirname(fn_out))
                inse.to_pickle(fn_out)

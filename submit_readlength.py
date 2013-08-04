#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=indent
'''
author:     Richard Neher, adapted from Fabio Zanini
date:       02/08/13
content:    Submit script for read length distribution and q control
'''
# Standard modules
import subprocess as sp
import os
import sys
import numpy as np
import glob

# Environment vars
JOBSUBMIT = os.path.realpath(__file__)      # This file
JOBDIR = os.path.dirname(JOBSUBMIT)+'/'
JOBSCRIPT = JOBDIR+'read_length_distribution.py'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'


# Cluster parameters
cluster_time = '0:59:59'
vmem = '8G'
VERBOSE = 1



# Script
if __name__ == '__main__':

    dir_list = glob.glob('/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/subsample/adapterID*')
    dir_list.extend(glob.glob('/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/subsample/unclass*'))
    for dirname in dir_list:
        sample_dir = dirname.split('/')[-1]
        # Create the path and the model (every run has a slightly
        # different model) and run an array of jobs with -t
        qsub_list = ['qsub','-cwd',
                     '-o',JOBLOGOUT,
                     '-e',JOBLOGERR,
                     '-l', 'h_rt='+cluster_time,
                     '-l', 'h_vmem='+vmem,
                     JOBSCRIPT, sample_dir]
        qsub_list = map(str,qsub_list)
        if VERBOSE:
            print ' '.join(qsub_list)
        sp.call(qsub_list)


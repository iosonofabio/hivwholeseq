#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/08/13
content:    Submit script for demultiplexing.
'''
# Standard modules
import subprocess as sp
import os
import sys
import numpy as np

# Environment vars
JOBSUBMIT = os.path.realpath(__file__)      # This file
JOBDIR = os.path.dirname(JOBSUBMIT)+'/'
JOBSCRIPT = JOBDIR+'demultiplex.py'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'


# Cluster parameters
cluster_time = '23:59:59'
vmem = '8G'
VERBOSE = 1



# Script
if __name__ == '__main__':

        # Create the path and the model (every run has a slightly
        # different model) and run an array of jobs with -t
        qsub_list = ['qsub','-cwd',
                     '-o',JOBLOGOUT,
                     '-e',JOBLOGERR,
                     '-l', 'h_rt='+cluster_time,
                     '-l', 'h_vmem='+vmem,
                     JOBSCRIPT]
        qsub_list = map(str,qsub_list)
        if VERBOSE:
            print ' '.join(qsub_list)
        sp.call(qsub_list)


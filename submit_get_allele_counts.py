#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       06/08/13
content:    Submit script for extracting allele and insert counts from the mapped
            reads.
'''
# Modules
import os
import sys
import numpy as np
import subprocess as sp
from map_HIV_HXB2 import load_adapter_table



# Globals
VERBOSE = 1
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'

# Submit vars
JOBSUBMIT = os.path.realpath(__file__)      # This file
JOBDIR = os.path.dirname(JOBSUBMIT)+'/'
JOBSCRIPT = JOBDIR+'get_allele_counts.py'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'

# Cluster parameters
cluster_time = '0:59:59'
vmem = '4G'
    


# Script
if __name__ == '__main__':
    
    # Iterate over adapter IDs
    adapter_table = load_adapter_table(data_folder)
    for adaID in adapter_table['ID']:

        # Call script
        qsub_list = ['qsub','-cwd',
                     '-o', JOBLOGOUT,
                     '-e', JOBLOGERR,
                     '-N', 'ac_'+'{:02d}'.format(adaID),
                     '-l', 'h_rt='+cluster_time,
                     '-l', 'h_vmem='+vmem,
                     JOBSCRIPT,
                     data_folder+'adapterID_'+'{:02d}'.format(adaID)+'/'
                    ]
        qsub_list = map(str,qsub_list)
        if VERBOSE:
            print ' '.join(qsub_list)
        sp.call(qsub_list)

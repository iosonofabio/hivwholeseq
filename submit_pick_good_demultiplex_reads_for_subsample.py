#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/08/13
content:    Take demultiplexed reads and pick a subsample to build a good consensus.
            (submit script)
'''
# Modules
import os
import sys
import numpy as np
import subprocess as sp



# Globals
VERBOSE = 1

# Submit vars
JOBSUBMIT = os.path.realpath(__file__)      # This file
JOBDIR = os.path.dirname(JOBSUBMIT)+'/'
JOBSCRIPT = JOBDIR+'pick_good_demultiplex_reads_for_subsample.py'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'

# Cluster parameters
cluster_time = '0:59:59'
vmem = '2G'
    


# Script
if __name__ == '__main__':
    
    # Iterate over adapter IDs
    adapter_table = load_adapter_table(data_folder)
    for adaID in adapter_table['ID']:

        # Call script
        qsub_list = ['qsub','-cwd',
                     '-o',JOBLOGOUT,
                     '-e',JOBLOGERR,
                     '-N', 'pick '+'{:02d}'.format(adaID),
                     '-l', 'h_rt='+cluster_time,
                     '-l', 'h_vmem='+vmem,
                     JOBSCRIPT,
                     adaID]
        qsub_list = map(str,qsub_list)
        if VERBOSE:
            print ' '.join(qsub_list)
        sp.call(qsub_list)

#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/08/13
content:    Map the reads to their own consensus, produced with a preliminary
            mapping to HXB2. This is the final mapping.
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
data_folder_subsample = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/subsample/'
adapters_table_file = 'adapters_table.dat'
stampy_bin = '/ebio/ag-neher/share/programs/bundles/stampy-1.0.22/stampy.py'
subsrate = '0.02'
use_filtered_reads = True

# Prepare I/O filenames
if use_filtered_reads: file_suffix = '_filtered_trimmed'
else: file_suffix = ''
consensus_file = 'consensus'+file_suffix+'.fasta'

# Submit vars
JOBSUBMIT = os.path.realpath(__file__)      # This file
JOBDIR = os.path.dirname(JOBSUBMIT)+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'

# Cluster parameters (short time for subsample only)
cluster_time = '23:59:59'
vmem = '8G'
    


# Script
if __name__ == '__main__':
  
    # Iterate over all adapter IDs
    # Note: the hash tables for mapping can be done on the cluster frontend,
    # because the HIV genome is small
    adapter_table = load_adapter_table(data_folder)
    for adaID in adapter_table['ID']:

        # Directory to read
        dirname = 'adapterID_'+'{:02d}'.format(adaID)+'/'

        # 1. Make genome index file
        if not os.path.isfile(data_folder_subsample+dirname+'consensus'+file_suffix+'.stidx'):
            sp.call([stampy_bin,
                     '--species=HIV',
                     '-G', data_folder_subsample+dirname+'consensus'+file_suffix,
                     data_folder_subsample+dirname+consensus_file,
                     ])
        
        # 2. Build a hash file
        if not os.path.isfile(data_folder_subsample+dirname+'consensus'+file_suffix+'.sthash'):
            sp.call([stampy_bin,
                     '-g', data_folder_subsample+dirname+'consensus'+file_suffix,
                     '-H', data_folder_subsample+dirname+'consensus'+file_suffix,
                     ])

        # 3. Map (using STAMPY)
        # Note: no --solexa option as of 2013 (illumina v1.8)
        qsub_list = ['qsub','-cwd',
                     '-o',JOBLOGOUT,
                     '-e',JOBLOGERR,
                     '-N', 'stampy_'+'{:02d}'.format(adaID),
                     '-l', 'h_rt='+cluster_time,
                     '-l', 'h_vmem='+vmem,
                     stampy_bin,
                     '-g', data_folder_subsample+dirname+'consensus'+file_suffix,
                     '-h', data_folder_subsample+dirname+'consensus'+file_suffix, 
                     '-o', data_folder+dirname+'mapped_to_self'+file_suffix+'.sam',
                     '--substitutionrate='+subsrate,
                     '-M',
                     data_folder+dirname+'read1'+file_suffix+'.fastq',
                     data_folder+dirname+'read2'+file_suffix+'.fastq']
        qsub_list = map(str,qsub_list)
        if VERBOSE:
            print ' '.join(qsub_list)
        sp.call(qsub_list)

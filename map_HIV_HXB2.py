#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/08/13
content:    Map the reads to HXB2 to start. The resulting consensus sequence is
            used for recursive mapping as a second step (other script).
'''
# Modules
import os
import sys
import numpy as np
import subprocess as sp



# Globals
VERBOSE = 1
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/subsample/'
adapters_table_file = 'adapters_table.dat'
HXB2_file = 'HXB2.fa.gz'
stampy_bin = '/ebio/ag-neher/share/programs/bundles/stampy-1.0.22/stampy.py'
subsrate = '0.05'

# Submit vars
JOBSUBMIT = os.path.realpath(__file__)      # This file
JOBDIR = os.path.dirname(JOBSUBMIT)+'/'
JOBSCRIPT = JOBSUBMIT
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'

# Cluster parameters
cluster_time = '23:59:59'
vmem = '8G'
    



# Functions
def load_adapter_table(data_folder):
    '''Load table of adapters and samples'''
    table = np.loadtxt(data_folder+adapters_table_file,
                       dtype=[('seq', 'S6'), ('ID', int), ('sample', 'S50')],
                       ndmin=1)
    return table



# Script
if __name__ == '__main__':

    # 1. Make genome index file (HXB2)
    if not os.path.isfile(data_folder+'HIV_HXB2.stidx'):
        sp.call([stampy_bin,
                 '--species=HIV',
                 '-G', data_folder+'HIV_HXB2',
                 HXB2_file,
                 ])
    
    # 2. Build a hash file
    if not os.path.isfile(data_folder+'HIV_HXB2.sthash'):
        sp.call([stampy_bin,
                 '-g', data_folder+'HIV_HXB2',
                 '-H', data_folder+'HIV_HXB2',
                 ])
    
    # 3. Map (using STAMPY)
    adapter_table = load_adapter_table(data_folder)
    for adaID in adapter_table['ID']:

        # FIXME
        if adaID == 2:
            continue

        # Directory to read
        dirname = 'adapterID_'+'{:02d}'.format(adaID)+'/'

        # Call stampy
        # Note: no --solexa option as of 2013 (illumina v1.8)
        qsub_list = ['qsub','-cwd',
                     '-o',JOBLOGOUT,
                     '-e',JOBLOGERR,
                     '-N', 'stampy_'+'{:02d}'.format(adaID),
                     '-l', 'h_rt='+cluster_time,
                     '-l', 'h_vmem='+vmem,
                     stampy_bin,
                     '-g', data_folder+'HIV_HXB2',
                     '-h', data_folder+'HIV_HXB2', 
                     '-o', data_folder+dirname+'mapped_to_HXB2.sam',
                     '--substitutionrate='+subsrate,
                     '-M', data_folder+dirname+'read1.fastq', data_folder+dirname+'read2.fastq']
        qsub_list = map(str,qsub_list)
        if VERBOSE:
            print ' '.join(qsub_list)
        sp.call(qsub_list)

#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/08/13
content:    Map the unclassified reads to phiX as a control.
'''
# Modules
import os
import subprocess as sp



# Globals
VERBOSE = 1
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'
phiX_file = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/phiX_genome.fasta'
stampy_bin = '/ebio/ag-neher/share/programs/bundles/stampy-1.0.22/stampy.py'

# Submit vars
JOBSUBMIT = os.path.realpath(__file__)      # This file
JOBDIR = os.path.dirname(JOBSUBMIT)+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'

# Cluster parameters
cluster_time = '23:59:59'
vmem = '8G'
    


# Script
if __name__ == '__main__':

    # 1. Make genome index file
    if not os.path.isfile(data_folder+'phiX.stidx'):
        sp.call([stampy_bin,
                 '--species=phiX',
                 '-G', data_folder+'phiX',
                 phiX_file,
                 ])
    
    # 2. Build a hash file
    if not os.path.isfile(data_folder+'phiX.sthash'):
        sp.call([stampy_bin,
                 '-g', data_folder+'phiX',
                 '-H', data_folder+'phiX',
                 ])
    
    # 3. Map (using STAMPY)
    # Directory to read
    dirname = 'unclassified_reads/'

    # Call stampy
    # Note: no --solexa option as of 2013 (illumina v1.8)
    qsub_list = ['qsub','-cwd',
                 '-o',JOBLOGOUT,
                 '-e',JOBLOGERR,
                 '-N', 'stampy_phiX',
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 stampy_bin,
                 '-g', data_folder+'phiX',
                 '-h', data_folder+'phiX', 
                 '-o', data_folder+dirname+'mapped_to_phiX_filtered_trimmed.sam',
                 '-M', data_folder+dirname+'read1_filtered_trimmed.fastq', data_folder+dirname+'read2_filtered_trimmed.fastq']
    qsub_list = map(str,qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)

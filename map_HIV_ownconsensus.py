#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/08/13
content:    Map the reads to their own consensus, produced with an iterative
            mapping of HIV onto itself (HXB2 at first). This is the final mapping.
'''
# Modules
import os
import sys
import numpy as np
import subprocess as sp
from adapter_info import load_adapter_table, foldername_adapter
from mapping.mapping_utils import stampy_bin, subsrate



# Globals
VERBOSE = 1
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']


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
    # Note: the hash tables for mapping might have been done by the recursive
    # map already, depending on convergence of that
    adapter_table = load_adapter_table(data_folder)
    for adaID in adapter_table['ID']:

        # Directory to read
        dirname = foldername_adapter(adaID)

        # Find the last (fragmented) consensus
        g = os.walk(data_folder+'subsample/'+dirname)
        fns = g.next()[2]
        if VERBOSE >= 3:
            print fns
        fns = filter(lambda x: ('consensus' in x) and ('fragmented.fasta' in x), fns)
        consensus_numbers = map(lambda x: int(x.split('_')[1]), fns)
        cons_max = max(consensus_numbers)
        if VERBOSE >= 2:
            print cons_max

        # Fragmented consensus filenames
        file_pattern = 'consensus_'+str(cons_max)+'_fragmented'
        ref_file = file_pattern+'.fasta'
        index_file = file_pattern+'.stidx'
        hash_file = file_pattern+'.sthash'

        # 1. Make genome index file
        if not os.path.isfile(data_folder+'subsample/'+dirname+index_file):
            sp.call([stampy_bin,
                     '--species=HIV',
                     '-G', data_folder+'subsample/'+dirname+file_pattern,
                     data_folder+'subsample/'+dirname+ref_file,
                     ])
        
        # 2. Build a hash file
        if not os.path.isfile(data_folder+'subsample/'+dirname+hash_file):
            sp.call([stampy_bin,
                     '-g', data_folder+'subsample/'+dirname+file_pattern,
                     '-H', data_folder+'subsample/'+dirname+file_pattern,
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
                     '-g', data_folder+'subsample/'+dirname+file_pattern,
                     '-h', data_folder+'subsample/'+dirname+file_pattern, 
                     '-o', data_folder+dirname+'mapped_to_consensus_'+str(cons_max)+'_fragmented.sam',
                     '--substitutionrate='+subsrate,
                     '-M',
                     data_folder+dirname+'read1_filtered_trimmed.fastq',
                     data_folder+dirname+'read2_filtered_trimmed.fastq']
        qsub_list = map(str,qsub_list)
        if VERBOSE:
            print ' '.join(qsub_list)
        sp.call(qsub_list)

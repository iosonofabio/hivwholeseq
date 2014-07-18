# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/05/14
content:    Make symlink to raw (demultiplexed) data to the Short Read Archive.
'''
# Modules
import os
import argparse
from itertools import izip
import numpy as np

from hivwholeseq.adapter_info import foldername_adapter
from hivwholeseq.samples import load_samples_sequenced, load_sequencing_run
from hivwholeseq.filenames import get_seqrun_foldername


# Functions
def make_output_folders(data_folder, adaIDs, VERBOSE=0):
    '''Make output folders for symlinking'''
    from hivwholeseq.generic_utils import mkdirs
    mkdirs(data_folder)
    if VERBOSE >= 1:
        print 'Folder created:', data_folder

    for adaID in adaIDs + [-1]:
        mkdirs(data_folder+foldername_adapter(adaID))
        if VERBOSE >= 1:
            print 'Folder created:', data_folder+foldername_adapter(adaID)


def make_symlinks(dataset, VERBOSE=0):
    '''Make symlinks for fastq.gz from the SRA'''
    seq_run = dataset.name
    samples = dataset.samples
    data_folder = get_seqrun_foldername(seq_run)
    raw_root_folder = dataset.loc['raw data']

    import re
    seq_run_int = int(re.findall(r'\d+$', seq_run)[0])

    # Unclassified reads
    unclass_fn = '/'.join(raw_root_folder.split('/')[:-2])+'/'
    for fn in os.listdir(unclass_fn):
        if ('illumina_M' in fn) and ('RunId'+'{:04d}'.format(seq_run_int) in fn):
            unclass_fn = unclass_fn+fn+'/LaneId1/'
            break
    fn1 = unclass_fn+[fn for fn in os.listdir(unclass_fn) if 'L001_R1' in fn][0]
    fn2 = unclass_fn+[fn for fn in os.listdir(unclass_fn) if 'L001_R2' in fn][0]

    dst_folder = data_folder+foldername_adapter(-1)
    dst_fn1 = dst_folder+'read1.fastq.gz'
    dst_fn2 = dst_folder+'read2.fastq.gz'
    if not os.path.isfile(dst_fn1):
        os.symlink(fn1, dst_fn1)
    elif VERBOSE:
            print dst_fn1, 'exists already'
    if not os.path.isfile(dst_fn2):
        os.symlink(fn2, dst_fn2)
    elif VERBOSE:
            print dst_fn2, 'exists already'
    if VERBOSE:
        print 'Unclassified reads symlinked'

    # Samples
    for sn, sample in samples.iterrows():
        if str(sample['raw name']) != 'nan':
            raw_fn = str(sample['raw name'])
        elif str(sample['patient sample']) != 'nan':
            raw_fn = str(sample['patient sample'])
        else:
            print 'ERROR: Sample '+sn+': could not find raw data folder. Please fill in the table'
            continue

        sample_fn = raw_root_folder+[fn for fn in os.listdir(raw_root_folder) 
                                     if raw_fn in fn][0]+'/'

        fn1 = sample_fn+[fn for fn in os.listdir(sample_fn) if 'L001_R1' in fn][0]
        fn2 = sample_fn+[fn for fn in os.listdir(sample_fn) if 'L001_R2' in fn][0]
        
        adaID = sample['adapter']
        dst_folder = data_folder+foldername_adapter(adaID)
        dst_fn1 = dst_folder+'read1.fastq.gz'
        dst_fn2 = dst_folder+'read2.fastq.gz'
        if not os.path.isfile(dst_fn1):
            os.symlink(fn1, dst_fn1)
        elif VERBOSE:
                print dst_fn1, 'exists already'
        if not os.path.isfile(dst_fn2):
            os.symlink(fn2, dst_fn2)
        elif VERBOSE:
                print dst_fn2, 'exists already'
        if VERBOSE:
            print sn+' '+adaID+' reads symlinked'



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Demultiplex HIV reads')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28, test_tiny)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    
    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose

    dataset = load_sequencing_run(seq_run)
    data_folder = dataset.folder

    adapters = dataset.samples.adapter.tolist()
    make_output_folders(data_folder, adapters, VERBOSE=VERBOSE)

    make_symlinks(dataset, VERBOSE=VERBOSE)

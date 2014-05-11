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

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.demultiplex import make_output_folders, get_adapters_designed
from hivwholeseq.adapter_info import foldername_adapter


# Functions
def make_symlinks(dataset, VERBOSE=0):
    '''Make symlinks for fastq.gz from the SRA'''
    import re
    seq_run = int(re.findall(r'\d+$', dataset['name'])[0])

    root_folder = dataset['raw_data']['SRA']
    data_folder = dataset['folder']

    # Unclassified reads
    unclass_fn = root_folder+'../'
    for fn in os.listdir(unclass_fn):
        if ('illumina_M' in fn) and ('RunId'+'{:04d}'.format(seq_run) in fn):
            unclass_fn = unclass_fn+fn+'/LaneId1/'
            break
    fn1 = unclass_fn+[fn for fn in os.listdir(unclass_fn) if 'L001_R1' in fn][0]
    fn2 = unclass_fn+[fn for fn in os.listdir(unclass_fn) if 'L001_R2' in fn][0]
    
    dst_folder = data_folder+foldername_adapter(-1)
    os.symlink(fn1, dst_folder+'read1.fastq.gz')
    os.symlink(fn2, dst_folder+'read2.fastq.gz')
    if VERBOSE:
        print 'Unclassified reads symlinked'

    # Samples
    for (samplename, adaID) in izip(dataset['samples'], dataset['adapters']):
        sn = samplename
        if (len(sn) > 5) and (sn[-5:] in ('_PCR1', '_PCR2')):
            sn = sn[:-5]
        sn = sn.replace('-', '')
        sample_fn = root_folder
        for fn in os.listdir(sample_fn):
            if sn in fn:
                sample_fn = sample_fn+fn+'/'
                break

        fn1 = sample_fn+[fn for fn in os.listdir(sample_fn) if 'L001_R1' in fn][0]
        fn2 = sample_fn+[fn for fn in os.listdir(sample_fn) if 'L001_R2' in fn][0]
        
        dst_folder = data_folder+foldername_adapter(adaID)
        os.symlink(fn1, dst_folder+'read1.fastq.gz')
        os.symlink(fn2, dst_folder+'read2.fastq.gz')
        if VERBOSE:
            print samplename+' '+adaID+' reads symlinked'



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

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']


    adapters_designed = get_adapters_designed(dataset, VERBOSE=VERBOSE, summary=False)

    make_output_folders(data_folder, adapters_designed, VERBOSE=VERBOSE,
                        summary=False)

    make_symlinks(dataset, VERBOSE=VERBOSE)

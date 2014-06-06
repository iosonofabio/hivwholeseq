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

    # FIXME
    return

    # Samples
    for (samplename, adaID) in izip(dataset['samples'], dataset['adapters']):
        sn = samplename
        sample_fn = root_folder
        sample_fnn = None

        # FIXME: get subfolder: the sample name might differ, only trust the adaID
        # FIXME: also not enough, because Stefan put all into the same folder... so
        # I need to write heuristic... !
        for fn in os.listdir(sample_fn):
            if sn in fn:
                sample_fnn = sample_fn+fn+'/'
                break

        if sample_fnn is None:
            if (len(sn) > 5) and (sn[-5:] in ('_PCR1', '_PCR2')):
                sn = sn[:-5]

            for fn in os.listdir(sample_fn):
                if sn in fn:
                    sample_fnn = sample_fn+fn+'/'
                    break

        if sample_fnn is None:
            if (len(sn) > 7) and (sn[-7:] in ('_PCR1-2', '_PCR2-2', '_PCR1-3', '_PCR2-3')):
                sn = sn[:-7]

            for fn in os.listdir(sample_fn):
                if sn in fn:
                    sample_fnn = sample_fn+fn+'/'
                    break

        if sample_fnn is None:
            sn = sn.replace('-', '')
            for fn in os.listdir(sample_fn):
                if sn in fn:
                    sample_fnn = sample_fn+fn+'/'
                    break

        if sample_fnn is None:
            sn = sn.replace('-', '')
            if (len(sn) > 5) and (sn[-5:] in ('_PCR1', '_PCR2')):
                sn = sn[:-5]

            for fn in os.listdir(sample_fn):
                if sn in fn:
                    sample_fnn = sample_fn+fn+'/'
                    break

        if sample_fnn is None:
            sn = sn.replace('-', '')
            if (len(sn) > 6) and (sn[-6:] in ('_PCR12', '_PCR22', '_PCR13', '_PCR23')):
                sn = sn[:-6]

            for fn in os.listdir(sample_fn):
                if sn in fn:
                    sample_fnn = sample_fn+fn+'/'
                    break

        if sample_fnn is None:
            print samplename, adaID
            import ipdb; ipdb.set_trace()
            continue

        fn1 = sample_fnn+[fn for fn in os.listdir(sample_fnn) if 'L001_R1' in fn][0]
        fn2 = sample_fnn+[fn for fn in os.listdir(sample_fnn) if 'L001_R2' in fn][0]
        
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

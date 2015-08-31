# vim: fdm=marker
'''
author:     Fabio Zanini
date:       04/02/14
content:    Remove temporary files (premapping, mapping) to save space.
'''
# Modules
import os
import glob
import subprocess as sp
import argparse


# Functions
def remove_premapped_tempfiles(data_folder, adaID, VERBOSE=0):
    '''Remove the part files of multi-threaded premapping'''
    from hivwholeseq.sequencing.filenames import get_premapped_filename

    dirname = os.path.dirname(get_premapped_filename(data_folder, adaID, type='bam', part=1))+'/'
    fns = glob.glob(dirname+'premapped_*part*') + \
          glob.glob(dirname+'premapped_*unsorted*')  
    fns.append(dirname+'premapped.sam')
    for fn in fns:
        os.remove(fn)
        if VERBOSE >= 3:
            print 'File removed:', fn


def remove_mapped_tempfiles(data_folder, adaID, fragment='F', VERBOSE=0, rescue=False):
    '''Remove the part files of multi-threaded mapping'''
    from hivwholeseq.sequencing.filenames import get_mapped_filename

    dirname = os.path.dirname(get_mapped_filename(data_folder, adaID, 'F1'))+'/'
    fns = glob.glob(dirname+fragment+'*_part*') + \
          glob.glob(dirname+fragment+'*_unsorted*') + \
          glob.glob(dirname+fragment+'*.00*.bam')
    fns.append(dirname+fragment+'.sam')
    if rescue:
        fns.append(dirname+fragment+'_rescue.sam')

    for fn in fns:
        os.remove(fn)
        if VERBOSE >= 3:
            print 'File removed:', fn


def remove_mapped_init_tempfiles(pname, samplename_pat, 
                                 samplename, fragment='F',
                                 PCR=1,
                                 VERBOSE=0, remove_chunks=False,
                                 only_chunk=None):
    '''Remove the part files of multi-threaded mapping to initial patient consensus'''
    from hivwholeseq.patients.filenames import get_mapped_to_initial_foldername, \
            get_mapped_to_initial_filename

    dirname = get_mapped_to_initial_foldername(pname, samplename_pat, PCR=PCR)
    prefix = samplename+'_'+fragment
    if only_chunk is None:
        fns = glob.glob(dirname+prefix+'*_part*') + \
              glob.glob(dirname+prefix+'*_unsorted*') + \
              glob.glob(dirname+prefix+'_part*') + \
              glob.glob(dirname+prefix+'_unsorted*') + \
              glob.glob(dirname+prefix+'.sam') + \
              glob.glob(dirname+prefix+'*.sam') + \
              glob.glob(dirname+prefix+'*.00*.bam')
    else:
        fns = glob.glob(dirname+prefix+'*_part*_chunk_'+str(only_chunk)+'*') + \
              glob.glob(dirname+prefix+'*_unsorted*_chunk_'+str(only_chunk)+'*') + \
              glob.glob(dirname+prefix+'_part*_chunk_'+str(only_chunk)+'*') + \
              glob.glob(dirname+prefix+'_unsorted*_chunk_'+str(only_chunk)+'*') + \
              glob.glob(dirname+prefix+'_chunk_'+str(only_chunk)+'*.sam') + \
              glob.glob(dirname+prefix+'_chunk_'+str(only_chunk)+'*.00*.bam')

    if remove_chunks:
        fns = fns + glob.glob(dirname+prefix+'*_chunk_*.bam')
    fns = frozenset(fns)
    for fn in fns:
        os.remove(fn)
        if VERBOSE >= 3:
            print 'File removed:', fn


def gzip_demultiplexed_reads(data_folder, adaID, VERBOSE=0):
    '''Gzip FastQ demultiplexed files'''
    from hivwholeseq.sequencing.filenames import get_read_filenames

    fns = get_read_filenames(data_folder, adaID)
    for fn in fns:
        if not os.path.isfile(fn):
            continue
        sp.call(['gzip', fn])
        if VERBOSE >= 2:
            print 'Gzipped:', fn


def gunzip_demultiplexed_reads(data_folder, adaID, VERBOSE=0):
    '''Gunzip FastQ.gz demultiplexed files'''
    from hivwholeseq.sequencing.filenames import get_read_filenames

    fns = get_read_filenames(data_folder, adaID, gzip=True)
    for fn in fns:
        if not os.path.isfile(fn):
            continue
        sp.call(['gunzip', fn])
        if VERBOSE >= 2:
            print 'Gunzipped:', fn



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--run', required=True,
                        help='MiSeq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    VERBOSE = args.verbose
    submit = args.submit

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = dataset['adapters']

    # Iterate over all adaIDs
    for adaID in adaIDs:

        # Submit to the cluster self if requested
        if submit:
            fork_self(seq_run, adaID, VERBOSE=VERBOSE, threads=threads,
                      reference=refname, summary=summary)
            continue

        remove_premapped_tempfiles(data_folder, adaID, VERBOSE=VERBOSE)

        remove_mapped_tempfiles(data_folder, adaID, VERBOSE=VERBOSE)

        gzip_demultiplexed_reads(data_folder, adaID, VERBOSE=VERBOSE)

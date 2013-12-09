#!/usr/bin/env python
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
import time
import argparse
import pysam
import numpy as np
import subprocess as sp

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.adapter_info import load_adapter_table, foldername_adapter
from hivwholeseq.mapping_utils import stampy_bin, subsrate, bwa_bin, convert_sam_to_bam, \
        convert_bam_to_sam, get_number_reads
from hivwholeseq.filenames import get_consensus_filename, get_mapped_filename,\
        get_read_filenames, get_divided_filename
from hivwholeseq.filter_mapped_reads import match_len_min, trim_bad_cigars
from hivwholeseq.filter_mapped_reads import filter_reads as filter_mapped_reads
from hivwholeseq.fork_cluster import fork_map_to_consensus as fork_self
from hivwholeseq.samples import samples



# Globals
# Stampy parameters
stampy_gapopen = 60	        # Default: 40
stampy_gapextend = 5 	    # Default: 3
stampy_sensitive = True     # Default: False


# Cluster submit
import hivwholeseq
JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
JOBSCRIPT = JOBDIR+'map_to_consensus.py'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
# Different times based on how many threads are used
cluster_time = ['23:59:59', '0:59:59']
vmem = '8G'


# Functions
def get_index_file(data_folder, adaID, fragment, ext=True):
    '''Get the index filename, with or w/o extension'''
    filename = 'consensus_'+fragment
    filename = 'hash/'+filename
    filename = foldername_adapter(adaID)+filename
    if ext:
        filename = filename+'.stidx'
    return data_folder+filename


def get_hash_file(data_folder, adaID, fragment, ext=True):
    '''Get the index filename, with or w/o extension'''
    filename = 'consensus_'+fragment
    filename = 'hash/'+filename
    filename = foldername_adapter(adaID)+filename
    if ext:
        filename = filename+'.sthash'
    return data_folder+filename


def make_output_folders(data_folder, adaID, VERBOSE=0, bwa=False):
    '''Make the output folders if necessary for hash and map'''
    hash_foldername = os.path.dirname(get_hash_file(data_folder, adaID, 'F0'))
    map_foldername = os.path.dirname(get_mapped_filename(data_folder, adaID, 'F0'))
    foldernames = [hash_foldername, map_foldername]

    # Add BWA folder if requested
    if bwa:
        bwa_foldername = os.path.dirname(get_mapped_filename(data_folder, adaID, 'F0',
                                         bwa=True))
        foldernames.append(bwa_foldername)

    # Make the folders
    for dirname in foldernames:
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        if VERBOSE:
            print 'Folder created:', dirname


def make_index_and_hash(data_folder, adaID, fragment, VERBOSE=0):
    '''Make index and hash files for consensus'''
    frag_gen = fragment[:2]

    # NOTE: we can use --overwrite here, because there is no concurrency (every
    # job has its own hash)
    # 1. Make genome index file
    if not os.path.isfile(get_index_file(data_folder, adaID, frag_gen, ext=True)):
        sp.call([stampy_bin,
                 '--species="HIV fragment '+frag_gen+'"',
                 '--overwrite',
                 '-G', get_index_file(data_folder, adaID, frag_gen, ext=False),
                 get_consensus_filename(data_folder, adaID, frag_gen, trim_primers=True),
                 ])
        if VERBOSE:
            print 'Built index: '+adaID+' '+frag_gen
    
    # 2. Build a hash file
    if not os.path.isfile(get_hash_file(data_folder, adaID, frag_gen, ext=True)):
        sp.call([stampy_bin,
                 '--overwrite',
                 '-g', get_index_file(data_folder, adaID, frag_gen, ext=False),
                 '-H', get_hash_file(data_folder, adaID, frag_gen, ext=False),
                 ])
        if VERBOSE:
            print 'Built hash: '+adaID+' '+frag_gen


def make_bwa_hash(data_folder, adaID, fragment, VERBOSE=0):
    '''Make hash tables for BWA'''
    frag_gen = fragment[:2]

    # Make folder if necessary
    dirname =  os.path.dirname(get_hash_file(data_folder, adaID, 'F0'))

    # Call bwa index
    cons_filename = get_consensus_filename(data_folder, adaID, frag_gen)
    if VERBOSE:
        print 'Build BWA index: '+adaID+' '+frag_gen
    sp.call([bwa_bin,
             'index',
             cons_filename,
             ])

    # Move the hashes into subfolder
    from glob import glob
    import shutil
    shutil.copy(cons_filename, dirname)
    hash_files = glob(cons_filename+'.*')
    for hash_file in hash_files:
        shutil.copy(hash_file, dirname)
        os.remove(hash_file)


def map_bwa(data_folder, adaID, fragment, VERBOSE=0):
    '''Map using BWA'''
    frag_gen = fragment[:2]
    index_prefix = get_hash_file(data_folder, adaID, frag_gen, ext=False)+'.fasta'

    if VERBOSE:
        print 'Map via BWA: '+adaID+' '+frag_gen
    output_filename = get_mapped_filename(data_folder, adaID, frag_gen,
                                          type='sam', bwa=True)

    # Map
    with open(output_filename, 'w') as f:
        call_list = [bwa_bin,
                     'mem',
                     index_prefix] +\
                    get_read_filenames(data_folder, adaID, premapped=True,
                                       fragment=fragment)
        if VERBOSE >=2:
            print ' '.join(call_list)
        sp.call(call_list, stdout=f)


def map_stampy(data_folder, adaID, fragment, VERBOSE=0, bwa=False, threads=1, n_pairs=-1):
    '''Map using stampy, either directly or after BWA'''
    # Get input filenames
    frag_gen = fragment[:2]
    if bwa:
        if VERBOSE:
            print 'Map via stampy after BWA: '+adaID+' '+frag_gen

        bwa_filename = get_mapped_filename(data_folder, adaID, frag_gen,
                                           type='bam', bwa=True)
        if not os.path.isfile(bwa_filename):
            if VERBOSE >= 2:
                print 'Converting SAM to BAM'
            convert_sam_to_bam(bwa_filename)
        input_filename = bwa_filename

    else:
        if VERBOSE:
            print 'Map via stampy: '+adaID+' '+frag_gen
        input_filename = get_divided_filename(data_folder, adaID, fragment, type='bam')

    # parallelize if requested
    if threads == 1:

        # Get output filename
        output_filename = get_mapped_filename(data_folder, adaID, frag_gen, type='sam')

        # Map
        call_list = [stampy_bin,
                     '-g', get_index_file(data_folder, adaID, frag_gen, ext=False),
                     '-h', get_hash_file(data_folder, adaID, frag_gen, ext=False), 
                     '-o', output_filename,
                     '--substitutionrate='+subsrate,
                     '--gapopen', stampy_gapopen,
                     '--gapextend', stampy_gapextend]
        if bwa:
            call_list.append('--bamkeepgoodreads')
        if stampy_sensitive:
            call_list.append('--sensitive')

        # Take only a (random) subsample: stampy uses the fraction of reads
        # intead of the number
        if n_pairs > 0:
            n_pairs_tot = get_number_reads(input_filename, 'bam') / 2
            frac_pairs = 1.0 * n_pairs / n_pairs_tot
            random_seed = np.random.randint(1e5)
            call_list.append('-s', frac_pairs + random_seed)

        call_list = call_list + ['-M', input_filename]
        call_list = map(str, call_list)
        if VERBOSE >=2:
            print ' '.join(call_list)
        sp.call(call_list)

    else:

        # Submit map script
        jobs_done = np.zeros(threads, bool)
        job_IDs = np.zeros(threads, 'S30')
        for j in xrange(threads):
    
            # Get output filename
            output_filename =  get_mapped_filename(data_folder, adaID, frag_gen,
                                               type='sam', part=(j+1))
            # Map
            call_list = ['qsub','-cwd',
                         '-b', 'y',
                         '-S', '/bin/bash',
                         '-o', JOBLOGOUT,
                         '-e', JOBLOGERR,
                         '-N', 'm '+adaID+frag_gen+' p'+str(j+1),
                         '-l', 'h_rt='+cluster_time[threads >= 10],
                         '-l', 'h_vmem='+vmem,
                         stampy_bin,
                         '-g', get_index_file(data_folder, adaID, frag_gen, ext=False),
                         '-h', get_hash_file(data_folder, adaID, frag_gen, ext=False), 
                         '-o', output_filename,
                         '--processpart='+str(j+1)+'/'+str(threads),
                         '--substitutionrate='+subsrate,
                         '--gapopen', stampy_gapopen,
                         '--gapextend', stampy_gapextend]
            if bwa:
                call_list.append('--bamkeepgoodreads')
            if stampy_sensitive:
                call_list.append('--sensitive')
            call_list = call_list + ['-M', input_filename]
            call_list = map(str, call_list)
            if VERBOSE >= 2:
                print ' '.join(call_list)
            job_ID = sp.check_output(call_list)
            job_ID = job_ID.split()[2]
            job_IDs[j] = job_ID

        # Monitor output
        output_file_parts = [get_mapped_filename(data_folder, adaID, frag_gen,
                                                 type='bam', part=(j+1))
                              for j in xrange(threads)]
        time_wait = 10 # secs
        while not jobs_done.all():

            # Sleep some time
            time.sleep(time_wait)

            # Get the output of qstat to check the status of jobs
            qstat_output = sp.check_output(['qstat'])
            qstat_output = qstat_output.split('\n')[:-1] # The last is an empty line
            if len(qstat_output) < 3:
                jobs_done[:] = True
                break
            else:
                qstat_output = [line.split()[0] for line in qstat_output[2:]]

            time_wait = 10 # secs
            for j in xrange(threads):
                if jobs_done[j]:
                    continue

                if job_IDs[j] not in qstat_output:
                    # Convert to BAM for merging
                    if VERBOSE >= 1:
                        print 'Convert mapped reads to BAM for merging: adaID '+\
                               adaID+', part '+str(j+1)+ ' of '+ \
                               str(threads)
                    convert_sam_to_bam(output_file_parts[j])
                    # We do not need to wait if we did the conversion (it takes
                    # longer than some secs)
                    time_wait = 0
                    jobs_done[j] = True

        # Concatenate output files
        output_filename = get_mapped_filename(data_folder, adaID, frag_gen, type='bam',
                                              unsorted=True)
        if VERBOSE >= 1:
            print 'Concatenate premapped reads: adaID '+adaID
        pysam.cat('-o', output_filename, *output_file_parts)

        # Sort the file by read names (to ensure the pair_generator)
        output_filename_sorted = get_mapped_filename(data_folder, adaID, frag_gen,
                                                     type='bam',
                                                     unsorted=False)
        # NOTE: we exclude the extension and the option -f because of a bug in samtools
        if VERBOSE >= 1:
            print 'Sort mapped reads: adaID '+adaID
        pysam.sort('-n', output_filename, output_filename_sorted[:-4])

        # Reheader the file without BAM -> SAM -> BAM
        if VERBOSE >= 1:
            print 'Reheader mapped reads: adaID '+adaID
        header_filename = get_mapped_filename(data_folder, adaID, frag_gen,
                                              type='sam', part=1)
        pysam.reheader(header_filename, output_filename_sorted)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('-n', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--bwa', action='store_true',
                        help='Use BWA for premapping?')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use for mapping')
    parser.add_argument('--filter', action='store_true',
                        help='Filter reads immediately after mapping')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    submit = args.submit
    bwa = args.bwa
    threads = args.threads
    n_pairs = args.n
    filter_reads = args.filter

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = dataset['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Iterate over all requested samples
    for adaID in adaIDs:

        # Make output folders if necessary
        make_output_folders(data_folder, adaID, VERBOSE=VERBOSE)

        # If the script is called with no fragment, iterate over all
        samplename = dataset['samples'][dataset['adapters'].index(adaID)]
        if not fragments:
            fragments_sample = samples[samplename]['fragments']
        else:
            from re import findall
            fragments_all = samples[samplename]['fragments']
            fragments_sample = []
            for fragment in fragments:
                frs = filter(lambda x: fragment in x, fragments_all)
                if len(frs):
                    fragments_sample.append(frs[0])

        if VERBOSE >= 3:
            print 'adaID '+adaID+': fragments '+' '.join(fragments_sample)

        # Iterate over fragments
        for fragment in fragments_sample:

            # Submit to the cluster self if requested
            if submit:
                fork_self(seq_run, adaID, fragment,
                          VERBOSE=VERBOSE, bwa=bwa,
                          threads=threads, n_pairs=n_pairs,
                          filter_reads=filter_reads)
                continue

            if bwa:
                # Make BWA hashes
                make_bwa_hash(data_folder, adaID, fragment,
                              VERBOSE=VERBOSE)
    
                # Map via BWA first
                map_bwa(data_folder, adaID, fragment,
                        VERBOSE=VERBOSE)
    
            # Make stampy hashes
            make_index_and_hash(data_folder, adaID, fragment, VERBOSE=VERBOSE)
            
            # Map via stampy afterwards
            map_stampy(data_folder, adaID, fragment,
                       VERBOSE=VERBOSE, bwa=bwa,
                       threads=threads, n_pairs=n_pairs)

            # Filter reads after mapping if requested
            if filter_reads:
                filter_mapped_reads(data_folder, adaID, fragment, VERBOSE=VERBOSE)

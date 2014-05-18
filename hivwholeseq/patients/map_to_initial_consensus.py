#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       25/10/13
content:    Map reads to the initial consensus of the patient (in order to have
            the whole patient data aligned to one reference).
'''
# Modules
import os
import time
import argparse
import subprocess as sp
import numpy as np
import pysam
import warnings

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.filenames import get_divided_filename
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.mapping_utils import stampy_bin, subsrate, \
        convert_sam_to_bam, convert_bam_to_sam, get_number_reads
from hivwholeseq.samples import samples as samples_seq
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.patients.filenames import get_initial_index_filename, \
        get_initial_hash_filename, get_initial_consensus_filename, \
        get_mapped_to_initial_filename, get_map_initial_summary_filename, \
        get_filter_mapped_init_summary_filename
from hivwholeseq.fork_cluster import fork_map_to_initial_consensus as fork_self
from hivwholeseq.clean_temp_files import remove_mapped_init_tempfiles
from hivwholeseq.patients.filter_mapped_reads import filter_mapped_reads



# Globals
stampy_gapopen = 60	        # Default: 40
stampy_gapextend = 5 	    # Default: 3
stampy_sensitive = True     # Default: False


# Functions
def make_output_folders(pname, sample, VERBOSE=0):
    '''Make the output folders if necessary for hash and map'''
    hash_foldername = os.path.dirname(get_initial_hash_filename(pname, 'F0'))
    map_foldername = os.path.dirname(get_mapped_to_initial_filename(pname, sample, 'F0'))
    foldernames = [hash_foldername, map_foldername]

    for dirname in foldernames:
        mkdirs(dirname)
        if VERBOSE:
            print 'Folder created:', dirname


def make_index_and_hash(pname, fragment, VERBOSE=0):
    '''Make index and hash files for consensus'''
    # 1. Make genome index file
    stdout = sp.check_output([stampy_bin,
                              '--overwrite',
                              '--species="HIV fragment '+fragment+'"',
                              '-G', get_initial_index_filename(pname, fragment, ext=False),
                              get_initial_consensus_filename(pname, fragment),
                              ],
                              stderr=sp.STDOUT)
    if VERBOSE:
        print 'Built index: '+pname+' '+fragment
    
    # 2. Build a hash file
    stdout = sp.check_output([stampy_bin,
                              '--overwrite',
                              '-g', get_initial_index_filename(pname, fragment, ext=False),
                              '-H', get_initial_hash_filename(pname, fragment, ext=False),
                              ],
                              stderr=sp.STDOUT)
    if VERBOSE:
        print 'Built hash: '+pname+' '+fragment


def map_stampy_singlethread(patient, samplename, fragment, VERBOSE=0, n_pairs=-1,
                            summary=True, only_chunk=None):
    '''Map using stampy, single thread (no cluster queueing race conditions)'''
    pname = patient.id
    sample = patient.sample_table.loc[samplename]
    seq_run = sample['run']
    data_folder = MiSeq_runs[seq_run]['folder']
    adaID = sample['adaID']

    if VERBOSE:
        print 'Map via stampy (single thread): '+pname+' '+samplename+' '+fragment

    if summary:
        summary_filename = get_map_initial_summary_filename(pname, samplename, fragment)

    # Specific fragment (e.g. F5 --> F5bi)
    frag_spec = filter(lambda x: fragment in x, sample['fragments'])
    if not len(frag_spec):
        raise ValueError(str(patient)+', '+samplename+': fragment '+fragment+' not found.')
    frag_spec = frag_spec[0]

    input_filename = get_divided_filename(data_folder, adaID, frag_spec, type='bam', chunk=only_chunk)

    # Extract subsample of reads if requested
    if n_pairs > 0:
        from hivwholeseq.mapping_utils import extract_mapped_reads_subsample
        input_filename_sub = get_mapped_to_initial_filename(pname, samplename,
                                                            fragment,
                                                            type='bam')[:-4]+\
                '_unmapped.bam'
        n_written = extract_mapped_reads_subsample(input_filename,
                                                   input_filename_sub,
                                                   n_pairs, VERBOSE=VERBOSE)

    # Get output filename
    output_filename = get_mapped_to_initial_filename(pname, samplename, fragment,
                                                     type='sam', only_chunk=only_chunk)

    # Map
    call_list = [stampy_bin,
                 '-g', get_initial_index_filename(pname, fragment, ext=False),
                 '-h', get_initial_hash_filename(pname, fragment, ext=False),
                 '-o', output_filename,
                 '--overwrite',
                 '--substitutionrate='+subsrate,
                 '--gapopen', stampy_gapopen,
                 '--gapextend', stampy_gapextend]
    if stampy_sensitive:
        call_list.append('--sensitive')

    if n_pairs > 0:
        call_list = call_list + ['-M', input_filename_sub]
    else:
        call_list = call_list + ['-M', input_filename]
    call_list = map(str, call_list)
    if VERBOSE >=2:
        print ' '.join(call_list)
    sp.call(call_list)

    output_filename_bam = get_mapped_to_initial_filename(pname, samplename,
                                                         fragment,
                                                         type='bam',
                                                         only_chunk=only_chunk)
    convert_sam_to_bam(output_filename_bam)

    if summary:
        with open(summary_filename, 'a') as f:
            f.write('Stampy mapped (single thread).\n')

    if only_chunk is None:
        if VERBOSE >= 1:
            print 'Remove temporary files: sample '+samplename
        remove_mapped_init_tempfiles(pname, samplename, fragment, VERBOSE=VERBOSE, only_chunk=only_chunk)

    if summary:
        with open(summary_filename, 'a') as f:
            f.write('Temp mapping files removed.\n')
            f.write('\n')

    if n_pairs > 0:
        os.remove(input_filename_sub)


def map_stampy_multithread(patient, samplename, fragment, VERBOSE=0, threads=2, summary=True):
    '''Map using stampy, multithread (via cluster requests, queueing race conditions possible)'''
    import hivwholeseq
    JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
    JOBLOGOUT = JOBDIR+'logout/'
    JOBLOGERR = JOBDIR+'logerr/'
    cluster_time = ['23:59:59', '0:59:59']
    vmem = '8G'

    pname = patient.id
    sample = patient.sample_table.loc[samplename]
    seq_run = sample['run']
    data_folder = MiSeq_runs[seq_run]['folder']
    adaID = sample['adaID']

    if VERBOSE:
        print 'Map via stampy: '+pname+' '+samplename+' '+fragment

    if summary:
        summary_filename = get_map_initial_summary_filename(pname, samplename, fragment)

    # Specific fragment (e.g. F5 --> F5bi)
    frag_spec = filter(lambda x: fragment in x, sample['fragments'])
    if not len(frag_spec):
        raise ValueError(str(patient)+', '+samplename+': fragment '+fragment+' not found.')
    frag_spec = frag_spec[0]

    input_filename = get_divided_filename(data_folder, adaID, frag_spec, type='bam')

    # Submit map scripts in parallel to the cluster
    jobs_done = np.zeros(threads, bool)
    job_IDs = np.zeros(threads, 'S30')
    for j in xrange(threads):
    
        output_filename = get_mapped_to_initial_filename(pname, samplename,
                                                         fragment,
                                                         type='sam', part=(j+1))
        # Map
        call_list = ['qsub','-cwd',
                     '-b', 'y',
                     '-S', '/bin/bash',
                     '-o', JOBLOGOUT,
                     '-e', JOBLOGERR,
                     '-N', 'm '+samplename+fragment+' p'+str(j+1),
                     '-l', 'h_rt='+cluster_time[threads >= 10],
                     '-l', 'h_vmem='+vmem,
                     stampy_bin,
                     '--overwrite',
                     '-g', get_initial_index_filename(pname, fragment, ext=False),
                     '-h', get_initial_hash_filename(pname, fragment, ext=False),
                     '-o', output_filename,
                     '--processpart='+str(j+1)+'/'+str(threads),
                     '--substitutionrate='+subsrate,
                     '--gapopen', stampy_gapopen,
                     '--gapextend', stampy_gapextend]
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
    output_file_parts = [get_mapped_to_initial_filename(pname, samplename,
                                                        fragment,
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
                    print 'Convert mapped reads to BAM for merging: sample '+\
                           samplename+', part '+str(j+1)+ ' of '+ \
                           str(threads)
                convert_sam_to_bam(output_file_parts[j])
                # We do not need to wait if we did the conversion (it takes
                # longer than some secs)
                time_wait = 0
                jobs_done[j] = True

    if summary:
        with open(summary_filename, 'a') as f:
            f.write('Stampy mapped ('+str(threads)+' threads).\n')

    # Concatenate output files
    output_filename = get_mapped_to_initial_filename(pname, samplename,
                                                     fragment,
                                                     type='bam', unsorted=True)
    if VERBOSE >= 1:
        print 'Concatenate premapped reads: sample '+samplename
    pysam.cat('-o', output_filename, *output_file_parts)
    if summary:
        with open(summary_filename, 'a') as f:
            f.write('BAM files concatenated (unsorted).\n')

    # Sort the file by read names (to ensure the pair_generator)
    output_filename_sorted = get_mapped_to_initial_filename(pname, samplename,
                                                            fragment,
                                                            type='bam')
    # NOTE: we exclude the extension and the option -f because of a bug in samtools
    if VERBOSE >= 1:
        print 'Sort mapped reads: sample '+samplename
    pysam.sort('-n', output_filename, output_filename_sorted[:-4])
    if summary:
        with open(summary_filename, 'a') as f:
            f.write('Joint BAM file sorted.\n')

    # Reheader the file without BAM -> SAM -> BAM
    if VERBOSE >= 1:
        print 'Reheader mapped reads: sample '+samplename
    header_filename = get_mapped_to_initial_filename(pname, samplename,
                                                     fragment,
                                                     type='sam', part=1)
    pysam.reheader(header_filename, output_filename_sorted)
    if summary:
        with open(summary_filename, 'a') as f:
            f.write('Joint BAM file reheaded.\n')

    if VERBOSE >= 1:
        print 'Remove temporary files: sample '+samplename
    remove_mapped_init_tempfiles(pname, samplename, fragment, VERBOSE=VERBOSE)
    if summary:
        with open(summary_filename, 'a') as f:
            f.write('Temp mapping files removed.\n')
            f.write('\n')


def map_stampy(patient, samplename, fragment, VERBOSE=0, threads=1, n_pairs=-1,
               summary=True, only_chunk=None):
    '''Map using stampy'''
    if threads == 1:
        map_stampy_singlethread(patient, samplename, fragment,
                                VERBOSE=VERBOSE, n_pairs=n_pairs,
                                summary=summary, only_chunk=only_chunk)

    else:
        if (n_pairs != -1):
            raise ValueError('Cannot multithread a subset - why would you do this?')
        elif only_chunk:
            raise ValueError('Cannot multithread a chunk - why would you do this?')

        map_stampy_multithread(patient, samplename, fragment, threads=threads,
                               VERBOSE=VERBOSE, summary=summary)


def get_number_chunks(pname, samplename, fragment, VERBOSE=0):
    '''Find the number of chunks the reads have been chopped down into'''
    sample = patient.sample_table.loc[samplename]
    seq_run = sample['run']
    data_folder = MiSeq_runs[seq_run]['folder']
    adaID = sample['adaID']
    # Has this fragment been sequenced in that sample?
    frag_spec = filter(lambda x: fragment in x,
                       sample['fragments'])
    if not len(frag_spec):
        warnings.warn(pname+', '+samplename+': fragment '+fragment+' not found.')
        return
    else:
        frag_spec = frag_spec[0]

    n_chunks = 0
    input_filename = get_divided_filename(data_folder, adaID, frag_spec, type='bam', chunk=n_chunks+1)
    while(os.path.isfile(input_filename)):
        n_chunks += 1
        input_filename = get_divided_filename(data_folder, adaID, frag_spec, type='bam', chunk=n_chunks+1)

    return n_chunks



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Map to initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--samples', nargs='+',
                        help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use for mapping')
    parser.add_argument('--filter', action='store_true',
                        help='Filter reads immediately after mapping')
    parser.add_argument('--skiphash', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')
    parser.add_argument('--chunks', type=int, nargs='+', default=[None],
                        help='Only map some chunks (cluster optimization): 0 for automatic detection')

    args = parser.parse_args()
    pname = args.patient
    samples = args.samples
    fragments = args.fragments
    submit = args.submit
    VERBOSE = args.verbose
    threads = args.threads
    n_pairs = args.maxreads
    filter_reads = args.filter
    skip_hash = args.skiphash
    summary = args.summary
    only_chunks = args.chunks

    patient = get_patient(pname)

    # If no samples are mentioned, use all sequenced ones
    if not samples:
        samples = patient.samples
    if VERBOSE >= 3:
        print 'samples', samples

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    sample_init = patient.initial_sample

    for sample in samples:
        make_output_folders(pname, sample, VERBOSE=VERBOSE)

    for fragment in fragments:
        if not skip_hash:
            make_index_and_hash(pname, fragment, VERBOSE=VERBOSE)
            
        for samplename in samples:

            if only_chunks == [0]:
                only_chunks_sample = range(1, get_number_chunks(pname, samplename, fragment, VERBOSE=VERBOSE) + 1)
            else:
                only_chunks_sample = only_chunks

            for only_chunk in only_chunks_sample:
    
                if submit:
                    fork_self(pname, samplename, fragment,
                              VERBOSE=VERBOSE, threads=threads,
                              n_pairs=n_pairs,
                              filter_reads=filter_reads,
                              summary=summary,
                              only_chunks=[only_chunk])
                    continue
    
                if summary:
                    sfn = get_map_initial_summary_filename(pname, samplename, fragment)
                    with open(sfn, 'w') as f:
                        f.write('Call: python map_to_initial_consensus.py'+\
                                ' --patient '+pname+\
                                ' --samples '+samplename+\
                                ' --fragments '+fragment+\
                                ' --threads '+str(threads)+\
                                ' --verbose '+str(VERBOSE))
                        if n_pairs != -1:
                            f.write(' --maxreads '+str(n_pairs))
                        if filter_reads:
                            f.write(' --filter')
                        f.write('\n')
    
    
                map_stampy(patient, samplename, fragment,
                           VERBOSE=VERBOSE, threads=threads, n_pairs=n_pairs,
                           summary=summary, only_chunk=only_chunk)
    
                if filter_reads:
                    if only_chunk is not None:
                        print 'Post-map filtering only AFTER the chunks are pasted together!'
                        continue

                    if summary:
                        sfn = get_filter_mapped_init_summary_filename(pname, samplename, fragment)
                        with open(sfn, 'w') as f:
                            f.write('Call: python filter_mapped_reads.py'+\
                                    ' --patient '+pname+\
                                    ' --samples '+samplename+\
                                    ' --fragments '+fragment+\
                                    ' --verbose '+str(VERBOSE))
                            if n_pairs != -1:
                                f.write(' --maxreads '+str(n_pairs))
                            f.write('\n')


                    filter_mapped_reads(pname, samplename, fragment,
                                        VERBOSE=VERBOSE, maxreads=n_pairs,
                                        summary=summary)




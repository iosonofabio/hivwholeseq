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
        get_mapped_to_initial_filename, get_map_initial_summary_filename
from hivwholeseq.fork_cluster import fork_map_to_initial_consensus as fork_self



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
    sp.call([stampy_bin,
             '--overwrite',
             '--species="HIV fragment '+fragment+'"',
             '-G', get_initial_index_filename(pname, fragment, ext=False),
             get_initial_consensus_filename(pname, fragment),
             ])
    if VERBOSE:
        print 'Built index: '+pname+' '+fragment
    
    # 2. Build a hash file
    sp.call([stampy_bin,
             '--overwrite',
             '-g', get_initial_index_filename(pname, fragment, ext=False),
             '-H', get_initial_hash_filename(pname, fragment, ext=False),
             ])
    if VERBOSE:
        print 'Built hash: '+pname+' '+fragment


def map_stampy(patient, samplename, fragment, VERBOSE=0, threads=1, n_pairs=-1,
               summary=True):
    '''Map using stampy'''
    # Get input filenames
    pname = patient.id
    sample = patient.sample_table.loc[samplename]
    seq_run = sample['run']
    data_folder = MiSeq_runs[seq_run]['folder']
    adaID = sample['adaID']

    if VERBOSE:
        print 'Map via stampy: '+pname+' '+samplename+' '+fragment

    if summary:
        summary_filename = get_map_initial_summary_filename(pname, samplename,
                                                            fragment)

    # Has this fragment been sequenced in that sample?
    frag_spec = filter(lambda x: fragment in x,
                       sample['fragments'])
    if not len(frag_spec):
        warnings.warn(str(patient)+', '+samplename+': fragment '+fragment+' not found.')
        return
    else:
        frag_spec = frag_spec[0]

    input_filename = get_divided_filename(data_folder, adaID, frag_spec, type='bam')

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


    # parallelize if requested
    if threads == 1:

        # Get output filename
        output_filename = get_mapped_to_initial_filename(pname, samplename, fragment,
                                                         type='sam')

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
                                                             type='bam')
        convert_sam_to_bam(output_filename_bam)

        if summary:
            with open(summary_filename, 'a') as f:
                f.write('Stampy mapped (single thread).\n')


    else:

        # Submit map script
        jobs_done = np.zeros(threads, bool)
        job_IDs = np.zeros(threads, 'S30')
        for j in xrange(threads):
    
            # Get output filename
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

            if n_pairs > 0:
                call_list = call_list + ['-M', input_filename_sub]
            else:
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
                               sample+', part '+str(j+1)+ ' of '+ \
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
            print 'Concatenate premapped reads: sample '+sample
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
            print 'Sort mapped reads: sample '+sample
        pysam.sort('-n', output_filename, output_filename_sorted[:-4])
        if summary:
            with open(summary_filename, 'a') as f:
                f.write('Joint BAM file sorted.\n')

        # Reheader the file without BAM -> SAM -> BAM
        if VERBOSE >= 1:
            print 'Reheader mapped reads: sample '+sample
        header_filename = get_mapped_to_initial_filename(pname, samplename,
                                                         fragment,
                                                         type='sam', part=1)
        pysam.reheader(header_filename, output_filename_sorted)
        if summary:
            with open(summary_filename, 'a') as f:
                f.write('Joint BAM file reheaded.\n')

    if n_pairs > 0:
        os.remove(input_filename_sub)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Map to initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--samples', nargs='*',
                        help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--fragments', nargs='*',
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

    # Get the patient
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

            if submit:
                fork_self(pname, samplename, fragment,
                          VERBOSE=VERBOSE, threads=threads,
                          n_pairs=n_pairs,
                          filter_reads=filter_reads,
                          summary=summary)
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
                        f.write(' -n '+str(n_pairs))
                    if filter_reads:
                        f.write(' --filter')
                    f.write('\n')


            map_stampy(patient, samplename, fragment,
                       VERBOSE=VERBOSE, threads=threads, n_pairs=n_pairs,
                       summary=summary)

            if filter_reads:
                warnings.warn('Filtering not implemented yet')
                #TODO
                #filter_mapped_reads(pname, sample, fragment, VERBOSE=VERBOSE)

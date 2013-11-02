#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/10/13
content:    Premap demultiplex reads to HXB2 (or other reference) 'as is',
            before any filtering. This is just a rough mapping, just to get
            the classification into fragments right.

            Note: stampy can do some kind of multithreading (see below).
'''
# Modules
import os
import time
import subprocess as sp
import argparse
import numpy as np
import pysam

from mapping.datasets import MiSeq_runs
from mapping.filenames import get_HXB2_entire, get_HXB2_index_file,\
        get_custom_reference_filename, get_custom_index_filename_fun, \
        get_custom_hash_filename_fun, \
        get_HXB2_hash_file, get_read_filenames, get_premapped_file
from mapping.mapping_utils import stampy_bin, subsrate, \
        convert_sam_to_bam, convert_bam_to_sam
from mapping.reference import load_HXB2, load_custom_reference



# Globals
# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'premap_to_reference.py'
# It is hard to tell whether 1h is sufficient, because the final sorting takes
# quite some time. So for now give up and require 2h.
cluster_time = ['23:59:59', '1:59:59']
vmem = '8G'



# Functions
def fork_self(miseq_run, adaID, VERBOSE=0, bwa=False, threads=1, report=0,
              reference='HXB2'):
    '''Fork self for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+'{:02d}'.format(adaID)

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'premap '+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time[threads >= 30],
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', miseq_run,
                 '--adaIDs', adaID,
                 '--verbose', VERBOSE,
                 '--threads', threads,
                 '--reference', reference,
                ]
    if report:
        qsub_list.append('--report')
    qsub_list = map(str, qsub_list)
    if VERBOSE >= 2:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def make_output_folders(data_folder, adaID, VERBOSE=0):
    '''Make output folders'''
    output_filename = get_premapped_file(data_folder, adaID)
    dirname = os.path.dirname(output_filename)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
        if VERBOSE:
            print 'Folder created:', dirname


def make_index_and_hash(data_folder, reference='HXB2', VERBOSE=0):
    '''Make index and hash files for reference or consensus'''
    if VERBOSE:
        print 'Making index and hash files'

    # Check whether or not the reference is standard
    # NOTE: this is needed for samples which are far away from HXB2
    if reference == 'HXB2':
        ref_filename = get_HXB2_entire(cropped=True)
        index_filename_fun = get_HXB2_index_file
        hash_filename_fun = get_HXB2_hash_file
    else:
        ref_filename = get_custom_reference_filename(reference)
        index_filename_fun = get_custom_index_filename_fun(reference)
        hash_filename_fun = get_custom_hash_filename_fun(reference)

    # Make folder if necessary
    dirname = os.path.dirname(hash_filename_fun())
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
        if VERBOSE:
            print 'Folder created:', dirname

    # 1. Make genome index file for reference
    if not os.path.isfile(index_filename_fun(ext=True)):
        sp.call([stampy_bin,
                 '--species="HIV"',
                 '-G', index_filename_fun(ext=False),
                 ref_filename,
                 ])
    
    # 2. Build a hash file for reference
    if not os.path.isfile(hash_filename_fun(ext=True)):
        sp.call([stampy_bin,
                 '-g', index_filename_fun(ext=False),
                 '-H', hash_filename_fun(ext=False),
                 ])


def premap_stampy(data_folder, adaID, reference='HXB2', VERBOSE=0, threads=1):
    '''Call stampy for actual mapping'''
    if VERBOSE:
        print 'Map: adaID '+'{:02d}'.format(adaID)

    # Check whether or not the reference is standard
    # NOTE: this is needed for samples which are far away from HXB2
    if reference == 'HXB2':
        ref_filename = get_HXB2_entire(cropped=True)
        index_filename_fun = get_HXB2_index_file
        hash_filename_fun = get_HXB2_hash_file
    else:
        ref_filename = get_custom_reference_filename(reference)
        index_filename_fun = get_custom_index_filename_fun(reference)
        hash_filename_fun = get_custom_hash_filename_fun(reference)

    # Get input filenames (raw reads)
    input_filenames = get_read_filenames(data_folder, adaID, filtered=False)

    # parallelize if requested
    if threads == 1:

        # Get output filename
        output_filename =  get_premapped_file(data_folder, adaID, type='sam')
        # Map
        call_list = [stampy_bin,
                     '-g', index_filename_fun(ext=False),
                     '-h', hash_filename_fun(ext=False), 
                     '-o', output_filename,
                     '--substitutionrate='+subsrate,
                     '-M'] + input_filenames
        call_list = map(str, call_list)
        if VERBOSE >= 2:
            print ' '.join(call_list)
        sp.call(call_list)

        # Convert to compressed BAM
        convert_sam_to_bam(get_premapped_file(data_folder, adaID, type='bam'))

    # Multithreading works as follows: call qsub + stampy, monitor the process
    # IDs with qstat at regular intervals, and finally merge results with pysam
    else:

        # Submit map script
        jobs_done = np.zeros(threads, bool)
        job_IDs = np.zeros(threads, 'S30')
        for j in xrange(threads):
    
            # Get output filename
            output_filename =  get_premapped_file(data_folder, adaID, type='sam',
                                                  part=(j+1))
            # Map
            call_list = ['qsub','-cwd',
                         '-b', 'y',
                         '-S', '/bin/bash',
                         '-o', JOBLOGOUT,
                         '-e', JOBLOGERR,
                         '-N', 'pre '+'{:02d}'.format(adaID)+' p'+str(j+1),
                         '-l', 'h_rt='+cluster_time[threads >= 30],
                         '-l', 'h_vmem='+vmem,
                         stampy_bin,
                         '-g', index_filename_fun(ext=False),
                         '-h', hash_filename_fun(ext=False), 
                         '-o', output_filename,
                         '--processpart='+str(j+1)+'/'+str(threads),
                         '--substitutionrate='+subsrate,
                         '-M'] + input_filenames
            call_list = map(str, call_list)
            if VERBOSE >= 2:
                print ' '.join(call_list)
            job_ID = sp.check_output(call_list)
            job_ID = job_ID.split()[2]
            job_IDs[j] = job_ID

        # Monitor output
        output_file_parts = [get_premapped_file(data_folder, adaID, type='bam',
                                                part=(j+1)) for j in xrange(threads)]
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
                        print 'Convert premapped reads to BAM for merging: adaID '+\
                               '{:02d}'.format(adaID)+', part '+str(j+1)+ ' of '+ \
                               str(threads)
                    convert_sam_to_bam(output_file_parts[j])
                    # We do not need to wait if we did the conversion (it takes
                    # longer than some secs)
                    time_wait = 0
                    jobs_done[j] = True

        # Concatenate output files
        output_filename = get_premapped_file(data_folder, adaID, type='bam',
                                             unsorted=True)
        if VERBOSE >= 1:
            print 'Concatenate premapped reads: adaID '+'{:02d}'.format(adaID)
        pysam.cat('-o', output_filename, *output_file_parts)

        # Sort the file by read names (to ensure the pair_generator)
        output_filename_sorted = get_premapped_file(data_folder, adaID, type='bam',
                                                    unsorted=False)
        # NOTE: we exclude the extension and the option -f because of a bug in samtools
        if VERBOSE >= 1:
            print 'Sort premapped reads: adaID '+'{:02d}'.format(adaID)
        pysam.sort('-n', output_filename, output_filename_sorted[:-4])

        # Reheader the file without BAM -> SAM -> BAM
        if VERBOSE >= 1:
            print 'Reheader premapped reads: adaID '+'{:02d}'.format(adaID)
        header_filename = get_premapped_file(data_folder, adaID, type='sam', part=1)
        pysam.reheader(header_filename, output_filename_sorted)


def report_insert_size_distribution(data_folder, adaID,  VERBOSE=0):
    '''Produce figures of the insert size distribution'''
    from mapping.insert_size_distribution import get_insert_size_distribution, \
            plot_cumulative_histogram, plot_histogram

    bins = np.linspace(0, 1000, 100)
    isz, h = get_insert_size_distribution(data_folder, adaID, 'premapped',
                                          bins=bins, maxreads=10000,
                                          VERBOSE=VERBOSE)
    plot_cumulative_histogram(miseq_run, adaID, 'premapped', isz, savefig=True,
                              lw=2, c='b')
    plot_histogram(miseq_run, adaID, 'premapped', h, savefig=True,
                   lw=2, color='b')


def report_coverage(data_folder, adaID, reference='HXB2', VERBOSE=0):
    '''Produce a report on rough coverage on reference (ignore inserts)'''
    # Get the reference sequence

    # Check whether or not the reference is standard
    # NOTE: this is needed for samples which are far away from HXB2
    if reference == 'HXB2':
        refseq = load_HXB2(cropped=True)
    else:
        refseq = load_custom_reference(reference)

    # Prepare data structures
    coverage = np.zeros(len(refseq), int)
    
    # Parse the BAM file
    output_filename = get_premapped_file(data_folder, adaID, type='bam')
    with pysam.Samfile(output_filename, 'rb') as bamfile:
        for read in bamfile:
            if read.is_unmapped or (not read.is_proper_pair):
                continue

            # Proceed along CIGARs
            ref_pos = read.pos
            for (bt, bl) in read.cigar:
                if bt not in (0, 2):
                    continue
                # Treat deletions as 'covered'
                coverage[ref_pos: ref_pos+bl] += 1
                ref_pos += bl

    # Save results
    from mapping.filenames import get_coverage_figure_filename
    import matplotlib.pyplot as plt
    plt.plot(np.arange(len(refseq)), coverage + 1, lw=2, c='b')
    plt.xlabel('Position in cropped HXB2')
    plt.ylabel('Coverage')
    plt.yscale('log')
    plt.title('adaID '+'{:02d}'.format(adaID)+', premapped to HXB2',
              fontsize=18)

    plt.tight_layout()
    plt.savefig(get_coverage_figure_filename(data_folder, adaID, 'premapped'))



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use for mapping')
    parser.add_argument('--report', action='store_true',
                        help='Perform quality checks and save into a report')
    parser.add_argument('--reference', default='HXB2',
                        help='Use alternative reference, e.g. chimeras (the file must exist)')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    VERBOSE = args.verbose
    submit = args.submit
    threads = args.threads
    report = args.report
    refname = args.reference

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = MiSeq_runs[miseq_run]['adapters']

    # Make index and hash files for reference (if not present already), using a
    # cropped reference to reduce mapping ambiguities at LTRs (F1 vs F6)
    make_index_and_hash(data_folder, reference=refname, VERBOSE=VERBOSE)

    # Iterate over all adaIDs
    for adaID in adaIDs:

        # Submit to the cluster self if requested
        if submit:
            fork_self(miseq_run, adaID, VERBOSE=VERBOSE, threads=threads,
                      report=report, reference=refname)
            continue

        # Make output folders if necessary
        make_output_folders(data_folder, adaID, VERBOSE=VERBOSE)

        # Map roughly to reference
        premap_stampy(data_folder, adaID, VERBOSE=VERBOSE, reference=refname,
                      threads=threads)

        # Check quality and report if requested
        if report:

            # Report distribution of insert sizes
            report_insert_size_distribution(data_folder, adaID, VERBOSE=VERBOSE)

            # Report rough coverage of reference
            report_coverage(data_folder, adaID, reference=refname, VERBOSE=VERBOSE)

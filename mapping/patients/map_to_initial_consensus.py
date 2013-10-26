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
import argparse
import numpy as np
import pysam

from mapping.generic_utils import mkdirs
from mapping.datasets import MiSeq_runs
from mapping.mapping_utils import stampy_bin, subsrate, \
        convert_sam_to_bam, convert_bam_to_sam
from mapping.patients.patients import get_patient, get_sequenced_samples, \
        get_initial_sequenced_sample
from mapping.samples import samples as samples_seq



# Globals
# Stampy parameters
stampy_gapopen = 60	        # Default: 40
stampy_gapextend = 5 	    # Default: 3
stampy_sensitive = True     # Default: False

# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBSCRIPT = JOBDIR+'patients/map_to_initial_consensus.py'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
# Different times based on how many threads are used
cluster_time = ['23:59:59', '0:59:59']
vmem = '8G'



# Functions
def fork_self(pname, sample, fragment, VERBOSE=3, threads=1,
              n_pairs=-1, filter_reads=False):
    '''Fork self for each sample and fragment'''
    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'm '+sample+' '+fragment,
                 '-l', 'h_rt='+cluster_time[threads >= 10],
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--patient', pname,
                 '--samples', sample,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                 '--threads', threads,
                 '-n', n_pairs
                ]
    if filter_reads:
        qsub.append('--filter')
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)



def make_output_folders(pname, sample, VERBOSE=0):
    '''Make the output folders if necessary for hash and map'''
    hash_foldername = os.path.dirname(get_initial_hash_filename(pname, 'F0'))
    map_foldername = os.path.dirname(get_mapped_to_initial_filename(pname, sample, 'F0'))
    foldernames = [hash_foldername, map_foldername]

    # Make the folders (using custom -p option -- poor Guido!)
    for dirname in foldernames:
        if not os.path.isdir(dirname):
            mkdirs(dirname)
        if VERBOSE:
            print 'Folder created:', dirname


def make_index_and_hash(pname, fragment, VERBOSE=0):
    '''Make index and hash files for consensus'''
    # 1. Make genome index file
    if not os.path.isfile(get_initial_index_filename(pname, fragment)):
        sp.call([stampy_bin,
                 '--species="HIV fragment '+fragment+'"',
                 '-G', get_initial_index_filename(pname, fragment, ext=False),
                 get_initial_consensus_filename(pname, fragment),
                 ])
        if VERBOSE:
            print 'Built index: '+pname+' '+fragment
    
    # 2. Build a hash file
    if not os.path.isfile(get_initial_hash_file(pname, fragment)):
        sp.call([stampy_bin,
                 '-g', get_initial_index_filename(pname, fragment, ext=False),
                 '-H', get_initial_hash_filename(pname, fragment, ext=False),
                 ])
        if VERBOSE:
            print 'Built hash: '+pname+' '+fragment


def map_stampy(pname, sample, fragment
               VERBOSE=VERBOSE, threads=threads, n_pairs=n_pairs):
    '''Map using stampy'''
    # Resolve ambiguous fragment primers
    if fragment in ['F5a', 'F5b']:
        raise ValueError('F5a/F5b harmony not implemented yet!')
    else:
        frag_gen = fragment 
    if VERBOSE:
        print 'Map via stampy: '+sample+' '+frag_gen

    # Get input filenames
    miseq_run = samples_seq[sample]['run']
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']
    adaID = dataset['adapters'][dataset['samples'].index(sample)]
    input_filename = get_divided_filenames(data_folder, adaID, [fragment], type='bam')[0]

    # parallelize if requested
    if threads == 1:

        # Get output filename
        output_filename = get_mapped_to_initial_filename(pname, sample, frag_gen,
                                                         type='sam')

        # Map
        call_list = [stampy_bin,
                     '-g', get_initial_index_filename(pname, fragment, ext=False),
                     '-h', get_initial_hash_filename(pname, fragment, ext=False),
                     '-o', output_filename,
                     '--substitutionrate='+subsrate,
                     '--gapopen', stampy_gapopen,
                     '--gapextend', stampy_gapextend]
        if stampy_sensitive:
            call_list.append('--sensitive')

        # Take only a (random) subsample: stampy uses the fraction of reads
        # intead of the number
        if n_pairs > 0:
            with pysam.Samfile(input_filename, 'rb') as in_bam:
                n_pairs_tot = sum(1 for read in in_bam) / 2
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
            output_filename = get_mapped_to_initial_filename(pname, sample, frag_gen,
                                                             type='sam', part=(j+1))
            # Map
            call_list = ['qsub','-cwd',
                         '-b', 'y',
                         '-S', '/bin/bash',
                         '-o', JOBLOGOUT,
                         '-e', JOBLOGERR,
                         '-N', 'm '+sample+frag_gen+' p'+str(j+1),
                         '-l', 'h_rt='+cluster_time[threads >= 10],
                         '-l', 'h_vmem='+vmem,
                         stampy_bin,
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
        output_file_parts = [get_mapped_to_initial_filename(pname, sample, frag_gen,
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

        # Concatenate output files
        output_filename = get_mapped_to_initial_filename(pname, sample, frag_gen,
                                                         type='bam', unsorted=True)
        if VERBOSE >= 1:
            print 'Concatenate premapped reads: sample '+sample
        pysam.cat('-o', output_filename, *output_file_parts)

        # Sort the file by read names (to ensure the pair_generator)
        output_filename_sorted = get_mapped_to_initial_filename(pname, sample, frag_gen,
                                                                type='bam')
        # NOTE: we exclude the extension and the option -f because of a bug in samtools
        if VERBOSE >= 1:
            print 'Sort mapped reads: sample '+sample
        pysam.sort('-n', output_filename, output_filename_sorted[:-4])

        # Reheader the file without BAM -> SAM -> BAM
        if VERBOSE >= 1:
            print 'Reheader mapped reads: sample '+sample
        header_filename = get_mapped_to_initial_filename(pname, sample, frag_gen,
                                                         type='sam', part=1)
        pysam.reheader(header_filename, output_filename_sorted)



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
    parser.add_argument('-n', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use for mapping')
    parser.add_argument('--filter', action='store_true',
                        help='Filter reads immediately after mapping')

    args = parser.parse_args()
    pname = args.patient
    samples = args.samples
    fragments = args.fragments
    submit = args.submit
    VERBOSE = args.verbose
    threads = args.threads
    n_pairs = args.n
    filter_reads = args.filter

    # Get the patient
    patient = get_patient(pname)

    # If no samples are mentioned, use all sequenced ones
    if not samples:
        samples = get_sequenced_samples(patient)
    if VERBOSE >= 3:
        print 'samples', samples

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Get the initial sample
    sample_init = get_initial_sequenced_sample(patient)

    # Make output folders if necessary
    for sample in samples:
        make_output_folders(pname, sample, VERBOSE=VERBOSE)

    # Iterate over samples and fragments
    for fragment in fragments:
        
        # Make stampy hashes
        make_index_and_hash(pname, fragment, VERBOSE=VERBOSE)
            
        for sample in samples:

            # Submit to the cluster self if requested
            if submit:
                fork_self(pname, sample, fragment,
                          VERBOSE=VERBOSE, threads=threads,
                          n_pairs=n_pairs,
                          filter_reads=filter_reads)
                continue

            # Fragment F5 has two sets of primers, so it's a mess in this case,
            # because the initial consensus might not have all the bases :-/
            # FIXME: make sure the initial consensus has ;-)
            if fragment == 'F5':
                raise ValueError('F5a/F5b harmony not implemented yet!')

            # Map via stampy
            map_stampy(pname, sample, fragment,
                       VERBOSE=VERBOSE, threads=threads, n_pairs=n_pairs)

            # Filter reads after mapping if requested
            if filter_reads:
                import warnings
                warnings.warn('Filtering not implemented yet')
                #filter_mapped_reads(pname, sample, fragment, VERBOSE=VERBOSE)

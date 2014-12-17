#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       25/10/13
content:    Map reads to the initial reference of the patient (in order to have
            the whole patient data aligned to one reference).

            This script must cope with a flexible sequencing schedule, e.g.
            repeated runs, some fragments in one run and some in another, etc.
            It accepts a list of either patients XOR samples. Samples can be
            either sequenced samples or patient samples, in which case we iterate
            over sequencing runs that contribute to them.
'''
# Modules
import sys
import os
import time
import argparse
import subprocess as sp
import numpy as np
import pandas as pd
import pysam
import warnings

from hivwholeseq.sequencing.samples import SampleSeq
from hivwholeseq.patients.patients import load_patients, load_patient, Patient
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.mapping_utils import stampy_bin, subsrate, \
        convert_sam_to_bam, convert_bam_to_sam, get_number_reads
from hivwholeseq.patients.filenames import get_initial_index_filename, \
        get_initial_hash_filename, get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_mapped_to_initial_foldername, \
        get_map_initial_summary_filename
from hivwholeseq.cluster.fork_cluster import fork_map_to_initial_reference as fork_self
from hivwholeseq.clean_temp_files import remove_mapped_init_tempfiles
from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.sequencing.samples import load_samples_sequenced as lss



# Globals
stampy_gapopen = 60	        # Default: 40
stampy_gapextend = 5 	    # Default: 3
stampy_sensitive = True     # Default: False


# Functions
def get_input_filename(data_folder, adaID, frag_spec, type='bam', only_chunk=None,
                       filtered=True):
    '''Get filename of input for mapping to initial reference'''
    # We should take reads filtered after mapping to the auto-consensus
    if filtered:
        from hivwholeseq.sequencing.filenames import get_mapped_filename
        frag_gen = frag_spec[:2]
        fn = get_mapped_filename(data_folder, adaID, frag_gen, type='bam', filtered=True)
    else:
        from hivwholeseq.sequencing.filenames import get_divided_filename
        fn = get_divided_filename(data_folder, adaID, frag_spec, type='bam', chunk=only_chunk)
    return fn


def make_output_folders(pname, samplename, PCR=1, VERBOSE=0):
    '''Make the output folders if necessary for hash and map'''
    hash_foldername = os.path.dirname(get_initial_hash_filename(pname, 'F0'))
    map_foldername = get_mapped_to_initial_foldername(pname, samplename, PCR=PCR)

    if not os.path.isdir(hash_foldername):
        mkdirs(hash_foldername)
        if VERBOSE:
            print 'Folder created:', hash_foldername

    mkdirs(map_foldername)
    if VERBOSE:
        print 'Folder created:', map_foldername


def make_index_and_hash(pname, fragment, VERBOSE=0):
    '''Make index and hash files for reference'''
    # 1. Make genome index file
    stdout = sp.check_output([stampy_bin,
                              '--overwrite',
                              '--species="HIV fragment '+fragment+'"',
                              '-G', get_initial_index_filename(pname, fragment, ext=False),
                              get_initial_reference_filename(pname, fragment),
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


def map_stampy_singlethread(sample, fragment, VERBOSE=0, n_pairs=-1,
                            summary=True, only_chunk=None, filtered=True):
    '''Map using stampy, single thread (no cluster queueing race conditions)'''
    pname = sample.patient
    samplename_pat = sample['patient sample']
    seq_run = sample['seq run']
    data_folder = sample.sequencing_run['folder']
    adaID = sample['adapter']
    PCR = int(sample.PCR)

    if VERBOSE:
        print 'Map via stampy (single thread): '+samplename+' '+fragment

    if summary:
        summary_filename = get_map_initial_summary_filename(pname, samplename_pat, 
                                                            samplename, fragment,
                                                            PCR=PCR)

    # Specific fragment (e.g. F5 --> F5bi)
    frag_spec = filter(lambda x: fragment in x, sample.regions_complete)
    if not len(frag_spec):
        if summary:
            with open(summary_filename, 'a') as f:
                f.write('Failed (specific fragment for '+fragment+'not found).\n')

        raise ValueError(samplename+': fragment '+fragment+' not found.')
    else:
        frag_spec = frag_spec[0]

    input_filename = get_input_filename(data_folder, adaID, frag_spec, type='bam',
                                        only_chunk=only_chunk, filtered=filtered)

    # NOTE: we introduced fragment nomenclature late, e.g. F3a. Check for that
    if not os.path.isfile(input_filename):
        if fragment == 'F3':
            input_filename = input_filename.replace('F3a', 'F3')

    # Check existance of input file, because stampy creates output anyway
    if not os.path.isfile(input_filename):
        if summary:
            with open(summary_filename, 'a') as f:
                f.write('Failed (input file for mapping not found).\n')

        raise ValueError(samplename+', fragment '+fragment+': input file not found.')

    # Extract subsample of reads if requested
    if n_pairs > 0:
        from hivwholeseq.mapping_utils import extract_mapped_reads_subsample
        input_filename_sub = get_mapped_to_initial_filename(pname, samplename_pat,
                                                            samplename, fragment,
                                                            PCR=PCR,
                                                            type='bam')[:-4]+\
                '_unmapped.bam'
        n_written = extract_mapped_reads_subsample(input_filename,
                                                   input_filename_sub,
                                                   n_pairs, VERBOSE=VERBOSE)

    # Get output filename
    output_filename = get_mapped_to_initial_filename(pname, samplename_pat, 
                                                     samplename, fragment,
                                                     PCR=PCR,
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

    output_filename_bam = get_mapped_to_initial_filename(pname, samplename_pat,
                                                         samplename, fragment,
                                                         type='bam',
                                                         PCR=PCR,
                                                         only_chunk=only_chunk)
    convert_sam_to_bam(output_filename_bam)

    if summary:
        with open(summary_filename, 'a') as f:
            f.write('Stampy mapped (single thread).\n')

    if only_chunk is None:
        if VERBOSE >= 1:
            print 'Remove temporary files: sample '+samplename
        remove_mapped_init_tempfiles(pname, samplename_pat,
                                     samplename, fragment,
                                     PCR=PCR,
                                     VERBOSE=VERBOSE, only_chunk=only_chunk)

    if summary:
        with open(summary_filename, 'a') as f:
            f.write('Temp mapping files removed.\n')
            f.write('\n')

    if n_pairs > 0:
        os.remove(input_filename_sub)


def map_stampy_multithread(sample, fragment, VERBOSE=0, threads=2, summary=True,
                           filtered=True):
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

    input_filename = get_input_filename(data_folder, adaID, frag_spec, type='bam')

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


def map_stampy(sample, fragment, VERBOSE=0, threads=1, n_pairs=-1,
               summary=True, only_chunk=None, filtered=True):
    '''Map using stampy'''
    if threads == 1:
        map_stampy_singlethread(sample, fragment,
                                VERBOSE=VERBOSE, n_pairs=n_pairs,
                                summary=summary, only_chunk=only_chunk,
                                filtered=filtered)

    else:
        if (n_pairs != -1):
            raise ValueError('Cannot multithread a subset - why would you do this?')
        elif only_chunk:
            raise ValueError('Cannot multithread a chunk - why would you do this?')

        map_stampy_multithread(sample, fragment, threads=threads,
                               VERBOSE=VERBOSE, summary=summary,
                               filtered=filtered)


def get_number_chunks(pname, samplename, fragment, VERBOSE=0):
    '''Find the number of chunks the reads have been chopped down into'''
    sample = patient.sample_table.loc[samplename]
    seq_run = sample['run']
    data_folder = MiSeq_runs[seq_run]['folder']
    adaID = sample['adaID']
    # Has this fragment been sequenced in that sample?
    frag_spec = filter(lambda x: fragment in x, sample['fragments'])
    if not len(frag_spec):
        warnings.warn(pname+', '+samplename+': fragment '+fragment+' not found.')
        return
    else:
        frag_spec = frag_spec[0]

    n_chunks = 0
    input_filename = get_input_filename(data_folder, adaID, frag_spec, type='bam', chunk=n_chunks+1)
    while(os.path.isfile(input_filename)):
        n_chunks += 1
        input_filename = get_input_filename(data_folder, adaID, frag_spec, type='bam', chunk=n_chunks+1)

    return n_chunks



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Map to initial reference',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
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
    parser.add_argument('--skiphash', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')
    parser.add_argument('--chunks', type=int, nargs='+', default=[None],
                        help='Only map some chunks (cluster optimization): 0 for automatic detection')
    parser.add_argument('--unfiltered', action='store_false', dest='filtered',
                        help='Map unfiltered reads (for quick checks only)')
    parser.add_argument('--include-contaminated', action='store_true',
                        help='Include majorly contaminated samples in the map (e.g. 12879 F4)')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    submit = args.submit
    VERBOSE = args.verbose
    threads = args.threads
    n_pairs = args.maxreads
    skip_hash = args.skiphash
    summary = args.summary
    only_chunks = args.chunks
    filtered = args.filtered
    use_contaminated = args.include_contaminated

    # Collect all sequenced samples from patients
    samples_pat = lssp()
    if pnames is not None:
        samples_seq = []
        for pname in pnames:
            patient = load_patient(pname)
            patient.discard_nonsequenced_samples()
            for samplename_pat, sample_pat in patient.samples.iterrows():
                sample_pat = SamplePat(sample_pat)
                samples_seq.append(sample_pat.samples_seq)
        samples_seq = pd.concat(samples_seq)

    else:
        samples_seq = lss()
        ind = samples_pat.index.isin(samplenames)
        if ind.sum():
            samplenames_pat = samples_pat.index[ind]
            samples_seq = samples_seq.loc[samples_seq['patient sample'].isin(samplenames_pat)]
        else:
            samples_seq = samples_seq.loc[samples_seq.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples_seq.index.tolist()
        
    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    for samplename, sample in samples_seq.iterrows():
        sample = SampleSeq(sample)

        samplename_pat = sample['patient sample']
        sample_pat = samples_pat.loc[samplename_pat] 
        sample['patient'] = pname = sample_pat.patient
        PCR = int(sample.PCR)
        fragments_sample = sorted(set(sample.regions_generic) & set(fragments))

        if VERBOSE:
            print samplename, samplename_pat, pname, PCR

        if not skip_hash:
            make_output_folders(pname, samplename_pat, PCR=PCR, VERBOSE=VERBOSE)

        for fragment in fragments_sample:
            if VERBOSE:
                print fragment

            # Check for contamination
            contstr = sample['suspected contamination']
            if (not use_contaminated) and pd.notnull(contstr) and (fragment in contstr):
                print 'WARNING: This sample has a suspected contamination! Skipping.'
                continue

            if not skip_hash:
                make_index_and_hash(pname, fragment, VERBOSE=VERBOSE)
    
            # only_chunks is used when we fragment the mapping into small aliquots
            # to use the short_queue on the cluster (it's a hack). [0] means all
            # chunks, automatically detected by the presence of the input files;
            # [None] means normal mapping, [1, 2, ..] means to map only those chunks.
            if only_chunks == [0]:
                only_chunks_sample = range(1, get_number_chunks(pname, samplename, fragment, VERBOSE=VERBOSE) + 1)
            else:
                only_chunks_sample = only_chunks

            for only_chunk in only_chunks_sample:
    
                # If the input file if missing, skip
                input_filename = get_input_filename(sample.seqrun_folder,
                                                    sample.adapter,
                                                    sample.convert_region(fragment),
                                                    type='bam',
                                                    only_chunk=only_chunk,
                                                    filtered=filtered)
                if not os.path.isfile(input_filename):
                    if VERBOSE:
                        print 'WARNING: input file not found'
                    continue

                if submit:
                    fork_self(samplename, fragment,
                              VERBOSE=VERBOSE, threads=threads,
                              n_pairs=n_pairs,
                              summary=summary,
                              only_chunks=[only_chunk],
                              filtered=filtered)
                    continue

                if summary:
                    sfn = get_map_initial_summary_filename(pname, samplename_pat,
                                                           samplename, fragment,
                                                           PCR=PCR, only_chunk=only_chunk)
                    with open(sfn, 'w') as f:
                        f.write('Call: python map_to_initial_reference.py'+\
                                ' --samples '+samplename+\
                                ' --fragments '+fragment+\
                                ' --threads '+str(threads)+\
                                ' --verbose '+str(VERBOSE))
                        if n_pairs != -1:
                            f.write(' --maxreads '+str(n_pairs))
                        if only_chunk is not None:
                            f.write('--chunks '+str(only_chunk))
                        if not filtered:
                            f.write(' --unfiltered')
                        f.write('\n')
    
    
                map_stampy(sample, fragment,
                           VERBOSE=VERBOSE, threads=threads, n_pairs=n_pairs,
                           summary=summary, only_chunk=only_chunk, filtered=filtered)

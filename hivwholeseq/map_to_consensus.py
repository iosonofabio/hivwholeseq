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
from itertools import izip
import pysam
import numpy as np
import subprocess as sp

from hivwholeseq.adapter_info import load_adapter_table, foldername_adapter
from hivwholeseq.mapping_utils import stampy_bin, subsrate, convert_sam_to_bam, \
        convert_bam_to_sam, get_number_reads
from hivwholeseq.filenames import get_consensus_filename, get_mapped_filename,\
        get_read_filenames, get_divided_filename, get_map_summary_filename, \
        get_reference_consensus_ali_filename        
from hivwholeseq.filter_mapped_reads import match_len_min, trim_bad_cigars
from hivwholeseq.filter_mapped_reads import filter_reads as filter_mapped_reads
from hivwholeseq.fork_cluster import fork_map_to_consensus as fork_self
from hivwholeseq.samples import load_sequencing_run, SampleSeq
from hivwholeseq.clean_temp_files import remove_mapped_tempfiles



# Globals
# Stampy parameters
stampy_sensitive = True     # Default: False


# Cluster submit
import hivwholeseq
JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr/'
JOBLOGOUT = JOBDIR+'logout/'
vmem = '8G'


# Functions
def check_consensus_length(data_folder, adaID, fragment, VERBOSE=0):
    '''Check consensus length, and if too short or absent complain'''
    from Bio import AlignIO
    ali_fn = get_reference_consensus_ali_filename(data_folder, adaID, fragment)
    if not os.path.isfile(ali_fn):
        if VERBOSE >= 2:
            print 'Consensus alignment to reference not found', adaID, fragment
        return False

    ali = AlignIO.read(ali_fn, 'fasta')
    len_ref = len(ali[0].seq.ungap('-'))
    len_cons = len(ali[1].seq.ungap('-'))
    if len_cons < len_ref - 200:
        if VERBOSE >= 2:
            print 'Consensus alignment to reference too short: ref', len_ref, 'cons:', len_cons
        return False
    elif len_cons > len_ref + 200:
        if VERBOSE >= 2:
            print 'Consensus alignment to reference too long: ref', len_ref, 'cons:', len_cons
        return False

    if VERBOSE >= 2:
        print 'Consensus checked, has approximately the right length: ref', len_ref, 'cons:', len_cons
    return True


def estimate_cluster_time(threads, filter_reads, VERBOSE=0):
    '''Estimate which queue to request on the cluster'''
    if (threads >= 15) and (not filter_reads):
        cluster_time = '0:59:59'
    else:
        cluster_time = '23:59:59'

    return cluster_time


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


def make_output_folders(data_folder, adaID, VERBOSE=0):
    '''Make the output folders if necessary for hash and map'''
    hash_foldername = os.path.dirname(get_hash_file(data_folder, adaID, 'F0'))
    map_foldername = os.path.dirname(get_mapped_filename(data_folder, adaID, 'F0'))
    foldernames = [hash_foldername, map_foldername]

    # Make the folders
    for dirname in foldernames:
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        if VERBOSE:
            print 'Folder created:', dirname


def make_index_and_hash(data_folder, adaID, fragment, VERBOSE=0, summary=True):
    '''Make index and hash files for consensus'''
    frag_gen = fragment[:2]

    # NOTE: we can use --overwrite here, because there is no concurrency (every
    # job has its own hash)
    # 1. Make genome index file
    sp.call([stampy_bin,
             '--species="HIV fragment '+frag_gen+'"',
             '--overwrite',
             '-G', get_index_file(data_folder, adaID, frag_gen, ext=False),
             get_consensus_filename(data_folder, adaID, frag_gen, trim_primers=True),
             ])
    if VERBOSE:
        print 'Built index: '+adaID+' '+frag_gen
    
    # 2. Build a hash file
    sp.call([stampy_bin,
             '--overwrite',
             '-g', get_index_file(data_folder, adaID, frag_gen, ext=False),
             '-H', get_hash_file(data_folder, adaID, frag_gen, ext=False),
             ])
    if VERBOSE:
        print 'Built hash: '+adaID+' '+frag_gen

    if summary:
        with open(get_map_summary_filename(data_folder, adaID, frag_gen), 'a') as f:
            f.write('\n')
            f.write('Stampy index and hash written.')
            f.write('\n')


def map_stampy(data_folder, adaID, fragment, VERBOSE=0, threads=1,
               cluster_time='23:59:59', maxreads=-1, summary=True,
               rescue=False, dry=False):
    '''Map using stampy'''
    frag_gen = fragment[:2]

    if summary:
        summary_filename = get_map_summary_filename(data_folder, adaID, frag_gen,
                                                    rescue=rescue)

    # Set mapping penalty scores: softer for rescues and F3 and F5
    global subsrate
    if rescue:
        subsrate = '0.2'
        stampy_gapopen = 5	    # Default: 40
        stampy_gapextend = 1 	    # Default: 3

    elif frag_gen not in ('F3', 'F5'):
        stampy_gapopen = 60	    # Default: 40
        stampy_gapextend = 5 	    # Default: 3
    else:
        stampy_gapopen = 30	    # Default: 40
        stampy_gapextend = 2 	    # Default: 3

    if VERBOSE:
        print 'Map via stampy: '+adaID+' '+frag_gen

    if not rescue: 
        input_filename = get_divided_filename(data_folder, adaID, fragment, type='bam')

        # NOTE: we introduced fragment nomenclature late, e.g. F3a. Check for that
        if not os.path.isfile(input_filename):
            if frag_gen == 'F3':
                input_filename = input_filename.replace('F3a', 'F3')

    else:
        input_filename = get_divided_filename(data_folder, adaID, 'unmapped', type='bam')

    # Check existance of input file, because stampy creates output anyway
    if not os.path.isfile(input_filename):
        if summary:
            with open(summary_filename, 'a') as f:
                f.write('Failed (input file for mapping not found).\n')

        raise ValueError(samplename+', fragment '+fragment+': input file not found.')

    # parallelize if requested
    if threads == 1:

        output_filename = get_mapped_filename(data_folder, adaID, frag_gen, type='sam',
                                              rescue=rescue)

        # Map
        call_list = [stampy_bin,
                     '-g', get_index_file(data_folder, adaID, frag_gen, ext=False),
                     '-h', get_hash_file(data_folder, adaID, frag_gen, ext=False), 
                     '-o', output_filename,
                     '--overwrite',
                     '--substitutionrate='+subsrate,
                     '--gapopen', stampy_gapopen,
                     '--gapextend', stampy_gapextend]
        if stampy_sensitive:
            call_list.append('--sensitive')

        # Take only a (random) subsample: stampy uses the fraction of reads
        # intead of the number
        if maxreads > 0:
            # FIXME: figure out the -s option and the --numrecords option
            call_list.extend(['--numrecords', maxreads])
            
            #n_pairs_tot = get_number_reads(input_filename, 'bam') / 2
            #frac_pairs = 1.0 * maxreads / n_pairs_tot
            #random_seed = np.random.randint(1e5)
            #call_list.extend(['-s', frac_pairs + random_seed])

        call_list = call_list + ['-M', input_filename]
        call_list = map(str, call_list)
        if VERBOSE >=2:
            print ' '.join(call_list)

        if not dry:
            sp.call(call_list)

            if summary:
                with open(summary_filename, 'a') as f:
                    f.write('Stampy mapped (single thread).\n')

            output_filename = get_mapped_filename(data_folder, adaID, frag_gen, type='bam',
                                                  rescue=rescue)
            convert_sam_to_bam(output_filename)
        else:
            if summary:
                with open(summary_filename, 'a') as f:
                    f.write('Dry run works (single thread).\n')

            if VERBOSE >= 1:
                print 'Dry run works (single thread)'

            return

    else:

        # Submit map script
        jobs_done = np.zeros(threads, bool)
        job_IDs = np.zeros(threads, 'S30')
        for j in xrange(threads):
        
            # Get output filename
            output_filename =  get_mapped_filename(data_folder, adaID, frag_gen,
                                               type='sam', part=(j+1), rescue=rescue)
            # Map
            call_list = ['qsub','-cwd',
                         '-b', 'y',
                         '-S', '/bin/bash',
                         '-o', JOBLOGOUT,
                         '-e', JOBLOGERR,
                         '-N', 'm'+adaID.replace('-', '')+frag_gen+str(j+1),
                         '-l', 'h_rt='+cluster_time,
                         '-l', 'h_vmem='+vmem,
                         stampy_bin,
                         '-g', get_index_file(data_folder, adaID, frag_gen, ext=False),
                         '-h', get_hash_file(data_folder, adaID, frag_gen, ext=False), 
                         '-o', output_filename,
                         '--overwrite',
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

            if not dry:
                job_ID = sp.check_output(call_list)
                job_ID = job_ID.split()[2]
                job_IDs[j] = job_ID

        if dry:
            if summary:
                with open(summary_filename, 'a') as f:
                    f.write('Dry run works (multi thread).\n')

            if VERBOSE >= 1:
                print 'Dry run works (multi thread)'
            return

        # Monitor output
        output_file_parts = [get_mapped_filename(data_folder, adaID, frag_gen,
                                                 type='bam', part=(j+1), rescue=rescue)
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
                               adaID+', fragment '+frag_gen+', part '+str(j+1)+ ' of '+ \
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
        output_filename = get_mapped_filename(data_folder, adaID, frag_gen, type='bam',
                                              unsorted=True, rescue=rescue)
        if VERBOSE >= 1:
            print 'Concatenate mapped reads: adaID '+adaID+', fragment '+frag_gen
        pysam.cat('-o', output_filename, *output_file_parts)
        if summary:
            with open(summary_filename, 'a') as f:
                f.write('BAM files concatenated (unsorted).\n')

        # Sort the file by read names (to ensure the pair_generator)
        output_filename_sorted = get_mapped_filename(data_folder, adaID, frag_gen,
                                                     type='bam',
                                                     unsorted=False,
                                                     rescue=rescue)
        # NOTE: we exclude the extension and the option -f because of a bug in samtools
        if VERBOSE >= 1:
            print 'Sort mapped reads: adaID '+adaID+', fragment '+frag_gen
        pysam.sort('-n', output_filename, output_filename_sorted[:-4])
        if summary:
            with open(summary_filename, 'a') as f:
                f.write('Joint BAM file sorted.\n')

        # Reheader the file without BAM -> SAM -> BAM
        if VERBOSE >= 1:
            print 'Reheader mapped reads: adaID '+adaID+', fragment '+frag_gen
        header_filename = get_mapped_filename(data_folder, adaID, frag_gen,
                                              type='sam', part=1, rescue=rescue)
        pysam.reheader(header_filename, output_filename_sorted)
        if summary:
            with open(summary_filename, 'a') as f:
                f.write('Joint BAM file reheaded.\n')

    # FIXME: check whether temp files are all deleted
    if VERBOSE >= 1:
        print 'Remove temporary files: adaID '+adaID+', fragment '+frag_gen
    remove_mapped_tempfiles(data_folder, adaID, frag_gen, VERBOSE=VERBOSE, rescue=rescue)
    if summary:
        with open(summary_filename, 'a') as f:
            f.write('Temp mapping files removed.\n')
            f.write('\n')



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Map reads to sample consensus',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='+',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use for mapping')
    parser.add_argument('--filter', action='store_true',
                        help='Filter reads immediately after mapping')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')
    parser.add_argument('--rescue', action='store_true',
                        help='Look at to-be-rescued reads (less stringent mapping)')
    parser.add_argument('--only-patient', action='store_true', dest='use_pats',
                        help='Map only patient samples')
    parser.add_argument('--dry', action='store_true',
                        help='Dry run (do everything except actual mapping)')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    submit = args.submit
    threads = args.threads
    maxreads = args.maxreads
    filter_reads = args.filter
    summary = args.summary
    use_rescue = args.rescue
    use_pats = args.use_pats
    use_dry = args.dry

    if submit and use_dry:
        raise ValueError('Won\'t submit a dry run!')

    dataset = load_sequencing_run(seq_run)
    data_folder = dataset.folder

    # If the script is called with no adaID, iterate over all
    samples = dataset.samples
    if adaIDs is not None:
        samples = samples.loc[samples.adapter.isin(adaIDs)]
    if VERBOSE >= 3:
        print 'adaIDs', samples.adapter

    for (samplename, sample) in samples.iterrows():
        sample = SampleSeq(sample)

        if str(sample.PCR) == 'nan':
            if VERBOSE:
                print samplename+': PCR type not found, skipping'
            continue

        if use_pats and (sample['patient sample'] == 'nan'):
            if VERBOSE:
                print samplename+': not a patient sample, skipping'
            continue

        adaID = sample.adapter

        # Make output folders if necessary
        make_output_folders(data_folder, adaID, VERBOSE=VERBOSE)

        # If the script is called with no fragment, iterate over all
        fragments_sample = sample.regions_complete
        fragments_sample_gen = sample.regions_generic
        if fragments is not None:
            fragments_sample = [fr for (fr, frg) in izip(fragments_sample,
                                                         fragments_sample_gen)
                                if frg in fragments]

        if VERBOSE >= 3:
            print 'adaID '+adaID+': fragments '+' '.join(fragments_sample)

        # Iterate over fragments
        for fragment in fragments_sample:

            frag_gen = fragment[:2]

            # Submit to the cluster self if requested
            if submit:
                fork_self(seq_run, adaID, frag_gen,
                          VERBOSE=VERBOSE,
                          threads=threads, maxreads=maxreads,
                          filter_reads=filter_reads,
                          summary=summary,
                          rescue=use_rescue)
                continue

            if summary:
                sfn = get_map_summary_filename(data_folder, adaID, frag_gen,
                                               rescue=use_rescue)
                with open(sfn, 'w') as f:
                    f.write('Call: python map_to_consensus.py'+\
                            ' --run '+seq_run+\
                            ' --adaIDs '+adaID+\
                            ' --fragments '+frag_gen+\
                            ' --threads '+str(threads)+\
                            ' --verbose '+str(VERBOSE))
                    if maxreads != -1:
                        f.write(' --maxreads '+str(maxreads))
                    if filter_reads:
                        f.write(' --filter')
                    f.write('\n')

            
            if not check_consensus_length(data_folder, adaID, fragment, VERBOSE=VERBOSE):
                if VERBOSE:
                    print 'Consensus not suitable for mapping:', adaID, fragment
                continue

            make_index_and_hash(data_folder, adaID, fragment, VERBOSE=VERBOSE,
                                summary=summary)

            cluster_time = estimate_cluster_time(threads, filter_reads,
                                                 VERBOSE=VERBOSE) 
            
            map_stampy(data_folder, adaID, fragment,
                       VERBOSE=VERBOSE,
                       threads=threads,
                       cluster_time=cluster_time,
                       maxreads=maxreads,
                       summary=summary,
                       rescue=use_rescue,
                       dry=use_dry)

            if filter_reads:
                if not rescue:
                    filter_mapped_reads(data_folder, adaID, fragment, VERBOSE=VERBOSE,
                                        summary=summary)

                else:
                    # TODO
                    print 'Filter for rescue not implemented yet!'

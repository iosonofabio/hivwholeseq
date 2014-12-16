#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/10/13
content:    Premap demultiplex reads to HXB2 (or other reference) 'as is',
            before any filtering. This is just a rough mapping, just to get
            the classification into fragments right.

            Note: stampy can do some kind of multithreading (see below).

            Operations in this script:
                1. taylor reference according to the primers used in the PCR
                2. make hash tables for mapper
                3. map to reference (including sorting and cleanup)
                4. generate report figures and summary file
'''
# Modules
import sys
import os
import time
import subprocess as sp
import argparse
import numpy as np
from Bio import SeqIO
import pysam

from hivwholeseq.sequencing.samples import load_sequencing_run
from hivwholeseq.sequencing.filenames import get_custom_reference_filename, \
        get_HXB2_hash_file, get_read_filenames, get_premapped_filename,\
        get_reference_premap_filename, get_premap_summary_filename, \
        get_reference_premap_index_filename, get_reference_premap_hash_filename,\
        get_coverage_figure_filename, get_insert_size_distribution_cumulative_filename,\
        get_insert_size_distribution_filename
from hivwholeseq.mapping_utils import stampy_bin, convert_sam_to_bam, convert_bam_to_sam
from hivwholeseq.fork_cluster import fork_premap as fork_self
from hivwholeseq.clean_temp_files import remove_premapped_tempfiles




# Functions
def make_output_folders(data_folder, adaID, VERBOSE=0, summary=True):
    '''Make output folders'''
    from hivwholeseq.generic_utils import mkdirs
    outfiles = [get_premapped_filename(data_folder, adaID)]
    if summary:
        outfiles.append(get_coverage_figure_filename(data_folder, adaID, 'premapped'))
    for outfile in outfiles:
        dirname = os.path.dirname(outfile)
        mkdirs(dirname)
        if VERBOSE:
            print 'Folder created:', dirname


def get_fragments(sample):
    '''Get fragments to trim the mapping reference to the relevant part'''
    fragments = sample.regions
    if str(fragments) in ('', 'nan'):
        if sample.PCR == 1:
            fragments = ('F1o', 'F6o')
        elif sample.PCR == 2:
            fragments = ('F1i', 'F6i')
        else:
            if VERBOSE:
                print 'Sample is neither PCR1 nor PCR2. Using whole reference.'
            fragments = None
    else:
        fragments = fragments.split(' ')
        if sample.PCR == 1:
            fragments = ['F'+fr+'o' for fr in fragments]
        elif sample.PCR == 2:
            fragments = ['F'+fr+'i' for fr in fragments]
        else:
            raise ValueError('Sample is neither PCR1 nor PCR2')

    print fragments

    return fragments


def make_reference(data_folder, adaID, fragments, refname, VERBOSE=0, summary=True):
    '''Make reference sequence trimmed to the necessary parts'''
    from hivwholeseq.reference import load_custom_reference
    seq = load_custom_reference(refname)

    output_filename = get_reference_premap_filename(data_folder, adaID)

    if fragments is None:
        seq_trim = seq
    else:
        # Look for the first fwd and the last rev primers to trim the reference
        # NOTE: this works even if F1 or F6 are missing (e.g. only F2-5 are seq-ed)!
        # If more than one primer is used for the first or last fragment, take the
        # longest reference
        from hivwholeseq.sequencing.primer_info import primers_PCR, primers_coordinates_HXB2
        if '+' in fragments[0]:
            fragment_subs = [fragments[0][:2]+fsub+fragments[0][-1]
                             for fsub in fragments[0][2:-1].split('+')]
            fr_pos_subs = [primers_coordinates_HXB2[fsub][0][0]
                           for fsub in fragment_subs]
            fragments[0] = fragment_subs[np.argmin(fr_pos_subs)]

        pr_fwd = primers_PCR[fragments[0]][0]

        if '+' in fragments[-1]:
            fragment_subs = [fragments[-1][:2]+fsub+fragments[-1][-1]
                             for fsub in fragments[-1][2:-1].split('+')]
            fr_pos_subs = [primers_coordinates_HXB2[fsub][1][1]
                           for fsub in fragment_subs]
            fragments[-1] = fragment_subs[np.argmax(fr_pos_subs)]

        pr_rev = primers_PCR[fragments[-1]][1]

        smat = np.array(seq)

        # Get all possible primers from ambiguous nucleotides and get the best match
        from hivwholeseq.sequence_utils import expand_ambiguous_seq as eas
        pr_fwd_mat = np.array(map(list, eas(pr_fwd)), 'S1')
        n_matches_fwd = [(smat[i: i + len(pr_fwd)] == pr_fwd_mat).sum(axis=1).max()
                         for i in xrange(len(seq) - len(pr_fwd))]
        pr_fwd_pos = np.argmax(n_matches_fwd)

        pr_rev_mat = np.array(map(list, eas(pr_rev)), 'S1')
        n_matches_rev = [(smat[i: i + len(pr_rev)] == pr_rev_mat).sum(axis=1).max()
                         for i in xrange(pr_fwd_pos + len(pr_fwd), len(seq) - len(pr_rev))]
        # Here you come from the right, i.e. look in the 3' LTR first
        pr_rev_pos = len(seq) - len(pr_rev) - 1 - np.argmax(n_matches_rev[::-1]) 

        output = [['Reference name:', refname]]
        output.append(['FWD primer:', fragments[0], str(pr_fwd_pos), pr_fwd])
        output.append(['REV primer:', fragments[-1], str(pr_rev_pos), pr_rev])
        output = '\n'.join(map(' '.join, output))
        if VERBOSE:
            print output

        if summary:
            with open(get_premap_summary_filename(data_folder, adaID), 'a') as f:
                f.write(output)
                f.write('\n')

        # The reference includes both the first fwd primer and the last rev one
        seq_trim = seq[pr_fwd_pos: pr_rev_pos + len(pr_rev)]
        seq_trim.id = '_'.join([seq_trim.id, str(pr_fwd_pos + 1),
                                str(pr_rev_pos + len(pr_rev))])
        seq_trim.name = '_'.join([seq_trim.name, str(pr_fwd_pos + 1),
                                  str(pr_rev_pos + len(pr_rev))])
        seq_trim.description = ' '.join([seq_trim.description,
                                         'from', str(pr_fwd_pos + 1),
                                         'to', str(pr_rev_pos + len(pr_rev)),
                                         '(indices from 1, extremes included)'])

    SeqIO.write(seq_trim, output_filename, 'fasta')

    if summary:
        with open(get_premap_summary_filename(data_folder, adaID), 'a') as f:
            f.write('Reference sequence written to: '+output_filename)
            f.write('\n')


def make_index_and_hash(data_folder, adaID, VERBOSE=0, summary=True):
    '''Make index and hash files for reference or consensus'''
    if VERBOSE:
        print 'Making index and hash files: adaID', adaID

    # 1. Make genome index file for reference
    if os.path.isfile(get_reference_premap_index_filename(data_folder, adaID, ext=True)):
        os.remove(get_reference_premap_index_filename(data_folder, adaID, ext=True))
    stdout = sp.check_output([stampy_bin,
                              '--species="HIV"',
                              '--overwrite',
                              '-G', get_reference_premap_index_filename(data_folder, adaID, ext=False),
                              get_reference_premap_filename(data_folder, adaID),
                              ],
                              stderr=sp.STDOUT)
    if VERBOSE:
        print 'Built index: '+adaID
    
    # 2. Build a hash file for reference
    if os.path.isfile(get_reference_premap_hash_filename(data_folder, adaID, ext=True)):
        os.remove(get_reference_premap_hash_filename(data_folder, adaID, ext=True))
    stdout = sp.check_output([stampy_bin,
                              '--overwrite',
                              '-g', get_reference_premap_index_filename(data_folder, adaID, ext=False),
                              '-H', get_reference_premap_hash_filename(data_folder, adaID, ext=False),
                              ],
                              stderr=sp.STDOUT)
    if VERBOSE:
        print 'Built hash: '+adaID

    if summary:
        with open(get_premap_summary_filename(data_folder, adaID), 'a') as f:
            f.write('\n')
            f.write('Stampy index and hash written.')
            f.write('\n')


def premap_stampy(data_folder, adaID, VERBOSE=0, threads=1, summary=True, maxreads=-1,
                  subsrate=0.05, gapopen=40, gapextend=3):
    '''Call stampy for actual mapping'''
    if VERBOSE:
        print 'Premapping: adaID ', adaID

    if summary:
        summary_filename = get_premap_summary_filename(data_folder, adaID)

    # Stampy can handle both gzipped and uncompressed fastq inputs
    input_filenames = get_read_filenames(data_folder, adaID, gzip=True)
    if not os.path.isfile(input_filenames[0]):
        input_filenames = get_read_filenames(data_folder, adaID, gzip=False)
    if not all(map(os.path.isfile, input_filenames)):
        raise OSError('Input files for mapping not found: '+input_filenames[0])

    # parallelize if requested
    if threads == 1:
        call_list = [stampy_bin,
                     '--overwrite',
                     '-g', get_reference_premap_index_filename(data_folder, adaID, ext=False),
                     '-h', get_reference_premap_hash_filename(data_folder, adaID, ext=False), 
                     '-o', get_premapped_filename(data_folder, adaID, type='sam'),
                     '--insertsize=450',
                     '--insertsd=100',
                     '--substitutionrate='+str(subsrate),
                     '--gapopen='+str(gapopen),
                     '--gapextend='+str(gapextend),
                    ]
        if maxreads > 0:
            call_list.append('--numrecords='+str(maxreads))
        call_list.extend(['-M'] + input_filenames)
        call_list = map(str, call_list)
        if VERBOSE >= 2:
            print ' '.join(call_list)
        sp.call(call_list)

        if summary:
            with open(get_premap_summary_filename(data_folder, adaID), 'a') as f:
                f.write('\nStampy premapped (single thread).\n')

        # Convert to compressed BAM
        convert_sam_to_bam(get_premapped_filename(data_folder, adaID, type='bam'))

        if summary:
            with open(summary_filename, 'a') as f:
                f.write('\nSAM file converted to compressed BAM: '+\
                        get_premapped_filename(data_folder, adaID, type='bam')+'\n')

    else:

        # Multithreading works as follows: call qsub + stampy, monitor the process
        # IDs with qstat at regular intervals, and finally merge results with pysam
        output_file_parts = [get_premapped_filename(data_folder, adaID, type='bam',
                                                part=(j+1)) for j in xrange(threads)]

        # Submit map script
        jobs_done = np.zeros(threads, bool)
        job_IDs = np.zeros(threads, 'S30')
        
        # Submit map call
        import hivwholeseq
        JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
        JOBLOGOUT = JOBDIR+'logout'
        JOBLOGERR = JOBDIR+'logerr'
        cluster_time = ['23:59:59', '1:59:59']
        vmem = '8G'
        for j in xrange(threads):
            call_list = ['qsub','-cwd',
                         '-b', 'y',
                         '-S', '/bin/bash',
                         '-o', JOBLOGOUT,
                         '-e', JOBLOGERR,
                         '-N', adaID+' p'+str(j+1),
                         '-l', 'h_rt='+cluster_time[threads >= 30],
                         '-l', 'h_vmem='+vmem,
                         stampy_bin,
                         '--overwrite',
                         '-g', get_reference_premap_index_filename(data_folder, adaID, ext=False),
                         '-h', get_reference_premap_hash_filename(data_folder, adaID, ext=False), 
                         '-o', get_premapped_filename(data_folder, adaID, type='sam', part=(j+1)),
                         '--processpart='+str(j+1)+'/'+str(threads),
                         '--insertsize=450',
                         '--insertsd=100',
                         '--substitutionrate='+str(subsrate),
                         '--gapopen='+str(gapopen),
                         '--gapextend='+str(gapextend),
                         '-M'] + input_filenames
            call_list = map(str, call_list)
            if VERBOSE >= 2:
                print ' '.join(call_list)
            job_ID = sp.check_output(call_list)
            job_ID = job_ID.split()[2]
            job_IDs[j] = job_ID

        # Monitor output
        time_wait = 10 # secs
        while not jobs_done.all():

            # Sleep some time
            time.sleep(time_wait)

            # Get the output of qstat to check the status of jobs
            qstat_output = sp.check_output(['qstat'])
            qstat_output = qstat_output.split('\n')[:-1] # The last is an empty line
            if VERBOSE >=3:
                print qstat_output
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
                               adaID+', part '+str(j+1)+ ' of '+ \
                               str(threads)
                    convert_sam_to_bam(output_file_parts[j])
                    # We do not need to wait if we did the conversion (it takes
                    # longer than some secs)
                    time_wait = 0
                    jobs_done[j] = True

        if summary:
            with open(summary_filename, 'a') as f:
                f.write('Stampy premapped ('+str(threads)+' threads).\n')

        # Concatenate output files
        if VERBOSE >= 1:
            print 'Concatenate premapped reads: adaID '+adaID+'...',
        output_filename = get_premapped_filename(data_folder, adaID, type='bam', unsorted=True)
        pysam.cat('-o', output_filename, *output_file_parts)
        if VERBOSE >= 1:
            print 'done.'
        if summary:
            with open(summary_filename, 'a') as f:
                f.write('BAM files concatenated (unsorted).\n')

        # Sort the file by read names (to ensure the pair_generator)
        # NOTE: we exclude the extension and the option -f because of a bug in samtools
        if VERBOSE >= 1:
            print 'Sort premapped reads: adaID '+adaID
        output_filename_sorted = get_premapped_filename(data_folder, adaID, type='bam', unsorted=False)
        pysam.sort('-n', output_filename, output_filename_sorted[:-4])
        if summary:
            with open(summary_filename, 'a') as f:
                f.write('Joint BAM file sorted.\n')

        # Reheader the file without BAM -> SAM -> BAM
        if VERBOSE >= 1:
            print 'Reheader premapped reads: adaID '+adaID
        header_filename = get_premapped_filename(data_folder, adaID, type='sam', part=1)
        pysam.reheader(header_filename, output_filename_sorted)
        if summary:
            with open(summary_filename, 'a') as f:
                f.write('Joint BAM file reheaded.\n')

    if VERBOSE >= 1:
        print 'Remove temporary files: adaID '+adaID
    remove_premapped_tempfiles(data_folder, adaID, VERBOSE=VERBOSE)
    if summary:
        with open(summary_filename, 'a') as f:
            f.write('Temp premapping files removed.\n')
            f.write('\n')


def report_insert_size(data_folder, adaID, seq_run, VERBOSE=0, summary=True):
    '''Produce figures of the insert size distribution'''
    from hivwholeseq.sequencing.check_insert_distribution import get_insert_size_distribution, \
            plot_cumulative_histogram, plot_histogram

    bins = np.linspace(0, 1000, 100)
    isz, h = get_insert_size_distribution(data_folder, adaID, 'premapped',
                                          bins=bins, maxreads=10000,
                                          VERBOSE=VERBOSE)

    plot_cumulative_histogram(data_folder, adaID, 'premapped', isz, savefig=True,
                              title='run '+str(seq_run)+', adaID '+str(adaID)+', premap',
                              lw=2, c='b')
    plot_histogram(data_folder, adaID, 'premapped', h, savefig=True,
                   title='run '+str(seq_run)+', adaID '+str(adaID)+', premap',
                   lw=2, color='b')

    if summary:
        with open(get_premap_summary_filename(data_folder, adaID), 'a') as f:
            f.write('\nInsert size distribution plotted:\n')
            f.write(get_insert_size_distribution_cumulative_filename(data_folder, adaID, 'premapped')+'\n')
            f.write(get_insert_size_distribution_filename(data_folder, adaID, 'premapped')+'\n')


def report_coverage(data_folder, adaID, VERBOSE=0, summary=True):
    '''Produce a report on rough coverage on reference (ignore inserts)'''
    ref_filename = get_reference_premap_filename(data_folder, adaID)
    refseq = SeqIO.read(ref_filename, 'fasta')

    # Prepare data structures
    coverage = np.zeros(len(refseq), int)
    
    # Parse the BAM file
    unmapped = 0
    mapped = 0
    bamfilename = get_premapped_filename(data_folder, adaID, type='bam')
    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for read in bamfile:
            if read.is_unmapped or (not read.is_proper_pair) or (not len(read.cigar)):
                unmapped += 1
                continue

            # Proceed along CIGARs
            ref_pos = read.pos
            for (bt, bl) in read.cigar:
                if bt not in (0, 2):
                    continue
                # Treat deletions as 'covered'
                coverage[ref_pos: ref_pos+bl] += 1
                ref_pos += bl
            mapped += 1

    # Save results
    from hivwholeseq.sequencing.filenames import get_coverage_figure_filename
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(13, 6))
    ax.plot(np.arange(len(refseq)), coverage + 1, lw=2, c='b')
    ax.set_xlabel('Position')
    ax.set_ylabel('Coverage')
    ax.set_yscale('log')
    ax.set_title('adaID '+adaID+', premapped', fontsize=18)
    ax.set_xlim(-20, len(refseq) + 20)
    plt.tight_layout()

    from hivwholeseq.generic_utils import mkdirs
    from hivwholeseq.sequencing.filenames import get_figure_folder
    mkdirs(get_figure_folder(data_folder, adaID))
    plt.savefig(get_coverage_figure_filename(data_folder, adaID, 'premapped'))
    plt.close(fig)

    if summary:
        with open(get_premap_summary_filename(data_folder, adaID), 'a') as f:
            f.write('\nPremapping results: '+\
                    str(mapped)+' read pairs mapped, '+str(unmapped)+' unmapped\n')
            f.write('\nCoverage plotted: '+\
                    get_coverage_figure_filename(data_folder, adaID, 'premapped')+'\n')



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Map against HIV reference strain',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--run', required=True,
                        help='MiSeq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use for mapping')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--reference', default='HXB2',
                        help='Use alternative reference (the file must exist)')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')
    parser.add_argument('--trimmed', action='store_true',
                        help='Use trimmed reads as input')
    parser.add_argument('--subsrate', type=float, default=0.05,
                        help='Baseline substitution rate')
    parser.add_argument('--gapopen', type=int, default=40,
                        help='Penality for gap opening')
    parser.add_argument('--gapextend', type=int, default=3,
                        help='Penality for gap extension')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    VERBOSE = args.verbose
    submit = args.submit
    threads = args.threads
    refname = args.reference
    summary = args.summary
    use_trimmed = args.trimmed
    maxreads = args.maxreads
    subsrate = args.subsrate
    gapopen = args.gapopen
    gapextend = args.gapextend

    # Specify the dataset
    dataset = load_sequencing_run(seq_run)
    data_folder = dataset.folder

    # If the script is called with no adaID, iterate over all
    samples = dataset.samples
    if adaIDs is not None:
        samples = samples.loc[samples.adapter.isin(adaIDs)]
    if VERBOSE >= 2:
        print samples.index.tolist()

    # Iterate over all adaIDs
    for samplename, sample in samples.iterrows():
        adaID = str(sample.adapter)

        # Submit to the cluster self if requested
        if submit:
            fork_self(seq_run, adaID, VERBOSE=VERBOSE, threads=threads,
                      reference=refname, summary=summary, trimmed=use_trimmed,
                      subsrate=subsrate, gapopen=gapopen, gapextend=gapextend,
                      maxreads=maxreads)
            continue

        make_output_folders(data_folder, adaID, VERBOSE=VERBOSE, summary=summary)

        if summary:
            with open(get_premap_summary_filename(data_folder, adaID), 'w') as f:
                outstr = 'Call: python premap_to_reference.py --run '+seq_run+\
                        ' --adaIDs '+adaID+\
                        ' --threads '+str(threads)+\
                        ' --reference '+refname+\
                        ' --subsrate '+str(subsrate)+\
                        ' --gapopen '+str(gapopen)+\
                        ' --gapextend '+str(gapextend)+\
                        ' --verbose '+str(VERBOSE)
                if maxreads != -1:
                    outstr = outstr + ' --maxreads '+str(maxreads)
                if use_trimmed:
                    outstr = outstr+' --trimmed'
                outstr = outstr+'\n'
                f.write(outstr)

        fragments = get_fragments(sample)

        make_reference(data_folder, adaID, fragments, refname, VERBOSE=VERBOSE, summary=summary)

        make_index_and_hash(data_folder, adaID, VERBOSE=VERBOSE, summary=summary)

        premap_stampy(data_folder, adaID, VERBOSE=VERBOSE, threads=threads, summary=summary,
                      subsrate=subsrate, gapopen=gapopen, gapextend=gapextend,
                      maxreads=maxreads)

        #FIXME: reactivate this
        #if summary:
        #    report_coverage(data_folder, adaID, VERBOSE=VERBOSE, summary=summary)

        #    report_insert_size(data_folder, adaID, seq_run,
        #                       VERBOSE=VERBOSE, summary=summary)

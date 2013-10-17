#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       28/08/13
content:    Classify the reads into the 6 fragments via a rough mapping before
            the final iterative mapping. This script reduces the chances of
            unpaired reads, unmapped, overhangs at fragments, and so on.

            This mapping is supposed to be fast and is a single run of stampy,
            so we run BWA before that to speed up.
'''
# Modules
import os
import time
import subprocess as sp
import argparse
import numpy as np
import pysam
from Bio.Seq import reverse_complement

from mapping.datasets import MiSeq_runs
from mapping.adapter_info import load_adapter_table
from mapping.filenames import get_HXB2_entire, get_HXB2_index_file,\
        get_HXB2_hash_file, get_read_filenames, get_premapped_file
from mapping.mapping_utils import stampy_bin, subsrate, bwa_bin,\
        convert_sam_to_bam, get_range_good_cigars, pair_generator, \
        convert_bam_to_sam
from mapping.primer_info import primers_coordinates_HXB2_inner as pci



# Globals
# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'divide_reads_into_fragments.py'
cluster_time = ['23:59:59', '0:59:59']
vmem = '8G'



# Functions
def fork_self(miseq_run, adaID, VERBOSE=0, subsample=False, bwa=False, threads=1):
    '''Fork self for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+'{:02d}'.format(adaID)

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'divide '+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time[subsample or (threads >= 30)],
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', miseq_run,
                 '--adaIDs', adaID,
                 '--verbose', VERBOSE,
                 '--threads', threads,
                ]
    if subsample:
        qsub_list.append('--subsample')
    if bwa:
        qsub_list.append('--bwa')
    qsub_list = map(str, qsub_list)
    if VERBOSE >= 2:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def make_output_folders(data_folder, adaID, VERBOSE=0):
    '''Make output folders'''
    output_filename = get_premapped_file(data_folder, adaID, subsample=subsample)
    dirname = os.path.dirname(output_filename)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
        if VERBOSE:
            print 'Folder created:', dirname


def make_index_and_hash(data_folder, cropped=True, VERBOSE=0):
    '''Make index and hash files for reference or consensus'''
    if VERBOSE:
        print 'Making index and hash files'

    # Make folder if necessary
    dirname = os.path.dirname(get_HXB2_hash_file())
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
        if VERBOSE:
            print 'Folder created:', dirname

    # 1. Make genome index file for HXB2
    if not os.path.isfile(get_HXB2_index_file(ext=True)):
        sp.call([stampy_bin,
                 '--species="HIV"',
                 '-G', get_HXB2_index_file(ext=False),
                 get_HXB2_entire(cropped=cropped),
                 ])
    
    # 2. Build a hash file for HXB2
    if not os.path.isfile(get_HXB2_hash_file(ext=True)):
        sp.call([stampy_bin,
                 '-g', get_HXB2_index_file(ext=False),
                 '-H', get_HXB2_hash_file(ext=False),
                 ])


def make_bwa_hash(data_folder, cropped=True, VERBOSE=0):
    '''Make hash tables for BWA'''
    if VERBOSE:
        print 'Build BWA index'

    # Make folder if necessary
    dirname =  os.path.dirname(get_HXB2_hash_file())
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # Call bwa index
    ref_filename = get_HXB2_entire(cropped=cropped)
    sp.call([bwa_bin,
             'index',
             ref_filename,
             ])

    # Move the hashes into subfolder
    from glob import glob
    import shutil
    shutil.copy(ref_filename, dirname)
    hash_files = glob(ref_filename+'.*')
    for hash_file in hash_files:
        shutil.copy(hash_file, dirname)
        os.remove(hash_file)


def premap_bwa(data_folder, adaID, VERBOSE=0, subsample=False, cropped=True):
    '''Premap using BWA'''
    if VERBOSE:
        print 'Map via BWA: '+'{:02d}'.format(adaID)

    # Get input filenames
    index_prefix = os.path.dirname(get_HXB2_hash_file(ext=False))+\
            '/'+os.path.basename(get_HXB2_entire(cropped=cropped))
            
    readfiles = get_read_filenames(data_folder, adaID,
                                   subsample=subsample,
                                   filtered=True)

    # Get output filename
    output_filename = get_premapped_file(data_folder, adaID,
                                         type='sam', subsample=subsample,
                                         bwa=True)
    # Map
    with open(output_filename, 'w') as f:
        call_list = [bwa_bin,
                     'mem',
                     index_prefix] + readfiles
        if VERBOSE >=2:
            print ' '.join(call_list)
        sp.call(call_list, stdout=f)


def premap_stampy(data_folder, adaID, VERBOSE=0, subsample=False, bwa=False,
                  threads=1):
    '''Call stampy for actual mapping'''
    if VERBOSE:
        print 'Map: adaID '+'{:02d}'.format(adaID)

    # Get input filenames
    if bwa:
        bwa_filename = get_premapped_file(data_folder, adaID,
                                        type='bam', subsample=subsample,
                                        bwa=True)
        if not os.path.isfile(bwa_filename):
            if VERBOSE >= 2:
                print 'Converting SAM to BAM'
            convert_sam_to_bam(bwa_filename)
        readfiles = [bwa_filename]

    else:
        readfiles = get_read_filenames(data_folder, adaID,
                                       subsample=subsample,
                                       filtered=True)

    # parallelize if requested
    if threads == 1:

        # Get output filename
        output_filename =  get_premapped_file(data_folder, adaID, type='sam',
                                              subsample=True)
        # Map
        call_list = [stampy_bin,
                     '-g', get_HXB2_index_file(ext=False),
                     '-h', get_HXB2_hash_file(ext=False), 
                     '-o', output_filename,
                     '--substitutionrate='+subsrate]
        if bwa:
            call_list.append('--bamkeepgoodreads')
        call_list = call_list + ['-M'] + readfiles
        call_list = map(str, call_list)
        if VERBOSE >= 2:
            print ' '.join(call_list)
        sp.call(call_list)

    else:

        jobs_done = np.zeros(threads, bool)
        job_IDs = np.zeros(threads, 'S30')

        # Submit map script
        for j in xrange(threads):
    
            # Get output filename
            output_filename =  get_premapped_file(data_folder, adaID, type='sam',
                                                  subsample=subsample, part=(j+1))
            # Map
            call_list = ['qsub','-cwd',
                         '-b', 'y',
                         '-S', '/bin/bash',
                         '-o', JOBLOGOUT,
                         '-e', JOBLOGERR,
                         '-N', 'div '+'{:02d}'.format(adaID)+' p'+str(j+1),
                         '-l', 'h_rt='+cluster_time[subsample or (threads >= 30)],
                         '-l', 'h_vmem='+vmem,
                         stampy_bin,
                         '-g', get_HXB2_index_file(ext=False),
                         '-h', get_HXB2_hash_file(ext=False), 
                         '-o', output_filename,
                         '--processpart='+str(j+1)+'/'+str(threads),
                         '--substitutionrate='+subsrate]
            if bwa:
                call_list.append('--bamkeepgoodreads')
            call_list = call_list + ['-M'] + readfiles
            call_list = map(str, call_list)
            if VERBOSE >= 2:
                print ' '.join(call_list)
            job_ID = sp.check_output(call_list)
            job_ID = job_ID.split()[2]
            job_IDs[j] = job_ID

        # Check output
        while not jobs_done.all():

            # Sleep 10 secs
            time.sleep(10)

            # Get the output of qstat to check the status of jobs
            qstat_output = sp.check_output(['qstat'])
            qstat_output = qstat_output.split('\n')[:-1] # The last is an empty line
            if len(qstat_output) < 3:
                jobs_done[:] = True
                break
            else:
                qstat_output = [line.split()[0] for line in qstat_output[2:]]

            for j in xrange(threads):
                if jobs_done[j]:
                    continue

                if job_IDs[j] not in qstat_output:
                    jobs_done[j] = True

        # Merge output files
        output_file_parts = [get_premapped_file(data_folder, adaID, type='bam',
                                                 subsample=subsample, part=(j+1))
                              for j in xrange(threads)]
        for output_file in output_file_parts:
            convert_sam_to_bam(output_file)
        output_filename = get_premapped_file(data_folder, adaID, type='bam',
                                             subsample=subsample, unsorted=True)
        pysam.merge(output_filename, *output_file_parts)

        # Sort the file by read names (to ensure the pair_generator)
        output_filename_sorted = get_premapped_file(data_folder, adaID, type='bam',
                                                    subsample=subsample,
                                                    unsorted=False)

        # Note: we exclude the extension and the option -f because of a bug in samtools
        pysam.sort('-n', output_filename, output_filename_sorted[:-4])

        # Make SAM out of the BAM for checking
        output_filename = get_premapped_file(data_folder, adaID, type='sam',
                                             subsample=subsample)
        convert_bam_to_sam(output_filename)


def write_read_pair(reads, ranges, fileouts):
    '''Write read pair into FASTQ files'''

    (read1, read2) = reads

    # Find orientation
    if not reads[0].is_reverse:
        fq1 = read1.seq[ranges[0][0]: ranges[0][1]]
        qual1 = read1.qual[ranges[0][0]: ranges[0][1]]
        fq2 = reverse_complement(read2.seq[ranges[1][0]: ranges[1][1]])
        qual2 = read2.qual[ranges[1][0]: ranges[1][1]][::-1]
    else:
        fq1 = reverse_complement(read1.seq[ranges[0][0]: ranges[0][1]])
        qual1 = read1.qual[ranges[0][0]: ranges[0][1]][::-1]
        fq2 = read2.seq[ranges[1][0]: ranges[1][1]]
        qual2 = read2.qual[ranges[1][0]: ranges[1][1]]

    # Note: we lose the x:N:0 part of the qname after the
    # whitespace, but that's used for barcoding only
    fileouts[0].write("@%s\n%s\n+\n%s\n" % (reads[0].qname, fq1, qual1))
    fileouts[1].write("@%s\n%s\n+\n%s\n" % (reads[1].qname, fq2, qual2))


def divide_into_fragments(data_folder, adaID, F5_type, n_cycles=500,
                          VERBOSE=0, subsample=False):
    '''Divide mapped reads into fastq files for the fragments'''
    if VERBOSE:
        print 'Divide into fragments: adaID '+'{:02d}'.format(adaID)

    # Extract fragments
    frags_sort = ['F1', 'F2', 'F3', 'F4', F5_type, 'F6']
    frags_pos = []
    for fragment in frags_sort:
        co = pci[fragment]
        frags_pos.append([co[0][0], co[1][1]])
    frags_pos = np.array(frags_pos, int).T
    # Since the reference is cropped, subtract from the positions F1 start
    frags_pos -= frags_pos.min()
    # Note: now we are in the reference of the CROPPED HXB2, and start from 0!

    # Get the premapped reads
    bamfilename = get_premapped_file(data_folder, adaID, type='bam',
                                     subsample=subsample)
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)

    # Output fastq files
    fileouts = get_read_filenames(data_folder, adaID,
                                  subsample=subsample,
                                  premapped=True)
    
    # Iterate over the mapped reads and assign fragments
    n_mapped = 0
    n_unmapped = 0
    n_ambiguous = 0
    with pysam.Samfile(bamfilename, 'rb') as bamfile,\
         open(fileouts[0][0], 'w') as fo1_F1,\
         open(fileouts[0][1], 'w') as fo1_F2,\
         open(fileouts[0][2], 'w') as fo1_F3,\
         open(fileouts[0][3], 'w') as fo1_F4,\
         open(fileouts[0][4], 'w') as fo1_F5,\
         open(fileouts[0][5], 'w') as fo1_F6,\
         open(fileouts[0][6], 'w') as fo1_am,\
         open(fileouts[0][7], 'w') as fo1_um,\
         open(fileouts[1][0], 'w') as fo2_F1,\
         open(fileouts[1][1], 'w') as fo2_F2,\
         open(fileouts[1][2], 'w') as fo2_F3,\
         open(fileouts[1][3], 'w') as fo2_F4,\
         open(fileouts[1][4], 'w') as fo2_F5,\
         open(fileouts[1][5], 'w') as fo2_F6,\
         open(fileouts[1][6], 'w') as fo2_am,\
         open(fileouts[1][7], 'w') as fo2_um:

        # Collect the file handles
        file_handles = ((fo1_F1, fo1_F2, fo1_F3, fo1_F4, fo1_F5, fo1_F6),
                        (fo2_F1, fo2_F2, fo2_F3, fo2_F4, fo2_F5, fo2_F6))

        for reads in pair_generator(bamfile):

            # Set read1 and read2
            (read1, read2) = reads

            # If it is unmapped or unpaired, keep in dumpyard
            if (read1.is_unmapped or read2.is_unmapped or
                (not read1.is_proper_pair) or (not read2.is_proper_pair)):
                n_unmapped += 1
                write_read_pair(reads,
                                ((0, read1.rlen), (0, read2.rlen)),
                                (fo1_um, fo2_um))
                continue

            # Get the extremes of the good CIGARs (we do some trimming to
            # avoid reading back into the adapters and stuff)
            ranges_read = []
            ranges_ref = []
            for read in reads:
                # Get the range of the good chunk, but with lenient criteria
                # (we want to keep the primers, but HXB2 might be quite distant
                # from our sequences hence generate lots of indels). These criteria
                # must include a bit of trimming, otherwise we read back into the
                # first 1-2 bases of the adapter (this happens only for short inserts).
                trim_both = 5 * (np.abs(read.isize) < n_cycles / 2 + 20)
                (rread, rref) = get_range_good_cigars(read.cigar, read.pos,
                                                      match_len_min=10,
                                                      trim_left=trim_both,
                                                      trim_right=trim_both)
                ranges_read.append(rread)
                ranges_ref.append(rref)

            # If no good CIGARs found in one of the reads, dump the pair
            if None in ranges_ref:
                n_unmapped += 1
                write_read_pair(reads,
                                ((0, read1.rlen), (0, read2.rlen)),
                                (fo1_um, fo2_um))
                continue

            # Find the fragment(s) compatible with each single read
            frags_reads = []
            for (start_ref, end_ref) in ranges_ref:
                ind = ((start_ref >= frags_pos[0]) &
                       (start_ref < frags_pos[1]) &
                       (end_ref > frags_pos[0]) &
                       (end_ref <= frags_pos[1]))
                frags_reads.append(ind.nonzero()[0])
            
            # Find the fragment(s) good for the whole pair
            frags_pair = np.intersect1d(*frags_reads)

            # A few cases can happen
            # 1. If the intersection is a single fragment, good
            if len(frags_pair) == 1:
                n_mapped += 1
                n_frag = frags_pair[0]
                write_read_pair(reads,
                                ranges_read,
                                (file_handles[0][n_frag],
                                 file_handles[1][n_frag]))

            # 2. If 2+ fragments are possible (tie), put into a special bucket
            # (essentially excluded, because we want two independent measurements
            # in the overlapping region, but we might want to recover them)
            elif len(frags_pair) > 1:
                n_ambiguous += 1
                write_read_pair(reads,
                                ranges_read,
                                (fo1_am, fo2_am))
            
            # 3. If no fragments are possible (e.g. one read crosses the
            # fragment boundary, they map to different fragments), dump it
            else:
                n_unmapped += 1
                write_read_pair(reads,
                                ((0, read1.rlen), (0, read2.rlen)),
                                (fo1_um, fo2_um))



    if VERBOSE >= 3:
        print 'Mapped: '+str(n_mapped)+', ambiguous: '+str(n_ambiguous)+\
                ', unmapped: '+str(n_unmapped)



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
    parser.add_argument('--subsample', action='store_true',
                        help='Apply only to a subsample of the reads')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--bwa', action='store_true',
                        help='Use BWA for premapping?')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use for mapping')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    VERBOSE = args.verbose
    subsample = args.subsample
    submit = args.submit
    bwa = args.bwa
    threads = args.threads

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # Set the number of cycles of the kit (for trimming adapters in short inserts)
    n_cycles = dataset['n_cycles']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']

    # Make index and hash files for HXB2, using a cropped reference to
    # reduce the problems of LTRs
    if bwa:
        make_bwa_hash(data_folder, cropped=True, VERBOSE=VERBOSE)
    make_index_and_hash(data_folder, cropped=True, VERBOSE=VERBOSE)

    # Iterate over all adaIDs
    for adaID in adaIDs:

        # Submit to the cluster self if requested
        if submit:
            fork_self(miseq_run, adaID, VERBOSE=VERBOSE, subsample=subsample,
                      bwa=bwa, threads=threads)
            continue
        
        # Make output folders
        make_output_folders(data_folder, adaID, VERBOSE=VERBOSE)

        # Map roughly to HXB2
        if bwa:
            premap_bwa(data_folder, adaID,
                       VERBOSE=VERBOSE, subsample=subsample)
        premap_stampy(data_folder, adaID,
                      VERBOSE=VERBOSE, subsample=subsample, bwa=bwa,
                      threads=threads)

        # Set the primers for fragment 5
        F5_type = dataset['primerF5'][dataset['adapters'].index(adaID)]

        # Divide into fragments and unclassified
        divide_into_fragments(data_folder, adaID, F5_type,
                              VERBOSE=VERBOSE, subsample=subsample,
                              n_cycles=n_cycles)

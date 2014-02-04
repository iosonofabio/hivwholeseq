#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/14
content:    Map PacBio reads to NL4-3 using stampy or custom pipeline.
'''
# Modules
import os
import argparse
from itertools import izip
from Bio import SeqIO
from Bio.Seq import reverse_complement as revcom
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#from hivwholeseq_pacbio.test_align import trim_reference, align_overlap_seqan
import seqanpy as sap
import hivwholeseq_pacbio.seqan_module.seqanpy as sap2
from hivwholeseq_pacbio.samples import samples
from hivwholeseq_pacbio.datasets import data_folder_dict
from hivwholeseq_pacbio.filenames import get_premapped_file, \
        get_reference_premap_filename

# Globals



# Functions
def fork_self(seq_run, sample, maxreads=-1, reference='NL4-3', VERBOSE=0):
    '''Submit map script to the cluster'''
    import subprocess as sp
    
    import hivwholeseq_pacbio
    JOBDIR = hivwholeseq_pacbio.__path__[0].rstrip('/')+'/'
    JOBLOGOUT = JOBDIR+'logout/'
    JOBLOGERR = JOBDIR+'logerr/'

    if VERBOSE:
        print 'Forking to the cluster'

    JOBSCRIPT = JOBDIR+'premap_to_reference.py'
    cluster_time = '23:59:59'
    vmem = '8G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'map pb '+sample,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--sample', sample,
                 '--maxreads', maxreads,
                 '--reference', reference,
                 '--verbose', VERBOSE,
                ]
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def make_output_folders(data_folder, samplename, VERBOSE=0):
    '''Make output folders'''
    from hivwholeseq.generic_utils import mkdirs
    outfiles = [get_premapped_file(data_folder, samplename)]
    for outfile in outfiles:
        dirname = os.path.dirname(outfile)
        mkdirs(dirname)
        if VERBOSE:
            print 'Folder created:', dirname


def trim_reference(refm, readm, band=100, VERBOSE=0):
    '''Trim the reference around a seed of the read'''

    seed_len = 30
    matches_min = 25
    rl = len(readm)

    # Try out a few seeds: middle of the read, a bit earlier, a bit later
    seed_starts = [max(0, rl / 2 - seed_len / 2),
                   max(0, rl / 2 - 5 * seed_len),
                   min(rl - seed_len, rl / 2 + 5 * seed_len)]

    for seed_start in seed_starts:
        seed = readm[seed_start: seed_start + seed_len]
        mismatches = 0

        # Scan the reference for incomplete matches of the seed
        ####################################################################
        ## NOTE: this is the algorithm that takes long - might be optimized
        #pos_seed = -1
        #score_old = -1
        #for i in xrange(len(refm) - (rl - seed_start) + band / 2):
        #    score = (seed == refm[i: i + seed_len]).sum()
        #    if score == seed_len:
        #        pos_seed = i
        #        break
        #    elif score > score_old:
        #        score_old = score
        #        pos_seed = i
        ####################################################################
        (pos_seed, score) = sap2.find_seed(refm.tostring(), seed.tostring())

        if score >= matches_min:
            if VERBOSE >= 3:
                print 'Position:', pos_seed, 'score:', score
        
            pos_trim_start = pos_seed - seed_start - band / 2
            pos_trim_end = pos_seed + (rl - seed_start) + band / 2
            ref_trim = refm[max(0, pos_trim_start): min(len(refm), pos_trim_end)]

            if pos_trim_end > len(ref_trim):
                ref_trim = np.concatenate([ref_trim, np.repeat('-', pos_trim_end - len(ref_trim))])
            if pos_trim_start < 0:
                ref_trim = np.concatenate([np.repeat('-', -pos_trim_start), ref_trim])

            return (pos_trim_start, ref_trim)
    
    raise ValueError('Seed not found in reference')


def map_read(refm, readm, qual, qname, band=100, VERBOSE=0):
    '''Make a mapped read out of an oerlap alignment with the reference'''

    if len(refm) > len(readm) + 2 * band:

        # Guess the strand by the A bias, or else try the other one
        if (readm == 'A').mean() > (readm == 'T').mean():
            is_reverse = False
        else:
            readm = np.fromstring(revcom(readm.tostring()), 'S1')
            qual = qual[::-1]
            is_reverse = True

        try:
            (pos_ref, refm) = trim_reference(refm, readm, band=band, VERBOSE=VERBOSE)

        except ValueError:
            try:
                readm = np.fromstring(revcom(readm.tostring()), 'S1')
                qual = qual[::-1]
                is_reverse = not is_reverse
                (pos_ref, refm) = trim_reference(refm, readm, band=band, VERBOSE=VERBOSE)

            except ValueError:

                # Unmapped: bring it back to its original strand
                if is_reverse:
                    readm = np.fromstring(revcom(readm.tostring()), 'S1')
                    qual = qual[::-1]

                readout = pysam.AlignedRead()
                readout.qname = qname
                readout.seq = readm.tostring()
                readout.flag = 0x0004
                readout.rname = 0
                readout.pos = 0
                readout.mapq = 0
                readout.cigar = []
                readout.mrnm = 0
                readout.mpos = 0
                readout.isize = 0
                readout.qual = qual
                readout.tags = ()

                return readout


    # Align pairwise ref and read
    (score, ref, read) = sap.align_overlap(refm.tostring(),
                                           readm.tostring(),
                                           band=band)

    # Trim alignment of side gaps (readm and qual are untouched)
    read = read.lstrip('-')
    offset = len(ref) - len(read)
    read = read.rstrip('-')
    ref = ref[offset: offset + len(read)]

    pos_ref += offset

    if VERBOSE >= 3:
        print ref[:50]
        print read[:50]
        print

    # Build cigar
    pos_ali = 0
    cigar = []
    lr = len(read)
    while pos_ali < lr:
        # Insertion
        if ref[pos_ali] == '-':
            bt = 1
            bl = lr - len(ref[pos_ali:].lstrip('-')) - pos_ali
    
        # Deletion
        elif read[pos_ali] == '-':
            bt = 2
            bl = lr - len(read[pos_ali:].lstrip('-')) - pos_ali

        # Match
        else:
            bt = 0

            # Find the end of the match block
            bl_ref = ref[pos_ali:].find('-')
            if bl_ref == -1:
                bl_ref = lr - pos_ali

            bl_read = read[pos_ali:].find('-')
            if bl_read == -1:
                bl_read = lr - pos_ali

            bl = min(bl_ref, bl_read)
        
        cigar.append((bt, bl))
        pos_ali += bl

    # A read starting with an insertion has the pos_ref shifted
    if len(cigar) and (cigar[0][0] == 1):
        pos_ref += cigar[0][1]

    if pos_ref < 0:
        import ipdb; ipdb.set_trace()

    readout = pysam.AlignedRead()
    readout.qname = qname
    readout.seq = readm.tostring()
    readout.qual = qual
    if is_reverse:
        readout.flag = 0x0010
    else:
        readout.flag = 0
    readout.rname = 0
    readout.pos = pos_ref
    readout.mapq = 99
    readout.cigar = cigar
    readout.mrnm = 0
    readout.mpos = 0
    readout.isize = sum(bl for (bt, bl) in cigar if bt in (0, 2))
    readout.tags = ()

    return readout



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--run', default='Upp23',
                        help='PacBio run to analyze (e.g. Upp23)')
    parser.add_argument('--sample', required=True,
                        help='Sample to analyze (e.g. S1)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to map')
    parser.add_argument('--reference', default='NL4-3',
                        help='Use alternative reference (the file must exist)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    seq_run = args.run
    samplename = args.sample
    VERBOSE = args.verbose
    maxreads = args.maxreads
    refname = args.reference
    submit = args.submit

    # Submit to the cluster if requested
    if submit:
        import sys
        fork_self(seq_run, samplename, maxreads=maxreads, reference=refname,
                  VERBOSE=VERBOSE)
        sys.exit(0)


    # Specify the dataset
    data_folder = data_folder_dict[seq_run]
    sample = samples.set_index('name').loc[samplename]

    # Parse raw circular consensus reads
    reads_iter = SeqIO.parse(data_folder+'ccs_reads/'+sample['filename']+\
                             '_ccs_reads.fastq.txt', 'fastq')
    
    # Get reference
    from hivwholeseq.reference import load_custom_reference
    refseq = load_custom_reference(refname)
    refm = np.array(refseq)

    # Make output folders
    make_output_folders(data_folder, samplename, VERBOSE=VERBOSE)

    # Copy mapping reference into folder
    SeqIO.write(refseq, get_reference_premap_filename(data_folder, samplename),
                'fasta')

    # Map reads
    reads_mapped = []
    bamfilename = get_premapped_file(data_folder, samplename)
    with pysam.Samfile(bamfilename, "wb",
                       referencenames=['HIV'],
                       referencelengths=[len(refm)]) as bamfile:
        for i, read in enumerate(reads_iter):
            if i == maxreads:
                break

            if VERBOSE >= 2:
                if not ((i + 1) % 10):
                    print (i + 1)

            readm = np.array(read)
            qual = (np.array(read.letter_annotations['phred_quality'],
                             np.int8) + ord('!')).tostring()
            
            read_mapped = map_read(refm, readm, qual, "read_"+str(i),
                                   VERBOSE=VERBOSE)
            
            bamfile.write(read_mapped)

            if VERBOSE >= 1:
                reads_mapped.append(read_mapped)

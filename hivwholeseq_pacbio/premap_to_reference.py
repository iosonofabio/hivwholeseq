#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/14
content:    Map PacBio reads to NL4-3 using stampy or custom pipeline.
'''
# Modules
import os
import gzip
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
from hivwholeseq_pacbio.samples import sample_table
from hivwholeseq_pacbio.datasets import PacBio_runs
from hivwholeseq_pacbio.filenames import get_premapped_file, \
        get_reference_premap_filename
from hivwholeseq.reference import load_custom_reference

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


def hash_set_sequence(seq, block_size):
    '''Hash the reference on the fly (note: __hash__ includes random stuff!)'''
    refs = ''.join(seq)
    h = {refs[i: i + block_size].__hash__(): i for i in xrange(len(refs) - block_size)}
    return h


def make_output_folders(data_folder, samplename, VERBOSE=0):
    '''Make output folders'''
    from hivwholeseq.generic_utils import mkdirs
    outfiles = [get_premapped_file(data_folder, samplename)]
    for outfile in outfiles:
        dirname = os.path.dirname(outfile)
        mkdirs(dirname)
        if VERBOSE:
            print 'Folder created:', dirname


def map_read(refm, ref_hash, readm, qual, qname, VERBOSE=0,
             trim_edges=True, band=100, hash_block_size=20, n_hash_matches=5,
             score_match=3, score_mismatch=-3, score_gapext=-1, score_gapopen=-7):
    '''Make a mapped read out of an oerlap alignment with the reference'''

    # Hash the read
    is_reverse = False
    read_hash = hash_set_sequence(readm, hash_block_size).items()
    np.random.shuffle(read_hash)
    n_matches = 0
    for h, pos_read in read_hash:
        if h in ref_hash:
            pos_ref = ref_hash[h]
            n_matches += 1
            if n_matches == n_hash_matches:
                break

    # if no hashes found, look in the opposite direction
    else:
        is_reverse = True
        readm = np.fromstring(revcom(readm.tostring()), 'S1')
        qual = qual[::-1]
        read_hash = hash_set_sequence(readm, hash_block_size).items()
        np.random.shuffle(read_hash)
        n_matches = 0
        for h, pos_read in read_hash:
            if h in ref_hash:
                pos_ref = ref_hash[h]
                n_matches += 1
                if n_matches == n_hash_matches:
                    break

        # If no hashes found, unmapped
        else:
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


    ref_start = max(0, pos_ref - pos_read - band)
    ref_end = min(len(refm), pos_ref + len(readm) - pos_read + band)
    ref_trimmed = refm[ref_start: ref_end]

    # Align pairwise ref and read
    (score, ref, read) = sap.align_overlap(ref_trimmed.tostring(),
                                           readm.tostring(),
                                           band=-1,
                                           score_match=score_match,
                                           score_mismatch=score_mismatch,
                                           score_gapext=score_gapext,
                                           score_gapopen=score_gapopen)

    # Trim alignment of side gaps (readm and qual are untouched)
    read = read.lstrip('-')
    offset = len(ref) - len(read)
    read = read.rstrip('-')
    ref = ref[offset: offset + len(read)]
    pos_ref = ref_start + offset

    if VERBOSE >= 3:
        print 'Alignment score:', score, 'of', score_match * len(ref),
        print 'pos: ['+str(pos_ref)+', '+str(pos_ref + len(ref.replace('-', '')))+']'
        print 'ref: ', ref[:40], '   ', ref[-40:]
        print 'read:', read[:40], '   ', read[-40:]
        print

    if trim_edges:
        if (pos_ref == 0) and (ref[0] == '-'):
            ref = ref.lstrip('-')
            offset = len(read) - len(ref)
            read = read[offset:]

        if (pos_ref + len(ref.replace('-', '')) == len(refm)) and (ref[-1] == '-'):
            ref = ref.rstrip('-')
            read = read[:len(ref)]

        if VERBOSE >= 3:
            print 'After edge trimming:'
            print 'ref: ', ref[:40], '   ', ref[-40:]
            print 'read:', read[:40], '   ', read[-40:]

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
    readout.tags = ()
    isize = sum(bl for (bt, bl) in cigar if bt in (0, 2))
    if is_reverse:
        readout.isize = -isize
    else:
        readout.isize = isize

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
    dataset = PacBio_runs[seq_run]
    data_folder = dataset['folder']
    sample = sample_table.set_index('name').loc[samplename]

    refseq = load_custom_reference(refname)
    refm = np.array(refseq)

    ref_hash = hash_set_sequence(refseq, block_size=20)

    make_output_folders(data_folder, samplename, VERBOSE=VERBOSE)

    fn_ref_out = get_reference_premap_filename(data_folder, samplename)
    SeqIO.write(refseq, fn_ref_out, 'fasta')

    fn_in = data_folder+'ccs_reads/'+sample['filename']+'_reads.fastq.gz'
    bamfilename = get_premapped_file(data_folder, samplename)
    with gzip.open(fn_in, 'rb') as f, \
         pysam.Samfile(bamfilename, "wb",
                       referencenames=['HIV'],
                       referencelengths=[len(refm)]) as bamfile:

        reads_iter = SeqIO.parse(f, 'fastq')
        reads_mapped = []

        for i, read in enumerate(reads_iter):
            if i == maxreads:
                break

            if VERBOSE >= 2:
                if not ((i + 1) % 10):
                    print (i + 1)

            readm = np.array(read)
            qual = (np.array(read.letter_annotations['phred_quality'],
                             np.int8) + ord('!')).tostring()
            
            read_mapped = map_read(refm, ref_hash, readm, qual, "read_"+str(i),
                                   VERBOSE=VERBOSE)
            
            bamfile.write(read_mapped)

            if VERBOSE >= 1:
                reads_mapped.append(read_mapped)

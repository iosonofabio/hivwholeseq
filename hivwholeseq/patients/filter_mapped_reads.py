#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       18/03/14
content:    Filter reads after mapping to the initial reference. If there are
            several sequencing runs contributing to one patient sample, merge them
            into the filtered reads.
'''
# Modules
import sys
import os
import argparse
import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO

from hivwholeseq.patients.patients import load_patient
from hivwholeseq.sequencing.filter_mapped_reads import plot_distance_histogram, \
        plot_distance_histogram_sliding_window, get_distance_from_consensus, \
        check_overhanging_reads, trim_bad_cigar
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_filter_mapped_init_summary_filename, \
        get_mapped_filtered_filename
from hivwholeseq.mapping_utils import convert_sam_to_bam, pair_generator
from hivwholeseq.fork_cluster import fork_filter_mapped_init as fork_self
from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.sequencing.samples import load_samples_sequenced as lss


# Functions
def filter_mapped_reads(sample, fragment,
                        PCR=1,
                        maxreads=-1,
                        VERBOSE=0,
                        n_cycles=600,
                        max_mismatches=100,
                        match_len_min=30,
                        trim_bad_cigars=3,
                        summary=True):
    '''Filter the reads to good chunks'''
    pname = sample.patient
    samplename_pat = sample.name
    samplenames_seq = sample['samples seq'].index.tolist()

    if VERBOSE >= 1:
        print 'Filtering reads:', pname, samplename_pat, fragment, PCR

    reffilename = get_initial_reference_filename(pname, fragment)
    refseq = SeqIO.read(reffilename, 'fasta')
    ref = np.array(refseq)

    outfilename = get_mapped_filtered_filename(pname, samplename_pat, fragment,
                                               type='bam', PCR=PCR)
    trashfilename = outfilename[:-4]+'_trashed.bam'

    infilenames = [get_mapped_to_initial_filename(pname, samplename_pat,
                                                 samplename, fragment,
                                                 type='bam', PCR=PCR)
                   for samplename in samplenames_seq]
    infilenames = filter(os.path.isfile, infilenames)
    if not len(infilenames):
        print ('WARNING: No mapped files found: '+', '.join([pname, samplename_pat,
                                                              fragment, str(PCR)]))
        return

    # Take reads evenly distributed across sequencing repetitions
    maxreads /= len(infilenames)

    if VERBOSE >= 2:
        print 'Input mapped filenames:',
        if len(infilenames) >= 2:
            print ''
        print '\n'.join(infilenames)

    # Use first file as template for the new bamfile
    infilename = infilenames[0]
    if not os.path.isfile(infilename):
        convert_sam_to_bam(infilename)
 
    with pysam.Samfile(infilename, 'rb') as bamfile:
        with pysam.Samfile(outfilename, 'wb', template=bamfile) as outfile,\
             pysam.Samfile(trashfilename, 'wb', template=bamfile) as trashfile:
 
            n_good = 0
            n_wrongname = 0
            n_unmapped = 0
            n_unpaired = 0
            n_mutator = 0
            n_badcigar = 0
            binsize = 200
            hist_distance_from_consensus = np.zeros(n_cycles + 1, int)
            hist_dist_along = np.zeros((len(ref) // binsize + 1, n_cycles + 1), int)

            # Iterate over input files, the first is already open
            for infilename in infilenames:

                if infilename != infilename[0]:
                    file_open = lambda: pysam.Samfile(infilename, 'rb')
                    file_close = lambda f: f.close()

                    if not os.path.isfile(infilename):
                        convert_sam_to_bam(infilename)

                else:
                    file_open = lambda: bamfile
                    file_close = lambda f: None

                try:
                    bamfile = file_open()
    
                    for irp, reads in enumerate(pair_generator(bamfile)):
                        if irp == maxreads:
                            break
                    
                        (read1, read2) = reads
                        i_fwd = reads[0].is_reverse
    
                        # Check a few things to make sure we are looking at paired reads
                        if read1.qname != read2.qname:
                            n_wrongname += 1
                            raise ValueError('Read pair '+str(irp)+': reads have different names!')
    
                        # Ignore unmapped reads
                        if read1.is_unmapped or read2.is_unmapped:
                            if VERBOSE >= 2:
                                print 'Read pair '+read1.qname+': unmapped'
                            n_unmapped += 1
                            map(trashfile.write, reads)
                            continue
                    
                        # Ignore not properly paired reads (this includes mates sitting on
                        # different fragments)
                        if (not read1.is_proper_pair) or (not read2.is_proper_pair):
                            if VERBOSE >= 2:
                                print 'Read pair '+read1.qname+': not properly paired'
                            n_unpaired += 1
                            map(trashfile.write, reads)
                            continue
                            
                        # Mismappings are often characterized by many mutations:
                        # check the number of mismatches and skip reads with too many
                        dc = get_distance_from_consensus(ref, reads, VERBOSE=VERBOSE)
                        hist_distance_from_consensus[dc.sum()] += 1
                        hbin = (reads[i_fwd].pos + reads[i_fwd].isize / 2) // binsize
                        hist_dist_along[hbin, dc.sum()] += 1
                        if (dc.sum() > max_mismatches):
                            if VERBOSE >= 2:
                                print 'Read pair '+read1.qname+': too many mismatches '+\
                                        '('+str(dc[0])+' + '+str(dc[1])+')'
                            n_mutator += 1
                            map(trashfile.write, reads)
                            continue
    
                        # Trim the bad CIGARs from the sides, if there are any good ones
                        skip = trim_bad_cigar(reads, match_len_min=match_len_min,
                                               trim_left=trim_bad_cigars,
                                               trim_right=trim_bad_cigars)
                        if skip:
                            n_badcigar += 1
                            map(trashfile.write, reads)
                            continue
    
                        # TODO: we might want to incorporate some more stringent
                        # criterion here, to avoid short reads, cross-overhang, etc.
    
                        # Write the output
                        n_good += 1
                        map(outfile.write, reads)

                finally:
                    file_close(bamfile)

    if VERBOSE >= 1:
        print 'Read pairs: '
        print 'Good:', n_good
        print 'Unmapped:', n_unmapped
        print 'Unpaired:', n_unpaired
        print 'Many-mutations:', n_mutator
        print 'Bad CIGARs:', n_badcigar
        print

    if summary:
        sfn = get_filter_mapped_init_summary_filename(pname, samplename_pat, fragment, PCR=PCR)
        with open(sfn, 'a') as f:
            f.write('Filter results: pname '+pname+', '+samplename_pat+', '+fragment+'\n')
            f.write('Total:\t\t\t'+str(irp + 1)+'\n')
            f.write('Good:\t\t\t'+str(n_good)+'\n')
            f.write('Unmapped:\t\t'+str(n_unmapped)+'\n')
            f.write('Unpaired:\t\t'+str(n_unpaired)+'\n')
            f.write('Many-mutations:\t\t'+str(n_mutator)+'\n')
            f.write('Bad CIGARs:\t\t'+str(n_badcigar)+'\n')

        # FIXME: implement or not?
        #plot_distance_histogram(data_folder, adaID, frag_gen,
        #                        hist_distance_from_consensus,
        #                        savefig=True)

        #plot_distance_histogram_sliding_window(data_folder, adaID, frag_gen,
        #                                       len(ref),
        #                                       hist_dist_along,
        #                                       binsize=binsize,
        #                                       savefig=True)


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Filter mapped reads',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')
    parser.add_argument('--PCR', default='all',
                        help='PCR to analyze (1, 2, or all)')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    submit = args.submit
    VERBOSE = args.verbose
    n_pairs = args.maxreads
    summary = args.summary
    PCR = args.PCR

    # Collect all sequenced samples from patients
    samples_pat = lssp()
    if pnames is not None:
        samples_seq = []
        for pname in pnames:
            patient = load_patient(pname)
            patient.discard_nonsequenced_samples()
            for samplename_pat, sample_pat in patient.samples.iterrows():
                samples_seq.append(sample_pat['samples seq'])
        samples_seq = pd.concat(samples_seq)

    else:
        samples_seq = lss()
        ind = samples_pat.index.isin(samplenames)
        samplenames_pat = samples_pat.index[ind]
        samples_seq = samples_seq.loc[samples_seq['patient sample'].isin(samplenames_pat)]

    if PCR != 'all':
        samples_seq = samples_seq.loc[np.array(samples_seq['PCR'], int) == int(PCR)]

    samples_groups = samples_seq.groupby(['patient sample', 'PCR'])

    if VERBOSE >= 2:
        print 'samples', samples_groups.groups

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    for (samplename_pat, PCR), samplenames_seq in samples_groups.groups.iteritems():
        sample_pat = samples_pat.loc[samplename_pat].copy()
        samples_seq_group = samples_seq.loc[samples_seq.index.isin(samplenames_seq)]
        sample_pat.set_value('samples seq', samples_seq_group)
        pname = sample_pat.patient
        PCR = int(PCR)

        for fragment in fragments:
            if submit:
                fork_self(samplename_pat, fragment,
                          VERBOSE=VERBOSE,
                          n_pairs=n_pairs,
                          summary=summary)
                continue

            if summary:
                sfn = get_filter_mapped_init_summary_filename(pname, samplename_pat, fragment, PCR=PCR)
                with open(sfn, 'w') as f:
                    f.write('Call: python filter_mapped_reads.py'+\
                            ' --samples '+samplename_pat+\
                            ' --fragments '+fragment+\
                            ' --verbose '+str(VERBOSE))
                    if n_pairs != -1:
                        f.write(' --maxreads '+str(n_pairs))
                    f.write('\n')

            filter_mapped_reads(sample_pat, fragment,
                                PCR=PCR,
                                VERBOSE=VERBOSE, maxreads=n_pairs,
                                summary=summary)



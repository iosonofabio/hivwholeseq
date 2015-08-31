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
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.sequencing.samples import load_samples_sequenced as lss
from hivwholeseq.sequencing.filter_mapped_reads import plot_distance_histogram, \
        plot_distance_histogram_sliding_window, get_distance_from_consensus, \
        check_overhanging_reads
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_filter_mapped_init_summary_filename, \
        get_mapped_filtered_filename
from hivwholeseq.utils.mapping import convert_sam_to_bam, pair_generator
from hivwholeseq.cluster.fork_cluster import fork_filter_mapped_init as fork_self


# Functions
def filter_read_pair(reads,
                     ref,
                     hist_distance_from_consensus=None,
                     hist_dist_along=None,
                     binsize=None,
                     max_mismatches=100,
                     match_len_min=30,
                     trim_bad_cigars=3,
                     VERBOSE=0):
    '''Filter read pair'''
    from hivwholeseq.utils.mapping import trim_short_cigars_pair

    (read1, read2) = reads

    # Check names to make sure we are looking at paired reads, this would
    # screw up the whole bamfile
    if read1.qname != read2.qname:
        n_wrongname += 1
        raise ValueError('Read pair '+str(irp)+': reads have different names!')

    # Ignore unmapped reads
    if read1.is_unmapped or read2.is_unmapped:
        if VERBOSE >= 2:
            print 'Read pair '+read1.qname+': unmapped'
        return 'unmapped'

    # Ignore not properly paired reads (this includes mates sitting on
    # different fragments)
    if (not read1.is_proper_pair) or (not read2.is_proper_pair):
        if VERBOSE >= 2:
            print 'Read pair '+read1.qname+': not properly paired'
        return 'unpaired'
        
    i_fwd = reads[0].is_reverse
    i_rev = not i_fwd
    readf = reads[i_fwd]
    readr = reads[i_rev]

    # Mismappings are often characterized by many mutations:
    # check the number of mismatches and skip reads with too many
    dc = get_distance_from_consensus(ref, reads, VERBOSE=VERBOSE)
    
    if hist_distance_from_consensus is not None:
        hist_distance_from_consensus[dc.sum()] += 1

    if hist_dist_along is not None:
        hbin = (readf.pos + readf.isize / 2) // binsize
        hist_dist_along[hbin, dc.sum()] += 1

    if (dc.sum() > max_mismatches):
        if VERBOSE >= 2:
            print 'Read pair '+read1.qname+': too many mismatches '+\
                    '('+str(dc[0])+' + '+str(dc[1])+')'
        return 'mutator'

    # Trim the bad CIGARs from the sides, if there are any good ones
    # FIXME: this must have a bug with leading insertions
    # NOTE: I rewrote the function, now simpler, it should work
    skip = trim_short_cigars_pair(reads, match_len_min=match_len_min,
                                  trim_pad=trim_bad_cigars, throw=False)
    if skip:
        return 'bad_cigar'

    # Check the reads are still long enough after trimming
    if (len(read1.seq) < 100):
        if VERBOSE >= 2:
            print 'Read too short:', read1.qname, len(read1.seq)
        return 'tiny'
    
    if (len(read2.seq) < 100):
        if VERBOSE >= 2:
            print 'Read too short:', read2.qname, len(read2.seq)
        return 'tiny'

    # NOTE: cross-overhang and similar stuff should never happen, because we
    # filter only insert sizes > 400 after premapping. Nonetheless...
    if readf.isize < 300:
        if VERBOSE >= 2:
            print 'Insert too small:', readf.isize
        return 'tiny'

    return 'good'


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
    samplenames_seq = sample.samples_seq.index.tolist()

    if VERBOSE >= 1:
        print 'Filtering reads:', pname, samplename_pat, fragment, PCR

    reffilename = get_initial_reference_filename(pname, fragment)
    refseq = SeqIO.read(reffilename, 'fasta')
    ref = np.array(refseq)

    outfilename = get_mapped_filtered_filename(pname, samplename_pat, fragment,
                                               type='bam', PCR=PCR,
                                               decontaminated=False)
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
            n_tiny = 0
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

                        pair_type = filter_read_pair(reads, ref,
                                                     hist_distance_from_consensus,
                                                     hist_dist_along,
                                                     binsize,
                                                     max_mismatches=max_mismatches,
                                                     match_len_min=match_len_min,
                                                     trim_bad_cigars=trim_bad_cigars,
                                                     VERBOSE=VERBOSE)
                    
                        if pair_type == 'unmapped':
                            n_unmapped += 1
                            map(trashfile.write, reads)

                        elif pair_type == 'unpaired':
                            n_unpaired += 1
                            map(trashfile.write, reads)

                        elif pair_type == 'mutator':
                            n_mutator += 1
                            map(trashfile.write, reads)

                        elif pair_type == 'bad_cigar':
                            n_badcigar += 1
                            map(trashfile.write, reads)

                        elif pair_type == 'tiny':
                            n_tiny += 1
                            map(trashfile.write, reads)

                        else:
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
        print 'Tiny:', n_tiny
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
            f.write('Tiny:\t\t\t'+str(n_tiny)+'\n')



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Filter mapped reads',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map')
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
    parser.add_argument('--PCR', default='1',
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
                sample_pat = SamplePat(sample_pat)
                samples_seq.append(sample_pat.samples_seq)
        samples_seq = pd.concat(samples_seq)

    elif samplenames is not None:
        samples_seq = lss()
        ind = samples_pat.index.isin(samplenames)
        samplenames_pat = samples_pat.index[ind]
        samples_seq = samples_seq.loc[samples_seq['patient sample'].isin(samplenames_pat)]

    else:
        samples_seq = lss()
        samples_seq = samples_seq.loc[samples_seq['patient sample'].isin(samples_pat.index)]


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
        sample_pat.samples_seq = samples_seq_group
        pname = sample_pat.patient
        PCR = int(PCR)

        for fragment in fragments:
            if submit:
                fork_self(samplename_pat, fragment,
                          VERBOSE=VERBOSE,
                          n_pairs=n_pairs,
                          PCR=PCR,
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



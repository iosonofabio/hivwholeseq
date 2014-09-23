#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/09/14
content:    Take the mapped reads and check for cross-contamination between the
            samples, including from reference sequence such as 38304.
'''
# Modules
import os
import sys
import argparse
import numpy as np
from Bio import SeqIO

from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.patients.patients import SamplePat
from hivwholeseq.reference import load_custom_reference
from hivwholeseq.sequence_utils import pretty_print_pairwise_ali
from hivwholeseq.patients.filenames import get_decontaminate_summary_filename
from hivwholeseq.fork_cluster import fork_decontaminate_reads_patient as fork_self



# Globals
refnames = ['38304', '38540', 'LAI-III']
maj_contnames = ['12879', '18798', '6154']



# Functions
def get_number_reads_summary(fn):
    '''Get the number of reads from the summary file'''
    with open(fn, 'r') as f:
        for line in f:
            if line[:4] == 'Good':
                n_reads_good = 2 * int(line.rstrip('\n').split()[1])
            elif line[:12] == 'Contaminated':
                n_reads_cont = 2 * int(line.rstrip('\n').split()[1])
    return (n_reads_good, n_reads_cont)


def check_already_decontaminated(sample, fragment, PCR):
    '''Check from the summary and output file whether decontamination has been done'''
    import os
    from hivwholeseq.mapping_utils import get_number_reads

    if (fragment == 'F4') and (sample.name in maj_contnames):
        return True

    fn_sum = get_decontaminate_summary_filename(sample.patient, sample.name, fragment, PCR=PCR)
    fn_in = sample.get_mapped_filtered_filename(fragment, PCR=PCR, decontaminated=False)
    fn_out = sample.get_mapped_filtered_filename(fragment, PCR=PCR, decontaminated=True)

    if not all(map(os.path.isfile, [fn_sum, fn_in, fn_out])):
        print sample.patient, sample.name, fragment, PCR
        return False

    n_reads_in = get_number_reads(fn_in)
    n_reads_out = get_number_reads(fn_out)
    (n_reads_good, n_reads_cont) = get_number_reads_summary(fn_sum)

    if n_reads_out != n_reads_good:
        print sample.patient, sample.name, fragment, PCR, n_reads_in, n_reads_out, n_reads_good, n_reads_cont
        return False

    if n_reads_good < 0.5 * n_reads_in:
        print sample.patient, sample.name, fragment, PCR, n_reads_in, n_reads_out, n_reads_good, n_reads_cont
        return False

    return True


def trim_align_overlap(ali):
    '''Trim wings of an overlap alignment'''
    (ali1, ali2) = ali
    start = len(ali2) - len(ali2.lstrip('-'))
    end = len(ali2.rstrip('-'))
    ali1 = ali1[start: end]
    ali2 = ali2[start: end]
    return (ali1, ali2)


def filter_contamination(bamfilename, bamfilename_out, contseqs, samplename, VERBOSE=0,
                         deltascore_max_pure=24, deltascore_max_cont=60,
                         maxreads=-1,
                         **kwargs):
    '''Fish contaminated reads from mapped reads

    The function checks for a maximal distance to the expected consensus, and only
    if it's more than that it checks all other samples.
    
    Args:
      deltascore_max_pure (int): the maximal delta in alignment score to the 
                                 consensus to be considered pure
      deltascore_max_cont (int): the maximal delta in alignment score to any other
                                 sample to be considered a contamination
      **kwargs: passed down to the pairwise alignment function
    '''
    import pysam
    from collections import defaultdict
    from operator import itemgetter
    from seqanpy import align_overlap

    from hivwholeseq.mapping_utils import pair_generator

    if 'score_match' in kwargs:
        score_match = kwargs['score_match']
    else:
        score_match = 3

    contseqs = contseqs.copy()
    consseq = contseqs.pop(samplename)

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        with pysam.Samfile(bamfilename_out, 'wb', template=bamfile) as bamfileout:
            n_good = 0
            n_cont = defaultdict(int)

            for irp, reads in enumerate(pair_generator(bamfile)):
                if irp == maxreads:
                    break

                if VERBOSE >= 2:
                    if not ((irp + 1) % 100):
                        if not ((irp + 1) == 100):
                            sys.stdout.write('\x1b[1A')
                        print irp + 1

                for read in reads:

                    # Look for distance to the own consensus, it that's zero move on
                    alignments_read = {}
                    deltas_read = {}
                    (score, ali1, ali2) = align_overlap(consseq, read.seq, **kwargs)
                    (ali1, ali2) = trim_align_overlap((ali1, ali2))
                    scoremax = len(ali1) * score_match
                    delta_read = scoremax - score
                    deltas_read[samplename] = delta_read
                    alignments_read[samplename] = (ali1, ali2)
                    # NOTE: here we might want to soften up a bit, e.g. less
                    # than 3 changes or so
                    if delta_read <= deltascore_max_pure:
                        if VERBOSE >= 4:
                            print 'Read is very close to its own consensus', scoremax, score, delta_read
                            pretty_print_pairwise_ali([ali1, ali2], width=90, name1='ref', name2='read')
                        continue

                    # Otherwise, move on to all other sequences
                    for contname, contseq in contseqs.iteritems():
                        (score, ali1, ali2) = align_overlap(contseq, read.seq, **kwargs)
                        (ali1, ali2) = trim_align_overlap((ali1, ali2))
                        scoremax = len(ali1) * score_match
                        delta_read = scoremax - score
                        deltas_read[contname] = delta_read
                        alignments_read[contname] = (ali1, ali2)

                    if VERBOSE >= 5:
                        print samplename
                        for key, d in deltas_read.iteritems():
                            print key, d
                        
                    (contname, delta_read) = min(deltas_read.iteritems(), key=itemgetter(1))

                    # Again, the correct consensus has precedence
                    if deltas_read[samplename] == delta_read:
                        contname = samplename

                    (ali1, ali2) = alignments_read[contname]

                    # The read may be closest to its own consensus, if not very close
                    if contname == samplename:
                        if VERBOSE >= 4:
                            print 'Read is closest to its consensus', scoremax, score, delta_read
                            pretty_print_pairwise_ali([ali1, ali2], width=90, name1='ref', name2='read')

                    # The read may come from another consensus, but it may be a bit
                    # distant from its own consensus
                    elif (delta_read <= deltascore_max_cont):
                        n_cont[contname] += 1

                        if VERBOSE >= 2:
                            print 'Good:', n_good, 'cont:', n_cont

                        if VERBOSE >= 3:
                            print 'Read is contaminated by', contname, scoremax, score, delta_read
                            pretty_print_pairwise_ali([ali1, ali2], width=90, name1='ref', name2='read')

                        break

                    # Finally, the read might be closest to somebody else, but not really
                    # close to anything... FIXME: for now accept those
                    else:
                        if VERBOSE >= 4:
                            print 'Read is close to nothing really', scoremax, score, delta_read
                            pretty_print_pairwise_ali([ali1, ali2], width=90, name1='ref', name2='read')

                else:
                    n_good += 1
                    bamfileout.write(reads[0])
                    bamfileout.write(reads[1])

    n_cont = dict(n_cont)

    return (n_good, n_cont)




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Decontaminate reads',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--PCR', type=int, default=None,
                        help='Analyze only reads from this PCR (e.g. 1)')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    submit = args.submit
    VERBOSE = args.verbose
    maxreads = args.maxreads
    summary = args.summary
    PCR = args.PCR

    samples = lssp()
    if pnames is not None:
        samples_focal = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples_focal = samples.loc[samples.index.isin(samplenames)]
    else:
        samples_focal = samples

    if VERBOSE >= 2:
        print 'samples', samples_focal.index.tolist()

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    if submit:
        for fragment in fragments:
            for samplename, sample in samples_focal.iterrows():
                # Skip the three F4 samples that are majorly contaminated
                if (fragment == 'F4') and (samplename in maj_contnames):
                    continue
                
                sample = SamplePat(sample)
                if PCR is None:
                    PCRs_sample = (1, 2)
                else:
                    PCRs_sample = [PCR]
                for PCR_sample in PCRs_sample:
                    bamfilename = sample.get_mapped_filtered_filename(fragment, PCR=PCR_sample)
                    if not os.path.isfile(bamfilename):
                        continue
                
                    if check_already_decontaminated(sample, fragment, PCR_sample):
                        continue

                    fork_self(samplename, fragment, VERBOSE=VERBOSE, maxreads=maxreads,
                              summary=summary, PCR=PCR_sample)

        sys.exit()

    for fragment in fragments:
        consensi = {refname: ''.join(load_custom_reference(refname+'_'+fragment))
                    for refname in refnames}
        for samplename, sample in samples.iterrows():
            sample = SamplePat(sample)
            try:
                consensi[samplename] = sample.get_consensus(fragment)
            except IOError:
                continue

        # Some consensi are bogus and must be deleted
        if fragment == 'F4':
            for consname in maj_contnames:
                del consensi[consname]

        for samplename, sample in samples_focal.iterrows():
            # Those three samples in F4 must be treated separately 
            if (fragment == 'F4') and (samplename in maj_contnames):
                if VERBOSE:
                    print samplename, fragment, 'majorly contaminated sample, skipping'
                continue

            sample = SamplePat(sample)
            pname = sample.patient

            if PCR is None:
                PCRs_sample = (1, 2)
            else:
                PCRs_sample = [PCR]

            for PCR_sample in PCRs_sample:
                bamfilename = sample.get_mapped_filtered_filename(fragment, PCR=PCR_sample)
                if not os.path.isfile(bamfilename):
                    continue

                bamfilename_out = sample.get_mapped_filtered_filename(fragment,
                                                                      decontaminated=True,
                                                                      PCR=PCR_sample)

                # Exclude the same patient as potential contaminants
                consensi_sample = consensi.copy()
                for contname in consensi:
                    # Keep the other references
                    if contname not in samples.index:
                        continue

                    # Keep the same sample
                    if contname == samplename:
                        continue

                    if samples.loc[contname].patient == pname:
                        del consensi_sample[contname]

                print samplename,
                if VERBOSE >= 2:
                    print ''
                (n_good, n_cont) = filter_contamination(bamfilename, bamfilename_out,
                                                        consensi_sample, samplename,
                                                        VERBOSE=VERBOSE,
                                                        deltascore_max_cont=60,
                                                        maxreads=maxreads)

                if VERBOSE:
                    print 'good:', n_good, 'contaminated:', n_cont

                if summary:
                    sfn = get_decontaminate_summary_filename(pname, samplename, fragment,
                                                             PCR=PCR_sample)
                    with open(sfn, 'w') as f:
                        f.write('Call: python decontaminate_reads.py'+\
                                ' --samples '+samplename+\
                                ' --fragments '+fragment+\
                                ' --verbose '+str(VERBOSE))
                        if maxreads != -1:
                            f.write(' --maxreads '+str(maxreads))
                        f.write('\n')
                        f.write('Good: '+str(n_good)+'\n')
                        f.write('Contaminated: '+str(sum(n_cont.itervalues()))+'\n')
                        f.write('Contamination sources:\n')
                        for contname, n_conti in n_cont.iteritems():
                            f.write('{:<20s}'.format(contname)+' '+'{:>7d}'.format(n_conti)+'\n')

# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/10/13
content:    Study PCR-mediated recombination from the plasmid mixes.

            This script is used both for:
                - MIX1 PCR2 probes (Tue28, adaID TS18)
                - MIX1 PCR1 probes (Tue42, adaID N1-S1)
                - MIX2 PCR2 probes (Tue28, adaID TS19)
                - MIX2 PCR1 probes (Tue42, adaID N3-S3)
'''
# Modules
import argparse
import numpy as np
from collections import Counter
from itertools import izip
from operator import itemgetter, attrgetter
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment as MSA
import pysam
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alphal
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.filenames import get_consensus_filename, get_mapped_filename, \
        get_allele_counts_filename, get_coverage_filename
from hivwholeseq.mapping_utils import align_muscle, pair_generator
from hivwholeseq.sequencing.minor_allele_frequency import filter_nus
from hivwholeseq.sequencing.coverage_tuples import get_coverage_tuples


# Globals
mix_adaIDs = {'Tue28': {'mix1': 'TS18', 'mix2': 'TS19'},
              'Tue42': {'mix1': 'N1-S1', 'mix2': 'N3-S3'}}

mix1_references_adaIDs = [(0.5, 'TS2'), (0.5, 'TS4')]
mix2_references_adaIDs = [(0.045, 'TS2'), (0.95, 'TS4'), (0.005, 'TS7')]



# Functions
def align_consensi_mix1(seq_run, fragment):
    '''Align consensi for mix1'''

    dataset_mix = MiSeq_runs[seq_run]
    dataset_pure = MiSeq_runs['Tue28']

    # Get the consensus of the mix1, and align it with the two consensi of
    # the pure strains
    adaID_mix1 = mix_adaIDs[seq_run]['mix1']
    adaIDs = [adaID_mix1] + map(itemgetter(1), mix1_references_adaIDs)
    datasets = [dataset_mix] + ([dataset_pure] * len(mix1_references_adaIDs))
    consensi = []
    for (adaID, dataset) in izip(adaIDs, datasets):
        data_folder = dataset['folder']
        i = dataset['adapters'].index(adaID)
        sample = dataset['samples'][i]
        seq = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment),
                         'fasta')
        consensi.append((sample, seq))

    alignment = align_muscle(*(map(itemgetter(1), consensi)))
    resorted = []
    for i, adaID in enumerate(adaIDs):
        for seq in alignment:
            if str(adaID)+'_F' in seq.name:
                seq.name = consensi[i][0]
                resorted.append(seq)
                break
    alignment = MSA(resorted)
    return alignment


def check_consensus_mix1(seq_run, fragment, alignment=None):
    '''Check the consensus switch for mix1'''
    if alignment is None:
        alignment = align_consensi_mix1(seq_run, fragment)
    ali = np.array(alignment)

    # Look for polymorphisms
    ind_likea = (ali[0] == ali[1]) & (ali[0] == ali[2])
    ind_like1 = (ali[0] == ali[1]) & (ali[0] != ali[2])
    ind_like2 = (ali[0] != ali[1]) & (ali[0] == ali[2])
    ind_liken = (ali[0] != ali[1]) & (ali[0] != ali[2])
    print 'Mix1,', fragment
    print 'conserved:', ind_likea.sum()
    print 'like', alignment[1].name+' only:', ind_like1.sum()
    print 'like', alignment[2].name+' only:', ind_like2.sum()
    print 'like none:', ind_liken.sum()
    print 

    # It jumps depending on the fragment, but it's constant within a single
    # fragment (i.e. there is amplification bias, but no strong recombination)


def align_consensi_mix2(seq_run, fragment):
    '''Align consensi for mix2'''

    dataset_mix = MiSeq_runs[seq_run]
    dataset_pure = MiSeq_runs['Tue28']

    # Get the consensus of the mix1, and align it with the two consensi of
    # the pure strains
    adaID_mix2 = mix_adaIDs[seq_run]['mix2']
    adaIDs = [adaID_mix2] + map(itemgetter(1), mix2_references_adaIDs)
    datasets = [dataset_mix] + ([dataset_pure] * len(mix2_references_adaIDs))

    consensi = []
    for (adaID, dataset) in izip(adaIDs, datasets):
        data_folder = dataset['folder']
        i = dataset['adapters'].index(adaID)
        sample = dataset['samples'][i]
        seq = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment),
                         'fasta')
        consensi.append((sample, seq))

    alignment = align_muscle(*(map(itemgetter(1), consensi)))
    resorted = []
    for i, adaID in enumerate(adaIDs):
        for seq in alignment:
            if str(adaID)+'_F' in seq.name:
                seq.name = consensi[i][0]
                resorted.append(seq)
                break
    alignment = MSA(resorted)
    return alignment


def check_consensus_mix2(seq_run, fragment, alignment=None):
    '''Check the consensus switch for mix1'''
    if alignment is None:
        alignment = align_consensi_mix2(seq_run, fragment)
    ali = np.array(alignment)

    # Look for polymorphisms
    ind_likea = (ali[0] == ali[1]) & (ali[0] == ali[2]) & (ali[0] == ali[3])
    ind_like1 = (ali[0] == ali[1]) & (ali[0] != ali[2]) & (ali[0] != ali[3])
    ind_like2 = (ali[0] != ali[1]) & (ali[0] == ali[2]) & (ali[0] != ali[3])
    ind_like3 = (ali[0] != ali[1]) & (ali[0] != ali[2]) & (ali[0] == ali[3])
    ind_liken = (ali[0] != ali[1]) & (ali[0] != ali[2]) & (ali[0] != ali[3])
    print 'Mix2,', fragment
    print 'conserved:', ind_likea.sum()
    print 'like', alignment[1].name+' only:', ind_like1.sum()
    print 'like', alignment[2].name+' only:', ind_like2.sum()
    print 'like', alignment[3].name+' only:', ind_like3.sum()
    print 'like none:', ind_liken.sum()
    print 

    # It jumps depending on the fragment, but it's constant within a single
    # fragment (i.e. there is amplification bias, but no strong recombination)
    # EXCEPTION: F6 (but there's something wrong in that consensus!)


def get_SNPs_mix1(seq_run, fragment):
    '''Get the SNPs of mix1'''

    # Get the positions of private alleles for each of the two references
    alignment = align_consensi_mix1(seq_run, fragment)
    ali = np.array(alignment)
    poss = (ali[1] != ali[2]).nonzero()[0]
    alls = ali[1:, poss]
    
    # Ignore private indels
    ind = (ali[:, poss] != '-').all(axis=0)
    poss = poss[ind]
    alls = alls[:, ind]

    # The reads are in the coordinare of their own consensus, not of the
    # alignment, hence construct a translation
    # Note: this algorithm works assuming no '-' natively in the reference (it
    # should not happen if our consensus building has worked properly)
    poss_ref = np.array([p - (ali[0, :p] == '-').sum() for p in poss], int)

    return poss, poss_ref, alls


def count_cross_reads_mix1(seq_run, fragment, maxreads=100, markers_min=4, VERBOSE=0):
    '''Count the number of reads jmping from one haplotype to the other
    
    This function calculates also the coverage normalization, and the switches
    between distant markers should not be double counted (we want to observe the
    saturation).
    
    '''
    if markers_min < 2:
        raise ValueError('To spot recombnation you need at least 2 markers!')

    # Get SNP positions and alleles
    poss, poss_ref, alls = get_SNPs_mix1(seq_run, fragment)

    # Go to the reads and ask how often you switch and where
    reads_identity = {'faithful': 0, 'cross': 0, 'borderline': 0}
    switch_counts = Counter()
    cocoverage = Counter()

    # Open BAM file
    adaID_mix1 = mix_adaIDs[seq_run]['mix1']
    bamfilename = get_mapped_filename(data_folder, adaID_mix1, fragment, type='bam',
                                      filtered=True)
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        # Proceed INSERT BY INSERT
        # (many recombination events are between read_fwd and read_rev)
        for irc, reads in enumerate(pair_generator(bamfile)):
            if irc == maxreads:
                if VERBOSE:
                    print 'Max reads reached:', maxreads
                break
        
            if VERBOSE >= 3:
                if not ((irc +1 ) % 10000):
                    print irc+1

            # Get the range of each read
            inds = []
            for read in reads:            
                start = read.pos
                end = start + sum(bl for (bt, bl) in read.cigar if bt in (0, 2))
                ind = ((poss_ref >= start) & (poss_ref < end)).nonzero()[0]
                inds.append(ind)

            # Check whether the read pair cover at least a few markers
            if len(np.unique(np.concatenate(inds))) < markers_min:
                continue

            # Assign binary identity to all covered markers
            poss_pair = [[], []]
            alls_pair = [[], []]
            for i, read in enumerate(reads):
                # If no markers, skip
                if len(inds[i]) == 0:
                    continue

                # Get markers covered by this read
                pos_SNP_read = poss_ref[inds[i]]
                alls_read = alls[:, inds[i]]

                # Prepare output structure
                all_read_bin = -np.ones_like(pos_SNP_read)

                # CIGARs are clean, iterate over them
                i_marker = 0	# Current focal marker
                seq = read.seq
                pos_ref = read.pos
                pos_read = 0
                for (bt, bl) in read.cigar:
                    # Ignore inserts
                    if bt == 1:
                        pos_read += bl
    
                    # Ignore deletions
                    elif bt == 2:
                        pos_ref += bl
    
                    # Only matches remain
                    else:
                        # Is any marker in this block? (iterate over all of them)
                        done = False
                        while not done:
                            if pos_ref + bl > pos_SNP_read[i_marker]:
                                al = seq[pos_read + pos_SNP_read[i_marker] - pos_ref]
                                if al == alls_read[0, i_marker]:
                                    all_read_bin[i_marker] = 1
                                elif al == alls_read[1, i_marker]:
                                    all_read_bin[i_marker] = 2
    
                                i_marker += 1
                                if i_marker == len(all_read_bin):
                                    done = 'all'

                            # if we encounter no marker in this CIGAR block, continue
                            else:
                                done = True

                        # once last marker found, skip the rest of the read
                        if done == 'all':
                            break

                        pos_ref += bl
                        pos_read += bl

                # Third allele and N are not relevant
                ind_cov_read = all_read_bin != -1

                poss_pair[i] = pos_SNP_read[ind_cov_read]
                alls_pair[i] = all_read_bin[ind_cov_read]

            # Merge positions and alleles for the whole insert
            poss_pair = map(list, poss_pair)
            poss_ins_tmp = np.sort(np.unique(np.concatenate(poss_pair)))
            poss_ins = []
            alls_ins_bin = []
            for pos in poss_ins_tmp:
                # If covered by one read only, take that allele
                if (pos in poss_pair[0]) and (pos not in poss_pair[1]):
                    poss_ins.append(pos)
                    alls_ins_bin.append(alls_pair[0][poss_pair[0].index(pos)])
                elif (pos not in poss_pair[0]) and (pos in poss_pair[1]):
                    poss_ins.append(pos)
                    alls_ins_bin.append(alls_pair[1][poss_pair[1].index(pos)])

                # If covered by both, take only if they agree
                else:
                    all0 = alls_pair[0][poss_pair[0].index(pos)]
                    all1 = alls_pair[1][poss_pair[1].index(pos)]
                    if all0 == all1:
                        poss_ins.append(pos)
                        alls_ins_bin.append(all0)

            # If too few markers are covered, continue
            if len(poss_ins) < markers_min:
                continue

            # Check whether the pair is crossing over
            n_all1 = alls_ins_bin.count(1)
            n_all2 = alls_ins_bin.count(2)
            if (n_all1 == 0) or (n_all2 == 0):
                reads_identity['faithful'] += 1

            elif (n_all1 == 1) or (n_all2 == 1):
                reads_identity['borderline'] += 1

            elif (n_all1 >= 2) and (n_all2 >= 2):
                reads_identity['cross'] += 1
                if VERBOSE >= 4:
                    print ''.join(map(str, alls_ins_bin))


            # Set the switch counts and coverage, for all pairs
            for i, posi in enumerate(poss_ins):
                alli = alls_ins_bin[i]
                for j, posj in enumerate(poss_ins[:i]):
                    allj = alls_ins_bin[j]
                    cocoverage[(posj, posi)] += 1
                    if alli != allj:
                        switch_counts[(posj, posi)] += 1    

    return (reads_identity, switch_counts, cocoverage)



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Errors in HIV pure strains')
    parser.add_argument('--run', required=True,
                        help='MiSeq run to analyze (e.g. Tue28)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('-n', type=int, default=100,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    seq_run = args.run
    fragments = args.fragments
    n_reads = args.n
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # MIX1
    fig, axs = plt.subplots(2, 3, figsize=(19, 12))
    for ax, fragment in izip(axs.ravel(), fragments):
        check_consensus_mix1(seq_run, fragment)
        (reads_identity,
         switch_counts,
         cocoverage) = count_cross_reads_mix1(seq_run, fragment,
                                              maxreads=n_reads, VERBOSE=VERBOSE)

        print reads_identity

        # Calculate and plot switches VS distance between markers (asinh?)
        lengths = []
        n_switch = []
        pos_frag = []
        for pair in cocoverage:
            lengths.append(pair[1] - pair[0])
            n_switch.append(1.0 * switch_counts[pair] / cocoverage[pair])
            pos_frag.append(0.5 * (pair[0] + pair[1]))

        # Linear fit
        m = np.dot(lengths, n_switch) / np.dot(lengths, lengths)
        q = 0

        # Plot
        colors = [cm.jet(int(255.0 * p / np.max(pos_frag))) for p in pos_frag]
        ax.scatter(lengths, n_switch, s=50, c=colors, label='data (color like starting\npos in fragment)')
        ax.plot([0, np.max(lengths)], [q, q + m * np.max(lengths)], lw=2,
                 ls='--', c='k', label='Fit: r = '+'{:1.0e}'.format(m)+' / base')
        ax.set_xlabel('dist [bases]')
        ax.set_ylabel('N of crossover')
        #ax.legend(loc=2, fontsize=18)
        ax.set_title('Mix1, '+seq_run+', '+fragment, fontsize=20)
        ax.set_ylim(-0.1, 0.5)

        # Follow along the fragment
        #TODO


    #plt.tight_layout()
    plt.ion()
    plt.show()

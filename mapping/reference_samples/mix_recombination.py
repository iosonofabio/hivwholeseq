# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/10/13
content:    Study PCR-mediated recombination from the plasmid mixes.
'''
# Modules
import argparse
import numpy as np
from collections import Counter
from operator import itemgetter, attrgetter
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment as MSA
import pysam

from mapping.miseq import alphal
from mapping.datasets import MiSeq_runs
from mapping.filenames import get_consensus_filename, get_mapped_filename, \
        get_allele_counts_filename, get_coverage_filename
from mapping.mapping_utils import align_muscle, pair_generator
from mapping.minor_allele_frequency import filter_nus
from mapping.coverage_tuples import get_coverage_tuples


# Globals
mix1_references_adaIDs = [(0.5, 2), (0.5, 4)]
mix2_references_adaIDs = [(0.045, 2), (0.95, 4), (0.005, 7)]



# Functions
def align_consensi_mix1(fragment):
    '''Align consensi for mix1'''
    # Specify the dataset
    miseq_run = 28
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # Get the consensus of the mix1, and align it with the two consensi of
    # the pure strains
    adaID_mix1 = 18
    adaIDs = [adaID_mix1] + map(itemgetter(1), mix1_references_adaIDs)
    consensi = []
    for adaID in adaIDs:
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


def check_consensus_mix1(fragment):
    '''Check the consensus switch for mix1'''
    alignment = align_consensi_mix1(fragment)
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


def count_cross_reads_mix1(fragment, maxreads=100, markers_min=4):
    '''Count the number of reads jmping from one haplotype to the other'''
    if markers_min < 2:
        raise ValueError('To spot recombnation you need at least 2 markers!')

    # Get the positions of private alleles for each of the two references
    alignment = align_consensi_mix1(fragment)
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

    # Go to the reads and ask how often you switch and where
    reads_identity = {'faithful': 0, 'cross': 0, 'borderline': 0}
    switch_counts = Counter()

    # Open BAM file
    adaID = 18
    bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                      filtered=True)
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        for irc, reads in enumerate(pair_generator(bamfile)):
            if irc == maxreads:
                if VERBOSE:
                    print 'Max reads reached:', maxreads
                break

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
            haplo_pair = [[], []]
            for i, read in enumerate(reads):
                # If no markers, skip
                if len(inds[i]) == 0:
                    continue

                # Get markers covered by this read
                pos_SNP_read = poss_ref[inds[i]]
                alls_read = alls[:, inds[i]]

                # Prepare output structure
                haplo_read = -np.ones_like(pos_SNP_read)

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
                                    haplo_read[i_marker] = 1
                                elif al == alls_read[1, i_marker]:
                                    haplo_read[i_marker] = 2
    
                                i_marker += 1
                                if i_marker == len(haplo_read):
                                    done = 'all'

                            # if we encounter no marker in this CIGAR block, continue
                            else:
                                done = True

                        # once last marker found, skip the rest of the read
                        if done == 'all':
                            break

                        pos_ref += bl
                        pos_read += bl

                poss_pair[i] = pos_SNP_read
                haplo_pair[i] = haplo_read

            # Check whether the pair is crossing over
            haplo_all = np.concatenate(haplo_pair)
            if (1 in haplo_all) != (2 in haplo_all):
                reads_identity['faithful'] += 1
            elif ((haplo_all == 1).sum() >= markers_min // 2) and \
                 ((haplo_all == 2).sum() >= markers_min // 2):
                reads_identity['cross'] += 1
                if VERBOSE >= 3:
                    print map(list, haplo_pair)

                # Check the crossover points
                poss_pairl = map(list, poss_pair)
                poss_ins = np.sort(np.unique(np.concatenate(poss_pair)))
                haplo_ins = -np.ones_like(poss_ins, int)
                for i, pos in enumerate(poss_ins):
                    # If that marker is covered by only one of the reads, fine
                    if (pos in poss_pair[0]) and (pos not in poss_pair[1]):
                        haplo_ins[i] = haplo_pair[0][poss_pairl[0].index(pos)]
                    elif (pos not in poss_pair[0]) and (pos in poss_pair[1]):
                        haplo_ins[i] = haplo_pair[1][poss_pairl[1].index(pos)]

                    # or else, check whether the two reads agree
                    else:
                        h0 = haplo_pair[0][poss_pairl[0].index(pos)]
                        h1 = haplo_pair[1][poss_pairl[1].index(pos)]
                        if h0 == h1:
                            haplo_ins[i] = h0
                        # If they do not agree, skip marker (it stays at -1)

                ind = haplo_ins != -1
                poss_ins = poss_ins[ind]
                haplo_ins = haplo_ins[ind]

                # Find switches
                for i in xrange(len(haplo_ins) - 1):
                    if haplo_ins[i + 1] != haplo_ins[i]:
                        switch_counts[(poss_ins[i], poss_ins[i+1])] += 1


            else:
                reads_identity['borderline'] +=1
    

    return (reads_identity, switch_counts)


def check_consensus_mix2(fragments):
    '''Check the consensus switch for mix1'''

    # Specify the dataset
    miseq_run = 28
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # Study the mix2
    adaID_mix2 = 19
    for fragment in fragments:

        # Get the consensus of the mix1, and align it with the two consensi of
        # the pure strains
        adaIDs = [adaID_mix2] + map(itemgetter(1), mix2_references_adaIDs)
        consensi = []
        for adaID in adaIDs:
            i = dataset['adapters'].index(adaID)
            sample = dataset['samples'][i]
            seq = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment),
                             'fasta')
            consensi.append((sample, seq))

        alignment = align_muscle(*(map(itemgetter(1), consensi)))
        resorted = []
        for adaID in adaIDs:
            for seq in alignment:
                if str(adaID)+'_F' in seq.name:
                    resorted.append(seq)
                    break
        ali = np.array(resorted)

        # Look for polymorphisms
        ind_likea = (ali[0] == ali[1]) & (ali[0] == ali[2]) & (ali[0] == ali[3])
        ind_like1 = (ali[0] == ali[1]) & (ali[0] != ali[2]) & (ali[0] != ali[3])
        ind_like2 = (ali[0] != ali[1]) & (ali[0] == ali[2]) & (ali[0] != ali[3])
        ind_like3 = (ali[0] != ali[1]) & (ali[0] != ali[2]) & (ali[0] == ali[3])
        ind_liken = (ali[0] != ali[1]) & (ali[0] != ali[2]) & (ali[0] != ali[3])
        print 'Mix2,', fragment
        print 'conserved:', ind_likea.sum()
        print 'like', consensi[1][0]+' only:', ind_like1.sum()
        print 'like', consensi[2][0]+' only:', ind_like2.sum()
        print 'like', consensi[3][0]+' only:', ind_like3.sum()
        print 'like none:', ind_liken.sum()
        print 

        # It jumps depending on the fragment, but it's constant within a single
        # fragment (i.e. there is amplification bias, but no strong recombination)
        # EXCEPTION: F6 (but there's something wrong in that consensus!)



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Errors in HIV pure strains')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    fragments = args.fragments
    VERBOSE = args.verbose

    # Specify the dataset
    miseq_run = 28
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # MIX1
    for fragment in fragments:
        check_consensus_mix1(fragment)
        maxreads = 1000
        (reads_identity, switch_counts) = count_cross_reads_mix1(fragment,
                                                                maxreads=maxreads)

        mtuples = switch_counts.keys()
        coverage_tuples = get_coverage_tuples(data_folder, 18, fragment,
                                              mtuples,
                                              maxreads=maxreads, VERBOSE=VERBOSE)

        # Get the switch probabilities
        nu_cross = {tup: 1.0 * switch_counts[tup] / coverage_tuples[i]
                    for i, tup in enumerate(mtuples)}

        # Get the switch rate per base
        r_cross = {tup: nu / (np.max(tup) - np.min(tup))
                   for tup, nu in nu_cross.iteritems()}

        # Get an idea about the crossover rates
        cross = np.array(r_cross.values())

        print reads_identity
        print 'log10 (crossover rate [per base]):', (np.log10(cross)).mean(), \
            '+-', (np.log10(cross)).std()


# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/10/13
content:    Study PCR-mediated recombination from the plasmid mixes.
'''
# Modules
import argparse
import numpy as np
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


def check_reads_mix1(fragment):
    '''Check reads jumping from one haplotype to the other'''
    alignment = align_consensi_mix1(fragment)
    return alignment

        #import ipdb; ipdb.set_trace()


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
        alignment = check_reads_mix1(fragment)

        # Position of private alleles
        ali = np.array(alignment)
        indi1 = ((ali[0] == ali[1]) & (ali[0] != ali[2])).nonzero()[0]
        indi2 = ((ali[0] != ali[1]) & (ali[0] == ali[2])).nonzero()[0]

        # Get the allele frequencies of the mix
        adaID = 18
        counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
        coverage = np.load(get_coverage_filename(data_folder, adaID, fragment))
        nu_filtered = filter_nus(counts, coverage)

        # Walk along the genome and tag the allele frequencies
        nus = []
        pos_af = 0
        for pos, alleles in enumerate(ali.T):
            if alleles[0] == '-':
                continue
            
            if (pos in indi1) or (pos in indi2):
                nus.append((pos,
                            pos_af,
                            alleles,
                            (nu_filtered[alphal.index(alleles[1]), pos_af],
                             nu_filtered[alphal.index(alleles[2]), pos_af])))

            pos_af += 1

        ## Plot the allele frequencies
        #import matplotlib.pyplot as plt
        #x = np.array(map(itemgetter(0), nus))
        #y = np.array(map(itemgetter(3), nus))
        #plt.plot(x, y + 1e-6, label=('NL4-3', 'SF162'))
        #plt.yscale('log')
        #plt.xlabel('position [bases]')
        #plt.ylabel(r'$\nu$')
        #plt.title('Mix1, '+fragment)

        #plt.legend()

        #plt.ion()
        #plt.show()

        # Vectorize the positions and alleles, forgetting indels
        poss = np.array(map(itemgetter(1), nus))
        alls = np.array(map(itemgetter(2), nus))
        ind_notindel = -(alls == '-').any(axis=1)
        poss = poss[ind_notindel]
        alls = alls[ind_notindel]

        # Go to the reads and ask how often you switch
        # Open BAM file
        # Note: the reads should already be filtered of unmapped stuff at this point
        maxreads = 1000
        reads_identity = {'faithful': 0, 'cross': 0}
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
                    ind = ((poss >= start) & (poss < end)).nonzero()[0]
                    inds.append(ind)

                # Check whether the reads cover at least a few markers
                if len(np.unique(np.concatenate(inds))) < 4:
                    continue

                haplo_identities = []
                for i, read in enumerate(reads):
                    if len(inds[i]) == 0:
                        haplo_identities.append([])
                        continue
                    pos_SNP_read = poss[inds[i]]
                    alls_read = alls[inds[i]]

                    haplo = -np.ones_like(pos_SNP_read)

                    # CIGARs are clean
                    i_marker = 0
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
                            # Is the marker in this block?
                            done = False
                            while not done:
                                if pos_ref + bl > pos_SNP_read[i_marker]:
                                    seqa = seq[pos_read + pos_SNP_read[i_marker] - pos_ref]
                                    if seqa == alls_read[i_marker, 1]:
                                        haplo[i_marker] = 1
                                    elif seqa == alls_read[i_marker, 2]:
                                        haplo[i_marker] = 2
    
                                    i_marker += 1
                                    if i_marker == len(haplo):
                                        done = 'all'
                                else:
                                    done = True
                            if done == 'all':
                                break

                            pos_ref += bl
                            pos_read += bl

                    haplo_identities.append(haplo)

                ## FIXME
                #import ipdb; ipdb.set_trace()

                haplo_all = np.concatenate(haplo_identities)
                if (1 in haplo_all) != (2 in haplo_all):
                    reads_identity['faithful'] += 1
                elif ((haplo_all == 1).sum() >= 2) and ((haplo_all == 2).sum() >= 2):
                    reads_identity['cross'] += 1
                    if VERBOSE >= 3:
                        print map(list, haplo_identities)
        

        print reads_identity


        #check_consensus_mix2(fragments)

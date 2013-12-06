# vim: fdm=marker
'''
author:     Fabio Zanini
date:       13/11/13
content:    Study the in vitro recombinaiton happening during library preparation
            and sequencing (rather than previous PCRs).
'''
# Modules
import argparse
import numpy as np
from collections import Counter
from operator import itemgetter, attrgetter
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment as MSA
import pysam
from matplotlib import cm
import matplotlib.pyplot as plt

from mapping.miseq import alpha, alphal
from mapping.datasets import MiSeq_runs
from mapping.filenames import get_consensus_filename, get_mapped_filename, \
        get_allele_counts_filename, get_coverage_filename, get_divided_filenames
from mapping.mapping_utils import pair_generator



# Functions
def span_fragments(reads, fragpri_pos):
    '''Assign read pair to fragments'''
    i_fwd = reads[0].is_reverse
    i_rev = not i_fwd

    # Insert coordinates
    start_fwd = reads[i_fwd].pos
    end_rev = reads[i_fwd].pos + reads[i_fwd].isize

    # If short inserts are not trimmed yet, we easily have the fwd read going
    # beyond reads[i_fwd].pos + reads[i_fwd].isize. Ignore that stuff, it's being
    # trimmed later on. For assignment, check only the insert AS A WHOLE.
    frag_start = ((start_fwd >= fragpri_pos[0]) &
                  (start_fwd < fragpri_pos[1])).nonzero()[0]
    frag_end = ((end_rev > fragpri_pos[0]) &
                (end_rev <= fragpri_pos[1])).nonzero()[0]
    frags_pair = np.arange(frag_start, frag_end + 1)

    return frags_pair





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


    ##########################################################################
    # Look for cross-fragment reads
    ##########################################################################
    # Focus on mix1
    adaID = 18

    ## We are going to have a special file for these, but for now filter the unmapped
    from mapping.primer_info import primers_coordinates_HXB2_inner as pcis
    from mapping.primer_info import primers_coordinates_HXB2_outer as pcos
    from mapping.primer_info import primers_inner, primers_outer
    from mapping.trim_and_divide import test_outer_primer
    from mapping.trim_and_divide import assign_to_fragment
    unmapped_filename = get_divided_filenames(data_folder, adaID, fragments)[-2]
    crossmapped_filename = get_divided_filenames(data_folder, adaID, fragments)[-3]
    fragments = ['F1', 'F2', 'F3', 'F4', 'F5b', 'F6']
    # This structure contains the fragment coordinates as of
    # - outer primers (row 0)
    # - inner primers (row 1)
    # - trimmed of all primers (row 2)
    frags_pos = np.zeros((3, 2, len(fragments)), int)
    for i, fragment in enumerate(fragments):
        pci = pcis[fragment]
        pco = pcos[fragment]
        frags_pos[0, :, i] = (pco[0][0], pco[1][1])
        frags_pos[1, :, i] = (pci[0][0], pci[1][1])
        frags_pos[2, :, i] = (pci[0][1], pci[1][0])
    # Since the reference is cropped, subtract from the positions F1 start
    # Note: now we are in the reference of the CROPPED HXB2, and start from 0!
    frags_pos -= frags_pos[1].min()

    ## Make primers with masks for ambiguous nucleotides
    #pr_outs = []
    #for fragment in fragments:
    #    ptmps = primers_outer[fragment]
    #    for i, ptmp in enumerate(ptmps):
    #        ptmp = np.ma.array(list(ptmp), mask=[p not in alpha[:4] for p in ptmp])
    #        ptmps[i] = ptmp
    #    pr_outs.append(ptmps)

    #with pysam.Samfile(unmapped_filename, 'rb') as input_file:
    #    with pysam.Samfile(crossmapped_filename, 'wb', template=input_file) as output_file:
    #        for reads in pair_generator(input_file):
    #            # Truly unmapped stuff we do not want
    #            if reads[0].is_unmapped or (not reads[0].is_proper_pair):
    #                if VERBOSE >= 3:
    #                    print 'Read pair unmapped/unpaired/tiny:', reads[0].qname
    #                continue

    #            # Stuff from the outer primers we do not want
    #            if test_outer_primer(reads, pr_outs, frags_pos[1, 1, -1]):
    #                if VERBOSE >= 3:
    #                    print 'Read pair from outer primer:', reads[0].qname
    #                continue

    #            # Check fragment assignment
    #            frags_pair = assign_to_fragment(reads, frags_pos[1])
    #            if len(frags_pair) > 0:
    #                continue

    #            # The other reads are good
    #            output_file.write(reads[0])
    #            output_file.write(reads[1])

    # Get reference
    from mapping.reference import load_HXB2
    refseq = load_HXB2(cropped=True)

    # Get alignment of the two pure strains with HXB2 to get sequence coordinates
    alignment = AlignIO.read(data_folder+'ali_references_mix1.fasta', 'fasta')
    ali = np.array(alignment)

    # Get positions and alleles of polymorphisms
    poss = (ali[0] != ali[1]).nonzero()[0]
    poss = poss[(ali[:, poss] != '-').all(axis=0)]
    poss_ref = np.array([pos - (ali[2, :pos] == '-').sum() for pos in poss])
    alls = ali[:, poss]

    crossmapped_filename = get_divided_filenames(data_folder, adaID, fragments)[-3]
    isizes = []
    alls_reads = []
    with pysam.Samfile(crossmapped_filename, 'rb') as input_file:
        for reads in pair_generator(input_file):
            is_fwd = reads[0].is_reverse
            is_rev = not is_fwd

            # Fill alleles matrix
            alls_read = np.ma.masked_all(alls.shape[1], 'S1')
            for read in reads:
                pos_ref = read.pos
                pos_read = 0
                for (bt, bl) in read.cigar:
                    if bt == 1:
                        pos_read += bl
                    elif bt == 2:
                        ind = (poss_ref >= pos_ref) & (poss_ref < pos_ref + bl)
                        alls_read[ind] = '-'
                        pos_ref += bl
                    elif bt == 0:
                        ind = (poss_ref >= pos_ref) & (poss_ref < pos_ref + bl)
                        poss_block = poss_ref[ind] - pos_ref
                        alls_read[ind] = np.fromstring(read.seq[pos_read: pos_read + bl], 'S1')[poss_block]
                        pos_read += bl
                        pos_ref += bl
                    else:
                        raise ValueError('CIGAR code not expected')

            # TODO: one should check separately that the two reads in the insert
            # agree if overlappinga and not just overwrite
            alls_reads.append(alls_read)

            # Increment insert sizes histogram
            isizes.append(reads[is_fwd].isize)

            # Find fragments and print info
            frags_pair = assign_to_fragment(reads, frags_pos[1])
            frags_span_pair = span_fragments(reads, frags_pos[1])
            print frags_pair, frags_span_pair,

            print reads[is_fwd].pos,
            print reads[is_fwd].pos + sum(bl for (bt, bl) in reads[is_fwd].cigar
                                          if bt in (0, 2)),
            print reads[is_rev].pos,
            print reads[is_rev].pos + sum(bl for (bt, bl) in reads[is_rev].cigar
                                          if bt in (0, 2)),
            print reads[is_fwd].isize,

            print len(reads[is_fwd].cigar), len(reads[is_rev].cigar)


    # Check switching illegitimate inserts
    alls_reads = np.ma.array(alls_reads)
    ## Markers covered by at least one read
    #ind = (-alls_reads.mask.all(axis=0)).nonzero()[0]
    #poss = poss[ind]
    #poss_ref = poss_ref[ind]
    #alls = alls[:, ind]
    #alls_reads = alls_reads[:, ind]

    # Find switches
    alls_reads_bin = -np.ones(alls_reads.shape, int)
    alls_reads_bin[alls_reads == alls[0]] = 1
    alls_reads_bin[alls_reads == alls[1]] = 2
    switch_counts = Counter()
    cocoverage = Counter()
    for a in alls_reads_bin:
        ind = (a != -1).nonzero()[0]
        for ii, i in enumerate(ind):
            alli = a[i]
            for j in ind[:ii]:
                allj = a[j]
                cocoverage[(poss_ref[j], poss_ref[i])] += 1
                if alli != allj:
                    switch_counts[(poss_ref[j], poss_ref[i])] += 1

    # Plot switch counts
    lengths = []
    n_switch = []
    pos_genome = []
    for (pair, cov) in cocoverage.most_common(200):
        lengths.append(pair[1] - pair[0])
        n_switch.append(1.0 * switch_counts[pair] / cov)
        pos_genome.append(0.5 * (pair[1] + pair[0]))


    # Linear fit
    m = np.dot(lengths, n_switch) / np.dot(lengths, lengths)
    q = 0
        
    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    colors = [cm.jet(int(255.0 * p / np.max(pos_genome))) for p in pos_genome]
    ax.scatter(lengths, n_switch, s=50, c=colors,
               label='data (color like starting\npos in genome)')
    ax.plot([0, np.max(lengths)], [q, q + m * np.max(lengths)], lw=2,
             ls='--', c='k', label='Fit: r = '+'{:1.0e}'.format(m)+' / base')
    ax.set_xlabel('dist [bases]')
    ax.set_ylabel('N of crossover')
    ax.legend(loc=2, fontsize=16)
    ax.set_title('Mix1, illigitimate reads', fontsize=20)
    ax.set_ylim(-0.1, 0.5)

    plt.tight_layout()

    # Plot along the genome (in HXB2)
    locations = poss_ref
    rates = np.ma.masked_all(len(locations) - 1)
    for i, posi in enumerate(locations[:-1]):
        j = i+1
        posj = locations[j]
        cov = cocoverage[(posi, posj)]
        if cov > 50:
            rates[i] = 1.0 * switch_counts[(posi, posj)] / cov / (posj - posi)

    fig2, ax2 = plt.subplots(1, 1, figsize=(14, 6))
    ax2.plot(0.5 * (locations[1:] + locations[:-1]), rates, lw=2)
    ax2.scatter(0.5 * (locations[1:] + locations[:-1]), rates, s=50, edgecolor='none')
    ax2.set_xlabel('Position [bases]')
    ax2.set_ylabel(r'$\rho$', fontsize=20)

    # Plot overlapping regions
    # Get fragments
    F5_primer = dataset['primerF5'][dataset['adapters'].index(adaID)]
    fragments = ['F1', 'F2', 'F3', 'F4', F5_primer, 'F6']
    # This structure contains the inner primers coordinates in cropped ref
    from mapping.primer_info import primers_coordinates_HXB2_inner as pcis_HXB2
    from mapping.primer_info import primers_coordinates_HXB2_outer as pcos_HXB2
    frags_pos = np.zeros((2, len(fragments)), int)
    frags_pos_out = np.zeros((2, len(fragments)), int)
    for i, fragment in enumerate(fragments):
        pci = pcis_HXB2[fragment]
        frags_pos[:, i] = (pci[0][0], pci[1][1])
        pco = pcos_HXB2[fragment]
        frags_pos_out[:, i] = (pco[0][0], pco[1][1])
    start = frags_pos.min()
    frags_pos -= start
    frags_pos_out -= start
    # Inner primers
    if frags_pos is not None:
        for i, frag_pos in enumerate(frags_pos.T):
            ax2.plot(frag_pos, 2 * [-0.01 - 0.01 * (i % 2)],
                     c=cm.jet(int(255.0 * i / len(frags_pos.T))), lw=2)

    # Outer primers
    if frags_pos_out is not None:
        for i, frag_pos in enumerate(frags_pos_out.T):
            ax2.plot(frag_pos, 2 * [-0.03 - 0.01 * (i % 2)],
                     c=cm.jet(int(255.0 * i / len(frags_pos_out.T))), lw=2)

    ax2.set_xlim(-500, 9500)
    ax2.set_title('PCR recombination rate across HIV genome\nin illegitimate reads')


    plt.tight_layout()
    plt.ion()
    plt.show()

    ##########################################################################
    

# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/10/13
content:    Study PCR-mediated recombination from the plasmid mixes.

            This script is used for:
                - PCR1 HiFi (Tue48, adaID N4-S1)
                - PCR1 Taq (Tue48, adaID N5-S1)
                - PCR2 HiFi (Tue48, adaID N6-S1)
                - PCR2 Taq (Tue48, adaID N1-S3)
'''
# Modules
import argparse
import numpy as np
from collections import Counter, defaultdict
from itertools import izip
from operator import itemgetter, attrgetter
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment as MSA
import pysam
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alphal, alpha
from hivwholeseq.sequencing.samples import load_sample_sequenced
from hivwholeseq.sequencing.filenames import get_consensus_filename, get_mapped_filename, \
        get_allele_counts_filename, get_coverage_filename
from hivwholeseq.reference import load_custom_reference
from hivwholeseq.utils.mapping import align_muscle, pair_generator
from hivwholeseq.sequencing.minor_allele_frequency import filter_nus
from hivwholeseq.utils.sequence import expand_ambiguous_seq
from hivwholeseq.sequencing.filenames import get_allele_frequencies_filename



# Globals
mix_RNA_strains = ['HXB2'] # The other one we have to build ourselves... see below
# In addition, hopefully by LAI III they mean HXB2 (!)
ref2_pol = load_custom_reference('38540_pol')

# Some markers are bad (polymorphic in the RNA strains)!
bad_markers = defaultdict(list)
bad_markers.update({'N4-S1 F2': [947, 653, 668, 1208, 1802],
               'N5-S1 F2': [947, 653, 668, 1208, 1802],
               'N6-S1 F2': [923, 629, 644, 1184],
               'N1-S3 F2': [923, 629, 644, 1184],
               'N4-S1 F4': [768, 392],
               'N5-S1 F4': [768, 392],
               'N4-S1 F4': [716, 340],
               'N1-S3 F4': [716, 340],
              })




# Functions
def get_samplename(PCRtype):
    '''Get samplename from input arg'''
    samplename = 'RNA_mix_'+PCRtype[:4]+'_Taq'
    if 'Taq' not in PCRtype:
        samplename = samplename + 'HiFi'
    return samplename


def guess_second_reference(seq_run, adaID, fragment, ref1,
                           threshold=0.05):
    '''Guess second reference 38540'''
    dataset = MiSeq_runs[seq_run]
    samplename = dataset['samples'][dataset['adapters'].index(adaID)]
    fragments_sample = samples[samplename]['fragments']
    frag_spec = [f for f in fragments_sample if fragment in f][0]

    # Get the consensus of the mix, and align it with NL4-3
    cons = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment), 'fasta')
    for feat in ref1.features:
        if frag_spec in feat.id:
            ref_start = feat.location.nofuzzy_start
            ref_end = feat.location.nofuzzy_end
            ref1_trim = ref1[ref_start: ref_end]
            break

    ali = align_muscle(cons, ref1_trim, sort=True)

    # Trim primers
    ali = ali[:, len(ali[0]) - len(str(ali[0].seq).lstrip('-')): \
                 len(str(ali[0].seq).rstrip('-'))]

    # Scan the genome: any position with only 1 allele > 0.01 is conserved,
    # the other are either NL4-3 or 38540: store 38540 by exclusion
    # we miss insertions this way, but this is not important for this analysis
    ref2 = np.ma.masked_all(len(cons), 'S1')
    afs = np.load(get_allele_frequencies_filename(data_folder, adaID, fragment))

    pos_cons = 0
    for i in xrange(len(ali[0])):
        bcon = ali[0, i]
        bref1 = ali[1, i] 

        # Ignore inserts to consensus for now
        if bcon == '-':
            continue

        # If ref1 has a gap, the major allele must be ref2
        elif bref1 == '-':
            ref2[pos_cons] = alpha[np.argmax(afs[:, pos_cons])]
            pos_cons += 1

        # If none has a gap, check how many alleles above 0.01;
        # if one, take it (conserved betweem ref1 and ref2)
        elif (afs[:, pos_cons] >= threshold).sum() == 1:
            ref2[pos_cons] = alpha[np.argmax(afs[:, pos_cons])]
            pos_cons += 1

        # if more than one, take the first non-ref1
        else:
            all_ind = np.argsort(afs[:, pos_cons])[::-1]
            if alpha[all_ind[0]] != bref1:
                ref2[pos_cons] = alpha[all_ind[0]]
            else:
                ref2[pos_cons] = alpha[all_ind[1]]

            pos_cons += 1

    ref2.fill_value = 'N'
    ref2s = ref2.tostring()
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet.IUPAC import ambiguous_dna
    ref2_rec = SeqRecord(Seq(ref2s, ambiguous_dna),
                         id='38540_guessed',
                         name='38540_guessed',
                         description='')

    return ref2_rec


def get_annotated_references(mix_RNA_strains):
    '''Annotate the references with fragments'''
    from hivwholeseq.annotate_genomewide_consensus import annotate_sequence
    refs = []
    for name in mix_RNA_strains:
        ref = load_custom_reference(name)
        annotate_sequence(ref, features=['PCR primers'])
        ref.id = name
        refs.append(ref)

    return refs


def align_mix_references(seq_run, adaID, fragment, refs, VERBOSE=0):
    '''Align sequences for mix: ref1, ref2, and mix consensus'''

    dataset = MiSeq_runs[seq_run]
    samplename = dataset['samples'][dataset['adapters'].index(adaID)]
    fragments_sample = samples[samplename]['fragments']
    frag_spec = [f for f in fragments_sample if fragment in f][0]

    # Get the consensus of the mix, and align it with the two references
    fn = get_consensus_filename(data_folder, adaID, fragment)
    cons = SeqIO.read(fn, 'fasta')

    # Trim references to fragment
    if VERBOSE >= 2:
        print 'Trimming references...',
    refs_trim = []
    for ref in refs:
        # Sometimes the reference is short already (e.g. pol)
        if len(ref) < len(cons):
            refs_trim.append(ref)
            continue

        for feat in ref.features:
            if frag_spec in feat.id:
                ref_start = feat.location.nofuzzy_start
                ref_end = feat.location.nofuzzy_end
                ref = ref[ref_start: ref_end]
                break
        refs_trim.append(ref)
    if VERBOSE >= 2:
        print 'done.\nAligning...',
    alignment = align_muscle(cons, refs_trim[0], refs_trim[1], sort=True)
    if VERBOSE >= 2:
        print 'done.',

    # Trim to actual fragment, excluding primers
    # The consensus cannot have external gaps, but the references still can
    # if they do not cover a whole fraghment, e.g. pol
    mix_start = len(alignment[0]) - len(str(alignment[0].seq).lstrip('-'))
    mix_end = len(str(alignment[0].seq).rstrip('-'))
    alignment = alignment[:, mix_start: mix_end]

    return alignment


def check_consensus_alleles_mix(adaID, fragment, alignment, VERBOSE=0):
    '''Check the consensus switch for mix1'''
    ali = np.array(alignment)

    # Sometimes the ref2 is trimmed (e.g. only pol)
    if ali[2, 0] == '-':
        ref2_start = (ali[2] != '-').nonzero()[0][0]
        ref2_end = (ali[2] != '-').nonzero()[0][-1] + 1
        ali = ali[:, ref2_start: ref2_end]
    else:
        ref2_start = 0
        ref2_end = ali.shape[1]

    # Look for polymorphisms
    ind_likea = (ali[0] == ali[1]) & (ali[0] == ali[2])
    ind_like1 = (ali[0] == ali[1]) & (ali[0] != ali[2])
    ind_like2 = (ali[0] != ali[1]) & (ali[0] == ali[2])
    ind_liken = (ali[0] != ali[1]) & (ali[0] != ali[2])
    print 'Mix,', adaID, fragment
    print 'Region covered by all references:', ref2_start, 'to', ref2_end, \
            '(0,', str(alignment.get_alignment_length())+')'
    print 'conserved:', ind_likea.sum()
    print 'like', alignment[1].name+' only:', ind_like1.sum()
    print 'like', alignment[2].name+' only:', ind_like2.sum()
    print 'like none:', ind_liken.sum()
    print 

    # If both are present, show the crossover points
    if VERBOSE >= 2:
        if ind_like1.sum() and ind_like2.sum():
            poss_SNPs = (ind_like1 | ind_like2 | ind_liken).nonzero()[0]
            alls_SNPs = []
            for pos in poss_SNPs:
                if ind_like1[pos]:
                    alls_SNPs.append('1')
                elif ind_like2[pos]:
                    alls_SNPs.append('2')
                else:
                    alls_SNPs.append('N')
            for i in xrange(len(poss_SNPs) / 10):
                poss_SNPs_t = poss_SNPs[i * 10: (i+1) * 10]
                alls_SNPs_t = alls_SNPs[i * 10: (i+1) * 10]
                print '\t'.join(map(str, poss_SNPs_t + ref2_start))
                print '\t'.join(alls_SNPs_t)
                print


def get_SNPs_mix(alignment):
    '''Get the SNPs of mix1'''

    # Get the positions of private alleles for each of the two references
    ali = np.array(alignment)

    # Sometimes the ref2 is trimmed (e.g. only pol)
    if ali[2, 0] == '-':
        ref2_start = (ali[2] != '-').nonzero()[0][0]
        ref2_end = (ali[2] != '-').nonzero()[0][-1] + 1
        ali = ali[:, ref2_start: ref2_end]
    else:
        ref2_start = 0
        ref2_end = ali.shape[1]

    # Get the SNP positions, excluding indels
    poss = ((ali[1] != ali[2]) & (ali != '-').all(axis=0)).nonzero()[0]
    
    # Ignore ambiguous bases in addition to indels
    ind = -(((ali[:, poss] != 'A') & (ali[:, poss] != 'C') & \
             (ali[:, poss] != 'G') & (ali[:, poss] != 'T')).any(axis=0))
    poss = poss[ind]
    alls = ali[:, poss]

    # FIXME: ignore all positions for which the mix is like the less frequently
    # represented reference, because the reference might actually be wrong there
    # because of the cell culture passages!
    j_refM = int((alls[0] == alls[2]).sum() > (alls[0] == alls[1]).sum()) + 1
    j_refm = (not (j_refM - 1)) + 1
    ind = -(alls[0] == alls[j_refm])
    poss = poss[ind]
    alls = alls[:, ind]

    # Take only alleles of the two references
    alls = alls[1:]

    # The reads are in the coordinare of their own consensus, not of the
    # alignment, hence construct a translation
    # Note: this algorithm works assuming no '-' natively in the reference (it
    # should not happen if our consensus building has worked properly)
    poss += ref2_start
    poss_ref = np.array([p - (ali[0, :p] == '-').sum() for p in poss], int)

    return poss, poss_ref, alls


def get_proportions_reference_alleles(seq_run, adaID, fragment, alignment):
    '''Get the allele proportions of the first/second reference along the fragment'''
    _, poss_ref, alls = get_SNPs_mix(alignment)
    afs = np.load(get_allele_frequencies_filename(data_folder, adaID, fragment))

    af_frac = np.zeros((2, len(poss_ref)))
    for i, (pos_ref, alls_pos) in enumerate(izip(poss_ref, alls.T)):
        for j in xrange(2):
            af_frac[j, i] = afs[alphal.index(alls_pos[j]), pos_ref]

    return af_frac


def count_cross_reads_mix(seq_run, adaID, fragment, poss_ref, alls,
                          maxreads=100, markers_min=4, VERBOSE=0):
    '''Count the number of reads jmping from one haplotype to the other
    
    This function calculates also the coverage normalization, and the switches
    between distant markers should not be double counted (we want to observe the
    saturation).
    
    '''
    if markers_min < 2:
        raise ValueError('To spot recombnation you need at least 2 markers!')

    # Go to the reads and ask how often you switch and where
    reads_identity = {'faithful': 0, 'cross': 0, 'borderline': 0, 'few markers': 0}
    switch_counts = Counter()
    cocoverage = Counter()

    # Open BAM file
    bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
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
                reads_identity['few markers'] += 1
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
                reads_identity['few markers'] += 1
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
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('-n', type=int, default=100,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--PCRtype', default='PCR1',
                        help='Mix to study (PCR1, PCR2, PCR1Taq, PCR2Taq)')

    args = parser.parse_args()
    fragments = args.fragments
    n_reads = args.n
    VERBOSE = args.verbose

    # Get the samplename
    samplename = get_samplename(args.PCRtype)

    sample = load_sample_sequenced(samplename)
    adaID = sample.adapter

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Annotate refs
    # Try other refs for ref2
    #mix_RNA_strains.append('AY901967')
    mix_RNA_strains.append('HIV1_CON_2002_subtype_C')
    refs = get_annotated_references(mix_RNA_strains)

    for fragment in fragments:
        ref2 = guess_second_reference(seq_run, adaID, fragment, refs[0])

        # FIXME: Use the actual sequenced pol for the time being (1000 bp)
        use_pol = False
        if use_pol:
            if fragment != 'F2':
                raise ValueError('Pol only in F2')

            ali_pol = align_muscle(ref2, ref2_pol, sort=True)
            ali_pol = ali_pol[:, len(ali_pol[0]) - len(str(ali_pol[1].seq).lstrip('-')):\
                                 len(str(ali_pol[1].seq).rstrip('-'))]

            if VERBOSE >= 3:
                for i in xrange(len(ali_pol[0]) / 50):
                    alit1 = str(ali_pol[0, i * 50: (i+1) * 50].seq)
                    alit2 = str(ali_pol[1, i * 50: (i+1) * 50].seq)
                    alit12 = []
                    for (a1, a2) in izip(alit1, alit2):
                        if a1 == a2:
                            alit12.append(' ')
                        else:
                            alit12.append('x')
                    alit12 = ''.join(alit12)
                    print 'GUESS', alit1
                    print '     ', alit12
                    print 'POL  ', alit2
                    print

            # The sequenced pol contains 28 ambiguous sites, too many to check all
            # seqs out. If the guessed seq contains a nucleotide that is in the seq
            # of possible ones, set it
            ref2_pol_correct = []
            for i in xrange(len(ali_pol[0])):
                a1 = ali_pol[0, i]
                a2 = ali_pol[1, i]
                if (a2 not in ['A', 'C', 'G', 'T']) and (a1 in expand_ambiguous_seq([a2])):
                    ref2_pol_correct.append(a1)
                else:
                    ref2_pol_correct.append(a2)
            ref2_pol_correct = ''.join(ref2_pol_correct)

            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            from Bio.Alphabet.IUPAC import ambiguous_dna
            refs = [refs[0],
                    SeqRecord(Seq(ref2_pol_correct, ambiguous_dna),
                              id='38540_pol', name='38540_pol', description='')]

        ali = align_mix_references(seq_run, adaID, fragment, refs)
        alim = np.array(ali)

        AlignIO.write(ali, data_folder+'ali_mix_'+adaID+'_'+fragment+'.fasta', 'fasta')

        check_consensus_alleles_mix(adaID, fragment, ali, VERBOSE=VERBOSE)

        poss, poss_ref, alls = get_SNPs_mix(ali)

        proportions = get_proportions_reference_alleles(seq_run, adaID, fragment, ali)

        if use_pol:
            # Check the alleles at 10%, what are the proportions there?
            afs = np.load(get_allele_frequencies_filename(data_folder, adaID, fragment))
            from hivwholeseq.utils.one_site_statistics import get_minor_allele_frequencies
            allm, num = get_minor_allele_frequencies(afs, alpha=alpha)
            poss_ref_8percent = (num > 0.08).nonzero()[0]
            if VERBOSE >= 2:
                print 'Positions of minor allele > 8%:', poss_ref_8percent
            # Alignment with ambiguous nucleotides in the Sangered pol
            ali_nc = align_mix_references(seq_run, adaID, fragment, [refs[0], ref2_pol])
            ali_ncm = np.array(ali_nc)
            poss_ref_ambi = ((ali_ncm[2] != 'A') & (ali_ncm[2] != 'C') & (ali_ncm[2] != 'G') & \
                             (ali_ncm[2] != 'T') & (ali_ncm[2] != '-')).nonzero()[0]

        # Get rid of strange high-recomb points by excluding a few "bad" markers
        # NOTE: these coordinates seems to wotk not only for the pol seq, but also
        # for the subtype C consensus
        if use_pol or ('subtype_C' in mix_RNA_strains[1]):
            ind = np.arange(len(poss_ref))
            for pos_bad in bad_markers[adaID+' '+fragment]:
                ind = ind[poss_ref[ind] != pos_bad]
            poss_ref = poss_ref[ind]
            alls = alls[:, ind]
            proportions = proportions[:, ind]



        # Plot proportions
        fig, ax = plt.subplots(1, 1)
        ax.plot(poss_ref, proportions[0] + 1e-5, lw=2, c='b', label='LAI III')
        ax.plot(poss_ref, proportions[1] + 1e-5, lw=2, c='g', label='38540')
        ax.set_xlabel('Position [bp]')
        ax.set_ylabel('Allele frequency')
        ax.set_yscale('log')
        ax.set_ylim(1e-5, 1.5)
        ax.set_title(samplename+', '+fragment, fontsize=20)

        print '\n'.join(map(''.join, np.array(ali)[:, poss_ref]))

        (reads_identity,
         switch_counts,
         cocoverage) = count_cross_reads_mix(seq_run, adaID, fragment, poss_ref, alls,
                                              maxreads=n_reads, VERBOSE=VERBOSE)

        print reads_identity

        # Calculate and plot switches VS distance between markers (asinh?)
        lengths = []
        n_switch = []
        pos_frag = []
        pairs = np.array(cocoverage.keys())
        for pair in pairs:
            pair = tuple(pair)
            lengths.append(pair[1] - pair[0])
            n_switch.append(1.0 * switch_counts[pair] / cocoverage[pair])
            pos_frag.append(0.5 * (pair[0] + pair[1]))
        lengths = np.array(lengths)
        pos_frag = np.array(pos_frag)
        n_switch = np.array(n_switch)
        cocoverage_list = np.array([cocoverage[tuple(pair)] for pair in pairs])
        # Proportions of ref1/ref2 at the markers, to normalize by the number of visible chances
        props = []
        for pair in pairs:
            prop = proportions[:, (poss_ref >= pair[0]) & (poss_ref <= pair[1])].mean(axis=1)
            props.append(prop)
        props = np.array(props).T

        # Find out about outliers at high switch numbers and short distances
        outlier_ind = ((lengths < 400) & (n_switch > lengths * 0.001))


        # Linear fit
        m = np.dot(lengths, n_switch) / np.dot(lengths, lengths)
        q = 0

        # Plot
        fig, ax = plt.subplots(1, 1)
        lenmax = np.array([pos_frag, lengths]).sum(axis=0).max()
        colors = [cm.jet(int(255.0 * (p + l /2) / lenmax)) for (p, l) in izip(pos_frag, lengths)]
        ax.scatter(lengths, n_switch, s=50, c=colors, label='data (color like starting\npos in fragment)')
        ax.plot([0, np.max(lengths)], [q, q + m * np.max(lengths)], lw=2,
                 ls='--', c='k', label='Fit: r = '+'{:1.0e}'.format(m)+' / base')
        ax.set_xlabel('dist [bases]')
        ax.set_ylabel('N of crossover')
        ax.set_title(samplename+', '+fragment, fontsize=20)
        ax.set_ylim(-0.1, 0.5)
        ax.text(10, 0.3, 'r = '+'{:1.1e}'.format(m)+' per base')

        # Plot normalized
        mn = np.dot(lengths, n_switch / 2.0 / props.prod(axis=0)) / np.dot(lengths, lengths)
        fig, ax = plt.subplots(1, 1)
        lenmax = np.array([pos_frag, lengths]).sum(axis=0).max()
        colors = [cm.jet(int(255.0 * (p + l /2) / lenmax)) for (p, l) in izip(pos_frag, lengths)]
        ax.scatter(lengths, n_switch / 2.0 / props.prod(axis=0), s=50, c=colors, label='data (color like starting\npos in fragment)')
        ax.plot([0, np.max(lengths)], [q, q + mn * np.max(lengths)], lw=2,
                 ls='--', c='k', label='Fit: r = '+'{:1.0e}'.format(mn)+' / base')
        ax.set_xlabel('dist [bases]')
        ax.set_ylabel('N of crossover / chance')
        ax.set_title(samplename+', '+fragment, fontsize=20)
        ax.set_ylim(-0.1, 10)
        ax.text(400, 4, 'r = '+'{:1.1e}'.format(mn)+' per base')



    plt.tight_layout()
    plt.ion()
    plt.show()

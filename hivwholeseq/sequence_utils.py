# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/12/13
content:    Support module with sequence utility functions which for some reason
            are missing from Biopython.
'''
# Modules
from hivwholeseq.genome_info import genes as genes_all



# Functions
def expand_ambiguous_seq(seq, seqtype='DNA'):
    '''Expand an ambiguous seq into all possible unambiguous ones'''
    if seqtype == 'DNA':
        from Bio.Data.IUPACData import ambiguous_dna_values as ttable
    elif seqtype == 'RNA':
        from Bio.Data.IUPACData import ambiguous_rna_values as ttable

    # Make a list of possibilities for each character
    amb = (ttable[c] for c in seq)

    # Generate all combinations
    from itertools import product
    seqs = map(''.join, product(*amb))

    return seqs


def pretty_print_pairwise_ali(ali, name1='', name2='', width=50, len_name=10):
    '''Pretty print function for pairwise alignments'''
    from itertools import izip
    ali = map(''.join, ali)

    for i in xrange(len(ali[0]) / width + 1):
        ali1_t = ali[0][i * width: (i+1) * width]
        ali2_t = ali[1][i * width: (i+1) * width]
        match_t = []
        for (a1, a2) in izip(ali1_t, ali2_t):
            if a1 == a2:
                match_t.append(' ')
            else:
                match_t.append('x')
        match_t = ''.join(match_t)

        name1 = name1[:len_name]
        name2 = name2[:len_name]
        print ('{:<'+str(len_name)+'}').format(name1)+':', ali1_t
        print (' ' * (len_name + 1)), match_t
        print ('{:<'+str(len_name)+'}').format(name2)+':', ali2_t
        print


def correct_genbank_features_load(record):
    '''Prepare genbank file after loading with the correct features'''
    for feat in record.features:
        try:
            feat.id = feat.qualifiers['note'][-1]
        except KeyError, IndexError:
            pass


def correct_genbank_features_save(record, molecule='DNA'):
    '''Prepare genbank file before saving with the correct features'''
    # BUG: feature id is lost during write, fake it with the 'note' qualifier
    for feat in record.features:
        if feat.id != '<unknown id>':
            if 'note' not in feat.qualifiers:
                feat.qualifiers['note'] = []
            # If already there, ignore (this is acting IN PLACE!)
            elif feat.id == feat.qualifiers['note'][-1]:
                continue
            feat.qualifiers['note'].append(feat.id)

    # BUG: gb allows only names up to 16 chars
    if len(record.name) > 16:
        if '|' in record.name:
            record.name = record.name.split('|')[-1]
        record.name = record.name[:16]

    # Specify an alphabet explicitely
    from Bio.Alphabet.IUPAC import ambiguous_dna, ambiguous_rna, extended_protein
    if molecule == 'DNA':
        record.seq.alphabet = ambiguous_dna
    elif molecule == 'RNA':
        record.seq.alphabet = ambiguous_rna
    else:
        record.seq.alphabet = extended_protein


def annotate_sequence_genes(seq_rec, fragment='genomewide', genes=genes_all,
                            max_end_slippage=10, VERBOSE=0):
    '''Annotate sequence with genes, checking their consistency'''
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

    from hivwholeseq.patients.build_initial_reference import check_genes_consensus as cgc
    (gene_seqs, genes_good, gene_poss) = cgc(''.join(seq_rec), fragment, genes=genes,
                                             max_end_slippage=max_end_slippage,
                                             VERBOSE=VERBOSE)

    for gene in gene_seqs:
        gene_pos = gene_poss[gene]
        if len(gene_pos) == 1:
            location = FeatureLocation(*(gene_pos[0]))
        else:
            locations = [FeatureLocation(*exon_pos) for exon_pos in gene_pos]
            location = CompoundLocation(locations)

        feature = SeqFeature(location, type='gene', id=gene, strand=1)
        seq_rec.features.append(feature)


def find_seed_imperfect(seq, seed, threshold=0.7, VERBOSE=0):
    '''Imperfect match of a seed to a sequence'''
    import numpy as np

    seed = ''.join(seed)
    seq = ''.join(seq)
    pos = seq.find(seed)
    if pos != -1:
        return pos

    seed = np.fromstring(seed, 'S1')
    seq = np.fromstring(seq, 'S1')
    sl = len(seed)
    seql = len(seq)
    n_match = [(seed == seq[i: i + sl]).sum() for i in xrange(seql - sl)]
    pos = np.argmax(n_match)
    if n_match[pos] >= threshold * sl:
        return pos

    raise ValueError('Seed not found at specified threshold ('+str(threshold)+')')


def rfind_seed_imperfect(seq, seed, threshold=0.7, VERBOSE=0):
    '''Imperfect match of a seed to a sequence, from the right'''
    import numpy as np

    seed = ''.join(seed)
    seq = ''.join(seq)
    pos = seq.rfind(seed)
    if pos != -1:
        return pos

    seed = np.fromstring(seed, 'S1')
    seq = np.fromstring(seq, 'S1')
    sl = len(seed)
    seql = len(seq)
    n_match = [(seed == seq[i: i + sl]).sum() for i in xrange(seql - sl)]
    pos = len(n_match) - 1 - np.argmax(n_match[::-1])
    if n_match[pos] >= threshold * sl:
        return pos

    raise ValueError('Seed not found at specified threshold ('+str(threshold)+')')


def merge_sequences(seqs, minimal_fraction_match=0.75, skip_initial=30, VERBOSE=0):
    '''Merge sequences from subsequent fragments into a genomewide one'''
    import numpy as np

    seqs = map(''.join, seqs)

    if VERBOSE:
        print 'Start with F1'
    seqs_all = [seqs[0]]
    for ifr, seq in enumerate(seqs[1:], 2):
        if VERBOSE:
            print 'merging in F'+str(ifr),

        # Avoid the first few bases because of low coverage
        pos_seed_new = skip_initial
        seed = seq[pos_seed_new: pos_seed_new + 30]
        # Find seed from the right and allow imperfect matches
        pos_seed = seqs_all[-1].rfind(seed)
        if pos_seed != -1:
            overlap_found = True
        else:
            refm = np.fromstring(seqs_all[-1], 'S1')
            seed = np.fromstring(seed, 'S1')
            sl = len(seed)
            n_match = np.array([(refm[i: i + sl] == seed).sum()
                                for i in xrange(len(refm) - sl)], int)
            pos_seed = len(n_match) - 1 - np.argmax(n_match[::-1])
            if n_match[pos_seed] >= minimal_fraction_match * sl:
                overlap_found = True
            else:
                overlap_found = False
                if VERBOSE:
                    print 'WARNING: Cannot merge consensus F'+str(ifr)+\
                            ' with previous ones'

        if overlap_found:
            seqs_all[-1] = seqs_all[-1][:pos_seed]
            seqs_all.append(seq[pos_seed_new:])
            if VERBOSE:
                print 'merged, total length:', len(''.join(seqs_all))
        else:
            seqs_all.append(('N' * 10) + seq)

    return ''.join(seqs_all)


def build_local_consensus(seqs, VERBOSE=0, store_allele_counts=False, full_cover=True):
    '''Build a local consensus from an MSA
    
    There is only ONE tricky point: what to do if some reads do not cover the whole
    block, e.g. at the end of a fragment because of low coverage?
    If full_cover == False, convert MSA gaps at the end of too short reads into N

    Args:
      seqs (list of SeqRecords): seqs to build consensus from
      store_allele_counts (bool): return also allele counts from the alignment
      full_cover (bool): if True, assume the reads fully cover the region (no gaps at edges)
    '''

    import numpy as np
    from hivwholeseq.miseq import alpha
    from hivwholeseq.mapping_utils import align_muscle

    ali = np.array(align_muscle(*seqs, sort=True), 'S1', ndmin=2)
    if full_cover:
        allele_counts = np.array([(ali == a).sum(axis=0) for a in alpha], int, ndmin=2)
    else:
        allele_counts = np.zeros((len(alpha), len(ali[0])),int)
        for i in xrange(len(seqs)):
            if ali[i, -1] == '-':
                first_finalgap = len(ali[i].tostring().rstrip('-'))
                ali[i, first_finalgap:] = 'X'
            for ai, a in enumerate(alpha):
                allele_counts[ai] += ali[i] == a

        cov = allele_counts.sum(axis=0)
        allele_counts = allele_counts[:, cov > 0]

    cons_local = []
    for counts in allele_counts.T:
        # Pick max count nucleotide, ignoring N
        maxinds = (counts[:-1] == counts.max()).nonzero()[0]
        if len(maxinds) < 1:
            cons_local.append('-')
            continue
        # Pick a random nucleotide in case of a tie
        elif len(maxinds) > 1:
            np.random.shuffle(maxinds)
        maxind = maxinds[0]
        cons_local.append(alpha[maxind])
    cons_local = np.array(cons_local, 'S1')

    ind_nongap = cons_local != '-'
    cons_local = ''.join(cons_local[ind_nongap])
    
    if store_allele_counts:
        allele_counts = allele_counts[:, ind_nongap]
        return (cons_local, allele_counts)

    return cons_local



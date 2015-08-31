# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/12/13
content:    Support module with sequence utility functions which for some reason
            are missing from Biopython.
'''
# Modules
from numpy import array
from hivwholeseq.utils.genome_info import genes as genes_all


# Globals
# Alphabet of nucleotides
alphas = 'ACGT-N'
alphal = list(alphas)
alpha = array(alphal, 'S1')

# Alphabet of amino acids
alphaas = 'ACDEFGHIKLMNPQRSTVWY*-X'
alphaal = list(alphaas)
alphaa = array(alphaal, 'S1')



# Functions
def align_pairwise(seq1, seq2, method='global', **kwargs):
    '''Align two sequences pairwise'''
    import numpy as np
    from seqanpy import align_global, align_overlap, align_ladder, align_local
    funcd = {'global': align_global,
             'overlap': align_overlap,
             'ladder': align_ladder,
             'local': align_local,
            }

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio.Alphabet import single_letter_alphabet

    if isinstance(seq1, np.ndarray):
        seq1 = ''.join(seq1)
    if isinstance(seq1, basestring):
        seq1 = SeqRecord(Seq(seq1, single_letter_alphabet), id='seq1', name='seq1',
                         description='')

    if isinstance(seq2, np.ndarray):
        seq2 = ''.join(seq2)
    if isinstance(seq2, basestring):
        seq2 = SeqRecord(Seq(seq2, single_letter_alphabet), id='seq2', name='seq2',
                         description='')

    score, ali1, ali2 = funcd[method](seq1, seq2, **kwargs)
    ali1 = SeqRecord(Seq(ali1, seq1.seq.alphabet), id=seq1.id, name=seq1.name,
                     description=seq1.description) 
    ali2 = SeqRecord(Seq(ali2, seq2.seq.alphabet), id=seq2.id, name=seq2.name,
                     description=seq2.description) 
    ali = MultipleSeqAlignment([ali1, ali2])
    return ali


def align_muscle(*seqs, **kwargs):
    '''Global alignment of sequences via MUSCLE'''
    import subprocess as sp
    from Bio import AlignIO, SeqIO
    from Bio.Align.Applications import MuscleCommandline
    
    if not len(seqs):
        return None

    # Convert to SeqRecord if required
    if isinstance(seqs[0], basestring):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import single_letter_alphabet
        seqs = [SeqRecord(Seq(s, single_letter_alphabet),
                          id='seq'+str(i+1),
                          name='seq'+str(i+1),
                          description='seq'+str(i+1))
                for i, s in enumerate(seqs)]

    muscle_cline = MuscleCommandline(diags=True, quiet=True)
    child = sp.Popen(str(muscle_cline),
                     stdin=sp.PIPE,
                     stdout=sp.PIPE,
                     stderr=sp.PIPE,
                     shell=True)
    SeqIO.write(seqs, child.stdin, "fasta")
    child.stdin.close()
    align = AlignIO.read(child.stdout, "fasta")
    child.stderr.close()
    child.stdout.close()

    if ('sort' in kwargs) and kwargs['sort']:
        from Bio.Align import MultipleSeqAlignment as MSA
        alisort = []
        for seq in seqs:
            for row in align:
                if row.id == seq.id:
                    alisort.append(row)
                    break
        align = MSA(alisort)

    return align


def align_codon_pairwise(seqstr, refstr, **kwargs):
    '''Pairwise alignment via codons
    
    Parameters:
       **kwargs: passed down to SeqAn alignment function
    '''
    from Bio.Seq import translate
    from seqanpy import align_global
    from itertools import izip

    if len(seqstr) % 3:
        raise ValueError('The length of the first sequence is not a multiple of 3')
    elif len(refstr) % 3:
        raise ValueError('The length of the second sequence is not a multiple of 3')

    seqpr = translate(seqstr)
    refpr = translate(refstr)
    (score, alis, alir) = align_global(seqpr, refpr, **kwargs)
    aliseq = []
    aliref = []
    poss = 0
    posr = 0
    for aas, aar in izip(alis, alir):
        if aas == '-':
            aliseq.append('---')
        else:
            aliseq.append(seqstr[poss: poss+3])
            poss += 3

        if aar == '-':
            aliref.append('---')
        else:
            aliref.append(refstr[posr: posr+3])
            posr += 3

    aliseq = ''.join(aliseq)
    aliref = ''.join(aliref)

    return (aliseq, aliref)


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


def reduce_ambiguous_seqs(seqs, seqtype='DNA'):
    '''Get the most specific ambiguous DNA/RNA sequence of a set'''
    if seqtype == 'DNA':
        from Bio.Data.IUPACData import ambiguous_dna_values as ttable
    elif seqtype == 'RNA':
        from Bio.Data.IUPACData import ambiguous_rna_values as ttable
    del ttable['X']

    ttable_back = {frozenset(value): key for (key, value) in ttable.iteritems()}

    from itertools import imap, izip
    seq = ''.join(imap(ttable_back.get, imap(frozenset, izip(*seqs))))

    return seq


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

    from hivwholeseq.store.store_initial_reference import check_genes_consensus as cgc
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
    from hivwholeseq.utils.mapping import align_muscle

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


def build_msa_haplotypes(haploc, VERBOSE=0, label=''):
    '''Build multiple sequence alignment from cluster of haplotypes'''
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet.IUPAC import ambiguous_dna
    
    seqs = [SeqRecord(Seq(seq, ambiguous_dna),
                      id=label+'count_'+str(count)+'_rank_'+str(i),
                      name=label+'count_'+str(count)+'_rank_'+str(i),
                      description='')
            for i, (seq, count) in enumerate(haploc.most_common())]

    from hivwholeseq.utils.mapping import align_muscle
    ali = align_muscle(*seqs, sort=True)

    return ali


def get_coordinates_genomic_region(ref, region):
    '''Get coordinates of region in an annotated sequence'''
    for feature in ref.features:
        if region == feature.id:
            return feature.location

    raise ValueError('Region not found')


def translate_with_gaps(seq):
    '''Translate sequence with gaps'''
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq, translate
    from Bio.Alphabet.IUPAC import protein

    L = len(seq)
    if L % 3:
        raise ValueError('The sequence length is not a multiple of 3')

    seqstr = ''.join(seq)
    prot = []
    for i in xrange(L // 3):
        codon = seqstr[3 * i: 3 * (i+1)]
        if codon == '---':
            prot.append('-')
        elif '-' in codon:
            raise ValueError('Non-aligned gaps found')
        else:
            prot.append(''.join(translate(codon)))
    prot = ''.join(prot)

    # Output in various formats
    if isinstance(seq, basestring):
        return prot
    elif isinstance(seq, Seq):
        return Seq(prot, protein)
    elif isinstance(seq, SeqRecord):
        return SeqRecord(Seq(prot, protein), id=seq.id, name=seq.name,
                         description=seq.description)
    else:
        import numpy as np
        return np.fromstring(prot, 'S1')


def translate_alignment(ali_sub, VERBOSE=0):
    '''Translate multiple sequence alignment'''
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment as MSA

    L = ali_sub.get_alignment_length()
    if L % 3:
        raise ValueError('The alignment length is not a multiple of 3')

    prots = []
    for seq in ali_sub:
        prot = SeqRecord(translate_with_gaps(seq.seq),
                         id=seq.id, name=seq.name, description=seq.description)
        prots.append(prot)

    return MSA(prots)


def get_subalignment(ali, ind):
    '''Get an alignment of only some sequences.
    
    NOTE: this is really a bugfix for Biopython, come on!
    '''
    from Bio.Align import MultipleSeqAlignment
    return MultipleSeqAlignment([ali._records[i] for i in ind], ali._alphabet)


def convert_alim_to_biopython(alim, seqtype='DNA'):
    '''Convert numpy matrix to biopython alignment'''
    from Bio.Align import MultipleSeqAlignment as MSA
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet.IUPAC import ambiguous_dna, ambiguous_rna, extended_protein

    if seqtype.upper() == 'DNA':
        alpha = ambiguous_dna
    elif seqtype.upper() == 'RNA':
        alpha = ambiguous_rna
    else:
        alpha = protein
        
    ali = MSA([SeqRecord(Seq(''.join(row), alpha),
                         id='seq'+str(i),
                         name='seq'+str(i),
                         description='')
               for i, row in enumerate(alim)])

    return ali


def get_degeneracy_dict():
    '''Get dictionary of degeneracies'''
    from collections import Counter
    from Bio.Data.CodonTable import standard_dna_table
    return dict(Counter(standard_dna_table.forward_table.itervalues()))


def get_codon_back_table():
    '''Get a complete back codon table'''
    from collections import defaultdict
    from Bio.Data.CodonTable import standard_dna_table
    
    table = defaultdict(list)
    for codon, aa in standard_dna_table.forward_table.iteritems():
        table[aa].append(codon)
    return dict(table)


def get_allele_frequencies_from_MSA(alim, alpha=alpha):
    '''Get allele frequencies from a multiple sequence alignment'''
    import numpy as np
    alim = np.asarray(alim)
    af = np.zeros((len(alpha), alim.shape[1]))
    for ia, nuc in enumerate(alpha):
        af[ia] = (alim == nuc).mean(axis=0)

    return af


def get_consensus_from_MSA(alim, alpha=alpha):
    '''Get consensus from multiple sequence alignment'''
    import numpy as np
    cons = np.zeros(alim.shape[1], 'S1')
    af = get_allele_frequencies_from_MSA(alim, alpha=alpha)
    icons = af.argmax(axis=0)
    cons = alpha[icons]
    return cons


def find_annotation(seq, annotation):
    '''Find an annotation in a SeqRecord'''
    for fea in seq.features:
        if fea.id == annotation:
            return fea
    raise ValueError(annotation+' not found')


def trim_to_refseq(seq, refseq):
    '''Trim sequence to a reference sequence'''
    from seqanpy import align_overlap

    (score, ali1, ali2) = align_overlap(seq, refseq, score_gapopen=-20)
    start = len(ali2) - len(ali2.lstrip('-'))
    end = len(ali2.rstrip('-'))

    return seq[start: end]


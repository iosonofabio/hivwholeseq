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


def pretty_print_pairwise_ali(ali, name1='', name2=''):
    '''Pretty print function for pairwise alignments'''
    from itertools import izip

    for i in xrange(len(ali[0]) / 51 + 1):
        ali1_t = str(ali[0, i * 51: (i+1) * 51].seq)
        ali2_t = str(ali[1, i * 51: (i+1) * 51].seq)
        match_t = []
        for (a1, a2) in izip(ali1_t, ali2_t):
            if a1 == a2:
                match_t.append(' ')
            else:
                match_t.append('x')
        match_t = ''.join(match_t)

        lh = min(max(map(len, [name1, name2])), 6)
        print name1[:lh]+(' ' * max(0, lh - len(name1)))+':', ali1_t
        print (' ' * (lh + 1)), match_t
        print name2[:lh]+(' ' * max(0, lh - len(name2)))+':', ali2_t
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




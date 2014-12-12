# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/12/13
content:    Information module on the HIV genome.
'''
# Modules
from Bio.Seq import reverse_complement as rc



# Globals
genes = ('gag', 'pol', 'env', 'vif', 'vpr', 'vpu', 'tat', 'rev', 'nef')
proteins = ('p17', 'p24', 'p2', 'p7', 'p1', 'p6', 'PR', 'RT', 'p15', 'IN', 'gp120', 'gp41')
RNA_structures = ('RRE', "LTR5'", "LTR3'")

# Relation protein-genes
protein_genes = {'p17': 'gag', 'p24': 'gag', 'p2': 'gag', 'p7': 'gag',
                 'p1': 'gag', 'p6': 'gag',
                 'PR': 'pol', 'RT': 'pol', 'p15': 'pol', 'IN': 'pol',
                 'gp120': 'env', 'gp41': 'env'}

# Edges of genes (fuzzy)
gag_edges = ['ATGGGTGCGAGAGCGTCAGTA', 'GACCCCTCGTCACAATAA']
pol_edges = ['TTTTTTAGGGAAAATTTG', 'ACAGGATGAGGATTAG']
env_edges = ['ATGAGAGTGANGGNGANNNNGANGA',
             'TAAGACAGGGCNNGGAAAGNNNTTTGCNATAA']
vif_edges = ['ATGGAAAACAGATGGCAGGTGA', 'ACAATGAATGGACACTAG']
vpr_edges = ['ATGGAACAAGCCCCAGAAGACCAG', 'AATGGATCCAGTAGATCCTAA']
vpu_edges = ['ATGCAACCTATACCAATAGTAGCAATAGTAGCATTAGTAGTAGCAATAATAATAGCAATAGTTGTGTGGTCC',             
             'GGGATGTTGATGATCTGTAG']
tat_edges = ['ATGGAGCCAGTAGATCCTAACC', 'TCATCAANNTNCTNTANCAAAGCA',
             'ACCCNNNNCCNNNNNCNNNAGGNACNNGACAGNCNCG', 'CAGAGACAGATCCATTCGATTAG']

rev_edges = ['ATGGCAGGAAGAAGC', 'TCATCAANNTNCTNTANCAAAGCA',
             'ACCCNNNNCCNNNNNCNNNAGGNACNNGACAGNCNCG', 'ACAGTATTGGAGTCAGGAACTAAAGAATAG']

nef_edges = ['ATGGGNNGCAANTGGTCAAAANGNA','GAGTNNTACAANNACTGCTGA']

gene_edges = {'gag': gag_edges,
              'pol': pol_edges,
              'env': env_edges,
              'vif': vif_edges,
              'vpr': vpr_edges,
              'vpu': vpu_edges,
              'tat': tat_edges,
              'rev': rev_edges,
              'nef': nef_edges}

# Edges of RNA structures
RRE_edges = ['AGGAGCTATGTTCCTTGGGT', 'ACCTAAGGGATACACAGCTCCT']
LTR5 = [None, 'CTCTAGCA']
LTR3 = ['TGGANGGGNTANTTNNNTC', None]

RNA_structure_edges = {'RRE': RRE_edges,
                       "LTR5'": LTR5,
                       "LTR3'": LTR3}


# Edges of other regions
env_peptide_edges = ['ATGAGAGTGAAGGAGAA', 'GCTCCTTGGGATGTTGATGATCTGTAGTGCT']
psi_element = ['CTCGGCTTGCT', 'AGCGGAGGCTAG']
V3_edges = ['ACAATGYACACATGGAATTARGCCA', rc('AGAAAAATTCYCCTCYACAATTAAA')]

other_edges = {'env peptide': env_peptide_edges,
               'psi': psi_element,
               'V3': V3_edges}


# All edges
all_edges = dict(gene_edges)
all_edges.update(RNA_structure_edges)
all_edges.update(other_edges)



# Functions
def find_region_edges(refm, edges, minimal_fraction_match=0.60):
    '''Find a region's edges in a sequence'''
    import numpy as np

    pos_edge = []

    # Start
    if edges[0] is None:
        start = -50
        pos_edge.append(None)
    else:
        seed = np.ma.array(np.fromstring(edges[0], 'S1'))
        seed[seed == 'N'] = np.ma.masked
        seed[seed == 'R'] = np.ma.masked
        seed[seed == 'W'] = np.ma.masked
        seed[seed == 'Y'] = np.ma.masked
        sl = len(seed)
        n_match = np.array([(refm[i: i + sl] == seed).sum()
                            for i in xrange(len(refm) - sl)], int)
        pos_seed = np.argmax(n_match)
        # Check whether a high fraction of the comparable (i.e. not masked)
        # sites match the seed
        if n_match[pos_seed] > minimal_fraction_match * (-seed.mask).sum():
            start = pos_seed
            pos_edge.append(start)
        else:
            start = -50
            pos_edge.append(None)

    # End
    if edges[1] is None:
        end = None
    else:
        seed = np.ma.array(np.fromstring(edges[1], 'S1'))
        seed[seed == 'N'] = np.ma.masked
        seed[seed == 'R'] = np.ma.masked
        seed[seed == 'W'] = np.ma.masked
        seed[seed == 'Y'] = np.ma.masked
        sl = len(seed)
        n_match = np.array([(refm[i: i + sl] == seed).sum()
                            for i in xrange(start + 50, len(refm) - sl)], int)
        pos_seed = np.argmax(n_match)
        if n_match[pos_seed] > minimal_fraction_match * (-seed.mask).sum():
            end = pos_seed + start + 50 + sl
        else:
            end = None
    pos_edge.append(end)

    return pos_edge


def find_region_edges_multiple(smat, edges, min_distance=2000):
    '''Find a multiple region (e.g. split gene)'''
    import numpy as np

    pos = 0
    pos_edges = []

    for i in xrange(len(edges) / 2):
        edges_i = edges[2 * i: 2 * (i + 1)]
        pos_edge = find_region_edges(smat[pos:], edges_i)
        pos_edge = [p + pos for p in pos_edge]
        pos_edges.append(pos_edge)
        pos = pos_edge[-1] + min_distance

    return pos_edges


# NOTE: duplicate of above, but more specific to genes
def locate_gene(refseq, gene, minimal_fraction_match='auto', VERBOSE=0,
                pairwise_alignment=False, output_compact=False):
    '''Locate a gene in a sequence
    
    Parameters:
        - minimal_fraction_match: if no location with at least e.g. 66% matches
          is found, the gene edge is not found.
    '''
    import numpy as np

    gene_edge = gene_edges[gene[:3]]
    # Deal with spliced genes: e.g. tat1 takes only the first exon, tat takes all
    if len(gene_edge) > 2:
        if len(gene) == 3:
            gene_edge = [gene_edge[0], gene_edge[-1]]
        else:
            i = 2 * (int(gene[3:]) - 1)
            gene_edge = gene_edge[i: i+2]

    # Automatic precision detection
    if minimal_fraction_match == 'auto':
        # The end of vif is variable
        if gene in ['vif', 'vpu']:
            minimal_fraction_match = 0.60
        else:
            minimal_fraction_match = 0.75

    # Make a string out of refseq, whatever you get in
    refseq = ''.join(refseq) 
    refm = np.fromstring(refseq, 'S1')

    # Try out pairwise local alignment
    if pairwise_alignment:
        from seqanpy import align_local as aol
        
        # Start
        seed = gene_edge[0]
        (score, ali_seq, ali_gen) = aol(refseq, seed)
        if score > minimal_fraction_match * len(seed.replace('N', '')):
            start_found = True
            start = refseq.find(ali_seq.replace('-', ''))
        else:
            start_found = False
            start = 0

        # End
        seed = gene_edge[1]
        (score, ali_seq, ali_gen) = aol(refseq[start:], seed, score_gapopen=-100)
        if score > minimal_fraction_match * len(seed.replace('N', '')):
            end_found = True
            end = start + refseq.find(ali_seq.replace('-', '')) + len(seed)
        else:
            end_found = False
            end = len(refseq)

        import ipdb; ipdb.set_trace()

    else:

        # Find start
        start_found = True
        start = refseq.find(gene_edge[0])
        # If perfect match does not work, try imperfect
        if start == -1:
            seed = np.ma.array(np.fromstring(gene_edge[0], 'S1'))
            seed[seed == 'N'] = np.ma.masked
            sl = len(seed)
            n_match = np.array([(refm[i: i + sl] == seed).sum()
                                for i in xrange(len(refm) - sl)], int)
            pos_seed = np.argmax(n_match)
            # Check whether a high fraction of the comparable (i.e. not masked)
            # sites match the seed
            if n_match[pos_seed] > minimal_fraction_match * (-seed.mask).sum():
                start = pos_seed
            else:
                start = 0
                start_found = False
        
        # Find end
        end_found = True
        end = refseq[start + 50:].find(gene_edge[1])
        if end != -1:
            end += len(gene_edge[1])
        else:
            seed = np.ma.array(np.fromstring(gene_edge[1], 'S1'))
            seed[seed == 'N'] = np.ma.masked
            sl = len(seed)
            n_match = np.array([(refm[i: i + sl] == seed).sum()
                                for i in xrange(start + 50, len(refm) - sl)], int)
            pos_seed = np.argmax(n_match)
            if n_match[pos_seed] > minimal_fraction_match * (-seed.mask).sum():
                end = pos_seed + sl
            else:
                end = len(refseq) - start - 50
                end_found = False
        end += start + 50


    if VERBOSE:
        print '{:<5s}'.format(gene),

        if start_found:
            print 'start: '+'{:>8d}'.format(start),
        else:
            print 'start not found',

        if end_found:
            print 'end: '+'{:>8d}'.format(end)
        else:
            print 'end not found'

    if (not start_found) and (not end_found):
        start = end = -1

    if not output_compact:
        return (start, end, start_found, end_found)
    else:
        if not start_found:
            start = None
        if not end_found:
            end = None
        return (start, end)

# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/12/13
content:    Information module on the HIV genome.
'''
# Globals
# Edges of genes (fuzzy)
gag_edges = ['ATGGGTGCGAGAGCGTCAGTA', 'GACCCCTCGTCACAATAA']
pol_edges = ['TTTTTTAGGGAAAATTTG', 'ACAGGATGAGGATTAG']
env_edges = ['ATGAGAGTGAGGGGGATGCAG', 'GAAAGGGCTTTGCTATAA']
vif_edges = ['ATGGAAAACAGATGGCAGGTGA', 'ACAATGAATGGACACTAG']
vpr_edges = ['ATGGAACAAGCCCCAGAAGACCAG', 'AGATCCTAA']
vpu_edges = ['ATGCAATCTTTAGAAA', 'AATTTGTAG']
tat_edges = ['ATGGAGCCAGTAGATCCTAACC', 'ATCAAAGCA',
             'ACCCAC', 'TTCGATTAG']
rev_edges = ['ATGGCAGGAAGAAGC', 'ATCAAAGCA',
             'ACCCAC', 'GGAACTAAAGAATAG']

gene_edges = {'gag': gag_edges,
              'pol': pol_edges,
              'env': env_edges,
              'vif': vif_edges,
              'vpr': vpr_edges,
              'vpu': vpu_edges,
              'tat': tat_edges,
              'rev': rev_edges}

# Edges of RNA structures
RRE_edges = ['AGGAGCTATGTTCCTTGGGT', 'ACCTAAGGGATACACAGCTCCT']
LTR5 = ['', 'CTCTAGCA']

RNA_structure_edges = {'RRE': RRE_edges,
                       "LTR5'": LTR5}


# Edges of other regions
env_peptide_edges = ['ATGAGAGTGAAGGAGAA', 'TGTAGTGCT']
psi_element = ['CTCGGCTTGCT', 'AGCGGAGGCTAG']

other_edges = {'env peptide': env_peptide_edges,
               'psi': psi_element}



# Functions
def find_region_edges(smat, edges):
    '''Find a region's edges in a sequence'''
    import numpy as np

    pos_edge = []

    # Gene start
    emat = np.fromstring(edges[0], 'S1')
    n_matches = [(emat == smat[pos: pos+len(emat)]).sum()
                 for pos in xrange(len(smat) - len(emat))]
    pos = np.argmax(n_matches)
    pos_edge.append(pos)

    # Gene end
    emat = np.fromstring(edges[1], 'S1')
    n_matches = [(emat == smat[pos: pos+len(emat)]).sum()
                 for pos in xrange(pos_edge[0], len(smat) - len(emat))]
    pos = np.argmax(n_matches) + pos_edge[0] + len(emat)
    pos_edge.append(pos)

    return pos_edge


def find_region_edges_multiple(smat, edges):
    '''Find a multiple region (e.g. split gene)'''
    import numpy as np

    pos = 0
    pos_edges = []

    for i in xrange(len(edges) / 2):
        edges_i = edges[2 * i: 2 * (i + 1)]
        pos_edge = find_region_edges(smat[pos:], edges_i)
        pos_edge = [p + pos for p in pos_edge]
        pos_edges.append(pos_edge)
        pos = pos_edge[-1]

    return pos_edges

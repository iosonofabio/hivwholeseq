# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/02/14
content:    Test script for syn/nonsyn allele frequency trajectories in a gene.
'''
import argparse
import numpy as np
from itertools import izip
from operator import attrgetter
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from hivwholeseq.miseq import alpha
from hivwholeseq.mapping_utils import align_muscle
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_allele_frequency_trajectories_filename, \
        get_consensi_alignment_filename


# Globals
gene_fragments = {'pol': 'F2',
                  'gag': 'F1',
                  'env': 'F5'}


# Functions
def locate_gene(consensus, gene, VERBOSE=0):
    '''Locate a gene in a single fragment consensus'''
    from hivwholeseq.genome_info import gene_edges
    gene_edge = gene_edges[gene]
    
    cons = np.array(consensus, 'S1')

    # Find start
    start_found = True
    start = consensus.seq.find(gene_edge[0])
    # If perfect match does not work, try imperfect
    if start == -1:
        seed = np.fromstring(gene_edge[0], 'S1')
        sl = len(seed)
        n_match = np.array([(cons[i: i + sl] == seed).sum()
                            for i in xrange(len(cons) - sl)], int)
        pos_seed = np.argmax(n_match)
        if n_match[pos_seed] > 2 * sl / 3:
            start = pos_seed
        else:
            start = 0
            start_found = False
    if VERBOSE:
        if start_found:
            print 'Gene start:', start
        else:
            print 'Gene start not found'

    # Find end
    end_found = True
    end = consensus.seq[start + 100:].find(gene_edge[1])
    if end == -1:
        seed = np.fromstring(gene_edge[1], 'S1')
        sl = len(seed)
        n_match = np.array([(cons[i: i + sl] == seed).sum()
                            for i in xrange(start + 100, len(cons) - sl)], int)
        pos_seed = np.argmax(n_match)
        if n_match[pos_seed] > 2 * sl / 3:
            end = pos_seed + sl
        else:
            end = len(consensus) - start - 100
            end_found = False
    end += start + 100
    if VERBOSE:
        if end_found:
            print 'Gene end:', end
        else:
            print 'Gene end not found'

    return (start, end, start_found, end_found)


def pretty_print_ali(ali, name1='', name2=''):
    '''Pretty print function for pairwise alignments'''
    for i in xrange(len(ali_prot[0]) / 50):
        ali1_t = str(ali[0, i * 50: (i+1) * 50].seq)
        ali2_t = str(ali[1, i * 50: (i+1) * 50].seq)
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




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--gene', required=True,
                        help='Gene to map (e.g. pol)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    gene = args.gene
    VERBOSE = args.verbose

    # Load the NL43 respective gene for comparison
    from hivwholeseq.reference import load_NL43
    ref = load_NL43()
    from hivwholeseq.annotate_genomewide_consensus import annotate_sequence
    annotate_sequence(ref, features=['gene'])
    ref_gene = ref.features[map(attrgetter('id'), ref.features).index(gene)].extract(ref)
    ref_protein = SeqRecord(ref_gene.seq.translate(),
                            id='NL43_'+gene+'_translated', name='NL43_'+gene+'_translated',
                            description='')

    # Get initial consensus and compare with reference (check for erroneous indels)
    # TODO: for now only genes within a single fragment and one exon (pol, not tat)
    fragment = gene_fragments[gene]        
    cons = SeqIO.read(get_initial_consensus_filename(pname, fragment), 'fasta')
    
    (gene_start, gene_end,
     start_found, end_found) = locate_gene(cons, gene, VERBOSE=VERBOSE)
    print 'Gene found in fragment', fragment, 'at coordinates:', gene_start, 'to', gene_end

    cons_gene = cons[gene_start: gene_end - gene_end % 3]
    cons_protein = SeqRecord(cons_gene.seq.translate(),
                             id=gene+'_translated', name=gene+'_translated',
                             description='')

    ali_nuc = align_muscle(ref_gene, cons_gene, sort=True)
    ali_prot = align_muscle(ref_protein, cons_protein, sort=True)

    print 'Protein'
    pretty_print_ali(ali_prot, 'NL43', 'CONS')
    print

    print 'DNA'
    pretty_print_ali(ali_nuc, 'NL43', 'CONS')


    # Align several time points
    ali_nuc_pat = AlignIO.read(get_consensi_alignment_filename(pname, fragment),
                               'fasta')
    (start_gene_ali, end_gene_ali,
     start_gene_found_ali, end_gene_found_ali) = locate_gene(ali_nuc_pat[0], gene, VERBOSE=VERBOSE)

    ali_nuc_pat = ali_nuc_pat[:, start_gene_ali: end_gene_ali]
    
    #FIXME
    ali_nuc_pat = ali_nuc_pat[1:]

    ali_prot_pat = align_muscle(*([SeqRecord(s.seq.ungap('-').translate(),
                                             id=s.id, name=s.name,
                                             description='')
                                   for s in ali_nuc_pat]),
                                sort=True)

    print ali_nuc_pat.format('phylip')
    print ali_prot_pat.format('phylip')

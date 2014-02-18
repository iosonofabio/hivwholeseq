# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/13
content:    Build the initial reference for the paitent.

            All reads are then mapped against this. It is irrelevant whether or not:
                - it covers the widest possible overlaps
                - it is actually the consensus of the first time point
            
            but it is relevant that:
                - it is a good reference (complete, no holes, ok with subtype)
                - it passes basic biological checks (genes in frame)

            This makes it much easier to work with the allele frequencies later
            on, without bothering too much about insertions as a first approx.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Alphabet.IUPAC import ambiguous_dna, unambiguous_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.miseq import alpha
from hivwholeseq.filenames import get_consensus_filename
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_foldername
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.samples import samples, date_to_integer
from hivwholeseq.primer_info import primers_inner
from hivwholeseq.primer_info import primers_coordinates_HXB2_inner as pci
from hivwholeseq.one_site_statistics import \
        build_consensus_from_allele_counts_insertions as build_consensus
from hivwholeseq.reference import load_custom_reference
from hivwholeseq.annotate_genomewide_consensus import extract_feature
from hivwholeseq.mapping_utils import align_muscle
from hivwholeseq.sequence_utils import pretty_print_pairwise_ali
from hivwholeseq.filenames import get_mapped_filename
from hivwholeseq.one_site_statistics import build_consensus_from_mapped_reads
from hivwholeseq.genome_info import genes, locate_gene, gene_edges
from hivwholeseq.primer_info import fragments_genes



# Functions
def check_seq_gene_consistency(gene_seq, gene,
                               check_start_M=True,
                               check_absence_stop=True,
                               VERBOSE=0):
    '''Check the consistency of a gene seq'''
    if len(gene_seq) % 3:
        if VERBOSE:
            print 'WARNING: gene', gene, 'has a length which is not '+\
                    'a multiple of 3.'
        return (False, 'has a length which is not multiple of 3')
    
    # Check the translation
    prot_seq = gene_seq.translate()

    # Some genes do not start with Met because of polygenic transcripts
    if check_start_M and (prot_seq[0] != 'M') and (gene not in ['pol']):
        if VERBOSE:
            print 'WARNING: gene', gene, 'does not start with methionine.'
        return (False, 'does not start with methionine')

    # Look at the protein
    protm = np.array(prot_seq)
    pos_stopcod = (protm == '*').nonzero()[0]

    # Some genes do not end with * because of polygenic transcripts
    if check_absence_stop and (len(pos_stopcod) == 0) and (gene not in ['gag']):
        if VERBOSE:
            print 'WARNING: gene', gene, 'has no stop codon.'
        return (False, 'has no stop codon')

    if len(pos_stopcod) and (pos_stopcod[0] != len(prot_seq) - 1):
        if VERBOSE:
            print 'WARNING: gene', gene, 'has floating stop codons.'
        return (False, 'has floating stop codons')

    return (True, '')



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose

    # Get the patient
    patient = get_patient(pname)

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Make dir for the patient if absent
    pfolder = get_foldername(pname)
    if not os.path.isdir(pfolder):
        os.mkdir(pfolder)
        if VERBOSE >= 1:
            print pname+': folder created.'
    
    # Get the first sequenced sample
    sample_init = patient.initial_sample
    seq_run = sample_init['run']
    adaID = sample_init['adaID']
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    for fragment in fragments:

        # Write output
        output_filename = get_initial_consensus_filename(pname, fragment)
        #conss = str(SeqIO.read(output_filename, 'fasta').seq)

        # Build consensus by assembly
        bamfilename = get_mapped_filename(data_folder, adaID, fragment,
                                          filtered=True)
        conss_ass = build_consensus_from_mapped_reads(bamfilename,
                                                      block_len=100,
                                                      VERBOSE=VERBOSE - 1)

        #consensi_built = {'counts+inserts': conss_ci,
        #                  'assembly': conss_ass}

        conss = conss_ass
                       
        # Locate genes in the consensus and check they are in frame
        gene_seqs = {}
        for gene in genes:

            # FIXME: only single-exon genes for now
            if len(gene_edges[gene]) > 2:
                continue

            if fragment not in fragments_genes[gene]:
                continue

            if VERBOSE:
                print gene,

            (gene_start, gene_end,
             start_found, end_found) = locate_gene(conss, gene, VERBOSE=0)

            if (not start_found) and (not end_found):
                if VERBOSE:
                    print 'not found!'
                continue

            elif start_found and end_found:
                if VERBOSE:
                    print gene_start, gene_end,

                gene_seq = Seq(conss[gene_start: gene_end], unambiguous_dna)
                (is_good, msg) = check_seq_gene_consistency(gene_seq, gene,
                                                            VERBOSE=0)

            # If the gene end is not present, trim to codon and ignore absence
            # of stop
            elif start_found and (not end_found):
                if VERBOSE:
                    print gene_start, 'no end',

                gene_end = gene_end - (gene_end - gene_start) % 3
                gene_seq = Seq(conss[gene_start: gene_end], unambiguous_dna)
                (is_good, msg) = check_seq_gene_consistency(gene_seq, gene,
                                                            check_absence_stop=False,
                                                            VERBOSE=0)


            # If the gene start is upstream, we might be out of frame
            elif (not start_found) and (end_found):
                if VERBOSE:
                    print 'no start', gene_end,
                for gene_start in [0, 1, 2]:
                    #gene_end = gene_end - (gene_end - gene_start) % 3
                    gene_seq = Seq(conss[gene_start: gene_end], unambiguous_dna)
                    (is_good, msg) = check_seq_gene_consistency(gene_seq, gene,
                                                         check_start_M=False,
                                                         VERBOSE=0)
                    if is_good or ('no stop codon' in msg):
                        break
                if not (is_good or ('no stop codon' in msg)):
                    msg = 'had some problems.'

            # Check for a stop codon a little bit downstream
            if 'no stop codon' in msg:
                gene_end_new = len(conss) - (len(conss) - gene_start) % 3
                gene_seq = Seq(conss[gene_start: gene_end_new], unambiguous_dna)
                prot_seq = gene_seq.translate()
                pos_first_stop = prot_seq.find('*')
                if pos_first_stop != -1:
                    gene_seq = gene_seq[:(pos_first_stop + 1) * 3]
                    msg = msg+', first stop codon found '+\
                            str(1 + pos_first_stop - (gene_end - gene_start) / 3)+\
                            ' codons downstream.'
                    gene_end = gene_end_new
                else:
                    msg = msg+' no stop codon found downstream at all.'

            if is_good:
                print '... in frame (OK)'
            else:
                print 'WARNING: gene', gene, msg

                ref = load_custom_reference('HXB2', 'gb')
                gene_ref = extract_feature(ref, gene)

                prot_seq = gene_seq.translate()
                prot_ref = gene_ref.seq.translate()

                ali = align_muscle(SeqRecord(prot_seq, id='cons'),
                                   SeqRecord(prot_ref, id='ref'), sort=True)
                pretty_print_pairwise_ali(ali, 'cons', 'HXB2')
                
            gene_seqs[gene] = gene_seq

        seq_in = SeqRecord(Seq(conss, unambiguous_dna),
                           id='cons_init_p'+pname,
                           name='cons_init_p'+pname,
                           description='Initial consensus of patient '+pname)

        ## Read the new consensus
        #input_filename = get_consensus_filename(data_folder, adaID, fragment)
        #seq_in = SeqIO.read(input_filename, 'fasta')

        # If absent, just copy the thing over
        if not os.path.isfile(output_filename):
            SeqIO.write(seq_in, output_filename, 'fasta')
            if VERBOSE >= 1:
                print pname+': initial consensus file created for sample', sample_init['name']

        # if present, check whether the sequences are the same (if so, no remapping
        # is needed!). Overwrite the file anyway, because single samples carry
        # their consensus (mapping reference) with them in the folder (not much
        # overhead and MUUUCH cleaner than otherwise).
        else:
            seq_out = SeqIO.read(output_filename, 'fasta')
            SeqIO.write(seq_in, output_filename, 'fasta')
            if str(seq_in.seq) != str(seq_out.seq):
                print 'NOTE: initial consensus updated to '+sample_init['name']+', remap!'
            

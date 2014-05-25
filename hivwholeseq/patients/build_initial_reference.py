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
from hivwholeseq.genome_info import locate_gene, gene_edges
from hivwholeseq.genome_info import genes as genes_all
from hivwholeseq.primer_info import fragments_genes



# Functions
def merge_initial_consensi(conss_frags, minimal_fraction_match=0.75, VERBOSE=0):
    '''Merge fragment consensi into a genome-wide one'''
    # Find overlaps
    if VERBOSE:
        print 'Starting with F1'
    conss_all = [conss_frags[0]]
    for ifr, conss in enumerate(conss_frags[1:], 2):
        if VERBOSE:
            print 'merging in F'+str(ifr)

        # Avoid the first few bases because of low coverage
        pos_seed_new = 30
        seed = conss[pos_seed_new: pos_seed_new + 30]
        pos_seed = conss_all[-1].find(seed)
        # Allow imperfect matches
        if pos_seed == -1:
            refm = np.fromstring(conss, 'S1')
            seed = np.fromstring(seed, 'S1')
            sl = len(seed)
            n_match = np.array([(refm[i: i + sl] == seed).sum()
                                for i in xrange(len(refm) - sl)], int)
            pos_seed = len(refm) - 1 - np.argmax(n_match[::-1])
            if n_match[pos_seed] < minimal_fraction_match * sl:
                overlap_found = False
                if VERBOSE:
                    print 'WARNING: Cannot merge consensus F'+str(ifr)+' with previous ones'
            else:
                overlap_found = True
        else:
            overlap_found = True

        if overlap_found:
            conss_all[-1] = conss_all[-1][:pos_seed]
            conss_all.append(conss[pos_seed_new:])
        else:
            conss_all.append(('N' * 10) + conss)

    return ''.join(conss_all)


def check_seq_gene_consistency(gene_seq, gene,
                               check_start_M=True,
                               check_absence_stop=True,
                               VERBOSE=0):
    '''Check the consistency of a gene seq'''
    if len(gene_seq) <= 0:
        if VERBOSE:
            print 'WARNING: gene', gene, 'is empty'
        return (False, 'is empty')

    if len(gene_seq) % 3:
        if VERBOSE and (gene not in ('tat1', 'tat2', 'rev1', 'rev2')):
            print 'WARNING: gene', gene, 'has a length which is not '+\
                    'a multiple of 3.'
        if gene not in ('tat1', 'tat2', 'rev1', 'rev2'):
            return (False, 'has a length which is not multiple of 3')
        else:
            return (True, 'has a length which is not multiple of 3 (OK)')
    
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

    return (True, 'in frame (OK)')


def check_gene_consensus(conss, fragment, gene, max_end_slippage=10,
                         start_search=0, VERBOSE=0):
    '''Check a single exon for consensus'''
    if VERBOSE:
        print gene,

    if (not start_found) and (not end_found):
        if VERBOSE:
            print 'not found!'
            genes_good[gene] = False
        return None

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
        gene_end_old = gene_end
        gene_end = len(conss) - (len(conss) - gene_start) % 3
        gene_seq = Seq(conss[gene_start: gene_end], unambiguous_dna)
        prot_seq = gene_seq.translate()
        pos_first_stop = prot_seq.find('*')
        if pos_first_stop != -1:
            gene_seq = gene_seq[:(pos_first_stop + 1) * 3]
            msg = msg+', first stop codon found '+\
                    str(1 + pos_first_stop - (gene_end_old - gene_start) / 3)+\
                    ' codons downstream'
            gene_end = gene_start + len(gene_seq)
            if (gene_end - gene_end_old) < max_end_slippage:
                is_good = True
                msg = msg + ' (OK)'
            else:
                msg = msg + ', this looks suspicious'
        else:
            msg = msg+' no stop codon found downstream at all.'

    if is_good:
        print msg
    else:
        print 'WARNING: gene', gene, msg

        # No reference alignment and crap for exons
        if len(gene) == 3:

            ref = load_custom_reference('HXB2', 'gb')
            gene_ref = extract_feature(ref, gene)
            prot_seq = gene_seq.translate()
            prot_ref = gene_ref.seq.translate()

            ali = align_muscle(SeqRecord(prot_seq, id='cons'),
                               SeqRecord(prot_ref, id='ref'), sort=True)
            pretty_print_pairwise_ali(ali, 'cons', 'HXB2')

    return (is_good, gene_seq, gene_start, gene_end)


def check_genes_consensus(conss, fragment, genes=genes_all, max_end_slippage=10, VERBOSE=0):
    '''Check gene consistency in a consensus'''
    # Locate genes in the consensus and check they are in frame
    genes_good = {}
    gene_seqs = {}
    gene_poss = {}
    for gene in genes:
        fragments_gene = fragments_genes[gene]

        # Single-exon genes
        if len(gene_edges[gene]) == 2:
            if (fragment != 'genomewide') and (fragment not in fragments_gene):
                continue

            if VERBOSE:
                print gene,

            # Locate gene
            (gene_start, gene_end,
             start_found, end_found) = locate_gene(conss, gene, VERBOSE=0)

            # If neither start nor end are found, skip
            if (not start_found) and (not end_found):
                if VERBOSE:
                    print 'not found!'
                    genes_good[gene] = False
                continue
        
            # If both start and end are found, ok
            elif start_found and end_found:
                if VERBOSE:
                    print gene_start, gene_end,
        
                gene_seq = Seq(conss[gene_start: gene_end], unambiguous_dna)
                gene_pos = (gene_start, gene_end)
                (is_good, msg) = check_seq_gene_consistency(gene_seq, gene,
                                                            VERBOSE=0)
        
            # If the gene end is not present, trim to codon and ignore if no stop
            elif start_found and (not end_found):
                if VERBOSE:
                    print gene_start, 'no end',
        
                gene_end = gene_end - (gene_end - gene_start) % 3
                gene_seq = Seq(conss[gene_start: gene_end], unambiguous_dna)
                gene_pos = (gene_start, gene_end)
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
                    gene_pos = (gene_start, gene_end)
                    (is_good, msg) = check_seq_gene_consistency(gene_seq, gene,
                                                         check_start_M=False,
                                                         VERBOSE=0)
                    if is_good or ('no stop codon' in msg):
                        break
                if not (is_good or ('no stop codon' in msg)):
                    msg = 'had some problems.'
        
            # Check for a stop codon a little bit downstream
            if 'no stop codon' in msg:
                gene_end_old = gene_end
                gene_end = len(conss) - (len(conss) - gene_start) % 3
                gene_seq = Seq(conss[gene_start: gene_end], unambiguous_dna)
                gene_pos = (gene_start, gene_end)
                prot_seq = gene_seq.translate()
                pos_first_stop = prot_seq.find('*')
                if pos_first_stop != -1:
                    gene_seq = gene_seq[:(pos_first_stop + 1) * 3]
                    msg = msg+', first stop codon found '+\
                            str(1 + pos_first_stop - (gene_end_old - gene_start) / 3)+\
                            ' codons downstream'
                    gene_end = gene_start + len(gene_seq)
                    gene_pos = (gene_start, gene_end)
                    if (gene_end - gene_end_old) < max_end_slippage:
                        is_good = True
                        msg = msg + ' (OK)'
                    else:
                        msg = msg + ', this looks suspicious'
                else:
                    msg = msg+' no stop codon found downstream at all.'
        
            if is_good:
                if VERBOSE:
                    print msg
            else:
                if VERBOSE:
                    print 'WARNING: gene', gene, msg,
        
                # No reference alignment and crap for exons
                if len(gene) == 3:
        
                    ref = load_custom_reference('HXB2', 'gb')
                    gene_ref = extract_feature(ref, gene)
                    if ('N' * 10) not in str(gene_seq):
                        print ''

                        prot_seq = gene_seq.translate()
                        prot_ref = gene_ref.seq.translate()
        
                        ali = align_muscle(SeqRecord(prot_seq, id='cons'),
                                           SeqRecord(prot_ref, id='ref'), sort=True)
                        pretty_print_pairwise_ali(ali, 'cons', 'HXB2')

                    else:
                        print 'gene not fully sequenced!'

                else:
                    print ''
    
            genes_good[gene] = is_good
            gene_seqs[gene] = gene_seq
            gene_poss[gene] = [gene_pos]

        # Multi-exon genes
        else:
            gene_seq = Seq('', ambiguous_dna)
            gene_pos = []
            start_search = 0
            for exon_num in xrange(len(gene_edges[gene]) // 2):
                exon = gene+str(exon_num + 1)
                fragments_exon = fragments_gene[exon_num]

                if (fragment != 'genomewide') and (fragment not in fragments_exon):
                    continue

                if VERBOSE and not len(gene_pos):
                    print gene,

                if VERBOSE and (fragment == 'genomewide'):
                    if not len(gene_pos):
                        print ''
                    print exon,


                # Locate exon
                if start_search >= len(conss):
                    if VERBOSE >= 2:
                        'exon not in consensus'
                    continue
                (exon_start, exon_end,
                 start_found, end_found) = locate_gene(conss[start_search:], exon,
                                                       VERBOSE=0)
                if start_found:
                    exon_start += start_search
                exon_end += start_search

                # Check exon
                if (not start_found) and (not end_found):
                    if VERBOSE >= 2:
                        print 'not found!'
                        is_good = False
                        break

                elif start_found and end_found:
                    if VERBOSE:
                        print exon_start, exon_end,
            
                    exon_seq = Seq(conss[exon_start: exon_end], unambiguous_dna)
                    exon_pos = (exon_start, exon_end)
                    (is_good, msg) = check_seq_gene_consistency(exon_seq, exon,
                                                                check_absence_stop=False,
                                                                check_start_M=False,
                                                                VERBOSE=0)

                elif start_found and (not end_found):
                    if VERBOSE >= 2:
                        print exon_start, 'no end',
        
                    exon_seq = Seq(conss[exon_start: exon_end], unambiguous_dna)
                    exon_pos = (exon_start, exon_end)
                    (is_good, msg) = check_seq_gene_consistency(exon_seq, exon,
                                                                check_absence_stop=False,
                                                                check_start_M=False,
                                                                VERBOSE=0)
        
                # If the gene start is upstream, we might be out of frame
                elif (not start_found) and (end_found):
                    if VERBOSE >= 2:
                        print 'no start', exon_end,
                    for exon_start in [0, 1, 2]:
                        exon_seq = Seq(conss[exon_start: exon_end], unambiguous_dna)
                        exon_pos = (exon_start, exon_end)
                        (is_good, msg) = check_seq_gene_consistency(gene_seq, exon,
                                                                    check_absence_stop=False,
                                                                    check_start_M=False,
                                                                    VERBOSE=0)
                        if is_good:
                            break
                    if not is_good:
                        msg = 'had some problems.'


                if VERBOSE >= 2:
                    if is_good:
                        print msg
                    else:
                        print 'WARNING: exon', exon, msg

                # Add up exons for the gene
                gene_seq = gene_seq + exon_seq
                gene_pos.append(exon_pos)
                gene_seqs[exon] = exon_seq
                gene_poss[exon] = [exon_pos]

                if start_found:
                    start_search = exon_end + 1000

            # Check whole gene
            if fragment == 'genomewide':
                (is_good, msg) = check_seq_gene_consistency(gene_seq, gene,
                                                            VERBOSE=0)

                if VERBOSE:
                    if is_good:
                        print msg
                    else:
                        print 'WARNING: gene', gene, msg

                gene_seqs[gene] = gene_seq
                genes_good[gene] = is_good
                gene_poss[gene] = gene_pos


    return (gene_seqs, genes_good, gene_poss)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose

    patient = get_patient(pname)

    # If the script is called with no fragment, iterate over all
    # NOTE: This assumes that the fragment has been sequenced. If not, we should
    # take a later time point or similia, but let's not worry about it yet.
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
        if fragment == 'genomewide':
            continue

        # Write output
        output_filename = get_initial_consensus_filename(pname, fragment)

        # Take consensus from folder, it was built by assisted assembly
        cons_rec = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment), 'fasta')
        conss = str(cons_rec.seq)

        ## Build consensus by assisted assembly
        #bamfilename = get_mapped_filename(data_folder, adaID, fragment,
        #                                  filtered=True)
        #conss = build_consensus_from_mapped_reads(bamfilename,
        #                                          maxreads=3000,
        #                                          block_len=100,
        #                                          VERBOSE=VERBOSE - 1)

        check_genes_consensus(conss, fragment, VERBOSE=VERBOSE)

        seq_in = SeqRecord(Seq(conss, unambiguous_dna),
                           id='cons_init_p'+pname,
                           name='cons_init_p'+pname,
                           description='Initial consensus of patient '+pname+', fragment '+fragment)

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
            
    # Merge all fragments if requested
    if 'genomewide' in fragments:
        conss_frags = []
        for fragment in ['F'+str(i) for i in xrange(1, 7)]:
            output_filename = get_initial_consensus_filename(pname, fragment)
            conss = str(SeqIO.read(output_filename, 'fasta').seq)
            conss_frags.append(conss)
        conss_genomewide = merge_initial_consensi(conss_frags, VERBOSE=VERBOSE,
                                                  minimal_fraction_match=0.75)

        ## Take consensus from folder, it was built by assisted assembly
        #cons_gw_rec = SeqIO.read(get_consensus_filename(data_folder, adaID, 'genomewide'), 'fasta')
        #conss_genomewide = str(cons_gw_rec.seq)
        
        gene_seqs, genes_good, gene_poss = check_genes_consensus(conss_genomewide, 'genomewide', VERBOSE=VERBOSE)
        if all(genes_good.values()):
            print 'Genomewide consensus approved: you can map the single fragments!'

        output_filename = get_initial_consensus_filename(pname, 'genomewide')
        seq_in = SeqRecord(Seq(conss_genomewide, unambiguous_dna),
                           id='cons_gw_init_p'+pname,
                           name='cons_gw_init_p'+pname,
                           description='Initial genomewide consensus of patient '+pname)
        SeqIO.write(seq_in, output_filename, 'fasta')

        




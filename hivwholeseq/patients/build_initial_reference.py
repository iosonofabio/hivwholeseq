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
import sys
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Alphabet.IUPAC import ambiguous_dna, unambiguous_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from hivwholeseq.samples import SampleSeq
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.miseq import alpha
from hivwholeseq.filenames import get_consensus_filename
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_foldername
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
from hivwholeseq.sequence_utils import merge_sequences



# Functions
def complement_consensus_PCR2(cons_rec, patient, fragment, samplen, VERBOSE=0):
    '''Complement consensus from PCR2 with wings from later PCR1 sample'''
    from hivwholeseq.sequence_utils import find_seed_imperfect, rfind_seed_imperfect

    found = False
    for _, sampletmp in patient.samples.iloc[samplen + 1:].iterrows():
        for _, sampleseqtmp in sampletmp['samples seq'].iterrows():
            sampleseqtmp = SampleSeq(sampleseqtmp)
            if int(sampleseqtmp.PCR) == 1:
                sampleseq_later = sampleseqtmp
                found = True
                break
        if found:
            break

    adaID_later = sampleseq_later['adapter']
    data_folder_later = sampleseq_later.sequencing_run.folder
    cons_rec_later = SeqIO.read(get_consensus_filename(data_folder_later, adaID_later, fragment), 'fasta')
    conss_later = str(cons_rec_later.seq)

    start = find_seed_imperfect(cons_rec_later, cons_rec[:20])
    end = rfind_seed_imperfect(cons_rec_later, cons_rec[-20:]) + 20

    if VERBOSE >= 1:
        print 'Complementing PCR2 consensus with later PCR1:',
        print sampleseq_later.name, sampleseq_later['seq run'], sampleseq_later.adapter

    frag_spec = sampleseq_later.regions_complete[sampleseq_later.regions_generic.index(fragment)]

    return (frag_spec, conss_later[:start]+cons_rec+conss_later[end:])



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Build initial reference',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 genomewide genomewide_new)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--samplenumber', type=int, default=0,
                        help='Index of the patient sample to use for the reference')
    parser.add_argument('--repnumber', type=int, default=0,
                        help='Index of the sequenced sample within that patient sample')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    repn = args.repnumber
    samplen = args.samplenumber

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Make dir for the patient if absent
    pfolder = patient.folder
    if not os.path.isdir(pfolder):
        os.mkdir(pfolder)
        if VERBOSE >= 1:
            print pname+': folder created.'
    
    # Get the first sequenced sample
    sample_init_pat = patient.samples.iloc[samplen]

    for fragment in fragments:
        if 'genomewide' in fragment:
            continue

        # Take any sequenced sample, for a consensus it should not matter
        if (patient.name in ('15241', '15319')) and (fragment in ('F4', 'F5', 'F6')):
            sample_init = SampleSeq(sample_init_pat['samples seq'].iloc[max(1, repn)])
        else:
            sample_init = SampleSeq(sample_init_pat['samples seq'].iloc[repn])

        seq_run = sample_init['seq run']
        adaID = sample_init['adapter']
        dataset = sample_init.sequencing_run
        data_folder = dataset.folder

        if VERBOSE:
            print 'Initial sample:', sample_init.name, sample_init['seq run'], sample_init.adapter

        cons_rec = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment), 'fasta')
        frag_spec = sample_init.regions_complete[sample_init.regions_generic.index(fragment)]

        # Complement PCR2 initial reference with tails from a later sample
        if int(sample_init.PCR) == 2:
            (frag_spec, cons_rec) = complement_consensus_PCR2(cons_rec, patient, fragment,
                                                              samplen, VERBOSE=VERBOSE)

        conss = str(cons_rec.seq)
        output_filename = get_initial_consensus_filename(pname, fragment)

        seq_in = SeqRecord(Seq(conss, unambiguous_dna),
                           id='cons_init_p'+pname+'_'+frag_spec,
                           name='cons_init_p'+pname+'_'+frag_spec,
                           description='Initial consensus of patient '+pname+', fragment '+frag_spec)

        # If absent, just copy the thing over
        if not os.path.isfile(output_filename):
            if VERBOSE >= 1:
                print pname+': initial consensus file created for sample', sample_init.name, 'fragment', fragment
            SeqIO.write(seq_in, output_filename, 'fasta')

        # if present, check whether the sequences are the same (if so, no remapping is needed)
        else:
            seq_out = SeqIO.read(output_filename, 'fasta')
            if str(seq_in.seq) != str(seq_out.seq):
                print 'NOTE: initial consensus changed for fragment '+fragment+': now '+sample_init.name+', remap!'
            SeqIO.write(seq_in, output_filename, 'fasta')
            
    # Copy genomewide consensus too if explicitely requested
    if ('genomewide' in fragments) or ('genomewide_new' in fragments):
        if 'genomewide' in fragments:
            from hivwholeseq.filenames import get_merged_consensus_filename

            sample_init = SampleSeq(sample_init_pat['samples seq'].iloc[repn])

            seq_run = sample_init['seq run']
            adaID = sample_init['adapter']
            dataset = sample_init.sequencing_run
            data_folder = dataset.folder

            # Take consensus from folder, it was built by assembling the single fragments
            input_filename = get_merged_consensus_filename(data_folder, adaID, ['F'+str(i+1) for i in xrange(6)])
            cons_rec = SeqIO.read(input_filename, 'fasta')
            conss = str(cons_rec.seq)

        else:
            consensi = [SeqIO.read(get_initial_consensus_filename(pname, 'F'+str(ifr+1)), 'fasta')
                           for ifr in xrange(6)]

            conss = merge_sequences(consensi, VERBOSE=VERBOSE)
            print len(conss)

        output_filename = get_initial_consensus_filename(pname, 'genomewide')
        seq_in = SeqRecord(Seq(conss, unambiguous_dna),
                           id='cons_init_p'+pname,
                           name='cons_init_p'+pname,
                           description='Initial consensus of patient '+pname+', genomewide')

        # If absent, just copy the thing over
        if not os.path.isfile(output_filename):
            if VERBOSE >= 1:
                if 'genomewide' in fragments:
                    print pname+': initial genomewide consensus file created for sample', sample_init.name
                else:
                    print pname+': initial genomewide consensus file created by de novo assembly'
            SeqIO.write(seq_in, output_filename, 'fasta')

        # if present, check whether the sequences are the same (if so, no remapping
        # is needed!). Overwrite the file anyway, because single samples carry
        # their consensus (mapping reference) with them in the folder (not much
        # overhead and MUUUCH cleaner than otherwise).
        else:
            seq_out = SeqIO.read(output_filename, 'fasta')
            if str(seq_in.seq) != str(seq_out.seq):
                if 'genomewide' in fragments:
                    print 'NOTE: initial genomewide consensus changed: now '+sample_init.name
                else:
                    print 'NOTE: initial genomewide consensus changed: de novo'
            SeqIO.write(seq_in, output_filename, 'fasta')


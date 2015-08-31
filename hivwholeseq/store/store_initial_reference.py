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

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.sequencing.samples import SampleSeq
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.samples import load_sample_sequenced, SamplePat
from hivwholeseq.utils.miseq import alpha
from hivwholeseq.sequencing.filenames import get_consensus_filename
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_foldername, get_initial_reference_foldername
from hivwholeseq.sequencing.primer_info import primers_inner
from hivwholeseq.sequencing.primer_info import primers_coordinates_HXB2_inner as pci
from hivwholeseq.utils.one_site_statistics import \
        build_consensus_from_allele_counts_insertions as build_consensus
from hivwholeseq.reference import load_custom_reference
from hivwholeseq.annotate_genomewide_consensus import extract_feature
from hivwholeseq.utils.mapping import align_muscle
from hivwholeseq.utils.sequence import pretty_print_pairwise_ali
from hivwholeseq.sequencing.filenames import get_mapped_filename
from hivwholeseq.utils.one_site_statistics import build_consensus_from_mapped_reads
from hivwholeseq.utils.genome_info import locate_gene, gene_edges
from hivwholeseq.utils.genome_info import genes as genes_all
from hivwholeseq.sequencing.primer_info import fragments_genes
from hivwholeseq.utils.sequence import merge_sequences



# Functions
def complement_consensus_PCR2(cons_rec, patient, fragment, samplen, VERBOSE=0):
    '''Complement consensus from PCR2 with wings from later PCR1 sample'''
    from hivwholeseq.utils.sequence import find_seed_imperfect, rfind_seed_imperfect

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
    parser.add_argument('--sample',
                        help='Use a specific sample (not the first time point) for the reference')
    parser.add_argument('--repnumber', type=int, default=0,
                        help='Index of the sequenced sample within that patient sample')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    repn = args.repnumber
    samplename = args.sample

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()

    mkdirs(get_initial_reference_foldername(pname))

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments
    
    if samplename is None:
        sample = SamplePat(patient.samples.iloc[samplen])
    else:
        sample = load_sample_sequenced(samplename)

    for fragment in fragments:
        sample_seq = SampleSeq(sample.samples_seq.iloc[repn])

        seq_run = sample_seq['seq run']
        adaID = sample_seq['adapter']
        dataset = sample_seq.sequencing_run
        data_folder = dataset.folder

        if VERBOSE:
            print 'Initial sample:', sample_seq.name, sample_seq['seq run'],
            print sample_seq.adapter

        cons_rec = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment),
                              'fasta')
        frag_spec = sample_seq.regions_complete[\
                            sample_seq.regions_generic.index(fragment)]

        # Complement PCR2 initial reference with tails from a later sample
        if int(sample_seq.PCR) == 2:
            (frag_spec, cons_rec) = complement_consensus_PCR2(cons_rec, patient,
                                                              fragment,
                                                              samplen,
                                                              VERBOSE=VERBOSE)

        conss = str(cons_rec.seq)
        output_filename = get_initial_reference_filename(pname, fragment)

        seq_in = SeqRecord(Seq(conss, unambiguous_dna),
                           id='cons_init_p'+pname+'_'+frag_spec,
                           name='cons_init_p'+pname+'_'+frag_spec,
                           description='Initial consensus of patient '+pname+\
                                       ', fragment '+frag_spec)

        # If absent, just copy the thing over
        if not os.path.isfile(output_filename):
            if VERBOSE >= 1:
                print pname+': initial consensus file created for sample', \
                        sample_seq.name, 'fragment', fragment
            SeqIO.write(seq_in, output_filename, 'fasta')

        # if present, check whether the sequences are the same (if so, no
        # remapping is needed)
        else:
            seq_out = SeqIO.read(output_filename, 'fasta')
            if str(seq_in.seq) != str(seq_out.seq):
                print 'NOTE: initial consensus changed for fragment '+fragment+\
                        ': now '+sample_seq.name+', remap!'
            SeqIO.write(seq_in, output_filename, 'fasta') 

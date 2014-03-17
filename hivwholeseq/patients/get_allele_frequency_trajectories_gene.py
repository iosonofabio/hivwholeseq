#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/13
content:    Collect the allele frequencies of all samples to the initial
            consensus and save them as a matrix into a single trajectory file.
'''
# Modules
import os
import argparse
from operator import attrgetter
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio import SeqIO

from hivwholeseq.miseq import alpha
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file, \
        get_allele_counts_insertions_from_file_unfiltered, \
        filter_nus
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories as plot_nus
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories_3d as plot_nus_3d
from hivwholeseq.genome_info import genes, locate_gene, gene_edges



# Globals



# Functions


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--gene', required=True,
                        help='Gene to analyze (e.g. gag)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele frequency trajectories to file')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the allele frequency trajectories')

    args = parser.parse_args()
    pname = args.patient
    gene = args.gene
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot

    # Get patient and the sequenced samples (and sort them by date)
    patient = get_patient(pname)

    # Get reference
    from hivwholeseq.annotate_genomewide_consensus import annotate_sequence
    cons_rec = SeqIO.read(get_initial_consensus_filename('20097', 'genomewide'), 'fasta')
    conss = str(cons_rec.seq)
        
    # Get allele freq trajectories
    aft_filename = get_allele_frequency_trajectories_filename(pname, 'genomewide')
    nus = np.load(aft_filename)

    # Extract gene from reference and allele freq trajectories
    from hivwholeseq.patients.build_initial_reference import check_genes_consensus
    gene_seqs, _, gene_poss = check_genes_consensus(conss, 'genomewide', genes=[gene],
                                                    VERBOSE=VERBOSE)
    conss_gene = gene_seqs[gene]
    gene_pos = gene_poss[gene]
    nus_gene = np.concatenate([nus[:, :, exon_pos[0]: exon_pos[1]]
                               for exon_pos in gene_pos], axis=2)

    # Plot
    if plot:
        import matplotlib.pyplot as plt

        times = patient.times()
        plot_nus(times, nus_gene, title='Patient '+pname+', '+gene, VERBOSE=VERBOSE)

        plot_nus_3d(times, nus_gene, title='Patient '+pname+', '+gene, VERBOSE=VERBOSE)

    if plot:
        plt.tight_layout()
        plt.ion()
        plt.show()

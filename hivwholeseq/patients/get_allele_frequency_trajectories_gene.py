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
from hivwholeseq.generic_utils import mkdirs


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
    parser.add_argument('--saveplot', action='store_true',
                        help='Save the plot to file')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--PCR1', action='store_true',
                        help='Show only PCR1 samples where possible (still computes all)')
    
    args = parser.parse_args()
    pname = args.patient
    gene = args.gene
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot
    saveplot = args.saveplot
    use_PCR1 = args.PCR1

    patient = get_patient(pname)
    times = patient.times()
    samplenames = patient.samples
    if use_PCR1:
        # Keep PCR2 only if PCR1 is absent
        ind = np.nonzero(map(lambda x: ('PCR1' in x[1]) or ((times == times[x[0]]).sum() == 1), enumerate(samplenames)))[0]
        times = times[ind]

    from hivwholeseq.annotate_genomewide_consensus import annotate_sequence
    cons_rec = SeqIO.read(get_initial_consensus_filename('20097', 'genomewide'), 'fasta')
    conss = str(cons_rec.seq)
        
    aft_filename = get_allele_frequency_trajectories_filename(pname, 'genomewide')
    aft = np.load(aft_filename)

    # Extract gene from reference and allele freq trajectories
    from hivwholeseq.patients.build_initial_reference import check_genes_consensus
    gene_seqs, _, gene_poss = check_genes_consensus(conss, 'genomewide', genes=[gene],
                                                    VERBOSE=VERBOSE)
    conss_gene = gene_seqs[gene]
    gene_pos = gene_poss[gene]
    aft_gene = np.concatenate([aft[:, :, exon_pos[0]: exon_pos[1]]
                               for exon_pos in gene_pos], axis=2)
    if save_to_file:
        aft_gene.dump(get_allele_frequency_trajectories_filename(pname, gene))

    if use_PCR1:
        aft = aft[ind]
        aft_gene = aft_gene[ind]

    # Plot
    if (plot is not None):
        import matplotlib.pyplot as plt
        if plot in ('2D', '2d', ''):
            plot_nus(times, aft_gene, title='Patient '+pname+', '+gene, VERBOSE=VERBOSE,
                     options=['syn-nonsyn'])
        elif plot in ('3D', '3d'):
            plot_nus_3d(times, aft_gene, title='Patient '+pname+', '+gene, VERBOSE=VERBOSE)

        plt.tight_layout()

    if (plot is not None) and (not saveplot):
        plt.ion()
        plt.show()

    if saveplot:
        fn = get_allele_frequency_trajectories_filename(pname, gene)
        mkdirs(os.path.dirname(fn))
        plt.savefig(fn)

# vim: fdm=marker
'''
author:     Fabio Zanini
date:       19/05/14
content:    Calculate and plot the propagator of allele frequencies.
'''
# Modules
import sys
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet.IUPAC import unambiguous_dna
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import patients as patients_all
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename



# Script
if __name__ == '__main__': 


    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele frequency trajectories')
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--genes', default=('gag', 'pol'), nargs='+',
                        help='Gene to analyze (e.g. pol)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the propagator to file')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the propagator')

    args = parser.parse_args()
    pnames = args.patients
    genes = map(lambda x: x.lower(), args.genes)
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot
    submit = args.submit

    if pnames is None:
        patients = [p for p in patients_all if (len(set(p.times())) >= 3) and (p.id != '15241')]
    else:
        patients = [p for p in patients_all if p.id in pnames]

    # Submit to the cluster self if requested (--save is assumed)
    if submit:
        fork_self(pnames, gene, VERBOSE=VERBOSE)
        sys.exit()

    fragments = ['F'+str(i) for i in xrange(1, 7)]

    # Prepare output data structures
    bins_to = np.insert(np.logspace(-2.5, 0, 10), 0, 0)
    bins_from = bins_to[(bins_to > 1e-3) & (bins_to < 0.5)]
    hist = np.zeros((len(bins_from) - 1, len(bins_to) - 1), float)
    hist_syn = np.zeros((len(bins_from) - 1, len(bins_to) - 1), float)
    hist_nonsyn = np.zeros((len(bins_from) - 1, len(bins_to) - 1), float)

    for patient in patients:
        pname = patient.id
        times = patient.times()
        samplenames = patient.samples

        # Keep PCR2 only if PCR1 is absent
        ind = np.nonzero(map(lambda x: ('PCR1' in x[1]) or ((times == times[x[0]]).sum() == 1),
                             enumerate(samplenames)))[0]
        times = times[ind]
        samplenames = [s for (i, s) in enumerate(samplenames) if i in ind] 
    
        for gene in genes:
            for fragment in fragments:
        
                # FIXME
                if gene == 'pol':
                    if fragment != 'F2':
                        continue
                
                elif gene == 'gag':
                    if fragment != 'F1':
                        continue
    
                else:
                    raise ValueError('Only pol/F2 and gag/F1 implemented!')

                if VERBOSE >= 1:
                    print pname, fragment
        
                act_filename = get_allele_count_trajectories_filename(pname, fragment)
                aft_filename = get_allele_frequency_trajectories_filename(pname, fragment)
        
                aft = np.load(aft_filename)
                act = np.load(act_filename)
        
                aft[np.isnan(aft)] = 0
                aft[(aft < 1e-5) | (aft > 1)] = 0
    
                # Extract gene from reference and allele freq trajectories
                from hivwholeseq.patients.build_initial_reference import check_genes_consensus
                cons_rec = SeqIO.read(get_initial_reference_filename(pname, fragment), 'fasta')
                conss = str(cons_rec.seq)
                gene_seqs, genes_good, gene_poss = check_genes_consensus(conss, fragment, genes=[gene],
                                                                         VERBOSE=VERBOSE)
                if not genes_good[gene]:
                    if VERBOSE >= 1:
                        print pname+':', gene, 'not found or with problems: skipping.'
                    continue
    
                gene_pos = gene_poss[gene]
                aft_gene = np.concatenate([aft[:, :, exon_pos[0]: exon_pos[1]]
                                   for exon_pos in gene_pos], axis=2)
                conss_gene = gene_seqs[gene]
                gene_len = len(conss_gene)
    
                # Collect counts
                for i in xrange(aft_gene.shape[0] - 1):
                    hist_frag = np.histogram2d(aft_gene[i].ravel(), aft_gene[i + 1].ravel(),
                                               bins=[bins_from, bins_to])
                    hist += hist_frag[0]
    
                    nu_syn = [[], []]
                    nu_nonsyn = [[], []]
                    cod_anc = MutableSeq('AAA', unambiguous_dna)
                    cod_new = MutableSeq('AAA', unambiguous_dna)
                    for j in xrange(gene_len // 3):
                        for jcod in xrange(3):
                            for ai in xrange(4):
                                cod_anc[:] = conss_gene[3 * j: 3 * (j+1)]
                                # Ancestral allele, skip (we only look at propagation of MINOR alleles)
                                if alpha[ai] == cod_anc[jcod]:
                                    continue
    
                                cod_new[:] = conss_gene[3 * j: 3 * (j+1)]
                                cod_new[jcod] = alpha[ai]
    
                                if str(cod_new.toseq().translate()) != str(cod_anc.toseq().translate()):
                                    nu_syn[0].append(aft_gene[i, ai, j + jcod])
                                    nu_syn[1].append(aft_gene[i + 1, ai, j + jcod])
                                else:
                                    nu_nonsyn[0].append(aft_gene[i, ai, j + jcod])
                                    nu_nonsyn[1].append(aft_gene[i + 1, ai, j + jcod])
    
                    if len(nu_syn[0]):
                        hist_syn += np.histogram2d(nu_syn[0], nu_syn[1], bins=[bins_from, bins_to])[0]
                    if len(nu_nonsyn[0]):
                        hist_nonsyn += np.histogram2d(nu_nonsyn[0], nu_nonsyn[1], bins=[bins_from, bins_to])[0]

    if plot:
        fig, ax = plt.subplots()
        z = hist[:, 1:-1]
        z  = (z.T / z.sum(axis=1)).T
        z = np.log10(z)
        im = ax.imshow(z.T, interpolation='nearest')
        ax.set_xlabel('Initial freq')
        ax.set_ylabel('Final freq')
        ax.set_xticks(np.arange(len(bins_from) - 1))
        ax.set_yticks(np.arange(len(bins_to) - 1))
        ax.set_xticklabels(map('{:1.1e}'.format, np.sqrt(bins_from[:-1] * bins_from[1:])), rotation=45, fontsize=10)
        ax.set_yticklabels(map('{:1.1e}'.format, np.sqrt(bins_to[1:-2] * bins_to[2:-1])), rotation=45, fontsize=10)
        ax.set_xlim(0.5, len(bins_from) - 1.5)
        ax.set_ylim(len(bins_to) - 3.5, -0.5)
        plt.colorbar(im)
        ax.set_title('Propagator for allele frequencies\nbetween consecutive time points, '+\
                     gene+'\n[log10 P(x1 | x0)]', fontsize=16)
        plt.tight_layout()

        fig_syn, ax = plt.subplots()
        z = hist_syn[:, 1:-1]
        z  = (z.T / z.sum(axis=1)).T
        z = np.log10(z)
        im = ax.imshow(z.T, interpolation='nearest')
        ax.set_xlabel('Initial freq')
        ax.set_ylabel('Final freq')
        ax.set_xticks(np.arange(len(bins_from) - 1))
        ax.set_yticks(np.arange(len(bins_to) - 1))
        ax.set_xticklabels(map('{:1.1e}'.format, np.sqrt(bins_from[:-1] * bins_from[1:])), rotation=45, fontsize=10)
        ax.set_yticklabels(map('{:1.1e}'.format, np.sqrt(bins_to[1:-2] * bins_to[2:-1])), rotation=45, fontsize=10)
        ax.set_xlim(0.5, len(bins_from) - 1.5)
        ax.set_ylim(len(bins_to) - 3.5, -0.5)
        plt.colorbar(im)
        ax.set_title('Propagator for allele frequencies\nbetween consecutive time points, '+\
                     '-'.join(genes)+'\n[log10 P(x1 | x0)], syn muts', fontsize=16)
        plt.tight_layout()

        fig_nonsyn, ax = plt.subplots()
        z = hist_nonsyn[:, 1:-1]
        z  = (z.T / z.sum(axis=1)).T
        z = np.log10(z)
        im = ax.imshow(z.T, interpolation='nearest')
        ax.set_xlabel('Initial freq')
        ax.set_ylabel('Final freq')
        ax.set_xticks(np.arange(len(bins_from) - 1))
        ax.set_yticks(np.arange(len(bins_to) - 1))
        ax.set_xticklabels(map('{:1.1e}'.format, np.sqrt(bins_from[:-1] * bins_from[1:])), rotation=45, fontsize=10)
        ax.set_yticklabels(map('{:1.1e}'.format, np.sqrt(bins_to[1:-2] * bins_to[2:-1])), rotation=45, fontsize=10)
        ax.set_xlim(0.5, len(bins_from) - 1.5)
        ax.set_ylim(len(bins_to) - 3.5, -0.5)
        plt.colorbar(im)
        ax.set_title('Propagator for allele frequencies\nbetween consecutive time points, '+\
                     '-'.join(genes)+'\n[log10 P(x1 | x0)], nonsyn muts', fontsize=16)
        plt.tight_layout()

        plt.ion()
        plt.show()

        # Save figs to tmp folder
        fig_syn.savefig('/ebio/ag-neher/home/fzanini/tmp/propagator_'+'-'.join(genes)+'_syn.png')
        fig_nonsyn.savefig('/ebio/ag-neher/home/fzanini/tmp/propagator_'+'-'.join(genes)+'_nonsyn.png')

# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/05/14
content:    Plot site frequency spectra for derived alleles.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet.IUPAC import unambiguous_dna

from hivwholeseq.miseq import alpha
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file, \
        get_allele_counts_insertions_from_file_unfiltered, \
        filter_nus
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories as plot_nus
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories_3d as plot_nus_3d
from hivwholeseq.patients.one_site_statistics import get_allele_frequency_trajectories



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele frequency trajectories')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--genes', default=('gag', 'pol'), nargs='+',
                        help='Gene to analyze (e.g. pol)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele frequency trajectories to file')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--PCR1', action='store_true',
                        help='Show only PCR1 samples where possible (still computes all)')

    args = parser.parse_args()
    pname = args.patient
    genes = map(lambda x: x.lower(), args.genes)
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot
    submit = args.submit
    use_PCR1 = args.PCR1

    patient = get_patient(pname)
    times = patient.times()
    samplenames = patient.samples
    # Keep PCR2 only if PCR1 is absent
    if use_PCR1:
        ind = np.nonzero(map(lambda x: ('PCR1' in x[1]) or ((times == times[x[0]]).sum()==1), enumerate(samplenames)))[0]
        times = times[ind]

    fragments = ['F'+str(i) for i in xrange(1, 7)]

    # Prepare output data structures
    bins = np.logspace(-2, -0.5, 11)
    binsc = np.sqrt(bins[1:] * bins[:-1])
    hist = np.zeros(len(bins) - 1, float)
    hist_syn = np.zeros_like(hist)
    hist_nonsyn = np.zeros_like(hist)

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

            # Get rid of gaps and low-coverage regions
            is_gap = ((aft.argmax(axis=1) == 6) | (act.sum(axis=1) < 100)).any(axis=0)
            if VERBOSE >= 2:
                print 'Fraction of gap sites (excluded):', is_gap.mean()

            if use_PCR1:
                aft = aft[ind]
                act = act[ind]

            # Get rid of ancestral alleles (not superfluous, because they might
            # become minority!) -- here we assume initial consensus as ancestral
            # for everybody later on, which is admittedly a bit so so.
            aft_der = aft.copy()
            aft_der[:, :, is_gap] = 0
            for i, ai in enumerate(aft[0].argmax(axis=0)):
                    aft_der[:, ai, i] = 0

            # Extract gene from reference and allele freq trajectories
            from hivwholeseq.store.store_initial_reference import check_genes_consensus
            cons_rec = SeqIO.read(get_initial_reference_filename(pname, fragment), 'fasta')
            conss = str(cons_rec.seq)
            gene_seqs, genes_good, gene_poss = check_genes_consensus(conss, fragment, genes=[gene],
                                                                     VERBOSE=VERBOSE)
            if not genes_good[gene]:
                if VERBOSE >= 1:
                    print pname+':', gene, 'not found or with problems: skipping.'
                continue
    
            gene_pos = gene_poss[gene]
            aft_der_gene = np.concatenate([aft_der[:, :, exon_pos[0]: exon_pos[1]]
                                           for exon_pos in gene_pos], axis=2)
            conss_gene = gene_seqs[gene]
            gene_len = len(conss_gene)

            hist += np.histogram(aft_der_gene.ravel(), bins=bins, density=False)[0]

            # Collect counts syn/nonsyn
            nu_syn = []
            nu_nonsyn = []
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
    
                        aftmp = aft_der_gene[:, ai, j + jcod]
                        aftmp = aftmp[(aftmp >= bins[0]) & (aftmp <= bins[-1])]
                        if not len(aftmp):
                            continue

                        if str(cod_new.toseq().translate()) != str(cod_anc.toseq().translate()):
                            nu_syn.extend(aftmp)
                        else:
                            nu_nonsyn.extend(aftmp)
    
            if len(nu_syn):
                hist_syn += np.histogram(nu_syn, bins=bins)[0]
            if len(nu_nonsyn):
                hist_nonsyn += np.histogram(nu_nonsyn, bins=bins)[0]


    # Normalize
    hist_norm = hist.copy()
    hist_norm /= hist_norm.sum()
    hist_norm /= bins[1:] - bins[:-1]

    hist_syn_norm = hist_syn.copy()
    hist_syn_norm /= hist_syn_norm.sum()
    hist_syn_norm /= bins[1:] - bins[:-1]

    hist_nonsyn_norm = hist_nonsyn.copy()
    hist_nonsyn_norm /= hist_nonsyn_norm.sum()
    hist_nonsyn_norm /= bins[1:] - bins[:-1]

    if plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot(binsc, hist_norm, lw=2, c='k', label='all')
        ax.plot(binsc, hist_syn_norm, lw=2, c='r', label='syn')
        ax.plot(binsc, hist_nonsyn_norm, lw=2, c='b', label='nonsyn')

        # Theoretical curves
        aind = -1
        alpha = hist_norm[aind]
        ax.plot([binsc[0], binsc[-1]], [alpha * ((binsc[aind] / binsc[0])**2),
                                        alpha * ((binsc[aind] / binsc[-1])**2)], lw=2, c='grey')
        ax.plot([binsc[0], binsc[-1]], [alpha * ((binsc[aind] / binsc[0])),
                                        alpha * ((binsc[aind] / binsc[-1]))], lw=2, c='grey')

        ax.set_xlabel('Freq')
        ax.set_ylabel('SFS [density = counts / sum / binsize]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(pname+', '+'-'.join(genes))
        ax.grid(True)
        ax.legend(loc=1, fontsize=12)

        plt.tight_layout()
        plt.ion()
        plt.show()



# vim: fdm=marker
'''
author:     Fabio Zanini
date:       04/11/13
content:    Study PCR amplification bias in the plasmid mixes.

            NOTE: we have found out already that, within a single fragment,
            one haplotype dominates.
'''
# Modules
import argparse
import numpy as np
from collections import Counter
from operator import itemgetter, attrgetter
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment as MSA
import pysam
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alphal
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.filenames import get_consensus_filename, get_mapped_filename, \
        get_allele_counts_filename, get_coverage_filename
from hivwholeseq.mapping_utils import align_muscle, pair_generator
from hivwholeseq.one_site_statistics import filter_nus
from hivwholeseq.coverage_tuples import get_coverage_tuples
from hivwholeseq.reference_samples.mix_recombination_PCR import align_consensi_mix1, \
        check_consensus_mix1, align_consensi_mix2, check_consensus_mix2


# Globals
mix_adaIDs = {'Tue28': {'mix1': 'TS18', 'mix2': 'TS19'},
              'Tue42': {'mix1': 'N1-S1', 'mix2': 'N3-S3'}}

mix1_references_adaIDs = [(0.5, 'TS2'), (0.5, 'TS4')]
mix2_references_adaIDs = [(0.045, 'TS2'), (0.95, 'TS4'), (0.005, 'TS7')]

colors = {'NL4-3': cm.jet(int(255.0 * 0 / 3)),
          'SF162': cm.jet(int(255.0 * 1 / 3)),
          'F10': cm.jet(int(255.0 * 2 / 3))}



# Functions
def get_ind_private_alleles_nogaps_mix1(ali):
    '''Get the positions of private alleles in the alignment'''
    poss = (((ali[0] == ali[1]) != (ali[0] == ali[2])) & \
            (ali != '-').all(axis=0)).nonzero()[0]
    return poss


def get_ind_private_alleles_nogaps_mix2(ali):
    '''Get the positions of private alleles in the alignment'''
    poss = ((ali[1] != ali[2]) & (ali[1] != ali[3]) & (ali[2] != ali[3]) &\
            ((ali[0] == ali[1:]).sum(axis=0) == 1) & \
            (ali != '-').all(axis=0)).nonzero()[0]
    return poss


def logit(y):
    return np.log10(y / (1-y))


def logit_rev(l):
    return 1.0 / (1 + 10**(-l))


def format_labels(ynu):
    labels = []
    for y in ynu:
        if 0.4 < y < 0.6:
            labels.append(r'$'+str(y)+r'$')
        elif y < 0.4:
            labels.append(r'$10^{'+'{:1.0f}'.format(np.log10(y))+r'}$')
        else:
            labels.append(r'$1 - 10^{'+'{:1.0f}'.format(np.log10(1 - y))+r'}$')

    return labels



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Errors in HIV pure strains')
    parser.add_argument('--run', required=True,
                        help='MiSeq run to analyze (e.g. Tue28)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    seq_run = args.run
    fragments = args.fragments
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    for fragment in fragments:

        # Prepare plots
        fig, axs = plt.subplots(2, 1, figsize=(7, 13))

        # MIX1
        ax = axs[0]
        adaID = mix_adaIDs[seq_run]['mix1']

        # Check consensus
        alignment = align_consensi_mix1(seq_run, fragment)
        ali = np.array(alignment)
        check_consensus_mix1(seq_run, fragment, alignment=alignment)

        # Check allele frequencies of the two strains (there is PCR-mediated
        # recombination - watch out!)
        counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
        coverage = np.load(get_coverage_filename(data_folder, adaID, fragment))
        nus = filter_nus(counts, coverage)
        pos_poly = get_ind_private_alleles_nogaps_mix1(ali)
        nus_ali = np.zeros((2, len(pos_poly)))
        for i, pos in enumerate(pos_poly):
            pos_self = (ali[0, :pos] != '-').sum()
            nus_ali[0, i] = nus[alphal.index(ali[1, pos]), pos_self]
            nus_ali[1, i] = nus[alphal.index(ali[2, pos]), pos_self]

        # Plot
        for j in xrange(2):
            c = colors[alignment[j+1].name]
            ax.plot(pos_poly, nus_ali[j], color=c,
                    label=alignment[j+1].name, lw=2)
            ax.scatter(pos_poly, nus_ali[j],
                       s=60, color=c, edgecolor='none')
        ax.set_xlabel('Position [bases]')
        ax.set_ylabel(r'$\nu$', fontsize=18)
        ax.set_ylim(-0.1, 1.1)
        ax.set_title('Mix1, '+str(fragment), fontsize=18)
        ax.legend(loc='upper center', fontsize=14)


        # MIX2
        ax = axs[1]
        adaID = mix_adaIDs[seq_run]['mix2']

        # Check consensus
        alignment = align_consensi_mix2(seq_run, fragment)
        ali = np.array(alignment)
        check_consensus_mix2(seq_run, fragment, alignment=alignment)

        # Check allele frequencies of the three strains (there is PCR-mediated
        # recombination - watch out!)
        counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
        coverage = np.load(get_coverage_filename(data_folder, adaID, fragment))
        nus = filter_nus(counts, coverage)
        pos_poly = get_ind_private_alleles_nogaps_mix2(ali)
        nus_ali = np.zeros((3, len(pos_poly)))
        for i, pos in enumerate(pos_poly):
            pos_self = (ali[0, :pos] != '-').sum()
            nus_ali[0, i] = nus[alphal.index(ali[1, pos]), pos_self]
            nus_ali[1, i] = nus[alphal.index(ali[2, pos]), pos_self]
            nus_ali[2, i] = nus[alphal.index(ali[3, pos]), pos_self]

        # Plot
        for j in xrange(3):
            y = nus_ali[j]
            y = np.abs(y - 1e-5)
            y = logit(y)
            c = colors[alignment[j+1].name]
            # Labels stuff
            adaID = alignment[j+1].id.split('_')[1]
            a = zip(*mix2_references_adaIDs)
            label = alignment[j+1].name+', '+str(100 * a[0][a[1].index(adaID)])+'%'
            ax.plot(pos_poly, y, color=c,
                    label=label, lw=2)
            ax.scatter(pos_poly, y, color=c, edgecolor='none',
                       s=80)
        ax.set_xlabel('Position [bases]')
        ax.set_ylabel(r'$\nu$', fontsize=18)
        ax.set_ylim(-6, +6)
        ynu = np.array([1e-6, 1e-4, 1e-2, 0.5, 1-1e-2, 1-1e-4, 1-1e-6])
        yticks = logit(ynu)
        yticklabels = format_labels(ynu)
        plt.yticks(yticks, yticklabels)
        ax.set_title('Mix2, '+str(fragment), fontsize=18)
        ax.legend(loc='upper center', fontsize=14)

        plt.tight_layout()
        plt.ion()
        plt.show()

        # NOTE
        savefig = False
        if savefig:
            fig.savefig('/ebio/ag-neher/home/fzanini/phd/sequencing/figures/'+\
                        'PCR_amplification_bias_'+str(fragment)+'.png')

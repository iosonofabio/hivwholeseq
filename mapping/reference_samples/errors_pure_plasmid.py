# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Check sequencing/PCR erros in pure plasmid HIV samples (they appear
            as minor variants).
'''
# Modules
import sys
import argparse
from collections import defaultdict
from itertools import izip
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

from mapping.datasets import MiSeq_runs
from mapping.miseq import alpha
from mapping.filenames import get_NL43_entire, get_NL43_fragmented, \
        get_F10_entire, get_F10_fragmented, \
        get_consensus_filename, get_allele_counts_filename, get_coverage_filename
from mapping.mapping_utils import align_muscle



# Globals
#             run adaID  sample  ref_function     fragmented_function
references = {(28, 2): ('NL4-3', get_NL43_entire, get_NL43_fragmented),
              (28, 7): ('F10', get_F10_entire, get_F10_fragmented)}
colors = {'fwd': 'b', 'rev': 'g'}
figures_folder = '/ebio/ag-neher/home/fzanini/phd/sequencing/figures/'



# Functions
def get_count(counts):
    '''Divide counts in fwd/rev'''
    count = {'fwd': counts[0] + counts[2], 'rev': counts[1] + counts[3]}
    return count


def get_consensus(count):
    '''Get fwd/rev consensi'''
    consensus = {key: alpha[[cou.argmax() for cou in count[key].T]] for key in count}
    return consensus


def get_minor_counts(count):
    '''Get fwd/rev minor counts'''
    minor_counts = {key: np.zeros(count[key].shape[-1], int) for key in count}
    for key in count:
        for i in xrange(len(minor_counts['fwd'])):
            cou = np.sort(count[key][:, i])[-2]
            minor_counts[key][i] = cou
    return minor_counts


def get_minor_nus(count):
    '''Get fwd/rev minor frequencies'''
    minor_nus = {key: np.zeros(count[key].shape[-1]) for key in count}
    for key in count:
        for i in xrange(len(minor_nus['fwd'])):
            cou = np.sort(count[key][:, i])[-2]
            minor_nus[key][i] = 1.0 * cou / (count[key][:, i].sum() + 1e-6)
    return minor_nus


def get_coverage(count):
    '''Get fwd/rev coverage'''
    coverage = {key: count[key].sum(axis=0) for key in count}
    return coverage


def consensus_vs_reference(miseq_run, adaID, fragment, VERBOSE=0):
    '''Check the consensus sequence VS the phiX reference'''
    # There are two seta of primer for F5
    if fragment == 'F5':
        frag_ref = MiSeq_runs[miseq_run]['primerF5'][MiSeq_runs[miseq_run]['adapters'].index(adaID)]
    else:
        frag_ref = fragment

    data_folder = MiSeq_runs[miseq_run]['folder']

    # Calculate consensus
    consensus_filename = get_consensus_filename(data_folder, adaID, fragment,
                                                trim_primers=True)
    consensus = SeqIO.read(consensus_filename, 'fasta')

    # Get reference
    sample, _, get_fragmented = references[(miseq_run, adaID)]
    refseq = SeqIO.read(get_fragmented(frag_ref, trim_primers=True), 'fasta')

    # Align them
    align = align_muscle(refseq, consensus)
    ali = np.array(align)
    
    errs = ali[0] != ali[1]
    if errs.any():
        print sample, 'fragment', fragment, 'ref len:', len(refseq)
        print 'pos mutation'
        print '---------------------'
        for pos in errs.nonzero()[0]:
            print pos, ali[0][pos], '->', ali[1][pos]


def minor_alleles_along_genome(miseq_run, adaID, fragment, VERBOSE=0, plot=False,
                               savefig=False):
    '''Show the minor alleles along the phiX genome'''
    # Get the counts
    sample = references[(miseq_run, adaID)][0]
    data_folder = MiSeq_runs[miseq_run]['folder']
    counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))

    # Study the minor alleles (sequencing errors)
    count = get_count(counts)
    minor_counts = get_minor_counts(count)
    minor_nus = get_minor_nus(count)
    cov = get_coverage(count)

    # Mean error frequency
    nu_err_mean = {}
    print 'Average error frequency:',
    for key in count:
        nue = minor_counts[key].mean() / cov[key].mean()
        nu_err_mean[key] = nue
        print key, '{:1.1e}'.format(nue),
    print ''

    # Plot
    if plot:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(2, 1, figsize=(7, 13))
        for key in count:
            axs[0].plot(minor_counts[key], lw=2, c=colors[key], alpha=0.5)
            axs[0].scatter(np.arange(len(minor_counts[key])), minor_counts[key],
                           s=20, c=colors[key], edgecolor='none',
                           label=key+': '+'{:1.1e}'.format(nu_err_mean[key]))
        #axs[0].set_xlabel('Position [bases]')
        axs[0].set_title('Error counts')
        axs[0].set_xlim(-100, len(minor_counts['fwd']) + 100)
        axs[0].set_yscale('log')
        axs[0].legend()
    
        for key in count:
            # Plot coverage first 
            axs[1].plot(1.0 / count[key].sum(axis=0), lw=0.5,
                        c=colors[key], ls='--', label='coverage '+key)
            axs[1].plot(minor_nus[key], lw=2, c=colors[key], alpha=0.5)
            axs[1].scatter(np.arange(len(minor_nus[key])), minor_nus[key],
                           s=20, c=colors[key], edgecolor='none')
        axs[1].set_xlabel('Position [bases]')
        axs[1].set_title('Error frequencies')
        axs[1].set_xlim(-100, len(minor_counts['fwd']) + 100)
        axs[1].set_yscale('log')
        axs[1].set_ylim(1e-5, 1e-1)
    
    
        plt.suptitle(sample+', fragment '+fragment, fontsize=20)
        plt.tight_layout(rect=(0, 0, 1, 0.95))
    
        #plt.ion()
        #plt.show()

        if savefig:
            fig.savefig(figures_folder+sample+'_minornu_'+fragment+'.png')


def spikes_motifs(miseq_run, adaID, fragment, VERBOSE=0, plot=False, savefig=False):
    '''Find the motifs around the spikes in error rates'''
    # Get the counts
    sample = references[(miseq_run, adaID)][0]
    data_folder = MiSeq_runs[miseq_run]['folder']
    counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))

    # Study the minor alleles (sequencing errors)
    count = get_count(counts)
    minor_counts = get_minor_counts(count)
    minor_nus = get_minor_nus(count)
    cov = get_coverage(count)
    consensus = get_consensus(count)

    # Pick the spikes, i.e. the positions with error rates in the tail of the
    # binomial, p < 0.001 or so
    from scipy.stats import binom
    spikes = defaultdict(list)
    for key in count:
        # Average error frequency
        p = 1.0 * minor_counts[key].mean() / cov[key].mean()
        k_crit = binom.ppf(0.999, cov[key], p)
        ind = (minor_counts[key] > k_crit).nonzero()[0]
        if len(ind):
            if VERBOSE:
                print 'Spikes: run', miseq_run, 'adaID', adaID, 'fragment', fragment, '('+key+')'
                print ' Pos WT Motifs'
                print '----------------------'
            for pos in ind:
                if VERBOSE:
                    print '{:4d}'.format(pos),
                    print '{:>2s}'.format(consensus[key][pos]),
                for le in xrange(1, 5):
                    if pos > le - 2:
                        mutm = list(consensus[key][pos - le + 1: pos + 1])
                        mutm[le - 1] = alpha[count[key][:, pos] == minor_counts[key][pos]][0]
                        if VERBOSE:
                            print ''.join(mutm),
                if VERBOSE:
                    print
                spikes[key].append((pos,
                                    (''.join(consensus[key][pos - le + 1: pos + 1]),
                                     ''.join(mutm)),
                                     1.0 * minor_counts[key][pos] / cov[key][pos]))

    # Order spikes by frequency
    from operator import itemgetter
    for key in spikes:
        spikes[key].sort(key=itemgetter(2), reverse=True)

    # Plot the error histogram and the binomial if requested
    if plot:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(2, 1, figsize=(7, 12))
        xf = lambda x: np.logspace(-5, np.log10(0.5), x)

        for i, key in enumerate(count):
            ax = axs[i]
            ax.hist(1.0 * minor_counts[key] / cov[key],
                    bins=xf(60), normed=True, label='Errors')
            for n, c in izip([2e4, 1e5], ('r', 'g')):
                yb = binom.pmf(np.floor(n * xf(200)), n, 1e-4)
                ax.plot(xf(200), yb / yb.max() * 0.9 * ax.get_ylim()[1], lw=2,
                        c=c,
                        label='Binomial ('+str(n)+')')
            ax.set_xlabel('Error frequency')
            ax.set_xscale('log')
            ax.set_xlim(1e-5, 0.5)
            ax.set_ylabel('Density (A.U.)')
            ax.legend(loc=1, fontsize=12)
            ax.set_title(key, fontsize=18)
    
            ym = 0.3 * ax.get_ylim()[1]
            for j, (pos, trans, nu_spike) in enumerate(spikes[key][:10]):
                ar = ax.arrow(nu_spike,
                              ym * 1.08**j, 0, -(ym * 1.08**j - 200),
                              edgecolor='k',
                              facecolor='k',
                              width=nu_spike / 100,
                              head_length=100,
                              overhang=0.4)
                txt = ax.text(nu_spike * 1.05, ym * 1.05 * 1.08**j,
                              str(pos)+': '+trans[0][-1]+' -> '+trans[1][-1],
                              fontsize=12)

        plt.suptitle(sample+', fragment '+fragment+', spikes', fontsize=20)
        plt.tight_layout(rect=(0, 0, 1, 0.97))
        plt.ion()
        plt.show()

        if savefig:
            fig.savefig(figures_folder+sample+'_spikes_'+fragment+'.png')

    return spikes



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Errors in HIV pure strains')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    for adaID in adaIDs:
        for fragment in fragments:
            #consensus_vs_reference(miseq_run, adaID, fragment, VERBOSE=VERBOSE)
            minor_alleles_along_genome(miseq_run, adaID, fragment,
                                       VERBOSE=VERBOSE, plot=True, savefig=False)
            #spikes_motifs(miseq_run, adaID, fragment,
            #              VERBOSE=VERBOSE, plot=True, savefig=False)

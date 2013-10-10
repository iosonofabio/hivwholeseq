# vim: fdm=marker
'''
author:     Fabio Zanini
date:       10/10/13
content:    Characterize errors in phiX.
'''
# Modules
import argparse
from collections import defaultdict, Counter
import numpy as np
from Bio import SeqIO

from mapping.datasets import MiSeq_runs
from mapping.miseq import alpha, read_types
from mapping.filenames import get_phix_filename, get_allele_counts_phix_filename



# Functions
def load_allele_counts(data_folder, VERBOSE=0):
    '''Load the precomputer allele counts'''
    allele_count_filename = get_allele_counts_phix_filename(data_folder)
    counts = np.load(allele_count_filename)
    return counts


def consensus_vs_reference(miseq_run, counts=None, VERBOSE=0):
    '''Check the consensus sequence VS the phiX reference'''
    if counts is None:
        counts = load_allele_counts(datasets[miseq_run]['folder'], VERBOSE=VERBOSE)

    # Calculate consensus
    count = counts.sum(axis=0)
    consensus = alpha[[cou.argmax() for cou in count.T]]

    if VERBOSE >= 1:
        # Get reference
        refseq = SeqIO.read(get_phix_filename(), 'fasta')
        ref = np.array(refseq)

        # Check differences
        ind = (consensus != ref).nonzero()[0]
        if len(ind):
            print 'PhiX run', miseq_run, 'differences from reference:'
            print 'Pos  Ref   Cons  Mutation'
            print '-------------------------'
            for i in ind:
                i0 = max(i-2, 0)
                i1 = min(i+3, len(refseq))
                print '{:4d}'.format(i), ''.join(ref[i0: i1]), ''.join(consensus[i0: i1]),\
                        ref[i], '->', consensus[i]
            print ''

    return ''.join(consensus)


def minor_alleles_along_genome(miseq_run, counts=None, VERBOSE=0, plot=False):
    '''Show the minor alleles along the phiX genome'''
    if counts is None:
        counts = load_allele_counts(datasets[miseq_run]['folder'], VERBOSE=VERBOSE)

    # Study the minor alleles (sequencing errors)
    count = counts.sum(axis=0)
    minor_counts = np.zeros(counts.shape[-1], int)
    minor_nus = np.zeros(counts.shape[-1])
    for i in xrange(len(minor_counts)):
        cou = np.sort(count[:, i])[-2]
        minor_counts[i] = cou
        minor_nus[i] = 1.0 * cou / (count[:, i].sum() + 1e-6)

    # Plot
    if plot:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(13, 7))
        axs[0].plot(minor_counts, lw=2)
        axs[0].scatter(np.arange(len(minor_counts)), minor_counts, s=20, c='b',
                       edgecolor='none')
        axs[0].set_xlabel('Position in phiX')
        axs[0].set_title('Sequencing errors')
        axs[0].set_xlim(-100, len(minor_counts) + 100)
        axs[0].set_yscale('log')
    
        # Plot coverage first
        axs[1].plot(1.0 / count.sum(axis=0), lw=1.5, c='grey', label='coverage')
    
        axs[1].plot(minor_nus, lw=2, c='b')
        axs[1].scatter(np.arange(len(minor_nus)), minor_nus, s=20, c='b',
                       edgecolor='none')
        axs[1].set_xlabel('Position in phiX')
        axs[1].set_title('Sequencing errors frequencies')
        axs[1].set_xlim(-100, len(minor_counts) + 100)
        axs[1].set_yscale('log')
    
    
        plt.suptitle('PhiX analysis: run '+str(miseq_run), fontsize=20)
        plt.tight_layout(rect=(0, 0, 1, 0.95))
    
        plt.ion()
        plt.show()

    return minor_counts, minor_nus


def characterize_motifs(miseq_run, counts=None, consensus=None, minor_counts=None, VERBOSE=0):
    '''Find the DNA motifs of the errors'''
    if counts is None:
        counts = load_allele_counts(datasets[miseq_run]['folder'], VERBOSE=VERBOSE)
    count = counts.sum(axis=0)

    if consensus is None:
        consensus = alpha[[cou.argmax() for cou in count.T]]
        consensus = ''.join(consensus)

    if minor_counts is None: 
        minor_counts = np.zeros(counts.shape[-1], int)
        for i in xrange(len(minor_counts)):
            minor_counts[i] = np.sort(count[:, i])[-2]

    pos_errors = (minor_counts > 0).nonzero()[0]
    if not len(pos_errors):
        return

    # Note: skip the first 300 positions, they are strange somehow (maybe mapping)
    pos_errors = pos_errors[pos_errors >= 300]

    motifs = [defaultdict(int) for l in xrange(1, 5)]
    for pos in pos_errors:
        muts = alpha[count[:, pos] == minor_counts[pos]]
        for mut in muts:
            # 1-mer motifs
            mot1 = (consensus[pos], mut)
            motifs[0][mot1] += 1

            # 2+ mer motifs
            for le in xrange(2, 5):
                if pos > le - 2:
                    wt = consensus[pos-le+1:pos+1]
                    mutm = list(wt)
                    mutm[le-1] = mut
                    mot = (wt, ''.join(mutm))
                motifs[le-1][mot] += 1

    # Normalize by nucleotide content in the consensus
    # Note: for longer motifs, we renormalize by all nucleotides in the consensus motif
    cons = np.array(list(consensus))
    abus = {a: (cons == a).mean() for a in alpha}
    for mots in motifs:
        for key in mots:
            for j in xrange(len(key[0])):
                mots[key] = 1.0 * mots[key] / (abus[key[0][j]] / 0.25)

    # Sort by abundance
    from operator import itemgetter
    motifs = [sorted(mots.iteritems(), key=itemgetter(1), reverse=True)
              for mots in motifs]

    if VERBOSE >= 1:
        print 'Motifs: run', miseq_run
        print 'Motif occurrences (rinormalized by abundance)'
        for i in xrange(len(motifs)):
            for (m, c) in motifs[i][:5]:
                print m[0], '->', m[1], '{:10f}'.format(c)

    return motifs


def spikes_motifs(miseq_run, counts=None, consensus=None, minor_counts=None, VERBOSE=0):
    '''Find the motifs around thes pikes in error rates'''
    if counts is None:
        counts = load_allele_counts(datasets[miseq_run]['folder'], VERBOSE=VERBOSE)
    count = counts.sum(axis=0)

    if consensus is None:
        consensus = alpha[[cou.argmax() for cou in count.T]]
        consensus = ''.join(consensus)

    if minor_counts is None: 
        minor_counts = np.zeros(counts.shape[-1], int)
        for i in xrange(len(minor_counts)):
            minor_counts[i] = np.sort(count[:, i])[-2]

    # Pick the spikes with p < 0.001 or so
    # Note: skip the first 300 positions, they are strange somehow (maybe mapping)
    n = count[:, 300:].sum(axis=0).mean()
    p = 1.0 * minor_counts[300:].mean() / n
    from scipy.stats import binom
    k_crit = binom.ppf(0.999, n, p)
    ind = (minor_counts[300:] > k_crit).nonzero()[0]
    spikes = []
    if len(ind):
        ind += 300
        print 'Spikes: run', miseq_run
        print ' Pos WT Motifs'
        print '----------------------'
        for pos in ind:
            print '{:4d}'.format(pos),
            print '{:>2s}'.format(consensus[pos]),
            for le in xrange(1, 5):
                if pos > le - 2:
                    mutm = list(consensus[pos - le + 1: pos + 1])
                    mutm[le - 1] = alpha[count[:, pos] == minor_counts[pos]][0]
                    print ''.join(mutm),
            print
            spikes.append((pos, (consensus[pos - le + 1: pos + 1], ''.join(mutm))))

    return spikes


def minor_nus_crossrun(data_runs, VERBOSE=0):
    '''Check errors across runs'''
    runs = sorted(data_runs.keys())
    minor_nus = {key: data_runs[key]['minor_nus'] for key in runs}
    
    # Take all pairs of runs and compare
    # Note: skip the first 300 positions, they are strange somehow (maybe mapping)
    import matplotlib.pyplot as plt
    for i, run1 in enumerate(runs[:-1]):
        nus1 = minor_nus[run1][300:]
        for run2 in runs[i+1:]:
            nus2 = minor_nus[run2][300:]

            fig, ax = plt.subplots(1, 1)
            ax.plot([1e-5, 0.5], [1e-5, 0.5], lw=1, c='grey')
            ax.scatter(nus1, nus2, s=20)
            ax.set_xlabel('run '+str(run1))
            ax.set_ylabel('run '+str(run2))
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(1e-5, 1e0)
            ax.set_ylim(1e-5, 1e0)
            ax.set_title('Minor allele frequencies in phiX')

    plt.tight_layout()
    plt.ion()
    plt.show()

def spikes_crossrun(data_runs, VERBOSE=0):
    '''Check the spikes across runs'''
    runs = sorted(data_runs.keys())
    spikes = {key: data_runs[key]['spikes'] for key in runs}

    # Take all pairs of runs and compare
    # Note: skip the first 300 positions, they are strange somehow (maybe mapping)
    for i, run1 in enumerate(runs[:-1]):
        spike1 = spikes[run1]
        for run2 in runs[i+1:]:
            spike2 = spikes[run2]

            # Track shared spikes
            shared = set(spike1) & set(spike2)
            if len(shared):
                print 'Shared spikes: runs ', run1, 'and', run2
                for share in shared:
                    print share,
                print





# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--runs', type=int, nargs='+', required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot results')

    args = parser.parse_args()
    miseq_runs = args.runs
    VERBOSE = args.verbose
    plot = args.plot

    # Crossrun data structures
    data_runs = defaultdict(dict)

    # Iterate over requested runs
    for miseq_run in miseq_runs:
    
        # Specify the dataset
        dataset = MiSeq_runs[miseq_run]
        data_folder = dataset['folder']
    
        # Get allele counts
        counts = load_allele_counts(data_folder, VERBOSE=VERBOSE)
        data_runs[miseq_run]['counts'] = counts
    
        # Check consensus
        consensus = consensus_vs_reference(miseq_run, counts, VERBOSE=0)
        data_runs[miseq_run]['consensus'] = consensus

        # Along genome
        minor_counts, minor_nus = minor_alleles_along_genome(miseq_run, counts, VERBOSE=VERBOSE,
                                                             plot=plot)
        data_runs[miseq_run]['minor_counts'] = minor_counts
        data_runs[miseq_run]['minor_nus'] = minor_nus

        # Error frequency
        print 'Average error frequency:', '{:1.1e}'.format(minor_counts[300:].mean() / counts.sum(axis=0).sum(axis=0).mean())

        # Motifs
        motifs = characterize_motifs(miseq_run, counts, minor_counts=minor_counts, VERBOSE=VERBOSE)
        data_runs[miseq_run]['motifs'] = motifs

        # Spikes
        spikes = spikes_motifs(miseq_run, counts, minor_counts=minor_counts, VERBOSE=VERBOSE)
        data_runs[miseq_run]['spikes'] = spikes

    # Check consistency across runs
    if len(data_runs) > 1:
        minor_nus_crossrun(data_runs)
        spikes_crossrun(data_runs)

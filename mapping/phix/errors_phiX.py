# vim: fdm=marker
'''
author:     Fabio Zanini
date:       10/10/13
content:    Characterize errors in phiX.

            Most can be done on allele frequencies without going down to the
            reads, but fwd and rev must be separated.
'''
# Modules
import os
import pysam
import argparse
from operator import itemgetter
from collections import defaultdict, Counter
import numpy as np
from Bio import SeqIO
from Bio.Seq import reverse_complement as rc

from mapping.datasets import MiSeq_runs
from mapping.miseq import alpha, read_types, alphal
from mapping.filenames import get_phix_filename, get_allele_counts_phix_filename, \
        get_mapped_phix_filename
from mapping.mapping_utils import get_ind_good_cigars



# Globals
match_len_min = 10
trim_bad_cigars = 3
rcf = {'fwd': lambda x: x, 'rev': rc}
colors = {'fwd': 'b', 'rev': 'g'}



# Functions
def get_count(counts):
    count = {'fwd': counts[0] + counts[2], 'rev': counts[1] + counts[3]}
    return count


def get_consensus(count):
    consensus = {key: alpha[[cou.argmax() for cou in count[key].T]] for key in count}
    return consensus


def get_minor_counts(count):
    minor_counts = {key: np.zeros(count[key].shape[-1], int) for key in count}
    for key in count:
        for i in xrange(len(minor_counts['fwd'])):
            cou = np.sort(count[key][:, i])[-2]
            minor_counts[key][i] = cou
    return minor_counts


def get_minor_nus(count):
    minor_nus = {key: np.zeros(counts.shape[-1]) for key in count}
    for key in count:
        for i in xrange(len(minor_nus['fwd'])):
            cou = np.sort(count[key][:, i])[-2]
            minor_nus[key][i] = 1.0 * cou / (count[key][:, i].sum() + 1e-6)
    return minor_nus


def get_coverage(count):
    coverage = {key: count[key].sum(axis=0) for key in count}
    return coverage


def load_allele_counts(data_folder, VERBOSE=0):
    '''Load the precomputer allele counts'''
    allele_count_filename = get_allele_counts_phix_filename(data_folder)
    counts = np.load(allele_count_filename)
    return counts


def consensus_vs_reference(miseq_run, counts=None, VERBOSE=0):
    '''Check the consensus sequence VS the phiX reference'''
    if counts is None:
        counts = load_allele_counts(MiSeq_runs[miseq_run]['folder'], VERBOSE=VERBOSE)

    # Calculate consensus
    count = get_count(counts)
    consensus = get_consensus(count)

    if VERBOSE >= 1:
        # Get reference
        refseq = SeqIO.read(get_phix_filename(), 'fasta')
        ref = np.array(refseq)

        # Check differences
        ind = {key: (con != ref).nonzero()[0] for key, con in consensus.iteritems()}
        for key, indo in ind.iteritems():
            if len(indo):
                print 'PhiX run', miseq_run, 'differences from reference ('+key+'):'
                print 'Pos  Ref   Cons  Mutation'
                print '-------------------------'
                for i in indo:
                    i0 = max(i-2, 0)
                    i1 = min(i+3, len(refseq))
                    print '{:4d}'.format(i), \
                            rcf[key](''.join(ref[i0: i1])), \
                            rcf[key](''.join(consensus[key][i0: i1])),\
                            rcf[key](ref[i]), '->', rcf[key](consensus[key][i])
                print ''

    return {key: ''.join(value) for key, value in consensus.iteritems()}


def minor_alleles_along_genome(miseq_run, counts=None, VERBOSE=0, plot=False):
    '''Show the minor alleles along the phiX genome'''
    if counts is None:
        counts = load_allele_counts(MiSeq_runs[miseq_run]['folder'], VERBOSE=VERBOSE)

    # Study the minor alleles (sequencing errors)
    count = get_count(counts)
    minor_counts = get_minor_counts(count)
    minor_nus = get_minor_nus(count)

    # Mean error frequency
    cov = get_coverage(count)
    print 'Average error frequency:',
    for key in count:
        print key, '{:1.1e}'.format(minor_counts[key][300:].mean() / cov[key][300:].mean()),
    print ''

    # Plot
    if plot:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(13, 7))
        for key in count:
            axs[0].plot(minor_counts[key], lw=2, c=colors[key], alpha=0.5)
            axs[0].scatter(np.arange(len(minor_counts[key])), minor_counts[key],
                           s=20, c=colors[key], edgecolor='none', label=key)
        axs[0].set_xlabel('Position in phiX')
        axs[0].set_title('Sequencing errors')
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
        axs[1].set_xlabel('Position in phiX')
        axs[1].set_title('Sequencing errors frequencies')
        axs[1].set_xlim(-100, len(minor_counts['fwd']) + 100)
        axs[1].set_yscale('log')
    
    
        plt.suptitle('PhiX analysis: run '+str(miseq_run), fontsize=20)
        plt.tight_layout(rect=(0, 0, 1, 0.95))
    
        plt.ion()
        plt.show()

    return minor_counts, minor_nus


def spikes_motifs(miseq_run, counts=None, consensus=None, minor_counts=None,
                  VERBOSE=0, plot=False):
    '''Find the motifs around the spikes in error rates'''
    if counts is None:
        counts = load_allele_counts(MiSeq_runs[miseq_run]['folder'], VERBOSE=VERBOSE)
    count = get_count(counts)

    if minor_counts is None: 
        minor_counts = get_minor_counts(count)

    # Get consensus
    if consensus is None:
        consensus = {key: ''.join(c) for key, c in get_consensus(count).iteritems()}

    # Pick the spikes, i.e. the positions with error rates in the tail of the
    # binomial, p < 0.001 or so
    # Note: skip the first 300 positions, they are strange somehow (maybe mapping)
    from scipy.stats import binom
    spikes = defaultdict(list)
    for key in count:
        n = count[key][:, 300:].sum(axis=0).mean()
        p = 1.0 * minor_counts[key][300:].mean() / n
        k_crit = binom.ppf(0.999, n, p)
        ind = (minor_counts[key][300:] > k_crit).nonzero()[0] + 300
        if len(ind):
            print 'Spikes: run', miseq_run, '('+key+')'
            print ' Pos WT Motifs'
            print '----------------------'
            for pos in ind:
                print '{:4d}'.format(pos),
                print '{:>2s}'.format(consensus[key][pos]),
                for le in xrange(1, 5):
                    if pos > le - 2:
                        mutm = list(consensus[key][pos - le + 1: pos + 1])
                        mutm[le - 1] = alpha[count[key][:, pos] == minor_counts[key][pos]][0]
                        print ''.join(mutm),
                print
                spikes[key].append((pos,
                                    (consensus[key][pos - le + 1: pos + 1], ''.join(mutm)),
                                    1.0 * minor_counts[key][pos] / count[key][:, pos].sum(axis=0)))

    # Order spikes by frequency
    from operator import itemgetter
    for key in spikes:
        spikes[key].sort(key=itemgetter(2), reverse=True)

    # Plot the error histogram and the binomial if requested
    if plot:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(11.1, 6.8))
        xf = lambda x: np.logspace(-5, np.log10(0.5), x)

        for i, key in enumerate(count):
            ax = axs[i]
            ax.hist(1.0 * minor_counts[key][300:] / count[key][:, 300:].sum(axis=0),
                    bins=xf(60), normed=True, label='Errors')
            yb = binom.pmf(np.floor(n * xf(200)), n, p)
            ax.plot(xf(200), yb / yb.max() * 0.9 * ax.get_ylim()[1], lw=2, c='r',
                    label='Binomial')
            ax.set_xlabel('Error frequency')
            ax.set_xscale('log')
            ax.set_xlim(1e-5, 0.5)
            ax.set_ylabel('Density (A.U.)')
            ax.legend(loc=1, fontsize=12)
            ax.set_title(key, fontsize=18)
    
            ym = 0.3 * ax.get_ylim()[1]
            for j, (pos, trans, nu_spike) in enumerate(spikes[key]):
                ar = ax.arrow(nu_spike,
                              ym * 1.08**j, 0, -(ym * 1.08**j - 500),
                              edgecolor='k',
                              facecolor='k',
                              width=nu_spike / 100,
                              head_length=400,
                              overhang=0.4)
                txt = ax.text(nu_spike * 1.05, ym * 1.05 * 1.08**j, str(pos), fontsize=12)

        fig.suptitle('run '+str(miseq_run), fontsize=18)
        plt.tight_layout()
        plt.ion()
        plt.show()

        #fig.savefig('/ebio/ag-neher/home/fzanini/phd/sequencing/figures/phix_run'+str(miseq_run)+'_spikes.png')

    return spikes


def characterize_motifs(miseq_run, spikes=None, counts=None, consensus=None, minor_counts=None,
                        VERBOSE=0, plot=False):
    '''Find the DNA motifs of the errors'''
    if counts is None:
        counts = load_allele_counts(MiSeq_runs[miseq_run]['folder'], VERBOSE=VERBOSE)
    count = get_count(counts)
    cov = get_coverage(count)

    if consensus is None:
        consensus = get_consensus(count)

    # Calculate the whole transition table for single mutations
    M1 = {key: np.zeros((4, 5)) for key in count}
    for ia1, a1 in enumerate(alpha[:4]):
        # FWD
        cons_a1 = (consensus['fwd'] == a1).nonzero()[0]
        # Exclude positions before 300
        cons_a1 = cons_a1[cons_a1 >= 300]

        # Exclude spikes if provided
        if spikes is not None:
            pos_spikes = map(itemgetter(0), spikes['fwd'])
            cons_a1 = np.array(list((set(cons_a1) - set(pos_spikes))), int)
    
        for ia2, a2 in enumerate(alpha[:5]):
            M1['fwd'][ia1, ia2] = 1.0 * count['fwd'][ia2][cons_a1].sum() / \
                    cov['fwd'][cons_a1].sum()

        # REV
        cons_a1 = (consensus['rev'] == a1).nonzero()[0]
        # Exclude positions before 300
        cons_a1 = cons_a1[cons_a1 >= 300]

        # Exclude spikes if provided
        if spikes is not None:
            pos_spikes = map(itemgetter(0), spikes['fwd'])
            cons_a1 = np.array(list((set(cons_a1) - set(pos_spikes))), int)

        ira1 = alphal.index(rc(a1))
        for ia2, a2 in enumerate(alpha[:5]):
            ira2 = alphal.index(rc(a2))
            M1['rev'][ira1, ira2] = 1.0 * count['rev'][ia2][cons_a1].sum() / \
                    cov['rev'][cons_a1].sum()

    # Plot
    if plot:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(14, 5))
        for i, key in enumerate(count):
            im = axs[i].imshow(np.log10(M1[key]), interpolation='nearest',
                               vmin=-6, vmax=-3)

            axs[i].set_xticks(np.arange(5))
            axs[i].set_yticks(np.arange(4))
            axs[i].set_xticklabels(alphal[:5])
            axs[i].set_yticklabels(alphal[:4])
            axs[i].set_title(key)


        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.2, 0.05, 0.6])
        fig.colorbar(im, cax=cbar_ax)
        fig.suptitle('Mutation matrix for MiSeq (phiX, run'+str(miseq_run)+')',
                     fontsize=18)
        plt.figtext(0.05, 0.5, 'from', rotation=90, fontsize=18)
        plt.figtext(0.45, 0.04, 'to', fontsize=18)
        plt.figtext(0.88, 0.85, '10^-x', fontsize=18)

        plt.ion()
        plt.show()

    # Calculate transition table for higher order motifs (assuming independence)
    M2 = {key: np.zeros((4, 4, 5)) for key in count}
    for ia11, a11 in enumerate(alpha[:4]):
        for ia12, a12 in enumerate(alpha[:4]):
            # FWD
            cons_a1 = ((consensus['fwd'][:-1] == a11) & \
                       (consensus['fwd'][1:] == a12)).nonzero()[0]
            # Exclude positions before 300
            cons_a1 = cons_a1[cons_a1 >= 300]
            # Exclude spikes if provided
            if spikes is not None:
                pos_spikes = map(itemgetter(0), spikes['fwd'])
                cons_a1 = np.array(list((set(cons_a1) - set(pos_spikes))), int)

            # Focus on the last allele
            cons_a1 += 1
            if VERBOSE >= 3:
                print a11+a12, len(cons_a1)

            for ia2, a2 in enumerate(alpha[:5]):
                M2['fwd'][ia11, ia12, ia2] = 1.0 * count['fwd'][ia2][cons_a1].sum() / \
                        cov['fwd'][cons_a1].sum()

            # REV
            cons_a1 = ((consensus['rev'][:-1] == a11) & \
                       (consensus['rev'][1:] == a12)).nonzero()[0]
            # Exclude positions before 300
            cons_a1 = cons_a1[cons_a1 >= 300]
            # Exclude spikes if provided
            if spikes is not None:
                pos_spikes = map(itemgetter(0), spikes['fwd'])
                cons_a1 = np.array(list((set(cons_a1) - set(pos_spikes))), int)

            # We already are at the last allele (see fwd)
            if VERBOSE >= 3:
                print a11+a12, len(cons_a1)

            ira11 = alphal.index(rc(a11))
            ira12 = alphal.index(rc(a12))
            for ia2, a2 in enumerate(alpha[:5]):
                ira2 = alphal.index(rc(a2))
                M2['rev'][ira12, ira11, ira2] = 1.0 * count['rev'][ia2][cons_a1].sum() / \
                        cov['rev'][cons_a1].sum()


    # Flatten first two dimensions
    Ma = np.array([[a1+a2 for a2 in alpha[:4]] for a1 in alpha[:4]]).T.ravel()
    M2f = {key: M.swapaxes(0,1).reshape((M.shape[0] * M.shape[1], M.shape[2]))
           for key, M in M2.iteritems()}

    # Plot
    if plot:
        import matplotlib.pyplot as plt
        fig2, axs = plt.subplots(1, 2, figsize=(10, 12))
        for i, key in enumerate(count):
            im = axs[i].imshow(np.log10(M2f[key]), interpolation='nearest',
                               vmin=-6, vmax=-3)

            axs[i].set_xticks(np.arange(5))
            axs[i].set_yticks(np.arange(4**2))
            axs[i].set_xticklabels(alphal[:5])
            axs[i].set_yticklabels(Ma)
            axs[i].set_title(key)


        fig2.subplots_adjust(right=0.8)
        cbar_ax = fig2.add_axes([0.85, 0.1, 0.05, 0.8])
        fig2.colorbar(im, cax=cbar_ax)
        fig2.suptitle('Double mutation matrix for MiSeq (phiX, run'+str(miseq_run)+')',
                     fontsize=18)
        plt.figtext(0.05, 0.5, 'from', rotation=90, fontsize=18)
        plt.figtext(0.45, 0.04, 'to', fontsize=18)
        plt.figtext(0.88, 0.92, '10^-x', fontsize=18)

        plt.ion()
        plt.show()

        fig.savefig('/ebio/ag-neher/home/fzanini/phd/sequencing/figures/phix_run'+str(miseq_run)+'_mutmatrix.png')
        fig2.savefig('/ebio/ag-neher/home/fzanini/phd/sequencing/figures/phix_run'+str(miseq_run)+'_mutmatrix2.png')


    return M1, M2


def errors_vs_quality(miseq_run, maxreads=1e3, VERBOSE=0, plot=False):
    '''Count the errors VS the phred score'''
    # In principle, the phred is q = 10 * log10(p) where p is the error rate,
    # so 30 -> 1e-3, 40 -> 1e-4 and so on.

    # Prepare output structure
    qual_errs = np.zeros(min(1e6, maxreads * 250), 'S1')
    qual_base = np.zeros(1e6, 'S1')
    
    # Calculate consensus
    counts = load_allele_counts(MiSeq_runs[miseq_run]['folder'], VERBOSE=VERBOSE)
    count = get_count(counts)
    consensus = get_consensus(count)

    # Open BAM file
    # Note: the reads should already be filtered of unmapped stuff at this point
    bamfilename = get_mapped_phix_filename(data_folder, filtered=False)
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        # Iterate over single reads (no linkage info needed)
        ierr = 0
        ibase = 0
        for i, read in enumerate(bamfile):

            if ierr >= len(qual_errs) - 1:
                break

            if i == maxreads:
                if VERBOSE >= 2:
                    print 'Max reads reached:', maxreads
                break
        
            # Print output
            if (VERBOSE >= 3) and (not ((i +1) % 10000)):
                print (i+1)
        
            # Exclude unmapped and unpaired
            if read.is_unmapped or (not read.is_proper_pair):
                continue

            # Divide by forward/reverse
            if read.is_reverse:
                key = 'rev'
            else:
                key = 'fwd'

            # Get appropriate consensus (they differ by a few bases)
            cons = consensus[key]
    
            # Read good CIGARs (to avoid mapping quirks)
            cigar = read.cigar
            len_cig = len(cigar)
            (good_cigars, first_good_cigar, last_good_cigar) = \
                    get_ind_good_cigars(cigar, match_len_min=match_len_min,
                                        full_output=True)
            
            # Sequence and position
            # Note: stampy takes the reverse complement already
            pos = read.pos
            seq = read.seq
            qual = read.qual

            # Iterate over CIGARs
            for ic, (block_type, block_len) in enumerate(cigar):

                # Inline block
                if block_type == 0:
                    # Exclude bad CIGARs
                    if good_cigars[ic]: 
            
                        # The first and last good CIGARs are matches:
                        # trim them (unless they end the read)
                        if (ic == first_good_cigar) and (ic != 0):
                            trim_left = trim_bad_cigars
                        else:
                            trim_left = 0
                        if (ic == last_good_cigar) and (ic != len_cig - 1):
                            trim_right = trim_bad_cigars
                        else:
                            trim_right = 0
            
                        seqb = np.array(list(seq[trim_left:block_len - trim_right]), 'S1')
                        qualb = np.array(list(qual[trim_left:block_len - trim_right]), 'S1')
                        # Check for errors
                        cloc = cons[pos + trim_left: pos + trim_left + len(seqb)]
                        errs = (cloc != seqb)
                        
                        # Reads with lots of mutations are mismapped
                        if errs.sum() > 20:
                            break
                        
                        elif errs.any():
                            qerrs = qualb[errs]
                            if len(qual_errs) - (ierr + len(qerrs)) <= 1:
                                qerrs = qerrs[:len(qual_errs) - ierr - 1]
                            
                            qual_errs[ierr: ierr + len(qerrs)] = qerrs
                            ierr += len(qerrs)

                            if VERBOSE >= 3:
                                print len(qerrs), qerrs, ierr

                        # Collect background
                        if len(qual_base) - ibase > len(qualb):
                            qual_base[ibase: ibase + len(qualb)] = qualb
                            ibase += len(qualb)

            
                    # Chop off this block
                    if ic != len_cig - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
                        pos += block_len
            
                # Deletion: the whole block is wrong, but has no quality!
                elif block_type == 2:
                    # Chop off pos, but not sequence
                    pos += block_len
            
                # Insertion
                # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
                # THEN the insert, FINALLY comes seq[391:]
                # The whole block is wrong, no matter what, but that might come
                # from mapping
                elif block_type == 1:
                    # Chop off seq, but not pos
                    if ic != len_cig - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
            
                # Other types of cigar?
                else:
                    raise ValueError('CIGAR type '+str(block_type)+' not recognized')

    # Convert and trim
    qual_errs = qual_errs[:ierr]
    qual_errs.sort()
    qual_errs = np.array(map(ord, qual_errs), int) - ord('!')
    qual_base = qual_base[:ibase]
    qual_base.sort()
    qual_base = np.array(map(ord, qual_base), int) - ord('!')

    if plot:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))

        # Cumulative
        ax = axs[0]
        x = qual_errs
        y = np.linspace(0, 1, len(qual_errs))
        ax.plot(x, y, lw=2, c='b', label='errs')
    
        xb = qual_base
        yb = np.linspace(0, 1, len(qual_base))
        ax.plot(xb, yb, lw=2, c='grey', label='baseline')

        ax.set_xlabel('Phred quality')
        ax.set_ylabel('Cumulative fraction of errors')
        ax.legend(loc=2)

        # Error probability
        # This goes via Bayes' theorem: P(e|q) = P(q|e) * P(e) / P(q)
        ax = axs[1]
        Pq = 1.0 * np.bincount(qual_base - 20) / len(qual_base)
        Pqe = 1.0 * np.bincount(qual_errs - 20) / len(qual_errs)
        Pe = 1.1e-4
        y = Pqe * Pe / (Pq + 1e-6)
        x = np.arange(20, 20 + len(y))
        ax.scatter(x[y > 0], y[y > 0], lw=2, c='k', label='Data')

        # Plot theoretical Phred curve
        ax.plot(x, 10**(-1.0 * x / 10), lw=2, c='grey', label='Theory')

        ax.set_xlabel('Phred quality')
        ax.set_ylabel('P(e | q)')
        ax.set_yscale('log')
        ax.set_ylim(1e-5, 1e-1)
        ax.legend(loc=1)

        fig.suptitle('run '+str(miseq_run), fontsize=18)
        plt.tight_layout(rect=(0, 0, 1, 0.96))
        plt.ion()
        plt.show()

    return (qual_errs, qual_base)


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
    
#        # Get allele counts
#        counts = load_allele_counts(data_folder, VERBOSE=VERBOSE)
#        data_runs[miseq_run]['counts'] = counts
#    
#        # Check consensus
#        consensus = consensus_vs_reference(miseq_run, counts, VERBOSE=VERBOSE)
#        data_runs[miseq_run]['consensus'] = consensus
#
#        # Along genome
#        minor_counts, \
#        minor_nus = minor_alleles_along_genome(miseq_run, counts,
#                                               VERBOSE=VERBOSE,
#                                               plot=plot)
#        data_runs[miseq_run]['minor_counts'] = minor_counts
#        data_runs[miseq_run]['minor_nus'] = minor_nus

#        # Spikes
#        spikes = spikes_motifs(miseq_run, VERBOSE=VERBOSE, plot=False)
#        data_runs[miseq_run]['spikes'] = spikes
#
#        # Motifs
#        motifs = characterize_motifs(miseq_run,
#                                     spikes=spikes,
#                                     VERBOSE=VERBOSE,
#                                     plot=plot)
#        data_runs[miseq_run]['motifs'] = motifs

        # Phred quality VS errors
        (qual_errs, qual_base) = errors_vs_quality(miseq_run, maxreads=1e6, VERBOSE=VERBOSE,
                                      plot=plot)
        data_runs[miseq_run]['qual_errs'] = (qual_errs, qual_base)

#    # Check consistency across runs
#    if len(data_runs) > 1:
#        minor_nus_crossrun(data_runs)
#        spikes_crossrun(data_runs)

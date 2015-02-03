# vim: fdm=marker
'''
author:     Fabio Zanini, Richard Neher
date:       19/01/15
content:    Estimate the number of template molecules to PCR, fragment by fragment.
'''
# Modules
import os
import sys
import argparse
from itertools import izip
from collections import defaultdict
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import hivwholeseq.utils.plot

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.samples import load_samples_sequenced, SamplePat
from hivwholeseq.patients.filenames import get_ntemplates_by_fragment_filename


# Globals
fragments = ['F'+str(i) for i in xrange(1, 7)]



# Functions
def align_fragments(c1, c2, VERBOSE=0):
    '''Align subsequence fragments'''
    import numpy as np
    from seqanpy import align_ladder
    from hivwholeseq.utils.sequence import pretty_print_pairwise_ali
    (score, a1, a2) = align_ladder(c1, c2, score_gapopen=-20)
    start2 = len(a2) - len(a2.lstrip('-'))
    end1 = len(a1.rstrip('-'))

    a1 = a1[start2: end1]
    a2 = a2[start2: end1]

    if VERBOSE >= 3:
        pretty_print_pairwise_ali((a1, a2), width=100,
                                  name1=fr1, name2=fr2)

    a1 = np.fromstring(a1, 'S1')
    a2 = np.fromstring(a2, 'S1')
    co1 = (a1 != '-').cumsum() - 1
    co2 = (a2 != '-').cumsum() - 1
    ind = (a1 != '-') & (a2 != '-')

    pos1 = co1[ind] + start2
    pos2 = co2[ind]

    return (pos1, pos2)


def get_chain_indices_nonnan(arr):
    '''Group valid numbers in an array into chains'''
    import numpy as np

    groups = []
    group = []
    for i, el in enumerate(arr):
        if np.isnan(el):
            # Empty groups need not get added
            if group:
                groups.append(group)
                group = []
        else:
            group.append(i)
    
    # Final group
    if group:
        groups.append(group)

    return map(set, groups)


def plot_template_estimate(data, samplename, figaxs=None):
    '''Plots for the estimate of template numbers'''

    if figaxs is not None:
        fig, axs = figaxs
    else:
        fig, axs = plt.subplots(1, 2, figsize=(13, 8))
    
    ax = axs[0]
    xmin = 1e-3
    xpl = np.linspace(xmin, 1 - xmin, 1000)

    # Plot diagonal
    ax.plot(xpl, xpl, lw=1.5, color='k')

    # Plot max discrepancy curve, x (1 - x) / var, given by two samples, when
    # x = (x1 + x2) / 2
    # var = ((x1 - x2) / 2)**2
    ypl = xpl * (1 - xpl) / ((xpl >= 0.5) - xpl)**2
    axs[1].plot(xpl, ypl, lw=1.5, color='k')

    for (samplename_p, fr1, fr2), (af1, af2) in data['af'].iteritems():
        if samplename_p != samplename:
            continue

        color = cm.jet(1.0 * int(fr1[1]) / 5)

        ax = axs[0]
        ax.scatter(af1, af2,
                   color=color,
                   label='-'.join((fr1, fr2)))

        # Plot the noise variance
        n = data['n'][(samplename, fr1, fr2)]
        std = np.sqrt(xpl * (1 - xpl) / n)

        y1 = xpl + std
        y2 = xpl - std
        ax.plot(xpl, y1, color=color, lw=1)
        ax.plot(xpl, y2, color=color, lw=1)

        # Plot the variances
        ax = axs[1]
        ax.scatter(data['mean'][(samplename, fr1, fr2)],
                   data['n_all'][(samplename, fr1, fr2)],
                   color=color)
        n = data['n'][(samplename, fr1, fr2)]
        if not np.isnan(n):
            ax.plot(xpl, [n] * len(xpl), lw=1.5, ls='-',
                    alpha=0.5,
                            color=color)

            ax = axs[0]
            ax.grid(True)
            ax.set_xscale('logit')
            ax.set_yscale('logit')
            ax.set_xlim(xmin, 1 - xmin)
            ax.set_ylim(xmin, 1 - xmin)
            ax.set_xlabel('leading fragment')
            ax.set_ylabel('trailing fragment')
            ax.set_title(samplename)

            ax = axs[1]
            ax.grid(True)
            ax.set_xlabel('Average frequency')
            ax.set_ylabel('x (1 - x) / var(x)')
            ax.set_xscale('logit')
            ax.set_yscale('log')
            ax.set_xlim(xmin, 1 - xmin)
            ax.set_ylim(ymin=1)

            plt.tight_layout()


def assign_templates(ns, samplenames):
    '''Assign template numbers to fragments from the overlaps'''
    from collections import defaultdict
    import numpy as np

    overlaps = zip(fragments[:-1], fragments[1:])

    ns = defaultdict(lambda: np.nan, ns)

    M = np.ma.masked_all((len(samplenames), len(fragments)), int)
    for isa, samplename in enumerate(samplenames):
        n_ov = [ns[(samplename, fr1, fr2)] for (fr1, fr2) in overlaps]
        n_ov = np.array(n_ov)

        # FIXME: for now take a hierarchical approach, but we
        # could exploit granularity of the SFS (10 templates don't produce
        # frequencies of 0.01 reliably)

        # FIXME: we should provide error bars, because we often rely on few sites

        # Group into chains (because we have some nans i.e. missing data)
        ind_grs = get_chain_indices_nonnan(n_ov)
        for ind_gr in ind_grs:
            
            # Singletons are not very useful: both fragments have only that
            # single overlap (either because on the other side there is a 
            # nan, or nothing ( <-F1 and F6-> ).
            if len(ind_gr) == 1:
                ind_gr = ind_gr.pop()
                M[isa, ind_gr: ind_gr + 2] = n_ov[ind_gr]
                continue

            # Chains are better. If a fragment has 1+ overlaps with high n,
            # it should have high n itself, so we resolve degeneracy this way:
            # the highest n propagates left and right, and the other fragments
            # are assigned what's left
            # NOTE: this is kind of a strict rule, maybe we'd like to average a bit
            i_frags_left = set(ind_gr) | set([max(ind_gr) + 1])
            i_max = max(ind_gr, key=n_ov.__getitem__)

            # First two fragments
            M[isa, i_max: i_max + 2] = n_ov[i_max]
            i_frags_left.remove(i_max)
            i_frags_left.remove(i_max + 1)

            # The other fragments
            while (i_frags_left):
                i_frag = i_frags_left.pop()
                if i_frag < i_max:
                    M[isa, i_frag] = n_ov[i_frag]
                else:
                    M[isa, i_frag] = n_ov[i_frag - 1]

    return M



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Estimate number of templates',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save alignment to file')
    parser.add_argument('--plot', action='store_true',
                        help='Plot frequencies in overlaps')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    VERBOSE = args.verbose
    use_save = args.save
    use_plot = args.plot

    samples = load_samples_sequenced()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()

    data = defaultdict(dict)
    for samplename, sample in samples.iterrows():
        sample = SamplePat(sample)
        if VERBOSE >= 1:
            print samplename

        for (fr1, fr2) in izip(fragments[:-1], fragments[1:]):
            try:
                ac1 = sample.get_allele_counts(fr1)
                ac2 = sample.get_allele_counts(fr2)
            except IOError:
                continue

            if VERBOSE >= 2:
                print fr1, fr2

            # Filter positions by coverage
            covmin = 100
            indcm1 = (ac1.sum(axis=0) >= covmin)
            indcm2 = (ac2.sum(axis=0) >= covmin)
            ac1 = ac1[:, indcm1]
            ac2 = ac2[:, indcm2]

            c1 = alpha[ac1.argmax(axis=0)]
            c2 = alpha[ac2.argmax(axis=0)]

            # Get only overlap
            (pos1, pos2) = align_fragments(c1, c2, VERBOSE=VERBOSE)
            ao1 = ac1[:, pos1]
            ao2 = ac2[:, pos2]

            # Get frequencies
            af1 = 1.0 * ao1 / ao1.sum(axis=0)
            af2 = 1.0 * ao2 / ao2.sum(axis=0)

            # Filter only polymorphic sites
            afmin = 3e-3
            indfm = (af1 >= afmin) & (af1 <= 1 - afmin) & (af2 >= afmin) & (af2 <= 1 - afmin)
            nsites = len(np.unique(indfm.nonzero()[1]))
            if VERBOSE >= 2:
                print 'Number of doubly polymorphic sites:', nsites

            # Estimate the template number
            mea = 0.5 * (af1[indfm] + af2[indfm])
            var = ((af1[indfm] - af2[indfm]) / 2)**2

            # In binomial sampling, the variance on k is var(k) = nx (1 - x), so
            # for the frequency var(k/n) = x (1 - x) / n
            n_all = mea * (1 - mea) / var
            
            # NOTE: pseudocounts that come from the F4 dilution estimate, so we
            # only listen to the data if there is enough data points to listen to
            len_pseudo = 1
            n_pseudo = sample.get_n_templates_dilutions()
            n_allp = np.concatenate([n_all, ([n_pseudo] * len_pseudo)])

            # NOTE: the estimate of n has a bad distribution because some points are
            # exactly on the diagonal, so we average the inverse (which is well
            # behaved) and also take the medians as alternatives
            n = 1.0 / (1.0 / n_allp).mean()
            ninv = n_allp.mean()
            nmed = np.median(n_allp)
            if VERBOSE >= 2:
                print fr1, fr2, n, ninv, nmed

            key = (samplename, fr1, fr2)
            data['af'][key] = (af1[indfm], af2[indfm])
            data['mean'][key] = mea
            data['var'][key] = var
            data['n_all'][key] = n_all
            data['n'][key] = n
            data['ninv'][key] = ninv
            data['nmed'][key] = nmed


    if use_plot:

        for samplename in samples.index:
            plot_template_estimate(data, samplename, figaxs=None)

    if use_plot:
        plt.ion()
        plt.show()


    if VERBOSE >= 2:
        print 'Assign template numbers overlaps -> fragments'
    M = assign_templates(data['n'], samples.index)

    if use_save:
        if VERBOSE >= 1:
            print 'Save results'

        fn_out = get_ntemplates_by_fragment_filename()
        with open(fn_out, 'w') as f:
            sep = '\t'
            f.write(sep.join(['Sample name', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6'])+'\n')
            for isa, (samplename, sample) in enumerate(samples.iterrows()): 
                row = [samplename]
                for n in M[isa]:
                    if np.ma.is_masked(n):
                        row.append('n.a.')
                    else:
                        row.append(str(n))
                f.write(sep.join(row)+'\n')












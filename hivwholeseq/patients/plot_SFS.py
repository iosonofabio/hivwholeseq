# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/05/14
content:    Plot site frequency spectra for derived alleles.
'''
# Modules
import os
import argparse
from operator import itemgetter
import numpy as np
from Bio import SeqIO

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_allele_counts_filename


# Functions
def get_allele_count(pname, samplename_pat, fragment, VERBOSE=0, use_PCR1=1):
    '''Get counts from a single patient sample'''

    fns = []
    samplenames_out = []
    fn1 = get_allele_counts_filename(pname, samplename_pat, fragment, PCR=1)
    fn2 = get_allele_counts_filename(pname, samplename_pat, fragment, PCR=2)
    if use_PCR1 == 0:
        for PCR, fn in enumerate((fn1, fn2), 1):
            if os.path.isfile(fn):
                fns.append(fn)
                samplenames_out.append((samplename_pat, PCR))
    elif use_PCR1 == 1:
        if os.path.isfile(fn1):
            fns.append(fn1)
            samplenames_out.append((samplename_pat, 1))
        elif os.path.isfile(fn2):
            fns.append(fn2)
            samplenames_out.append((samplename_pat, 2))
    elif use_PCR1 == 2:
        if os.path.isfile(fn1):
            fns.append(fn1)
            samplenames_out.append((samplename_pat, 1))

    acs = []
    for i, fn in enumerate(fns):
        ac = np.load(fn)
        if i == 0:
            ls = ac.shape[2]
            acs = np.zeros((len(fns), len(alpha), ls), int)

        # Average directly over read types?
        acs[i] = np.load(fn).sum(axis=0)

    return (samplenames_out, acs)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele frequency trajectories',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele frequency trajectories to file')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot
    submit = args.submit
    use_PCR1 = args.PCR1
    use_logit = args.logit

    # Bins for histogram
    itouch = 0
    bins = np.logspace(-2.8, -0.0001, 30)
    binsc = np.sqrt(bins[1:] * bins[:-1])
    binsw = bins[1:] - bins[:-1]
    hist = None

    samples = lssp()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    # Iterate over samples and fragments
    for fragment in fragments:
        if VERBOSE >= 1:
            print fragment

        for samplename_pat, sample in samples.iterrows():
            if VERBOSE >= 1:
                print samplename_pat,

            pname = sample.patient
            ref_fn = get_initial_reference_filename(pname, fragment)
            ref = SeqIO.read(ref_fn, 'fasta')
            refm = np.array(ref, 'S1')

            # FIXME: determine better PCR errors?
            ntemplate = sample['n templates']
            if ntemplate > 0.1:
                depthmax = np.maximum(1e-3, 1.0 / ntemplate)
            else:
                depthmax = 1e-3

            (sns, acs) = get_allele_count(pname, samplename_pat, fragment, VERBOSE=VERBOSE, use_PCR1=use_PCR1)

            # We might be getting back one count table, or two in case of both PCR types
            for ac in acs:
                cov = ac.sum(axis=0)

                if cov.mean() < 1000:
                    continue

                af = 1.0 * ac / cov
                af[np.isnan(af)] = np.ma.masked
                af[af < depthmax] = 0
                af[af > 1.0 - depthmax] = 1

                # Get rid of gaps and low-coverage regions
                is_gap = (af.argmax(axis=0) == 6) | (ac.sum(axis=0) < 100)
                if VERBOSE >= 2:
                    print 'Fraction of gap sites (excluded):', is_gap.mean()

                # Get rid of ancestral alleles
                af_der = af.copy()
                af_der[:, is_gap] = 0
                for i, nuc in enumerate(refm):
                    ai = alphal.index(nuc)
                    af_der[ai, i] = np.ma.masked

                histtmp = np.histogram(af_der, bins=bins, density=False)[0]
                if hist is None:
                    hist = histtmp
                else:
                    hist += histtmp

                if VERBOSE >= 1:
                    print 'ok',

            print ''

    # Get density from histogram counts
    histd = 1.0 * hist / hist.sum() / binsw

    if plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        x = binsc
        if use_logit:
            trfun = lambda xx: np.log10(xx / (1 - xx))
            trfuni = lambda y: 10**y / (1 + 10**y)
            x = trfun(x)
        ax.plot(x, histd, lw=2, c='k')
        alpha = histd[itouch]
        if use_logit:
            xext = trfun(np.array([binsc[0], binsc[-1]]))
            xmod = np.linspace(xext[0], xext[1], 100)
            xmodi = trfuni(xmod)
            ymod1 = alpha * ((binsc[itouch] / xmodi))
            ymod2 = alpha * ((binsc[itouch] / xmodi)**2)
        else:
            xmod = np.array([binsc[0], binsc[-1]])
            ymod1 = alpha * ((binsc[itouch] / xmod))
            ymod2 = alpha * ((binsc[itouch] / xmod)**2)
        ax.plot(xmod, ymod1, lw=2, c='g')
        ax.plot(xmod, ymod2, lw=2, c='r')

        ax.set_xlabel('Freq')
        ax.set_ylabel('SFS [density = counts / sum / binsize]')
        if use_logit:
            tickloc = np.array([0.0001, 0.01, 0.5, 0.99, 0.9999])
            ax.set_xticks(trfun(tickloc))
            ax.set_xticklabels(map(str, tickloc))
            from matplotlib.ticker import FixedLocator
            ticklocminor = np.concatenate([[10**po * x for x in xrange(2 , 10)] for po in xrange(-4, -1)] + \
                                          [[0.1 * x for x in xrange(2 , 9)]] + \
                                          [[1 - 10**po * (10 - x) for x in xrange(2, 10)] for po in xrange(-2, -5, -1)])
            ax.xaxis.set_minor_locator(FixedLocator(trfun(ticklocminor)))
            ax.set_xlim(-2.8, 1.1)

        else:
            ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title('SFS, '+str(fragments))
        ax.grid(True)

        plt.tight_layout()
        plt.ion()
        plt.show()



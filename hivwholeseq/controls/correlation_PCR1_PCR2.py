# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/03/14
content:    Test the correlation coefficient between PCR1 and PCR2 patient samples,
            e.g. in allele frequencies and coallele frequencies.
'''
# Modules
import sys
import os
import argparse
import numpy as np
from itertools import izip
import pysam
from Bio import SeqIO
from scipy.stats import pearsonr, spearmanr
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.sequencing.filter_mapped_reads import plot_distance_histogram, \
        plot_distance_histogram_sliding_window, get_distance_from_consensus, \
        check_overhanging_reads, trim_bad_cigar
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_filter_mapped_init_summary_filename, \
        get_allele_cocounts_filename, get_allele_frequency_trajectories_filename
from hivwholeseq.utils.mapping import convert_sam_to_bam, pair_generator
from hivwholeseq.utils.two_site_statistics import get_coallele_counts_from_file


# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Correlation PCR1/PCR2')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--samples', nargs='*',
                        help='Samples to analyze')
    parser.add_argument('--threshold', type=float, default=-3,
                        help='Log10 of the minimal allele frequency to trust')
    parser.add_argument('--allele-frequencies', action='store_true',
                        dest='allele_frequencies',
                        help='Look at allele frequencies')
    parser.add_argument('--allele-cofrequencies', action='store_true',
                        dest='allele_cofrequencies',
                        help='Look at allele cofrequencies')
    parser.add_argument('--saveplot', action='store_true',
                        help='Save the plot to file')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    samplenames = args.samples
    threshold = args.threshold
    use_allele_frequencies = args.allele_frequencies
    use_allele_cofrequencies = args.allele_cofrequencies
    saveplot = args.saveplot

    patient = get_patient(pname)

    # If the script is called with no fragment, iterate over all
    # NOTE: This assumes that the fragment has been sequenced. If not, we should
    # take a later time point or similia, but let's not worry about it yet.
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    for fragment in fragments:
        for samplename in samplenames:
            samplenames_both = [samplename+'_PCR1', samplename+'_PCR2']
            index_both = map(patient.samples.index, samplenames_both)

            # Check allele frequencies
            aft = np.load(get_allele_frequency_trajectories_filename(pname, fragment))
            aft = np.ma.array(aft[index_both])
            aft[:, (aft < 10**threshold).any(axis=0)] = np.ma.masked
            aftlog = np.log10(aft)

            if use_allele_frequencies:
                print 'Allele frequencies'
                print 'Difference, mean:', np.diff(aftlog[:, (aftlog > threshold).all(axis=0)]).mean()
                print 'Difference, std :', np.diff(aftlog[:, (aftlog > threshold).all(axis=0)]).std()
                import matplotlib.pyplot as plt
                colors = cm.jet(1.0 * np.arange(aft.shape[-1]) / aft.shape[-1])
                shapes = ('o', 'v', 's', '*', '^', 'p')
                plt.figure()
                for ai in xrange(aftlog.shape[1]):
                    plt.scatter(aftlog[0, ai].data,
                                aftlog[1, ai].data,
                                facecolor=colors[aftlog[0, ai].mask],
                                marker=shapes[ai],
                                s=40,
                                edgecolor='none', alpha=0.5,
                                label=alpha[ai])
                x = np.linspace(threshold, 0)
                plt.plot(x, x + 0.5, lw=2, c='r')
                plt.plot(x, x - 0.5, lw=2, c='r')
                plt.xlabel('log10 allele freq, PCR1')
                plt.ylabel('log10 allele freq, PCR2')
                plt.xlim(threshold - 0.1, 0)
                plt.ylim(threshold - 0.1, 0)
                plt.title(', '.join([pname, fragment, samplename]))
                plt.legend(loc=4, fontsize=10)

                if saveplot:
                    from hivwholeseq.patients.filenames import get_correlation_PCR1_PCR2_aft_figure_filename
                    fn = get_correlation_PCR1_PCR2_aft_figure_filename(pname, fragment, samplename)
                    plt.savefig(fn)

            if not use_allele_cofrequencies:
                continue

            # Check coallele frequencies
            coaft = np.ma.empty((2, len(alpha), len(alpha), aft.shape[-1], aft.shape[-1]), float)
            if VERBOSE >= 1:
                print 'Collecting cocounts and normalizing'
            for i, samplename in enumerate(samplenames_both):
                if VERBOSE >= 2:
                    print pname, fragment, samplename
                cocounts = np.ma.array(np.load(get_allele_cocounts_filename(pname, samplename, fragment)), int)
                cocoverage = cocounts.sum(axis=0).sum(axis=0)
                cocounts[:, :, cocoverage < 100] = np.ma.masked
                coaft[i] = 1.0 * cocounts / cocoverage
                del cocounts, cocoverage
            if VERBOSE >= 1:
                print 'Broadcasting allele frequencies, calculating LD, taking logs'
            aftbroad = np.einsum('mij...,m...kl->mikjl', aft, np.ones_like(aft))
            aftbroadT = aftbroad.swapaxes(1, 2).swapaxes(3, 4)
            coaft[:, (aftbroad < 10**threshold) | (aftbroadT < 10**threshold)] = np.ma.masked
            coaft[:, (coaft < 1e-3).any(axis=0)] = np.ma.masked
            LD = coaft - aftbroad * aftbroadT
            coaft = np.log10(coaft)

            #print 'Allele cofrequencies'
            #print 'Difference, mean:', (coaft[1] - coaft[0]).mean()
            #print 'Difference, std :', (coaft[1] - coaft[0]).std()

            #print 'LD'
            #print 'Difference, mean:', (LD[1] - LD[0]).mean()
            #print 'Difference, std :', (LD[1] - LD[0]).std()


            # Pick a subsample only
            if VERBOSE >= 1:
                print 'Getting random subsample'
            ind = np.random.randint(coaft[0].ravel().shape[0], size=1000000); ind.sort()
            # Getting distances between sites in pairs
            distind = np.unravel_index(ind, coaft[0].shape)[2:]
            distind = np.abs(distind[0] - distind[1])
            colorsind = cm.jet(np.array(255.0 * distind / coaft.shape[-1], int))

            # Plot cofrequencies
            if VERBOSE >= 1:
                print 'Plotting'
            plt.figure()
            plt.scatter(coaft[0].ravel()[ind],
                        coaft[1].ravel()[ind],
                        s=40, facecolor=colorsind,
                        edgecolor='none',
                        alpha=0.5)
            x = np.linspace(2 * threshold, 0)
            plt.plot(x, x + 0.5, lw=2, c='r')
            plt.plot(x, x - 0.5, lw=2, c='r')
            plt.xlabel('log10 allele cofreq, PCR1')
            plt.ylabel('log10 allele cofreq, PCR2')
            plt.xlim(1.5 * threshold - 0.1, 0.05)
            plt.ylim(1.5 * threshold - 0.1, 0.05)
            (r, Pv) = pearsonr(coaft[0].ravel()[ind], coaft[1].ravel()[ind])
            plt.text(1.47 * threshold, -0.25, 'r = '+'{:2.2f}'.format(r))
            plt.title(', '.join([pname, fragment, samplename]))
            plt.tight_layout()

            ## Make a correlation plot between distance and logdifference in cofreqs
            #logdiff = np.abs(np.log10(coaft[0].ravel()[ind] / coaft[1].ravel()[ind]))
            #ind2 = -logdiff.mask
            #distuni = np.arange(500)
            #logdiffuni = np.array([np.median(logdiff[(distind == d) & ind2].data) for d in distuni])
            #logdiffunim = np.array([np.mean(logdiff[(distind == d) & ind2].data) for d in distuni])
            #plt.figure()
            #plt.scatter(distuni, logdiffuni, color='b', label='Median')
            #plt.scatter(distuni, logdiffunim, color='r', label='Mean')
            #plt.xlabel('Distance between sites in the pair[bp]')
            #plt.ylabel('Median/mean of |logdiff| cofreqs')
            #plt.title(', '.join([pname, fragment, samplename]))
            #plt.legend(loc=2)
            #plt.tight_layout()

            ## Plot cumulative hists
            #distssample = [10, 100, 200, 300, 400, 450]
            #logdiffsamples = [(logdiff[(distind == d) & ind2].data) for d in distssample]
            #map(lambda x: x.sort(), logdiffsamples)
            #plt.figure()
            #for i in xrange(len(distssample)):
            #    plt.plot(logdiffsamples[i], np.linspace(0, 1, len(logdiffsamples[i]))[::-1],
            #             lw=2,
            #             color=cm.jet(int(255.0 * i / len(distssample))),
            #             label=str(distssample[i]))
            #plt.xlabel('|logdiff| cofreqs')
            #plt.ylabel('Fraction of pairs above this x')
            #plt.legend(loc=1)
            #plt.title(', '.join([pname, fragment, samplename]))
            #plt.tight_layout()

            # Repeat the LD, median, and cumulative distr plots for LD
            plt.figure()
            plt.scatter(LD[0].ravel()[ind],
                        LD[1].ravel()[ind],
                        s=40, facecolor=colorsind,
                        edgecolor='none',
                        alpha=0.5)
            plt.xlabel('LD, PCR1')
            plt.ylabel('LD, PCR2')
            x = np.linspace(-1, 1)
            plt.plot(x, x, lw=2, c='r')
            plt.plot([0] * len(x), x, lw=2, c='k', alpha=0.3)
            plt.plot(x, [0] * len(x), lw=2, c='k', alpha=0.3)
            plt.xlim(-0.03, 0.03)
            plt.ylim(-0.03, 0.03)
            (r, Pv) = pearsonr(LD[0].ravel()[ind], LD[1].ravel()[ind])
            plt.text(-0.028, 0.025, 'r = '+'{:2.2f}'.format(r))
            plt.title(', '.join([pname, fragment, samplename]))
            plt.tight_layout()

            # Plot stratified by distance
            distssample = [10, 30, 60, 100, 200, 400]
            LDsamples = []
            for d in distssample:
                LDd = []
                for i in xrange(2):
                    LDd.append(LD[i].ravel()[:, (distind >= 0.79 * d) & (distind <= 1.21 * d)])
                LDd = np.ma.array(LDd)
                LDd = LDd[:, -LDd.mask.any(axis=0)].data
                LDsamples.append(LDd)
            fig, axs = plt.subplots(2, len(distssample) // 2, figsize=(5 + 3 * (len(distssample) // 2), 10))
            slopes = []
            for i, (d, LDd, ax) in enumerate(izip(distssample, LDsamples, axs.ravel())):
                indsm = (np.abs(LDd) < 0.010).all(axis=0).nonzero()[0]
                np.random.shuffle(indsm)
                indsm = indsm[:300]
                indsm.sort()
                ax.scatter(LDd[0][indsm], LDd[1][indsm], edgecolor='none',
                           s=40,
                           facecolor=cm.jet(int(255.0 * i / len(distssample))),
                           label=str(d),
                           alpha=0.3)
                x = np.linspace(-0.05, 0.05)
                ax.set_xlabel('LD, PCR1')
                ax.set_ylabel('LD, PCR2')
                #ax.plot(x, x, lw=1, c='k', alpha=0.4)
                ax.plot([0] * len(x), x, lw=2, c='k', alpha=0.3)
                ax.plot(x, [0] * len(x), lw=2, c='k', alpha=0.3)
                ax.set_xlim(-0.03, 0.03)
                ax.set_ylim(-0.03, 0.03)
                (r, Pv) = pearsonr(LDd[0][indsm], LDd[1][indsm])
                (R, PvR) = spearmanr(LDd[0][indsm], LDd[1][indsm])
                alpha = np.inner(LDd[0][indsm], LDd[1][indsm]) / \
                        np.inner(LDd[0][indsm], LDd[0][indsm])
                ax.plot(x, alpha * x, lw=2, c='r')
                #(alpha, y0) = np.polyfit(LDd[0], LDd[1], deg=1)
                #ax.plot(x, y0 + alpha * x, lw=2, c='r')
                
                ax.text(-0.028, 0.01,
                        'r = '+'{:2.2f}'.format(r)+\
                        '\nR = '+'{:2.2f}'.format(R)+\
                        ',\n$\\alpha = $'+'{:2.2f}'.format(alpha))
                ax.set_title(str(d))
            fig.suptitle(', '.join([pname, fragment, samplename]))
            plt.tight_layout(rect=(0, 0, 1, 0.96))

            ## Make a correlation plot between distance and logdifference in cofreqs
            #LDdiff = np.log(LD[1].ravel()[ind] / LD[0].ravel()[ind])
            #LDdiff[np.isnan(LDdiff) | np.isinf(LDdiff)] = np.ma.masked
            #ind2 = -LDdiff.mask
            #distuni = np.arange(500)
            #LDdiffuni = np.array([np.median(LDdiff[(distind == d) & ind2].data) for d in distuni])
            #LDdiffunim = np.array([np.mean(LDdiff[(distind == d) & ind2].data) for d in distuni])
            #plt.figure()
            #plt.scatter(distuni, LDdiffuni, color='b', alpha=0.5, label='Median')
            #plt.scatter(distuni, LDdiffunim, color='r', alpha=0.5, label='Mean')
            #plt.xlabel('Distance between sites in the pair[bp]')
            #plt.ylabel('Median/mean of LDdiff')
            #plt.title(', '.join([pname, fragment, samplename]))
            #plt.legend(loc=2)
            #plt.tight_layout()

            ## Plot cumulative hists
            #distssample = [10, 100, 200, 300, 400, 450]
            #LDdiffsamples = [(LDdiff[(distind < d + 10) & (distind > d - 10) & ind2].data) for d in distssample]
            #map(lambda x: x.sort(), LDdiffsamples)
            #plt.figure()
            #for i in xrange(len(distssample)):
            #    plt.plot(LDdiffsamples[i], np.linspace(0, 1, len(LDdiffsamples[i]))[::-1],
            #             lw=2,
            #             color=cm.jet(int(255.0 * i / len(distssample))),
            #             label=str(distssample[i]))
            #plt.plot([0, 0], [0, 1], lw=1.5, ls='--', c='k', alpha=0.3)
            #plt.xlabel('log(LD(PCR1) / LD(PCR2))')
            #plt.ylabel('Fraction of pairs above this x')
            #plt.legend(loc=1)
            #plt.title(', '.join([pname, fragment, samplename]))
            #plt.tight_layout()


            plt.ion()
            plt.show()

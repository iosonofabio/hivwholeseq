# vim: fdm=marker
'''
author:     Fabio Zanini
date:       04/12/14
content:    Check the allele frequencies in overlaps of subsequent fragments.
'''
# Modules
import os
import sys
import argparse
from itertools import izip
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

import hivwholeseq.plot_utils
from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.patients.samples import SamplePat



# Functions
def get_map_overlap(sample, fr1, fr2):
    '''Get a coordinate map of the overlap between the two fragments'''
    import numpy as np
    from seqanpy import align_ladder

    seq1 = sample.get_reference(fr1)
    seq2 = sample.get_reference(fr2)
    (score, ali1, ali2) = align_ladder(seq1, seq2, score_gapopen=-20)
    start2 = len(ali2) - len(ali2.lstrip('-'))
    end1 = len(ali1.rstrip('-'))
    
    mapco = []
    pos1 = start2
    pos2 = 0
    for i in xrange(start2, end1):
        if (ali1[i] != '-') and (ali2[i] != '-'):
            mapco.append((pos1, pos2))

        if ali1[i] != '-':
            pos1 += 1

        if ali2[i] != '-':
            pos2 += 1

    return np.array(mapco, int)


def get_allele_frequency_overlap(samples, overlaps, cov_min=1000,
                                  VERBOSE=0, qual_min=30):
    '''Get allele frequency in the overlaps'''
    data = [] 
    for io, overlap in enumerate(overlaps):
        fr1 = overlap[:2]
        fr2 = 'F'+overlap[-1]

        for samplename, sample in samples.iterrows():
            sample = SamplePat(sample)

            if VERBOSE >= 1:
                print overlap, samplename

            # FIXME: actually use frequencies
            try:
                ac1 = sample.get_allele_counts(fr1, qual_min=qual_min)
                ac2 = sample.get_allele_counts(fr2, qual_min=qual_min)
            except IOError:
                continue

            coord_map = get_map_overlap(sample, fr1, fr2)

            acjoint = np.zeros((2, ac1.shape[0], coord_map.shape[0]), int)
            acjoint[0] = ac1[:, coord_map[:, 0]]
            acjoint[1] = ac2[:, coord_map[:, 1]]

            # Convert to frequencies
            afjoint = np.ma.array(acjoint)
            cov = acjoint.sum(axis=1)
            mask = (cov < cov_min).any(axis=0)
            mask = np.tile(mask, (afjoint.shape[0], afjoint.shape[1], 1))
            afjoint.mask = mask
            afjoint = (1.0 * afjoint.swapaxes(0, 1) / afjoint.sum(axis=1)).swapaxes(0, 1)

            data.append({'af': afjoint,
                         'samplename': samplename,
                         'overlap': overlap,
                         'io': io,
                         'n_templates': sample.get_n_templates_dilutions(),
                         'coverage': cov})

    return data



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele counts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--overlaps', nargs='*',
                        help='Fragments to analyze (e.g. F1-2)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')
    parser.add_argument('--qualmin', type=int, default=30,
                        help='Minimal quality of base to call')
    parser.add_argument('--covmin', type=int, default=100,
                        help='Minimal coverage: lower -> more, higher -> better')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    overlaps = args.overlaps
    VERBOSE = args.verbose
    qual_min = args.qualmin
    cov_min = args.covmin
    use_plot = args.plot
    use_logit = args.logit

    samples = lssp()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()

    if not overlaps:
        overlaps = ['F'+str(i)+'-'+str(i+1) for i in xrange(1, 6)]
    if VERBOSE >= 3:
        print 'ovarlaps', overlaps

    data = get_allele_frequency_overlap(samples, overlaps, cov_min=cov_min,
                                        VERBOSE=VERBOSE, qual_min=qual_min)

    if use_plot:

        # Plot allele frequencies
        fig, ax = plt.subplots()
        for ida, datum in enumerate(data):
            afjoint = datum['af']
            color = cm.jet(1.0 * ida / len(data))
            ax.scatter(afjoint[0].ravel(), afjoint[1].ravel(),
                       s=50,
                       color=color,
                       alpha=0.7,
                       label=datum['samplename']+', '+datum['overlap'],
                       edgecolor='none')
        
        xmin = 1e-4
        ax.plot([xmin, 1 - xmin], [xmin, 1 - xmin], lw=2, color='k', alpha=0.5)

        ax.set_xlabel('allele frequency leading fragment')
        ax.set_ylabel('allele frequency trailing fragment')
        ax.grid(True)
        #ax.legend(loc=2, fontsize=10, title='Overlap:')
        if (samplenames is not None) and len(samplenames) == 1:
            ax.set_title(samplenames[0])
        elif (pnames is not None) and (len(pnames) == 1):
            ax.set_title(pnames[0])

        if use_logit:
            ax.set_xscale('logit')
            ax.set_yscale('logit')
            ax.set_xlim(xmin, 1 - xmin)
            ax.set_ylim(xmin, 1 - xmin)

        # Plot stddev of a certain number of molecules in Poisson sampling
        n = 300
        x = np.linspace(-4, 0, 1000)
        x = 1.0 / (1 + 10**(-x))
        y = x - np.sqrt(x / n)
        ax.plot(np.concatenate([x, 1 - y[::-1]]), np.concatenate([y, 1 - x[::-1]]), lw=4, c='black', alpha=0.7)
        ax.plot(np.concatenate([y, 1 - x[::-1]]), np.concatenate([x, 1 - y[::-1]]), lw=4, c='black', alpha=0.7,
                label='Poisson noise, n = '+str(n))

        plt.tight_layout()

        ## Plot correlation between template number and y scatter
        ## Bin by allele frequency, as there is more scatter close to 0/1
        #bins = np.array([0.01, 0.05, 0.25, 0.75, 0.95, 0.99])
        #binsc = np.array([0.02, 0.12, 0.5, 0.88, 0.98])
        #dcorr = [[] for b in bins[:-1]]
        #for datum in data:
        #    afjoint = datum['af']
        #    n_templates = datum['n_templates']
        #    if np.isnan(n_templates):
        #        continue
        #    tmp = np.vstack([np.mean(afjoint, axis=0).ravel(),
        #                     np.abs(np.diff(afjoint, axis=0)[0].ravel())])
        #    for ib in xrange(len(bins) - 1):
        #        ind = (tmp[0] >= bins[ib]) & (tmp[0] <= bins[ib + 1]) & \
        #              (-tmp.mask.any(axis=0))
        #        dcorr[ib].extend(zip(tmp[1, ind].data,
        #                             np.repeat(n_templates, ind.sum())))

        #fig, ax = plt.subplots()
        #ax.set_xlabel('# templates from dilutions')
        #ax.set_ylabel('|af1 - af2|')
        #ax.grid(True)
        #ax.set_ylim(1e-8, 1)
        #ax.set_yscale('log')
        #if (samplenames is not None) and len(samplenames) == 1:
        #    ax.set_title(samplename)
        #elif (pnames is not None) and (len(pnames) == 1):
        #    ax.set_title(pnames[0])
 
        #for ip, (bc, dc) in enumerate(izip(binsc, dcorr)):
        #    color = cm.jet(1.0 * ip / len(binsc))
        #    (y, x) = zip(*dc)
        #    r = np.corrcoef(x, y)[0, 1]
        #    ax.scatter(x, y, s=50, c=color, alpha=0.5,
        #               label=str(bc)+', r = '+'{:1.2f}'.format(r))

        #ax.legend(loc=3, fontsize=10)

        plt.tight_layout()
        plt.ion()
        plt.show()



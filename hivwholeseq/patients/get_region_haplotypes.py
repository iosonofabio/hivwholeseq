# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/11/14
content:    Get local haplotypes from single read pairs, including insertions
            and deletions. This includes aggressive clustering to keep the
            multiple sequence alignments efficient.
'''
# Modules
import os
import argparse
from operator import itemgetter
from itertools import izip
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.patients.patients import load_patient
from hivwholeseq.utils.argparse import RoiAction
from hivwholeseq.utils.sequence import build_msa_haplotypes as build_msa



# Functions
def plot_haplotype_frequencies(times, hft, figax=None, title='',
                               picker=None):
    '''Plot haplotype frequencies'''
    import hivwholeseq.utils.plot
    from matplotlib import cm
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style('darkgrid')
    fs = 16

    if figax is None:
        fig, ax = plt.subplots(figsize=(12, 7))
    else:
        fig, ax = figax

    # TODO: The hard part is finding an ordering
    hft_cum = hft.cumsum(axis=1)

    # Randomize colors to make them more visible
    colors = cm.jet(1.0 * np.arange(hft.shape[1]) / hft.shape[1])
    np.random.shuffle(colors)

    # Use fake zero/one for logit plots
    freqmin = 1e-6

    # Plot first line
    ax.fill_between(times, hft_cum[:, 0], freqmin + np.zeros(hft.shape[0]), color=colors[0],
                    label=str(0),
                    picker=picker)
    for i in xrange(1, hft.shape[1]):
        ax.fill_between(times, hft_cum[:, i],
                        np.minimum(1-freqmin, hft_cum[:, i - 1]), color=colors[i],
                        label=str(i),
                        picker=picker)

    ax.set_xlabel('Time from infection [days]', fontsize=fs)
    ax.set_ylabel('Haplotype frequency', fontsize=fs)
    ax.set_ylim(1e-4, 1 - 1e-4)
    ax.set_xlim(times[0], times[-1])
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)

    if title:
        ax.set_title(title, fontsize=fs)

    return (fig, ax)
    


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get local haplotypes',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--region', required=True,
                        help='Genomic region (e.g. V3 p17)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot local haplotype trajectories')

    args = parser.parse_args()
    pname = args.patient
    region = args.region
    VERBOSE = args.verbose
    use_plot = args.plot

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()

    if VERBOSE >= 1:
        print patient.name, region


    hct, ind, seqs = patient.get_haplotype_count_trajectory(region,
                                                            aligned=True)
    hft = (1.0 * hct.T / hct.sum(axis=1)).T


    if use_plot:
        times = patient.times[ind]
        plot_haplotype_frequencies(times, hft, title=patient.code+', '+region)

        plt.tight_layout()
        plt.ion()
        plt.show()

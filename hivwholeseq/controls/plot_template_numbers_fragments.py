# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/02/15
content:    Plot the result of the estimate of fragment-specific template numbers.
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

from hivwholeseq.utils.miseq import alpha
from hivwholeseq.patients.samples import load_samples_sequenced, SamplePat



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Plot estimated number of templates',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    VERBOSE = args.verbose

    samples = load_samples_sequenced()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()


    # Plot distributions
    fig, ax = plt.subplots()
    for ifr in xrange(1, 7):
        fragment = 'F'+str(ifr)
        if VERBOSE >= 2:
            print fragment

        color = cm.jet(1.0 * ifr / 6.0)

        n = np.array(samples.loc[:, fragment+'q'])
        n = n[-np.isnan(n)]
        ax.plot(np.sort(n), 1.0 - np.linspace(0, 1, len(n)),
                color=color, lw=2,
                label=fragment)

    ax.set_xlabel('# templates')
    ax.set_ylabel('fraction of samples > x')
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlim(0.1, 1e5)
    ax.set_xscale('log')
    ax.grid(True)
    ax.legend(loc='lower left', title='Fragment:')
    ax.set_title('Number of templates to RT-PCR')

    plt.tight_layout()

    plt.ion()
    plt.show()

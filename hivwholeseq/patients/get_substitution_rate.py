# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/09/14
content:    Get substitution rate from divergence.
'''
# Modules
import os
import argparse
from operator import itemgetter
import numpy as np
from Bio import SeqIO
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import load_patients, Patient



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get divergence and diversity',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. F1 V3)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    pnames = patients.index.tolist()

    data = []

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        for ifr, region in enumerate(regions):
            if VERBOSE >= 1:
                print pname, region

            try:
                dg, ind = patient.get_divergence(region, cov_min=10)
            except ValueError:
                continue
            times = patient.times[ind]

            data.append({'pname': pname, 'region': region, 'dg': dg, 't': times})

    if VERBOSE >= 1:
        print 'Fit slopes'

    for d in data:
        pname = d['pname']
        region = d['region']
        dg = d['dg']
        times = d['t']

        r = np.linalg.lstsq(times[:, np.newaxis], dg)[0][0]
        d['r'] = r


    if plot:
        if VERBOSE >= 1:
            print 'Plot'

        ## Plot divergence and the fit
        #fig, axs = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={'width_ratios':[4, 1]})
        #ax = axs[0]
        #ax.set_xlabel('Time from transmission [days]')
        #ax.set_ylabel('Divergence')
        ##ax.set_yscale('log')

        #for i, d in enumerate(data):
        #    pname = d['pname']
        #    region = d['region']
        #    dg = d['dg']
        #    times = d['t']
        #    r = d['r']

        #    ax.plot(times, dg, lw=2, label=', '.join([pname, region]),
        #            color=cm.jet(1.0 * i / len(data)))

        #    ax.plot(times, times * r, lw=2, ls='--',
        #            color=cm.jet(1.0 * i / len(data)))

        #    axs[1].scatter(0, r * 365.24, color=cm.jet(1.0 * i / len(data)),
        #                   label=', '.join([pname, region]))

        #ax.legend(loc=4, fontsize=12)
        #ax.grid(True)

        #ax = axs[1]
        #ax.xaxis.set_ticklabels('')
        #ax.set_ylabel('Substitution rate [changes/site/yr]')
        #ax.set_ylim(0.96 * 365.24 * min(data, key=itemgetter('r'))['r'],
        #            1.04 * 365.24 * max(data, key=itemgetter('r'))['r'])
        #ax.set_xlim(-0.5, 0.5)
        #ax.grid(True)
        #plt.tight_layout()

        # Plot just the slope
        nrows = 1 + (len(regions) - 1) // 4
        ncols = min(len(regions), 4)
        fig, axs = plt.subplots(nrows, ncols, figsize=(2.5 * ncols, 4 * nrows))
        axs = axs.ravel()
        
        # Clean up empty axes
        for ir in xrange(0, len(axs) - len(regions)):
            axs[-1 - ir].xaxis.set_visible(False)
            axs[-1 - ir].yaxis.set_visible(False)


        ylim = (0.96 * 365.24 * min(data, key=itemgetter('r'))['r'],
                1.04 * 365.24 * max(data, key=itemgetter('r'))['r'])

        for ir, region in enumerate(regions):
            datar = filter(lambda x: x['region'] == region, data)
            ax = axs[ir]
            ax.set_xlim(-0.5, 0.5)
            ax.set_ylim(*ylim)
            ax.set_yscale('log')
            ax.grid(True)
            ax.xaxis.set_ticklabels('')
            if ir % ncols:
                ax.yaxis.set_ticklabels('')
            else:
                ax.set_ylabel('Substitution rate [changes/site/yr]')
            ax.set_title(region)

            for i, d in enumerate(datar):
                pname = d['pname']
                r = d['r']
                if len(pnames) == 1:
                    x = 0
                else:
                    x = -0.2 + 0.4 * pnames.index(pname) / (len(pnames) - 1)
                ax.scatter(x, r * 365.24, color=cm.jet(1.0 * i / len(datar)),
                           label=', '.join([pname, region]))

        plt.tight_layout()

        plt.ion()
        plt.show()


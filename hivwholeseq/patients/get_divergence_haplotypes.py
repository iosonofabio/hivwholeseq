# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/12/14
content:    Study the divergence from founder sequence for haplotypes.
'''
# Modules
import os
import argparse
from itertools import izip
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.sequence_utils import align_muscle



# Functions
def get_distance_distributions(hct, alim, VERBOSE=0):
    '''Get haplotype distance distribution from initial major allele'''
    iseq0 = hct[0].argmax()
    seq0 = alim[iseq0]

    ds = []
    for it in xrange(hct.shape[0]):
        d = []
        for iseq, ali in enumerate(alim):
            c = hct[it, iseq]
            if c:
                d.extend([(ali != seq0).mean()] * c)

        d.sort()
        ds.append(d)

    return ds



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get divergence of local haplotypes',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--regions', required=True, nargs='+',
                        help='Genomic region (e.g. V3)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot local haplotype trajectories')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    use_plot = args.plot

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    data = []
    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        for region in regions:
            if VERBOSE >= 1:
                print patient.name, region

            if VERBOSE >= 2:
                print 'Get haplotype counts'
            (hct, ind, seqs) = patient.get_region_count_trajectories(region,
                                                                     VERBOSE=VERBOSE)
            
            times = patient.times[ind]

            if VERBOSE >= 2:
                print 'Align sequences'
            ali = align_muscle(*seqs, sort=True)
            alim = np.array(ali)

            if VERBOSE >= 2:
                print 'Get distributions'
            ds = get_distance_distributions(hct, alim, VERBOSE=VERBOSE)

            data.append({'pname': pname, 'ds': ds, 't': times,
                         'region': region})
        

    if use_plot:
        for pname, patient in patients.iterrows():
            patient = Patient(patient)

            nrows = 1 + (len(regions) - 1) // 3
            ncols = min(len(regions), 3)
            fig, axs = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows),
                                    sharey=True)
            if len(regions) > 1:
                axs = axs.ravel()
            else:
                axs = [axs]
            fig.suptitle(pname)
            
            # Clean up empty axes
            for ir in xrange(0, len(axs) - len(regions)):
                axs[-1 - ir].xaxis.set_visible(False)
                axs[-1 - ir].yaxis.set_visible(False)

            for ir, region in enumerate(regions):
                datum = filter(lambda x: (x['pname'] == pname) and
                                         (x['region'] == region),
                               data)[0]
                ax = axs[ir]
                ax.set_ylim(-0.05, 1.05)
                ax.set_xlim(0, 0.06)

                ax.grid(True)
                ax.set_title(region)

                ax.set_xlabel('Genetic distance [changes/bp]')

                if ir % ncols:
                    ax.yaxis.set_ticklabels('')
                else:
                    ax.set_ylabel('Fraction of haplotypes > x')

                times = datum['t']
                ds = datum['ds']

                for it, (t, d) in enumerate(izip(times, ds)):
                    color = cm.jet(1.0 * (t - patient.times[0]) / \
                                   (patient.times[-1] - patient.times[0]))
                    ax.plot(np.sort(d), 1.0 - np.linspace(0, 1, len(d)),
                            lw=2, color=color,
                            label=str(int(t))+' days',
                           )

            plt.tight_layout(rect=(0, 0, 1, 0.96))

            plt.savefig('/ebio/ag-neher/share/users/fzanini/phd/sequencing/figures/haplotype_divergence_'+pname+'.png')
            plt.ion()
            plt.show()

# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/14
content:    Study whether HIV in a patient tends to explore the mutational space
            in a way that comes closer to a subtype average (entropy profile).
'''
# Modules
import os, sys
import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alphal, alpha
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.reference import load_custom_reference, load_custom_alignment
from hivwholeseq.cross_sectional.get_subtype_allele_frequencies import (
    get_subtype_reference_alignment_allele_frequencies)


# Globals
refname = 'HXB2'
Sbins = [0, 0.05, 0.2, 0.5, 2]



# Functions
def trim_comap(comap, ref1, ref2, VERBOSE=0):
    '''Trim comap to relevant region'''
    from seqanpy import align_overlap

    (score, ali1, ali2) = align_overlap(ref1, ref2, score_gapopen=-20)
    # Start/end of subtype alignment
    start = len(ali2) - len(ali2.lstrip('-'))
    end = len(ali2.rstrip('-'))

    if VERBOSE >= 3:
        from hivwholeseq.utils.sequence import pretty_print_pairwise_ali
        pretty_print_pairwise_ali((ali1, ali2), width=100, name1='pat', name2='ali')

    #FIXME: is this correct?
    ind = (comap[:, 0] >= start) & (comap[:, 0] < end)
    comap = comap[ind]
    comap[:, 0] -= start

    return comap



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Mutations away from/towards subtype',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Genomic regions (e.g. V3 IN)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', type=int, default=0,
                        help='Plot results')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    use_plot = args.plot

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    

    data = []
    for region in regions:
        if VERBOSE >= 1:
            print region

        # Load HXB2
        refseq = load_custom_reference(refname, 'gb')
        refseq.name = refname
        refm = np.array(refseq)

        # Load subtype alignment
        afs = get_subtype_reference_alignment_allele_frequencies(region)


        for pname, patient in patients.iterrows():
            patient = Patient(patient)

            if VERBOSE >= 1:
                print region, pname

            comap = patient.get_map_coordinates_reference(region,
                                                          refname=(refname, region))

            aft, ind = patient.get_allele_frequency_trajectories(region)
            
            # Sometimes the alignment is trimmed
            comap = trim_comap(comap, alpha[aft[0].argmax(axis=0)], alpha[afs.argmax(axis=0)],
                               VERBOSE=VERBOSE)

            afsp = afs[:, comap[:, 0]]
            aft = aft[:, :, comap[:, 1]]

            # NOTE: ignore gaps for now
            afsp = afsp[:4]
            aft = aft[:, :4]

            # First, check the fraction of sites for which initial consensus = subtype consensus
            consi = aft[0].argmax(axis=0)
            conss = afsp.argmax(axis=0)
            print 'Fraction of sites at which patient initial consensus agrees with subtype:',
            print '{:2.0%}'.format((consi == conss).mean())

            # Check at the end
            consf = aft[-1].argmax(axis=0)
            print 'Fraction of sites at which patient final consensus agrees with subtype:',
            print '{:2.0%}'.format((consf == conss).mean())

            # Classify mutations in away/towards
            mclass = np.ma.masked_all(afsp.shape, 'S5')
            for ia in xrange(mclass.shape[0]):
                # From B
                ipos = (consi != ia) & (consi == conss)
                mclass[ia, ipos] = 'fromB'

                # To B
                ipos = (consi != ia) & (conss == ia)
                mclass[ia, ipos] = 'toB'

            aft_fromB = aft[:, mclass == 'fromB']
            aft_toB = aft[:, mclass == 'toB']

            if use_plot >= 2:
                fig, ax = plt.subplots()

                ax.plot(patient.times[ind], aft_fromB.mean(axis=1),
                        lw=2,
                        label='fromB')
                ax.plot(patient.times[ind], aft_toB.mean(axis=1),
                        lw=2,
                        label='toB')

                ax.set_xlabel('Time [days from infection]')
                ax.set_ylabel('Allele frequency')

                ax.grid(True)
                ax.legend(loc=2, fontsize=12)
                ax.set_title(pname+', '+region)
                
                plt.tight_layout()

            
            # Classify by entropy in the subtype
            from hivwholeseq.one_site_statistics import get_entropy
            S = get_entropy(afsp)

            # Add to global data structure
            for key in ('fromB', 'toB'):
                for ibin in xrange(len(Sbins) - 1):
                    pos = (S >= Sbins[ibin]) & (S < Sbins[ibin + 1])
                    posM = np.tile(pos, (afsp.shape[0], 1))
                    posM &= mclass == key

                    for it, t in enumerate(patient.times[ind]):
                        for ia, pos in zip(*(posM.nonzero())):

                            # Skip masked alleles
                            if aft.mask[it, ia, pos]:
                                continue

                            af = aft[it, ia, pos]
                            data.append((region, pname, key, ibin, ia, pos, t, af))



    data = pd.DataFrame.from_records(data=data,
                        columns=('region', 'patient', 'mclass', 'ibin', 'inuc', 'pos', 'time', 'af'))


    if use_plot:
        fig, ax = plt.subplots(figsize=(10, 7))
        marker = {'fromB': 'o', 'toB': 's'}
        ls = {'fromB': '--', 'toB': '-'}
        for key in ('toB', 'fromB'):
            for ibin in xrange(len(Sbins) - 1):
                color = cm.jet(1.0 * ibin / (len(Sbins) - 1))
                datum = (data.loc[(data['mclass'] == key) & (data['ibin'] == ibin)]
                         .loc[:, ['time', 'af']]
                         .groupby('time')
                         .mean())['af']

                x = np.array(datum.index)
                y = np.array(datum)

                # Multiply the fromB by 3, because entropy mostly refer to a 2 state model
                if key == 'fromB':
                    y *= 3

                ax.scatter(x, y,
                           s=10,
                           lw=2,
                           marker=marker[key],
                           color=color,
                           label=key+', S e ['+str(Sbins[ibin])+', '+str(Sbins[ibin+1])+']',
                          )

                xfit = np.linspace(0, x.max(), 1000)

                # Fit linear increase
                m = np.dot(x, y) / np.dot(x, x)
                yfit = m * xfit
                
                # Fit exponential saturation
                from scipy.optimize import curve_fit
                fun = lambda x, l, u: l * (1 - np.exp(- u/l * x))
                try:
                    l, u = curve_fit(fun, x, y, p0=(y[-1], 1e-5))[0]
                except RuntimeError:
                    continue

                # Constrain the fit to a maximal frequency of 1
                if l > 1:
                    l = 1
                    fun_red = lambda x, u: 1 - np.exp(- u * x)
                    try:
                        u = curve_fit(fun_red, x, y, p0=(1e-5,))[0][0]
                    except RuntimeError:
                        continue

                yfit = fun(xfit, l, u)
                ax.plot(xfit, yfit,
                        ls=ls[key],
                        lw=2,
                        alpha=0.5,
                        label=('l = '+'{:2.0e}'.format(l)+', '+
                               'u = '+'{:2.0e}'.format(u)+''),
                        color=color)


        ax.set_xlabel('Time [days from infection]')
        ax.set_ylabel('Allele frequency')

        ax.grid(True)
        ax.set_yscale('log')
        ax.set_ylim(1e-4, 1)

        if len(regions) == 1:
            ax.set_title(region)
        
        ax.set_title(',\n'.join(map(' '.join, (np.sort(np.array(patients.code)), regions))), fontsize=12)

        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', fontsize=12, bbox_to_anchor=(1, 0.5))

        plt.tight_layout(rect=(0, 0, 0.75, 1))

        plt.ion()
        plt.show()


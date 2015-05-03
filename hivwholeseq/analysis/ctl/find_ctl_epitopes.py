# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/05/15
content:    Find CTL epitopes in our patients based on the LANL table.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from hivwholeseq.patients.patients import load_patients, iterpatient
from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.cross_sectional.ctl_epitope_map import get_ctl_epitope_map



# Globals
regions = ['p17', 'p24', 'PR', 'RT', 'IN']



# Functions
def collect_data(pnames, regions, VERBOSE=0):
    '''Collect data for CTL epitope map'''
    import pandas as pd

    patients = load_patients(pnames=pnames)
    ctl_table = get_ctl_epitope_map(species='human')

    data = []
    for region in regions:
        for pname, patient in iterpatient(patients):
            if VERBOSE >= 1:
                print region, patient.code

            seq = patient.get_reference(region)
            prot = seq.seq.translate()

            ind = [i for i, epi in enumerate(ctl_table['Epitope'])
                   if epi in prot]

            datum = ctl_table.loc[ind].copy()
            datum['pcode'] = patient.code
            datum['region'] = region
            data.append(datum)
    
    data = pd.concat(data)

    return data


def get_sharing(data):
    '''Print info on sharing of epitopes'''
    share = (data.loc[:, ['Epitope', 'pcode', 'region']]
               .groupby(['region', 'Epitope'])
               .size())
    return share


def plot_sharing(share, regions=None):
    '''Plot the sharing of epitopes'''

    if regions is not None:
        share = share.loc[regions]

    # Every epitope is present in at least one patient
    y = np.bincount(share)[1:]
    x_left = np.arange(len(y)) + 0.5
    width = 1.0

    fig, ax = plt.subplots()
    ax.bar(x_left, y, width=width, color='steelblue')
    ax.set_xlabel('# patient with epitope')
    ax.set_ylabel('# epitopes')
    ax.set_xlim(0.5, len(y) + 0.5)
    ax.grid(True)

    regions = share.index.get_level_values('region').unique()
    ax.set_title(', '.join(regions))

    plt.tight_layout()




# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Study accumulation of minor alleles for different kinds of mutations',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', action=PatientsAction,
                        help='Patient to analyze')
    parser.add_argument('--regions', default=regions, nargs='+',
                        help='Regions to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot results')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot

    data = collect_data(pnames, regions, VERBOSE=VERBOSE)

    share = get_sharing(data)
    if VERBOSE >= 1:
        print share

    if plot:
        plot_sharing(share)

        plt.ion()
        plt.show()

# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/11/14
content:    Get the number of molecules to PCR (templates) using the dilution
            series.
'''
# Modules
import os
import datetime
import argparse
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.patients.patients import load_patients, Patient


# Functions
def get_dilution(dilstr):
    '''Get dilution numbers from string'''
    if not isinstance(dilstr, basestring):
        return dilstr

    else:
        dil_factor = float(dilstr.split()[0].split(':')[1])
        dil_tries = map(int, dilstr.split()[1][1:-1].split('/'))
        return [dil_factor] + dil_tries


def estimate_ntemplates_Poisson(dilution):
    '''Estimate number of templates from dilution info'''
    if np.isscalar(dilution):
        return np.nan
    
    #TODO: finish this!
    if (dilution[1] == 1) and (dilution[2] == 2):
        return np.log(2) * dilution[0]

    else:
        if (dilution[1] == 0):
            dil_good = dilution[0] / 10
            dil_bad = dilution[0]

        elif (dilution[1] == dilution[2]):
            dil_good = dilution[0]
            dil_bad = dilution[0] * 10

        else:
            raise ValueError('Dilution info not understood')

        return dil_good * np.log(1 + dil_bad / dil_good)


def get_template_number(dilstr):
    '''Get the Poisson estimate for template numbers from dilutions'''
    dilution = get_dilution(dilstr)
    # NOTE: for each fragment, we have TWO PCR reactions in parallel, but the
    # dilution series is based on ONE only (i.e. 1/12 of the 400 ul)
    n_templates =  2 * estimate_ntemplates_Poisson(dilution)

    return n_templates



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot template number info')

    args = parser.parse_args()
    pnames = args.patients
    VERBOSE = args.verbose
    use_plot = args.plot

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    data = [] 
    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        if VERBOSE:
            print pname, patient.code


        samples = patient.samples
        n_approx = samples['templates approx']
        dils = [get_dilution(x) for x in samples['dilutions']]
        n_dils = [2 * estimate_ntemplates_Poisson(x) for x in dils]

        # Attach sample date info
        age = np.array((datetime.datetime.now() - samples.date)) / 86400e9

        data.append({'n_approx': n_approx, 'n_dil': n_dils, 'age': age})

    if use_plot:
        fig, axs = plt.subplots(1, 2, figsize=(13, 8))
        fig.suptitle('Correlation between viral load and PCR dilutions', fontsize=18)
        ax = axs[0]
        agemax = np.concatenate([data_pat['age'] for data_pat in data]).max()
        nmax = 1e5

        for i, data_pat in enumerate(data):
            color = cm.jet(1.0 * data_pat['age'] / agemax)
            ax.scatter(data_pat['n_approx'],
                       np.array(data_pat['n_dil']) * (1.05 * i / len(data)),
                       s=30, color=color)
        
        ax.set_xlabel('# templates (from viral load)')
        ax.set_ylabel('# templates (from dilutions)')
        ax.set_ylim((1, nmax))
        ax.set_xlim((1, nmax))

        ax.plot([1, nmax], [1, nmax], lw=2)

        ax.grid(True)

        ax.set_xscale('log')
        ax.set_yscale('log')

        ax = axs[1]
        x = np.concatenate([data_pat['n_approx'] for data_pat in data])
        y = np.concatenate([data_pat['n_dil'] for data_pat in data])
        ind = -(np.isnan(x) | np.isnan(y) | (x == 0) | (y == 0))
        x = x[ind]
        y = y[ind]
        age = np.concatenate([data_pat['age'] for data_pat in data])[ind]

        h = np.histogram(y / x, bins=np.logspace(-3, 0, 8), density=False)
        ax.bar(h[1][:-1], h[0], width=np.diff(h[1]), color='grey')
        ax.set_xlabel('PCR efficiency')
        ax.set_ylabel('# of samples')
        ax.grid(True)

        ax.set_xscale('log')

        plt.tight_layout(rect=(0, 0, 1, 0.95))   


        plt.ion()
        plt.show()


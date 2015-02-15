# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/09/14
content:    Get substitution rate from divergence.
'''
# Modules
import os
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient

from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.paper_figures.plots import plot_substitution_rate




# Globals
pnames = ['20097', '15363', '15823', '15376', '20529', '9669', '15241', '15319']
regions = ['p17', 'PR', 'RT', 'vif', 'V3', 'RRE', 'nef']



# Functions
def plot_divergence_fit(data, VERBOSE=0):
    from operator import itemgetter

    fig, axs = plt.subplots(1, 2, figsize=(10, 5),
                            gridspec_kw={'width_ratios':[4, 1]})
    ax = axs[0]
    ax.set_xlabel('Time from transmission [days]')
    ax.set_ylabel('Divergence')
    #ax.set_yscale('log')

    for i, d in enumerate(data):
        pcode = d['pcode']
        region = d['region']
        dg = d['dg']
        times = d['t']
        r = d['r']

        ax.plot(times, dg, lw=2, label=', '.join([pcode, region]),
                color=cm.jet(1.0 * i / len(data)))

        ax.plot(times, times * r, lw=2, ls='--',
                color=cm.jet(1.0 * i / len(data)))

        axs[1].scatter(0, r * 365.24, color=cm.jet(1.0 * i / len(data)),
                       label=', '.join([pcode, region]))

    ax.legend(loc=4, fontsize=12)
    ax.grid(True)

    ax = axs[1]
    ax.xaxis.set_ticklabels('')
    ax.set_ylabel('Substitution rate [changes/site/yr]')
    ax.set_ylim(0.96 * 365.24 * min(data, key=itemgetter('r'))['r'],
                1.04 * 365.24 * max(data, key=itemgetter('r'))['r'])
    ax.set_xlim(-0.5, 0.5)
    ax.grid(True)
    plt.tight_layout()



# Script
if __name__ == '__main__':

    VERBOSE = 2
    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'substitution_rates.pickle'


    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    pcodes = patients.code.tolist()

    data = []
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        for ifr, region in enumerate(regions):
            if VERBOSE >= 1:
                print pname, patient.code, region

            try:
                dg, ind = patient.get_divergence(region, cov_min=10)
            except ValueError:
                continue
            times = patient.times[ind]

            data.append({'pcode': patient.code, 'region': region, 'dg': dg, 't': times})

    if VERBOSE >= 1:
        print 'Fit slopes'

    for d in data:
        pcode = d['pcode']
        region = d['region']
        dg = d['dg']
        times = d['t']

        r = np.linalg.lstsq(times[:, np.newaxis], dg)[0][0]
        d['r'] = r


    if VERBOSE >= 1:
        print 'Prepare data for plot'

    datap = pd.DataFrame([{'pcode': d['pcode'],
                           'region': d['region'],
                           'rate': d['r'] * 365.24}
                          for d in data])

    if VERBOSE >= 1:
        print 'Save data for plot to pickle'
    datap.to_pickle(fn_data)

    if VERBOSE >= 1:
        print 'Plot'
    # Plot divergence and the fit
    plot_divergence_fit(data, VERBOSE=VERBOSE)

    # Plot just the slope
    filename = foldername+'substitution_rates'
    for ext in ['png', 'pdf', 'svg']:
        plot_substitution_rate(datap,
                               VERBOSE=VERBOSE,
                               savefig=filename+'.'+ext)

    plot_substitution_rate(datap, VERBOSE=VERBOSE)

    plt.ion()
    plt.show()


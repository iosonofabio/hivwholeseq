# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/09/14
content:    Get substitution rate from divergence.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient




# Globals
pnames = ['20097', '15363', '15823', '15376', '20529', '9669', '15241', '15319']
regions = ['p17', 'p24', 'PR', 'RT', 'IN', 'vif', 'vpr', 'vpu', 'V3', 'RRE', 'gp41', 'nef']
window_size = 300



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


def plot_region_boxes(regions, ax, refname='HXB2', VERBOSE=0):
    '''Plot boxes for genomic regions of HXB2'''
    from matplotlib.patches import Rectangle

    from hivwholeseq.reference import load_custom_reference
    refseq = load_custom_reference(refname, format='gb')

    data = []
    
    y1 = 5
    height = 5
    pad = 2
    for feature in refseq.features:
        if feature.id in regions:
            x = [feature.location.nofuzzy_start, feature.location.nofuzzy_end]
            data.append({'name': feature.id,
                         'x1': x[0], 'x2': x[1], 'width': x[1] - x[0]})

    data = pd.DataFrame(data)
    data.sort('x1', inplace=True)
    data.index = np.arange(len(data))
    data['height'] = height
    data['parity'] = ((1 + np.arange(len(data))) % 2)
    data['row'] = 'bottom'; data.loc[np.array(data['parity'] == 1), 'row'] = 'top'
    data['y1'] = y1 + (height + pad) * data['parity']
    data['y2'] = data['y1'] + data['height']

    for _, datum in data.iterrows():
        r = Rectangle((datum['x1'], datum['y1']), datum['width'], datum['height'],
                      facecolor=[0.9] * 3,
                      edgecolor='k',
                      label=datum['name'])
        
        xt = datum['x1'] + 0.5 * datum['width']
        yt = datum['y1'] + 0.5 * datum['height']
        if datum['row'] == 'top':
            yt +=  height + 0
        else:
            yt -= height + 3.5

        ax.add_patch(r)
        ax.text(xt, yt,
                datum['name'],
                color='k', 
                fontsize=16,
                ha='center')


def plot_substitution_rate(data, regions, title='', VERBOSE=0):
    '''Plot the substitution rates'''

    if VERBOSE:
        print 'Plot substitution rates'
    
    # Sort patient codes
    pcodes = sorted(set(data['pcode']), key=lambda x: int(x[1:]))

    cfun = lambda pcode: cm.jet(1.0 * pcodes.index(pcode) / len(pcodes))

    fig, axs = plt.subplots(2, 1,
                            sharex=True,
                            figsize=(10, 7),
                            gridspec_kw={'height_ratios':[4, 1]})

    ax = axs[0]
    for pcode in pcodes:
        datum = data.loc[data['pcode'] == pcode].iloc[0]

        y = datum['rate']
        x = datum['x']
        color = cfun(pcode)

        ax.plot(x, y, color=color, lw=2, label=pcode)

    ax.set_ylim(1e-4, 2e-1)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=14)
    ax.grid(True)
    ax.set_yscale('log')

    ax.set_ylabel('Substitution rate\n[changes / year / site]', labelpad=10, fontsize=14)
    ax.legend(loc='upper center',
              title='Patients',
              ncol=4,
              fontsize=14)

    if title:
        ax.set_title(title)

    ax = axs[1]
    ax.set_ylim(-4, 26)
    ax.set_xlabel('Position in HXB2 [bp]', fontsize=14)
    ax.set_yticks([])
    ax.grid(True, axis='x')

    plot_region_boxes(regions, ax, VERBOSE=VERBOSE)

    plt.tight_layout(rect=(0, 0, 0.98, 1), h_pad=0.001)



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Measure substitution rates in a sliding window',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=pnames,
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+', default=regions,
                        help='Regions to highlight in the plot (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot results')
    parser.add_argument('--window', type=int, default=window_size,
                        help='Size of sliding window')


    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot
    window_size = args.window

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    pcodes = patients.code.tolist()

    data = []
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        if VERBOSE >= 1:
            print pname, patient.code

        aft, ind = patient.get_allele_frequency_trajectories('genomewide',
                                                             cov_min=100,
                                                             depth_min=100)
        times = patient.times[ind]

        from hivwholeseq.patients.get_divergence_diversity_local import (
        get_divergence_diversity_sliding)
        x, dg, _ = get_divergence_diversity_sliding(aft, window_size, VERBOSE=VERBOSE)

        coomap = patient.get_map_coordinates_reference('genomewide')
        coomapd = pd.Series(coomap[:, 0], index=coomap[:, 1])
        x_sub = np.array(coomapd.loc[x.astype(int)])

        # Exclude missing positions
        ind_pos = -np.isnan(x_sub)
        x_sub = x_sub[ind_pos]
        dg = dg[:, ind_pos]

        # Pad a bit left and right for low coverage
        pad = 50
        x_sub = x_sub[pad: -pad]
        dg = dg[:, pad: -pad]

        data.append({'pcode': patient.code,
                     'x': x_sub,  # integer positions simplify a few things
                     'dg': dg,
                     't': times})

    if VERBOSE >= 1:
        print 'Fit slopes'

    for d in data:
        pcode = d['pcode']
        dg = d['dg']
        times = d['t']

        r = np.ma.masked_all(dg.shape[1])
        for pos, dg_pos in enumerate(dg.T):
            ind = -dg_pos.mask
            if ind.sum() < 3:
                continue

            y = dg_pos[ind].data
            x = times[ind]
            # Linear fit
            r[pos] = np.dot(y, x) / np.dot(x, x)

        d['r'] = r

    if VERBOSE >= 1:
        print 'Prepare data for plot'

    datap = pd.DataFrame([{'pcode': d['pcode'],
                           'x': d['x'],
                           'rate': d['r'] * 365.24}
                          for d in data])

    if plot:
        if VERBOSE >= 1:
            print 'Plot'

        #plot_divergence_fit(data, VERBOSE=VERBOSE)

        plot_substitution_rate(datap, regions, VERBOSE=VERBOSE)

        plt.ion()
        plt.show()

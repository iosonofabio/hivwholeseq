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
from hivwholeseq.paper_figures.plots import plot_substitution_rate_sliding




# Globals
pnames = ['20097', '15823', '15376', '20529', '9669', '15241', '15319']
regions = ["LTR5'", 'p17', 'p24', 'PR', 'RT', 'IN', 'vif', 'vpr', 'vpu',
           'V1', 'V2', 'V3', 'V4', 'V5',
           'RRE', 'gp41', 'nef',
           "LTR3'"]



# Functions
def get_subs_matrix(datap):
    from collections import Counter

    cou = Counter()
    for _, row in datap.iterrows():
        cou += Counter(row['x'].astype(int))

    poss = np.sort([pos for pos, c in cou.iteritems() if c == len(datap)])
    M = np.zeros((len(datap), len(poss)))
    for i, (_, row) in enumerate(datap.iterrows()):
        ind = np.array(pd.Series(row['x'].astype(int)).isin(poss))
        M[i] = row['rate'][ind]

    return M


def calc_fold_change(M):
    Mm = np.maximum(-4, np.log10(M))
    Mmean = Mm.mean(axis=0)
    avg = (Mm - Mmean).std(axis=0).mean() * np.log2(10)
    std = (Mm - Mmean).std(axis=0).std() * np.log2(10)
    print 'log2(fold change) = {:.1G} +- {:.1G}'.format(avg, std)




# Script
if __name__ == '__main__':

    VERBOSE = 2
    window_size = 300

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'substitution_rates_sliding.pickle'


    if not os.path.isfile(fn_data):
        from hivwholeseq.analysis.substitution_rate.rate_sliding_window import (
            collect_data, fit_substitution_rate)

        if VERBOSE >= 1:
            print 'Collect data'
        data = collect_data(pnames, VERBOSE=VERBOSE)

        if VERBOSE >= 1:
            print 'Fit substitution rates'
        datap = fit_substitution_rate(data, VERBOSE=VERBOSE)

        if VERBOSE >= 1:
            print 'Save data for plot to pickle'
        datap.to_pickle(fn_data)

    else:
        datap = pd.read_pickle(fn_data)


    # Print fold changes excluding p9
    M = get_subs_matrix(datap.loc[np.array(-datap['pcode'].isin(['p9']))])
    calc_fold_change(M)

    filename = foldername+'substitution_rates_sliding'
    for ext in ['png', 'pdf', 'svg']:
        plot_substitution_rate_sliding(datap,
                                       regions,
                                       VERBOSE=VERBOSE,
                                       savefig=filename+'.'+ext)

    plot_substitution_rate_sliding(datap, regions, VERBOSE=VERBOSE)
    plt.ion()
    plt.show()


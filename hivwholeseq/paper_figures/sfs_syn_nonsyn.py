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
from hivwholeseq.paper_figures.plots import plot_sfs_syn_nonsyn




# Globals
pnames = ['20097', '15363', '15823', '9669', '20529', '15241', '15376', '15319']
regions = ['p17', 'p24', 'PR', 'RT', 'IN', 'vif', 'vpu', 'nef', 'gp41']



# Functions
def store_data(data, fn):
    '''Store data to file for the plots'''
    import cPickle as pickle
    with open(fn, 'wb') as f:
        pickle.dump(data, f, protocol=-1)


def load_data(fn):
    '''Load the data for the plots'''
    import cPickle as pickle
    with open(fn, 'rb') as f:
        return pickle.load(f)


def test_significance(datap):
    '''Test statistically the difference between syn and nonsyn'''
    # Fisher's exact on 0.5
    M = np.zeros((2, 2), int)
    bin05 = (datap['afbin_center'] < 0.5).nonzero()[0][-1] + 1
    M[0] = datap.iloc[:bin05].loc[:, ['syn', 'nonsyn']].sum(axis=0)
    M[1] = datap.iloc[bin05 + 1:].loc[:, ['syn', 'nonsyn']].sum(axis=0)
    from scipy.stats import fisher_exact
    print "Fisher's exact test on 0.5: ", fisher_exact(M)



# Script
if __name__ == '__main__':

    VERBOSE = 2

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'sfs_syn_nonsyn.pickle'


    if not os.path.isfile(fn_data):
        from hivwholeseq.analysis.sfs.syn_nonsyn import (
            collect_data, get_sfs, get_number_fixed)

        if VERBOSE >= 1:
            print 'Collect data'
        data = collect_data(pnames, regions, VERBOSE=VERBOSE)

        if VERBOSE >= 1:
            print 'Calculate SFS'
        datap = get_sfs(data,
                        bins_af='linear',
                        normalize='fraction',
                        attrnames=['mutclass'],
                        VERBOSE=VERBOSE)

        data_fixed = get_number_fixed(data)

        if VERBOSE >= 1:
            print 'Save data for plot to pickle'
        store_data({'plot': datap, 'fixed': data_fixed}, fn_data)

    else:
        d = pd.read_pickle(fn_data)
        datap = d['plot']
        data_fixed = d['fixed']

    test_significance(datap)

    #filename = foldername+'sfs_synnonsyn_loglinear'
    #for ext in ['png', 'pdf', 'svg']:
    #    plot_sfs_syn_nonsyn(datap,
    #                        VERBOSE=VERBOSE,
    #                        savefig=filename+'.'+ext)

    #plot_sfs_syn_nonsyn(datap, VERBOSE=VERBOSE)
    #plt.ion()
    #plt.show()

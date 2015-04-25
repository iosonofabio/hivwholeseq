# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/09/14
content:    Get Linkage Disequilibrium (LD).
'''
# Modules
import os
import numpy as np
import pandas as pd

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient

from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.paper_figures.plots import plot_LD




# Globals
pnames = ['20097', '15363', '15823', '9669', '20529', '15241', '15376', '15319']
fragments = ['F'+str(i) for i in xrange(1, 7)]



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



# Script
if __name__ == '__main__':

    VERBOSE = 2

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'LD.pickle'

    if not os.path.isfile(fn_data):
        from hivwholeseq.analysis.recombination.LD_simple import (
            collect_data)

        if VERBOSE >= 1:
            print 'Collect data'
        data = collect_data(pnames, fragments=fragments, VERBOSE=VERBOSE)

        if VERBOSE >= 1:
            print 'Save data for plot to pickle'
        store_data(data, fn_data)

    else:
        data = pd.read_pickle(fn_data)

    filename = foldername+'LD_Dp'
    for ext in ['png', 'pdf', 'svg']:
        plot_LD(data,
                VERBOSE=VERBOSE,
                savefig=filename+'.'+ext)

    plot_LD(data, VERBOSE=VERBOSE)

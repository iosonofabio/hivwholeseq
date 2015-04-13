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
pnames = ['20097', '15363', '15823', '15376', '20529', '9669', '15241', '15319']
regions = ["LTR5'", 'p17', 'p24', 'PR', 'RT', 'IN', 'vif', 'vpr', 'vpu', 'V3',
           'RRE', 'gp41', 'nef',
           "LTR3'"]



# Functions



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

    # Plot just the slope
    filename = foldername+'substitution_rates_sliding'
    #for ext in ['png', 'pdf', 'svg']:
    #    plot_substitution_rate_sliding(datap,
    #                                   regions,
    #                                   VERBOSE=VERBOSE,
    #                                   savefig=filename+'.'+ext)

    plot_substitution_rate_sliding(datap, regions, VERBOSE=VERBOSE)

    plt.ion()
    plt.show()


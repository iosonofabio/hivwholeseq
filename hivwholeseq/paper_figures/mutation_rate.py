# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/02/15
content:    Estimate the mutation rate matrix of HIV in vivo.
'''
# Modules
import os
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient

from hivwholeseq.analysis.mutation_rate.explore_divergence_synonymous import (
    collect_data_mutation_rate, fit_mutation_rate)
from hivwholeseq.analysis.mutation_rate.comparison_Abram import comparison_Abram2010

from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.paper_figures.plots import plot_mutation_rate




# Globals
pnames = ['20097', '15363', '15823', '15376', '20529', '9669', '15241', '15319']
regions = ['p17', 'p24', 'PR', 'RT', 'vif', 'vpu', 'nef']



# Script
if __name__ == '__main__':

    VERBOSE = 2
    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_fits = fn_data + 'mutation_rate.pickle'
    fn_comp = fn_data + 'mutation_rate_comparison_Abram2010.pickle'
    fn_data = fn_data + 'mutation_rate_data.pickle'

    if (not os.path.isfile(fn_data)) or (not os.path.isfile(fn_fits)):
        data = collect_data_mutation_rate(regions, pnames, VERBOSE=VERBOSE)
        fits = fit_mutation_rate(data, VERBOSE=VERBOSE, Smin=0.1)
        comp = comparison_Abram2010(fits)


        if VERBOSE >= 1:
            print 'Save data for plot to pickle'
        data.to_pickle(fn_data)
        fits.to_pickle(fn_fits)
        comp.to_pickle(fn_comp)

    else:
        data = pd.read_pickle(fn_data)
        fits = pd.read_pickle(fn_fits)
        comp = pd.read_pickle(fn_comp)

        
    datap = {'data': data, 'fits': fits, 'comp': comp}

    filenames = [foldername+'mutation_rate', foldername+'mutation_rate_fits']
    for ext in ['png', 'pdf', 'svg']:
        plot_mutation_rate(datap,
                           VERBOSE=VERBOSE,
                           savefig=[fn+'.'+ext for fn in filenames])

    plot_mutation_rate(datap, VERBOSE=VERBOSE)

    plt.ion()
    plt.show()


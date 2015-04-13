# vim: fdm=indent
'''
author:     Fabio Zanini
date:       08/04/15
content:    Plot divergence in terms of consensus and population.
'''
# Modules
import os
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.utils.generic import mkdirs

from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.paper_figures.plots import plot_divergence_cons_pop_diversity


# Globals
pnames = ['20097', '15363', '15823', '9669', '20529', '15241', '15376', '15319']
regions = ['p17', 'p24', 'p6', 'p7', 'PR', 'RT', 'p15', 'IN', 'vif', 'vpu', 'vpr', 'nef']



# Script
if __name__ == '__main__':

    VERBOSE = 2
    window_size = 300

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'divergence_cons_pop_diversity.pickle'

    if not os.path.isfile(fn_data):
        from hivwholeseq.analysis.divdiv.divergence_diversity import (
            collect_data)

        if VERBOSE >= 1:
            print 'Collect data'
        data = collect_data(pnames, regions, VERBOSE=VERBOSE)

        if VERBOSE >= 1:
            print 'Save data for plot to pickle'
        data.to_pickle(fn_data)

    else:
        data = pd.read_pickle(fn_data)

    # Plot
    filename = foldername+'divergence_consensus_population_diversity'
    for ext in ['png', 'pdf', 'svg']:
        plot_divergence_cons_pop_diversity(data,
                                           VERBOSE=VERBOSE,
                                           savefig=filename+'.'+ext)

    plot_divergence_cons_pop_diversity(data, VERBOSE=VERBOSE)



    plt.ion()
    plt.show()


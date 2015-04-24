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
            collect_data, get_sfs)

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

        if VERBOSE >= 1:
            print 'Save data for plot to pickle'
        datap.to_pickle(fn_data)

    else:
        datap = pd.read_pickle(fn_data)

    #filename = foldername+'sfs_syn_nonsyn'
    #for ext in ['png', 'pdf', 'svg']:
    #    plot_sfs_syn_nonsyn(datap,
    #                        regions,
    #                        VERBOSE=VERBOSE,
    #                        savefig=filename+'.'+ext)

    plot_sfs_syn_nonsyn(datap, VERBOSE=VERBOSE)
    plt.ion()
    plt.show()

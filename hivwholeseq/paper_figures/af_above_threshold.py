# vim: fdm=indent
'''
author:     Fabio Zanini
date:       24/03/15
content:    Make a coverage plot for the methods.
'''
# Modules
import os
import numpy as np
import pandas as pd

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.paper_figures.plots import plot_af_above_threshold
from hivwholeseq.utils.sequence import find_annotation



# Globals
regions = ['p17', 'p24', 'PR', 'RT', 'p15', 'IN', 'gp41']
pnames = ['20097', '15363', '15823', '9669', '20529', '15241', '15376', '15319']



# Functions



# Script
if __name__ == '__main__':

    VERBOSE = 2
    threshold = 0.01

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'af_above_threshold.pickle'

    if not os.path.isfile(fn_data):
        from hivwholeseq.analysis.exploration.af_above_threshold_entropy \
                import collect_data, calculate_fraction_above_threshold

        data = collect_data(pnames, regions, VERBOSE=VERBOSE)

        datap = calculate_fraction_above_threshold(data)

        datap.to_pickle(fn_data)

    else:
        datap = pd.read_pickle(fn_data)
        
    filename = (foldername+'fraction_alleles_above_'+
                str(threshold).replace('.', '')+
                '_entropy')
    for ext in ['png', 'pdf', 'svg']:
        plot_af_above_threshold(datap,
                                VERBOSE=VERBOSE,
                                savefig=filename+'.'+ext)

    plot_af_above_threshold(datap,
                            VERBOSE=VERBOSE)


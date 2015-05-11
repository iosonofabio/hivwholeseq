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
from hivwholeseq.paper_figures.plots import plot_entropy_correlation
from hivwholeseq.utils.sequence import find_annotation



# Globals
# NOTE: we exclude env because of the whole antubody thing
# The correlation in gp41 is lower and saturates fast
regions = ['p17', 'p24', 'PR', 'RT', 'p15', 'IN', 'vif', 'vpu', 'nef']



# Functions



# Script
if __name__ == '__main__':

    VERBOSE = 2

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'entropy_correlation_subtype_patient.pickle'

    if not os.path.isfile(fn_data):
        from hivwholeseq.analysis.exploration.correlate_diversity_patient_subtype \
                import collect_data, calculate_correlation

        data = collect_data(None, regions, VERBOSE=VERBOSE)

        datap = calculate_correlation(data)

        datap.to_pickle(fn_data)

    else:
        datap = pd.read_pickle(fn_data)
        
    filename = foldername+'entropy_correlation_patient_subtype'
    for ext in ['png', 'pdf', 'svg']:
        plot_entropy_correlation(datap,
                                 VERBOSE=VERBOSE,
                                 savefig=filename+'.'+ext)

    plot_entropy_correlation(datap,
                             VERBOSE=VERBOSE)



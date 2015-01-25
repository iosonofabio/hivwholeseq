# vim: fdm=marker
'''
author:     Fabio Zanini
date:       10/12/14
content:    Plot allele frequencies in overlap of consecutive fragments.
'''
# Modules
import os
import numpy as np

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.samples import load_samples_sequenced
from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.controls.check_allele_frequency_overlap import get_allele_frequency_overlap

from hivwholeseq.paper_figures.plots import plot_allele_frequency_overlap



# Globals
pnames = ['15319']
overlaps = ['F1-2', 'F2-3', 'F3-4', 'F4-5', 'F5-6']
cov_min = 1000
qual_min = 30




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
    foldername = get_figure_folder(username, 'controls')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'allele_frequency_overlap.pickle'

    if not os.path.isfile(fn_data):
        samples = load_samples_sequenced()
        if pnames is not None:
            samples = samples.loc[samples.patient.isin(pnames)]
        data = get_allele_frequency_overlap(samples, overlaps, cov_min=cov_min,
                                            VERBOSE=VERBOSE, qual_min=qual_min)


        #data = compress_data(data)
        store_data(data, fn_data)
    else:
        data = load_data(fn_data)
        
    filename = foldername+'allele_frequency_overlap'
    for ext in ['png', 'pdf', 'svg']:
        plot_allele_frequency_overlap(data,
                                      VERBOSE=VERBOSE,
                                      savefig=filename+'.'+ext)

    plot_allele_frequency_overlap(data, VERBOSE=VERBOSE)

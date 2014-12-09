# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/12/14
content:    Plot the quality along the read pair.
'''
# Modules
import os

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.patients import load_patients
from hivwholeseq.paper_figures.filenames import get_figure_folder

from hivwholeseq.patients.get_template_number import get_template_numbers
from hivwholeseq.paper_figures.plots import plot_template_numbers


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
    fn_data = fn_data + 'n_templates.pickle'

    if not os.path.isfile(fn_data):

        patients = load_patients()

        data = get_template_numbers(patients)

        store_data(data, fn_data)
    else:
        data = load_data(fn_data)
        
    filename = foldername+'n_templates'
    for ext in ['png', 'pdf', 'svg']:
        plot_template_numbers(data,
                              VERBOSE=VERBOSE,
                              savefig=filename+'.'+ext)

    plot_template_numbers(data,
                          VERBOSE=VERBOSE)

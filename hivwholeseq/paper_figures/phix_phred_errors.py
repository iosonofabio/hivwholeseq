# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/12/14
content:    Plot the quality along the read pair.
'''
# Modules
import os

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.sequencing.samples import load_sequencing_run
from hivwholeseq.paper_figures.filenames import get_figure_folder

from hivwholeseq.phix.phred_errors_phix import phred_errors_position
from hivwholeseq.paper_figures.plots import plot_phred_errors


# Functions
def compress_data(counts, seq_run):
    '''Compress data for plots, discarding useless info'''
    import numpy as np
    from scipy.stats import percentileofscore as pof

    data = []
    datum = {'counts': counts,
             'axes': ('read1/2', 'position in read', 'phred', 'correct/error'),
             'label': seq_run}
    data.append(datum)

    return data


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
    fn_data = fn_data + 'phix_phred.pickle'

    if not os.path.isfile(fn_data):
        seq_run = 'Tue42'
        maxreads = 10000
        # Part of PhiX looks weird, they gotta be illumina controls
        posrange = [2300, 5000]

        dataset = load_sequencing_run(seq_run)
        data_folder = dataset['folder']
        n_cycles = dataset['cycles']

        counts = phred_errors_position(data_folder, n_cycles,
                                       posrange=posrange,
                                       maxreads=maxreads, VERBOSE=VERBOSE,
                                       use_consensus=False)

        data = compress_data(counts, seq_run)
        store_data(data, fn_data)
    else:
        data = load_data(fn_data)

    filename = foldername+'phix_phred_errors'
    for ext in ['png', 'pdf', 'svg']:
        plot_phred_errors(data,
                          VERBOSE=VERBOSE,
                          savefig=filename+'.'+ext)


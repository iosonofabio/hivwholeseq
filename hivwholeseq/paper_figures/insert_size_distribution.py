# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/12/14
content:    Plot the quality along the read pair.
'''
# Modules
import os
import numpy as np

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.sequencing.samples import load_sequencing_run, SampleSeq
from hivwholeseq.paper_figures.filenames import get_figure_folder

from hivwholeseq.paper_figures.plots import plot_insert_size_distribution


# Functions
def compress_data(datar):
    '''Compress data for plots, discarding useless info'''
    data = []
    for (h, label) in datar:
        datum = {'bins': h[1],
                 'counts': h[0],
                 'label': label,
                }
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
    fn_data = fn_data + 'insert_size_distribution.pickle'

    if not os.path.isfile(fn_data):
        data = []
        bins = np.linspace(10, 1000, 50)

        # The original Nextera XT
        seq_run = 'LinaNexteraXT'
        adaID = '01'
        label = 'Standard Nextera XT'
        dataset = load_sequencing_run(seq_run)
        sample = SampleSeq(dataset.samples.iloc[0])
        _, h = sample.get_insert_size_distribution('premapped', VERBOSE=VERBOSE,
                                                   bins=bins,
                                                   density=False)
        data.append((h, label))

        # A good sample
        seq_run = 'Tue48'
        adaID = 'N2-S2'
        label = 'Optimized'
        dataset = load_sequencing_run(seq_run)
        sample = SampleSeq(dataset.samples.loc[dataset.samples.adapter == adaID].iloc[0])
        _, h = sample.get_insert_size_distribution('premapped', VERBOSE=VERBOSE,
                                                   bins=bins,
                                                   density=False)
        data.append((h, label))

        data = compress_data(data)
        store_data(data, fn_data)
    else:
        data = load_data(fn_data)
        
    filename = foldername+'insert_size_distribution'
    for ext in ['png', 'pdf', 'svg']:
        plot_insert_size_distribution(data,
                                      VERBOSE=VERBOSE,
                                      savefig=filename+'.'+ext)

    plot_insert_size_distribution(data,
                                  VERBOSE=VERBOSE)

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
from hivwholeseq.sequencing.filenames import get_read_filenames
from hivwholeseq.sequencing.check_quality_along_read import quality_score_along_reads
from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.paper_figures.plots import plot_cuts_quality_along_reads


# Functions
def compress_data(quality, thresholds=(10, 20, 30, 35)):
    '''Compress data for plots, discarding useless info'''
    import numpy as np
    from scipy.stats import percentileofscore as pof

    data = [[], []]
    for i, qual in enumerate(quality):
        cycles = len(qual)
        for qthresh in thresholds:
            x = np.arange(len(qual))
            y = np.array([100 - pof(qual[k], qthresh) for k in xrange(cycles)])
            data[i].append({'threshold': qthresh,
                            'x': x,
                            'y': y})
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
    fn_data = fn_data + 'quality_along_read.pickle'

    if not os.path.isfile(fn_data):
        seq_run = 'Tuen47'
        adaID = 'N5-S2'
        maxreads = 20000

        dataset = load_sequencing_run(seq_run)
        data_folder = dataset.folder
        read_len = dataset.cycles // 2
        reads_filenames = get_read_filenames(data_folder, adaID, gzip=True)

        quality = quality_score_along_reads(read_len, reads_filenames,
                                            randomreads=(maxreads >= 1),
                                            maxreads=maxreads, VERBOSE=VERBOSE)

        data = compress_data(quality)
        store_data(data, fn_data)
    else:
        data = load_data(fn_data)

    filename = foldername+'quality_along_read'
    for ext in ['png', 'pdf', 'svg']:
        plot_cuts_quality_along_reads(data,
                                      VERBOSE=VERBOSE,
                                      savefig=filename+'.'+ext)


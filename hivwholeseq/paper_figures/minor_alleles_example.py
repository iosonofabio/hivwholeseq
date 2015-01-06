# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/12/14
content:    Plot the quality along the read pair.
'''
# Modules
import os

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.samples import load_sample_sequenced
from hivwholeseq.paper_figures.filenames import get_figure_folder

from hivwholeseq.paper_figures.plots import plot_minor_allele_example


# Functions
def compress_data(counts, samplename, fragment):
    '''Compress data for plots, discarding useless info'''
    def get_minor_freqs(counts, cov_min=1000):
        '''Get minor freqs from counts'''
        import numpy as np
        af = np.zeros(counts.shape[-1])
        for pos, count in enumerate(counts.T):
            cov = count.sum()
            if cov >= cov_min:
                af_pos = 1.0 * count / cov
                af[pos] = np.sort(af_pos)[-2]
        return af

    data = []
    datum = {'freq_minor': get_minor_freqs(counts),
             'samplename': samplename,
             'fragment': fragment}
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

    # Patient data
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'minor_alleles_example.pickle'

    # Control
    foldername_ref = get_figure_folder(username, 'controls')
    fn_data_ref = foldername_ref+'data/'
    fn_data_ref = fn_data_ref + 'ref_minor_alleles.pickle'

    if not os.path.isfile(fn_data):
        samplename = '27134'
        fragment = 'F1'

        sample = load_sample_sequenced(samplename)

        counts = sample.get_allele_counts(fragment, merge_read_types=True)

        data = compress_data(counts, samplename, fragment)
        store_data(data, fn_data)
    else:
        data = load_data(fn_data)

    # Add reference
    data = load_data(fn_data_ref) + data
        
    filename = foldername+'freq_minor_alleles_example'
    for ext in ['png', 'pdf', 'svg']:
        plot_minor_allele_example(data,
                                  VERBOSE=VERBOSE,
                                  savefig=filename+'.'+ext)

    plot_minor_allele_example(data,
                              VERBOSE=VERBOSE)

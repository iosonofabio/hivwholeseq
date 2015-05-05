# vim: fdm=marker
'''
author:     Fabio Zanini
date:       10/12/14
content:    Plot allele frequencies in overlap of consecutive fragments.
'''
# Modules
import os
import numpy as np

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.controls.check_allele_frequency_overlap import get_allele_frequency_overlap

from hivwholeseq.paper_figures.plots import plot_allele_frequency_overlap



# Globals
overlaps = ['F1-2', 'F2-3', 'F3-4', 'F4-5', 'F5-6']
cov_min = 1000
qual_min = 30




# Functions
def estimate_templates_overlaps(sample, data):
    '''Estimate templates for the overlaps'''
    for datum in data:
        af1, af2 = datum['af']

        # Filter only polymorphic sites
        afmin = 3e-3
        indfm = (af1 >= afmin) & (af1 <= 1 - afmin) & (af2 >= afmin) & (af2 <= 1 - afmin)
        nsites = len(np.unique(indfm.nonzero()[1]))

        # Estimate the template number
        mea = 0.5 * (af1[indfm] + af2[indfm])
        var = ((af1[indfm] - af2[indfm]) / 2)**2

        # In binomial sampling, the variance on k is var(k) = nx (1 - x), so
        # for the frequency var(k/n) = x (1 - x) / n
        n_all = mea * (1 - mea) / var
        
        # NOTE: pseudocounts that come from the F4 dilution estimate, so we
        # only listen to the data if there is enough data points to listen to
        len_pseudo = 1
        n_pseudo = sample.get_n_templates_dilutions()
        n_allp = np.concatenate([n_all, ([n_pseudo] * len_pseudo)])

        if VERBOSE >= 2:
            print datum['overlap'], 'Number of doubly polymorphic sites:', nsites, 'n_pseudo:', n_pseudo

        # NOTE: the estimate of n has a bad distribution because some points are
        # exactly on the diagonal, so we average the inverse (which is well
        # behaved) and also take the medians as alternatives
        n = 1.0 / (1.0 / n_allp).mean()

        datum['n'] = n


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
    pname = 'p11'
    n_time = 4

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'controls')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'allele_frequency_overlap.pickle'

    if not os.path.isfile(fn_data):
        patient = load_patient(pname)
        samples = patient.samples.iloc[[n_time]]
        sample = SamplePat(samples.iloc[0])
        data = get_allele_frequency_overlap(samples, overlaps, cov_min=cov_min,
                                            VERBOSE=VERBOSE, qual_min=qual_min)

        estimate_templates_overlaps(sample, data)


        store_data(data, fn_data)
    else:
        data = load_data(fn_data)
        
    filename = foldername+'allele_frequency_overlap'
    for ext in ['png', 'pdf', 'svg']:
        plot_allele_frequency_overlap(data,
                                      VERBOSE=VERBOSE,
                                      savefig=filename+'.'+ext)

    plot_allele_frequency_overlap(data, VERBOSE=VERBOSE)

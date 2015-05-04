# vim: fdm=indent
'''
author:     Fabio Zanini
date:       24/03/15
content:    Make a coverage plot for the methods.
'''
# Modules
import os
import numpy as np

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.paper_figures.plots import plot_coverage_example
from hivwholeseq.utils.sequence import find_annotation



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
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'coverage_example.pickle'

    if not os.path.isfile(fn_data):
        pcode = 'p1'
        n_time = 3

        patient = load_patient(pcode)
        ref = patient.get_reference('genomewide', 'gb')

        sample = SamplePat(patient.samples.iloc[n_time])

        cov = {}
        for fragment in ['F'+str(i) for i in xrange(1, 7)]:
            pos = find_annotation(ref, fragment).location.nofuzzy_start
            co = sample.get_coverage(fragment)
            cov[fragment] = {'pos': pos, 'cov': co}
        co = sample.get_coverage('genomewide')
        cov['genomewide'] = {'pos': pos, 'cov': co}

        data = {'sample': sample.name,
                'n_templates': patient.n_templates[n_time],
                'cov': cov,
               }

        store_data(data, fn_data)
    else:
        data = load_data(fn_data)
        
    samplename = data['sample']
    filename = foldername+'_'.join(['coverage_example', samplename])
    for ext in ['png', 'pdf', 'svg']:
        plot_coverage_example(data,
                              VERBOSE=VERBOSE,
                              savefig=filename+'.'+ext)

    plot_coverage_example(data,
                          VERBOSE=VERBOSE)


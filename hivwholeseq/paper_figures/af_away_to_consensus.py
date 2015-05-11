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
from hivwholeseq.paper_figures.plots import plot_n_muts_awayto, plot_af_entropy_awayto
from hivwholeseq.utils.sequence import find_annotation



# Globals
regions = ['p17', 'p24', 'PR', 'RT', 'IN', 'vif', 'vpu', 'nef', 'gp41']



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
    fn_data = fn_data + 'n_muts_af_entropy_awayto.pickle'

    if not os.path.isfile(fn_data):
        from hivwholeseq.analysis.to_from_subtype_consensus.number_mutations_time_entropy_af \
                import collect_data, bin_data, get_n_mutations_patientwise, \
                average_n_muts_patients, get_allele_frequency_entropy, \
                get_fraction_to_baseline

        data = collect_data(None, regions, VERBOSE=VERBOSE)
        bins = bin_data(data, ['time', 'time_rough', 'entropy', 'af'], VERBOSE=VERBOSE)
        datap = get_n_mutations_patientwise(data, attrnames=['time', 'awayto'])

        n_muts = average_n_muts_patients_total(datap, n_bootstrap=100, VERBOSE=3)
        frac = average_n_muts_patients_fraction(datap, n_bootstrap=100, VERBOSE=3)
        af_avg = get_allele_frequency_entropy(data, n_bootstrap=100, VERBOSE=3)
        frac0 = get_fraction_to_baseline(data, bins)

        datap = {'bins': bins,
                 'n_muts': n_muts,
                 'fraction': frac,
                 'fraction0': frac0,
                 'af_avg': af_avg}

        store_data(datap, fn_data)


    else:
        datap = load_data(fn_data)
        
    # Plot average number of mutations
    filename = foldername+'n_mutations_awayto_ratio'
    for ext in ['png', 'pdf', 'svg']:
        plot_n_muts_awayto({'bins': datap['bins'],
                            'n_muts': datap['n_muts'],
                            'fraction': datap['fraction'],
                            'fraction0': datap['fraction0'],
                           },
                           VERBOSE=VERBOSE,
                           savefig=filename+'.'+ext)

    plot_n_muts_awayto({'bins': datap['bins'],
                        'n_muts': datap['n_muts'],
                        'fraction': datap['fraction'],
                        'fraction0': datap['fraction0'],
                       },
                       VERBOSE=VERBOSE)

    # Plot allele frequency by entropy
    filename = foldername+'allele_freq_avg_away_to'
    for ext in ['png', 'pdf', 'svg']:
        plot_af_entropy_awayto({'bins': datap['bins'],
                                'af_avg': datap['af_avg']},
                               savefig=filename+'.'+ext)

    plot_af_entropy_awayto({'bins': datap['bins'],
                            'af_avg': datap['af_avg']},
                           VERBOSE=VERBOSE)


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
from hivwholeseq.sequencing.samples import load_sample_sequenced
from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.reference_samples.mix_recombination_RNA import (
   get_samplename, get_alignment, get_good_polymorphic_sites, get_cocounts,
   get_switches_conditional)

from hivwholeseq.paper_figures.plots import plot_rna_recombination



# Globals
strains = ['LAI-III', '38540']
fragment = 'F2'
maxreads = 10000
freqmin = 6e-3



# Functions
def compress_data(datar):
    '''Compress data for plots, discarding useless info'''
    data = []
    for (al_polyd, counts, cocounts, switchesn, refnames, PCRtype) in datar:
        datum = {'aldict': al_polyd,
                 'counts': dict(counts),
                 'cocounts': dict(cocounts),
                 'switchesn': dict(switchesn),
                 'refnames': refnames,
                 'PCRtype': PCRtype,
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
    fn_data = fn_data + 'mix_rna_recombination.pickle'

    if not os.path.isfile(fn_data):

        PCRtype = 'PCR2'
        samplename = get_samplename(PCRtype)
        sample = load_sample_sequenced(samplename)
        samplerefs = [load_sample_sequenced(strain) for strain in strains]
        ali = get_alignment((sample, samplerefs[0], samplerefs[1]),
                            fragment,
                            VERBOSE=VERBOSE)
        alim = np.array(ali)
        al_polyd = get_good_polymorphic_sites(alim, samplerefs, fragment, VERBOSE=VERBOSE)
        bamfilename = sample.get_mapped_filename(fragment, filtered=True)
        counts, cocounts = get_cocounts(bamfilename, np.sort(al_polyd.keys()),
                                        maxreads=maxreads,
                                        VERBOSE=VERBOSE)
        switchesn = get_switches_conditional(cocounts, al_polyd, VERBOSE=VERBOSE,
                                             freqmin=freqmin)

        data = compress_data([(al_polyd, counts, cocounts, switchesn, strains,
                               PCRtype)])
        store_data(data, fn_data)
    else:
        data = load_data(fn_data)
        
    filename = foldername+'mix_rna_recombination'
    for ext in ['png', 'pdf', 'svg']:
        plot_rna_recombination(data,
                               VERBOSE=VERBOSE,
                               savefig=filename+'.'+ext)

    plot_rna_recombination(data,
                           VERBOSE=VERBOSE)

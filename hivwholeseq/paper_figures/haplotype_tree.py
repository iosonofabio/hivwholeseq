# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/12/14
content:    Plot the quality along the read pair.
'''
# Modules
import os
import numpy as np

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.utils.sequence import align_muscle
from hivwholeseq.utils.tree import build_tree_fasttree, filter_rare_leaves
from hivwholeseq.paper_figures.plots import plot_haplotype_tree_example



# Functions
def compress_data(tree, pname, region):
    '''Compress data for plots, discarding useless info'''
    data = []
    datum = {'tree': tree,
             'pname': pname,
             'region': region}
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
    regions = ['V3', 'RT3']

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fd_data = foldername+'data/'
    mkdirs(fd_data)

    for region in regions:
        print region
        fn_data = fd_data + 'haplotype_tree_example_'+region+'.pickle'

        if not os.path.isfile(fn_data):
            pname = '15823'
            cutoff = 0.04

            patient = load_patient(pname)
            tree = patient.get_local_tree(region)

            filter_rare_leaves(tree, cutoff)

            data = compress_data(tree, pname, region)
            store_data(data, fn_data)
        else:
            data = load_data(fn_data)
            
        filename = foldername+'haplotype_tree_example_'+region
        for ext in ['png', 'pdf', 'svg']:
            plot_haplotype_tree_example(data,
                                        VERBOSE=VERBOSE,
                                        savefig=filename+'.'+ext)

        plot_haplotype_tree_example(data,
                                    VERBOSE=VERBOSE)

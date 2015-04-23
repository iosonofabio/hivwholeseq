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
def compress_data(tree, pcode, region, title):
    '''Compress data for plots, discarding useless info'''
    data = []
    datum = {'tree': tree,
             'pcode': pcode,
             'region': region,
             'title': title,
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
    regions = [{'name': 'p17', 'title': 'p17'},
               {'name': 'RT1', 'title': 'RT start'},
               {'name': 'IN2', 'title': 'integrase middle'},
               {'name': 'V3', 'title': 'V3'},
              ]
    pcode = 'p3'
    cutoff = 0.02

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fd_data = foldername+'data/'
    mkdirs(fd_data)

    fn_data = fd_data + 'haplotype_tree_example_'+pcode+'.pickle'
    if not os.path.isfile(fn_data):
        patient = load_patient(pcode)
        data = []
        for regdata in regions:
            region = regdata['name']
            title = regdata['title']
            print region, title
            tree = patient.get_local_tree(region)
            filter_rare_leaves(tree, cutoff)
            data.extend(compress_data(tree, patient.code, region, title))
        store_data(data, fn_data)
    else:
        data = load_data(fn_data)
            

    #filename = foldername+'_'.join(['haplotype_tree_example', pcode])
    #for ext in ['png', 'pdf', 'svg']:
    #    plot_haplotype_tree_example(data,
    #                                VERBOSE=VERBOSE,
    #                                legend='V3',
    #                                savefig=filename+'.'+ext)

    plot_haplotype_tree_example(data,
                                legend='V3',
                                VERBOSE=VERBOSE)


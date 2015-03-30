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
from hivwholeseq.utils.sequence import align_muscle
from hivwholeseq.utils.tree import build_tree_fasttree, filter_rare_leaves
from hivwholeseq.paper_figures.plots import plot_haplotype_tree_example

from hivwholeseq.website.filenames import get_tree_figure_filename


# Globals
pnames = ['20097', '15363', '15823', '15376', '9669', '15107', '15241', '15034', '15319']
region = 'V3'
cutoff = 0.04




# Functions
def compress_data(tree, pname, region):
    '''Compress data for plots, discarding useless info'''
    data = []
    datum = {'tree': tree,
             'pname': pname,
             'region': region}
    data.append(datum)

    return data



# Script
if __name__ == '__main__':

    VERBOSE = 2

    for pname in pnames:

        patient = load_patient(pname)
        tree = patient.get_local_tree(region)

        filter_rare_leaves(tree, cutoff)

        data = compress_data(tree, patient.code, region)
            
        filename = get_tree_figure_filename(patient.code, region, format='svg')

        plot_haplotype_tree_example(data,
                                    VERBOSE=VERBOSE,
                                    savefig=filename)

        plot_haplotype_tree_example(data,
                                    VERBOSE=VERBOSE,
                                    )


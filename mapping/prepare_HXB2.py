# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/08/13
content:    Chop HXB2 into the 6 fragments, which are used as "chromosomes" for
            mapping.
'''
# Modules
import os
import sys
import re
from itertools import izip
import numpy as np
import Bio.SeqIO as SeqIO

from mapping.reference import load_HXB2
from mapping.primer_info import primers_coordinates_HXB2_inner as pci
from mapping.filenames import get_HXB2_fragmented, get_HXB2_entire



# Globals
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']



# Script
if __name__ == '__main__':

    # Get the annotated sequence
    seq = load_HXB2(data_folder)

    # Extract fragments coordinates and sequences
    frag_coos = {f: [co[0][0], co[1][1]] for f, co in pci.iteritems()}
    frag_seqs = {f: seq[co[0]: co[1]] for f, co in frag_coos.iteritems()}

    # Relabel sequences
    for (f, fs) in frag_seqs.iteritems():
        fs.name = 'HXB2_reference_'+f
        fs.id = fs.name
        fs.description = re.sub('_', ' ', fs.id)

    # Make a cropped sequence from F1 to F6, to reduce LTR degeneracy problems
    start_F1 = frag_coos['F1'][0]
    end_F6 = frag_coos['F6'][1]
    seq_cropped = seq[start_F1: end_F6]
    seq_cropped.id = seq_cropped.name = seq.name+'_cropped_F1_F6'
    seq_cropped.description = seq.description+' (cropped to fragments F1-F6)'

    # Make output folder if necessary
    dirname = os.path.dirname(get_HXB2_entire(data_folder))
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # Write to file both the entire, the cropped, and the fragmented (in separated files)
    SeqIO.write(seq, get_HXB2_entire(data_folder), 'fasta')
    SeqIO.write(seq_cropped, get_HXB2_entire(data_folder, cropped=True), 'fasta')
    for (f, fs) in frag_seqs.iteritems():
        SeqIO.write(fs, get_HXB2_fragmented(data_folder, f), 'fasta')
        print f

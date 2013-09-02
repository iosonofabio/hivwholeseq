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

from mapping.sequence_utils.annotate_HXB2 import load_HXB2
from mapping.filenames import get_HXB2_fragmented, get_HXB2_entire



# Globals
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']



# Script
if __name__ == '__main__':

    # Get the annotated sequence
    seq = load_HXB2()

    # Extract fragments
    frags = [s for s in seq.features if s.type[:8] == 'fragment']
    frag_seqs = [s.extract(seq) for s in frags]

    # Relabel sequences
    for (fs, f) in izip(frag_seqs, frags):
        fs.name = re.sub(' ', '_', f.type)+' '+fs.name
        fs.id = fs.name
        fs.description = re.sub('complete genome', f.type, fs.description).split(';')[0]

    # Make a cropped sequence from F1 to F6, to reduce LTR degeneracy problems
    start_F1 = [f for f in seq.features
                if f.type == 'fragment F1'][0].location.nofuzzy_start
    end_F6 = [f for f in seq.features
              if f.type == 'fragment F6'][0].location.nofuzzy_end
    seq_cropped = seq[start_F1: end_F6]
    seq_cropped.id = seq_cropped.name = seq.name+'_cropped_F1_F6'
    seq_cropped.description = seq.description+' (cropped to fragments F1-F6)'

    # Write to file both the entire, the cropped, and the fragmented (in separated files)
    SeqIO.write(seq, get_HXB2_entire(data_folder), 'fasta')
    SeqIO.write(seq_cropped, get_HXB2_entire(data_folder, cropped=True), 'fasta')
    for frag_seq, fragment in izip(frag_seqs, ('F'+str(i) for i in xrange(1, 7))):
        SeqIO.write(frag_seq, get_HXB2_fragmented(data_folder, fragment), 'fasta')


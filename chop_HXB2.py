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


# Env var
module_group = '/ebio/ag-neher/share/programs/modules'
if module_group not in sys.path:
    sys.path.append(module_group)
from annotate_HXB2 import load_HXB2



# Globals
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'
data_folder = 'data/'



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

    # Write to file
    SeqIO.write(frag_seqs,
                data_folder+'HXB2_fragmented.fasta',
                'fasta')


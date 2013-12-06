# vim: fdm=indent
'''
author:     Fabio Zanini
date:       20/08/13
content:    Details about the MiSeq, used in the downstream analysis scripts.
'''
# Modules
from numpy import array


# Globals
# Alphabet of output nucleotides
alphas = 'ACGT-N'
alphal = list(alphas)
alpha = array(alphal, 'S1')

# Read types
read_types = ['read1 f', 'read1 r', 'read2 f', 'read2 r']
read_pair_types = ['read1 f read2 r', 'read1 r read2 f']

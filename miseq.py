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


# Generator that gives pairs of reads at a time (if the total number of reads is
# odd, skip the last one)
def pair_generator(iterable):
    '''Generator for pairs (the last item is lost if odd)'''
    it = iter(iterable)
    while True:
        try:
            a = it.next()
            b = it.next()
            yield (a, b)
        except StopIteration:
            raise

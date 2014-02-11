# vim: fdm=marker
'''
author:     Fabio Zanini
date:       04/02/14
content:    Alphabet and other settings of the PacBio RS II instrument.
'''
# Modules
from numpy import array


# Globals
# Alphabet of output nucleotides
alphas = 'ACGT-N'
alphal = list(alphas)
alpha = array(alphal, 'S1')

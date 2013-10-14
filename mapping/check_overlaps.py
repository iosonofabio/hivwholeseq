# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/10/13
content:    Check overlapping regions for PCR amplification biases.

            The simplest way of doing this is via allele frequencies. More
            thorough analyses on the reads could be performed.
'''
# Modules

from mapping.check_consensus import find_overlap
from mapping.mapping_utils import align_muscle

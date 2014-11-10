# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/11/14
content:    Store coverage for the website.
'''
# Modules
import os
import sys
import shutil

from hivwholeseq.website.filenames import get_propagator_filename as get_fn_out
from hivwholeseq.patients.filenames import get_propagator_filename



# Script
if __name__ == '__main__':

    fragments = ['F'+str(i) for i in xrange(1, 7)]

    # TODO: generalize this :-)
    dts = [[500, 1000]]

    for dt in dts:
        # NOTE: everything should be anonymized already
        fn_in = get_propagator_filename(['all'], fragments, dt)
        fn_out = get_fn_out(['allpats'], fragments, dt)
        shutil.copy(fn_in, fn_out)

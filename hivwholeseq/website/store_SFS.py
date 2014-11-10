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

from hivwholeseq.website.filenames import get_SFS_filename as get_fn_out
from hivwholeseq.patients.filenames import get_SFS_filename



# Script
if __name__ == '__main__':

    fragments = ['F'+str(i) for i in xrange(1, 7)]

    # NOTE: everything should be anonymized for SFS already
    fn_in = get_SFS_filename(['all'], fragments)
    fn_out = get_fn_out(['allpats'], fragments)
    shutil.copy(fn_in, fn_out)

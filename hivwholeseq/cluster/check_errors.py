# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/12/14
content:    Check errors on the cluster output (in files).
'''
# Modules
import os, glob, re

from hivwholeseq.cluster import JOBLOGERR, JOBLOGOUT
from hivwholeseq.cluster.fork_cluster import nothing



# Functions
def print_stdout(fn):
    '''Print stdout of filename'''

    # Accepts index or filename (we are sloppy, we should use a closure)
    if isinstance(fn, int):
        fn = fn_outs[fn]

    with open(JOBLOGOUT+fn, 'r') as f:
        print f.read()


def print_stderr(fn):
    '''Print stdout of filename'''

    # Accepts index or filename (we are sloppy, we should use a closure)
    if isinstance(fn, int):
        fn = fn_errs[fn]

    with open(JOBLOGERR+fn, 'r') as f:
        print f.read()





# Script
if __name__ == '__main__':


    fn_errs = []
    for fn in os.listdir(JOBLOGERR):
        # Skip hiden files (including .gitignore)
        if fn[0] == '.':
            continue

        with open(JOBLOGERR+fn, 'r') as f:
            for line in f:
                line = line.rstrip('\n').strip()
                if '.Xauthority' in line:
                    continue

                if line[:7] == 'stampy:':
                    continue

                if len(line):
                    fn_errs.append(fn)
                    break

    print 'Filenames with errors:'
    print '\n'.join(fn_errs)

    fn_outs = [re.sub(r'\.e', r'.o', fn) for fn in fn_errs]


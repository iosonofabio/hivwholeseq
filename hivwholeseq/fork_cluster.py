# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/12/13
content:    Module with all submit functions for the cluster. With this we can
            keep all cluster-specific code in one place.
'''
# Globals
import subprocess as sp

import hivwholeseq
JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
JOBLOGOUT = JOBDIR+'logout'
JOBLOGERR = JOBDIR+'logerr'



# Functions
def empty_log_folders():
    '''Empty log folders of old files'''
    import shutil, os

    shutil.rmtree(JOBLOGOUT)
    os.mkdir(JOBLOGOUT)

    shutil.rmtree(JOBLOGERR)
    os.mkdir(JOBLOGERR)


def fork_demultiplex(miseq_run, VERBOSE=0, summary=True):
    '''Submit demultiplex script to the cluster'''
    JOBSCRIPT = JOBDIR+'demultiplex.py'
    cluster_time = '23:59:59'
    vmem = '8G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'demux HIV',
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', miseq_run,
                 '--verbose', VERBOSE,
                ]
    if not summary:
        call_list.append('--no-summary')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)

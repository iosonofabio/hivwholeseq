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


def fork_demultiplex(seq_run, VERBOSE=0, summary=True):
    '''Submit demultiplex script to the cluster'''
    if VERBOSE:
        print 'Forking to the cluster'

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
                 '--run', seq_run,
                 '--verbose', VERBOSE,
                ]
    if not summary:
        call_list.append('--no-summary')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_premap(seq_run, adaID, VERBOSE=0, bwa=False, threads=1, report=0,
                reference='HXB2', summary=True):
    '''Submit premap script to the cluster for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+'{:02d}'.format(adaID)

    JOBSCRIPT = JOBDIR+'premap_to_reference.py'
    # It is hard to tell whether 1h is sufficient, because the final sorting takes
    # quite some time. So for now give up and require 2h.
    cluster_time = ['23:59:59', '1:59:59']
    vmem = '8G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'premap '+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time[threads >= 30],
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--adaIDs', adaID,
                 '--verbose', VERBOSE,
                 '--threads', threads,
                 '--reference', reference,
                ]
    if report:
        call_list.append('--report')
    if not summary:
        call_list.append('--no-summary')
    call_list = map(str, call_list)
    if VERBOSE >= 2:
        print ' '.join(call_list)
    return sp.check_output(call_list)




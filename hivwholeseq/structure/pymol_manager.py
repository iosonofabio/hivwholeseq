#!/usr/bin/env python2
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       27/01/15
content:    Control the PyMol server from here.
'''
# Modules
import argparse



# Functions
def start_pymol(VERBOSE=0, gui=False):
    '''Start PyMol server as an RPCbind'''
    import subprocess as sp
    # Args: -R = server mode, -q = quiet, -c = command-line, -K = keep alive
    if gui:
        args = '-Rq'
    else:
        args = '-RqcK'

    call_list = ['pymol', args]
    if VERBOSE >= 1:
        print ' '.join(call_list)

    # Note: pymol starts children processes, so it's quite useless to return
    # the PID of this one, unless we find a way to surf the process hierarchy
    # (not eager to do that, honestly)
    # FIXME: make options for stdout, now it's polluting the whole thing...
    sp.Popen(' '.join(call_list), shell=True)


def check_pymol(VERBOSE=0):
    '''Check for a PyMol server'''
    import psutil
    pids = []
    for pid in psutil.pids():
        p = psutil.Process(pid)
        cmdline = ' '.join(p.cmdline())
        if ('/bin/pymol' in cmdline) or ('pymol/__init__.py' in cmdline):
            pids.append(pid)
    return pids


def stop_pymol(VERBOSE=0):
    '''Stop PyMol server'''
    import psutil
    pids = check_pymol(VERBOSE=VERBOSE)

    for pid in pids:
        if VERBOSE >= 1:
            print pid
        p = psutil.Process(pid)
        p.kill()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Manage the pymol RPC server',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('action',
                        choices=['start', 'startgui', 'check', 'stop'],
                        help='start|startgui|check|stop')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level')

    args = parser.parse_args()
    action = args.action
    VERBOSE = args.verbose

    if action == 'start':
        start_pymol(VERBOSE=VERBOSE, gui=False)

    elif action == 'startgui':
        start_pymol(VERBOSE=VERBOSE, gui=True)

    elif action == 'check':
        output = '\n'.join(map(str, check_pymol(VERBOSE=VERBOSE)))
        if output:
            print output

    else:
        stop_pymol(VERBOSE=VERBOSE)


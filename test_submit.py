#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       19/08/13
content:    Test for submit scripts (module import and such).

            Note: for this and similar scripts to work, two things are required:
            1. to call qsub with the -b y option
            2. to have the parent folder of this one in your PYTHONPATH. This
            usually involves putting something like the following in .bash_profile:

                export PYTHONPATH=$PYTHONPATH:"<parent_folder>"

            Note also that the script assumes that your folder is called "mapping".
            More elaborate schemes can be devised for this thing, but none of them
            is clean. The clean way is to make a pure Python module out of this
            folder and install it using distutils, but that is not handy for
            rapidly-changing projects.
'''
# Modules
import sys
import os
import subprocess as sp


# Environment vars
JOBSUBMIT = os.path.realpath(__file__)      # This file
JOBDIR = os.path.dirname(JOBSUBMIT)+'/'
JOBSCRIPT = JOBSUBMIT
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'


# Cluster parameters
cluster_time = '00:00:59'
vmem = '1G'
VERBOSE = 1



# Script
if __name__ == '__main__':

    # Input args
    args = sys.argv

    # No args --> submit script
    if len(args) < 2:

        # Submit self
        qsub_list = ['qsub','-cwd',
                     '-o',JOBLOGOUT,
                     '-e',JOBLOGERR,
                     '-l', 'h_rt='+cluster_time,
                     '-l', 'h_vmem='+vmem,
                     '-b', 'y',
                     JOBSCRIPT,
                     'run']
        qsub_list = map(str,qsub_list)
        if VERBOSE:
            print ' '.join(qsub_list)
        sp.call(qsub_list)

    # Args --> run script
    else:
        print os.getenv('PYTHONPATH')
        from mapping.test import identity
        print identity(5)

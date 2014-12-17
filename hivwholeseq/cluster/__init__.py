# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/12/14
content:    Cluster utilities.
'''
# Modules
import os.path


# Globals
JOBDIR = os.path.dirname(__path__[0].rstrip('/')).rstrip('/')+'/'
JOBLOGOUT = JOBDIR+'cluster/logout/'
JOBLOGERR = JOBDIR+'cluster/logerr/'

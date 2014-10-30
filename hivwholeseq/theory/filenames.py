# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/10/14
content:    Support module with standard filenames for theoretical calculations.
'''
# Modules
from hivwholeseq.filenames import theory_folder



# Functions
def get_sfs_betatree_filename(N, alpha):
    '''Get the filename of the precomputed SF from beta trees'''
    filename = 'betatree_sfs_N_'+str(N)+'_alpha_'+str(float(alpha))+'.dat'
    return theory_folder+filename

# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/15
content:    Support module for file names.
'''
# Modules
from hivwholeseq.sequencing.filenames import reference_folder
structure_folder = reference_folder+'structures/'



# Globals
structure_names = {'PR': '1HSG'}



# Functions
def get_PDB_filename(region, seqid=None, VERBOSE=0, gzip=False):
    '''Get the filename of a PDB structure'''
    if seqid is None:
        seqid = structure_names[region]
    fn = region+'_'+seqid+'.pdb'
    if gzip:
        fn = fn+'.gz'
    fn = structure_folder+fn
    return fn


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
structures = {'PR': '1HSG',
              'p17': '1TAM',
              'p24': 'hex_3MGE',
              'RT': '2HMI',
             }

chains = {'1HSG': [0, 1],
          '1TAM': [0],
          'hex_3MGE': range(6),
          '2HMI': [2, 3],
         }



# Functions
def get_PDB_filename_and_chains(region, seqid=None, VERBOSE=0, gzip=False):
    '''Get the filename of a PDB structure'''
    if seqid is None:
        seqid = structures[region]
    fn = region+'_'+seqid+'.pdb'
    if gzip:
        fn = fn+'.gz'
    fn = structure_folder+fn
    return fn, chains[seqid]


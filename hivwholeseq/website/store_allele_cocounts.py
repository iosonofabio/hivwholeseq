#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/11/14
content:    Store allele cocounts for the website. We need two formats for
            downloads, a dense one and a sparse one.
'''
# Modules
import os
import sys
import argparse
import numpy as np

from hivwholeseq.miseq import alphal
from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, iterpatient
from hivwholeseq.website.filenames import get_allele_cocounts_filename as get_fn_out
from hivwholeseq.website.filenames import get_coverage_filename
from hivwholeseq.cluster.fork_cluster import fork_store_cocounts_website as fork_self


# Globals
s = \
'''# The data is stored in the following format. The matrix of cocounts has originally
# the following dimensions: 6 x 6 x L x L where 6 is the length of the alphabet:
# ACGT-N and the last two dimensions refer to the positions of the two alleles.
# The matrix is flattened in this file with C-style order. To reconstruct the
# original matrix, reshaping is needed, e.g. in Python (numpy):
# arr = np.loadtxt(<filename>)
# L = int(np.sqrt(len(arr) // 36))
# M = arr.reshape((6, 6, L, L))
'''


# Functions
def save_for_download(filename, use_tempfile=True, **data):
    '''Save data for download (anotations, etc)
    
    Using a temporary file requires more space on disk but less RAM.
    '''
    from itertools import izip
    import StringIO
    import zipfile

    sep = '\t'
    fmt = lambda x: sep.join(x)+'\n'

    filename_uncompressed = filename[:-3]+'tsv'

    if not use_tempfile:
        with zipfile.ZipFile(filename, 'w',
                             compression=zipfile.ZIP_DEFLATED) as zf:

            # Descriptor
            f = StringIO.StringIO()
            f.write(s)
            f.write('# days since infection: '+str(data['days since infection'])+'\n')

            np.savetxt(f, data['cocounts'].ravel(), fmt='%d')

            zf.writestr(filename_uncompressed, f.getvalue())

    else:
        with open(filename_uncompressed, 'w') as f:
            f.write(s)
            f.write('# days since infection: '+str(data['days since infection'])+'\n')
            np.savetxt(f, data['cocounts'].ravel(), fmt='%d')

        import subprocess as sp
        sp.call(['zip', '-j9', filename, filename_uncompressed])
        os.remove(filename_uncompressed)




# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Store haplotypes and trees',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--samplenumbers', nargs='+', type=int,
                        default=None,
                        help='Sample numbers to store')
    parser.add_argument('--fragments', nargs='+',
                        default=['F'+str(j) for j in xrange(1, 7)],
                        help='Fragments to store (e.g. F1 F3)')
    parser.add_argument('--submit', action='store_true',
                        help='Submit as a parallel job to the cluster')

    args = parser.parse_args()
    pnames = args.patients
    samplenumbers = args.samplenumbers
    fragments = args.fragments
    submit = args.submit


    patients = load_patients(pnames=pnames)
    for pname, patient in iterpatient(patients):

        # Sample by sample
        for i, sample in enumerate(patient.itersamples()):
            samplename = patient.code+'_sample_'+str(i+1)

            if (samplenumbers is not None) and (i not in samplenumbers):
                continue

            for region in fragments:
                print patient.code, patient.name, i+1, region,

                fn_out = get_fn_out(samplename, region, format='zip', type='nuc')
                if os.path.isfile(sample.get_allele_cocounts_filename(region)):
                    print 'input file found'
                else:
                    print 'input file not found'
                    continue

                if submit:
                    fork_self(pname, i, region, VERBOSE=1)
                    continue
                
                cocounts = sample.get_allele_cocounts(region)
                data = {'days since infection': sample['days since infection'],
                        'cocounts': cocounts,
                       }

                save_for_download(fn_out, use_tempfile=True, **data)

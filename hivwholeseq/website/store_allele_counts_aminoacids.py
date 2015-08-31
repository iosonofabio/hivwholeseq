# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/11/14
content:    Store allele counts for the website. We need two formats, one binary
            fast-access for the plots, and one annotated and standard for
            downloads.

            We store both allele count trajectories for whole patients and sample
            by sample.
'''
# Modules
import os
import sys
import numpy as np

from hivwholeseq.miseq import alphal
from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.website.filenames import get_allele_count_trajectories_filename as get_fn_out
from hivwholeseq.website.filenames import get_allele_counts_filename as get_fn2_out
from hivwholeseq.website.filenames import get_coverage_filename
from hivwholeseq.patients.filenames import get_allele_count_trajectories_filename


# Globals
s = \
'''The data is stored in the following format. Each file is in Tab Separated
Value (TSV) format. Each file contains a separate time point. The estimated
number of days from infection is in the filename and again in the first line
of the file.

Each row is a position in the gene of interest, in amino acids. The coordinates
refer to the same gene in the reference sequence for each patient. The genome
reference can be downloaded in Genbank (GB) format from the website, and the
nucleotide sequence of the gene itself can be extracted from the annotations.

Each column is a possible residue, with a 22 letters alphabet ending in -N. The
- (dash) indicates deletions (with respect to the patient reference sequence),
N an ambiguous nucleotide. The latter indicates ambiguous residues, should be
rare and can be ignored for most purposes.

Insertions are not included.
'''

proteins = ['RT', 'PR', 'p15', 'IN', 'p17', 'p24', 'nef', 'vpu', 'vpr', 'vif']



# Functions
def save_for_download(filename, **data):
    '''Save data for download (anotations, etc)'''
    from itertools import izip
    import StringIO
    import zipfile

    sep = '\t'
    fmt = lambda x: sep.join(x)+'\n'

    with zipfile.ZipFile(filename, 'w',
                         compression=zipfile.ZIP_DEFLATED) as zf:

        # Descriptor
        f = StringIO.StringIO()
        f.write(s)
        zf.writestr('README.txt', f.getvalue())

        for ac, time in izip(data['act'], data['times']):
            f = StringIO.StringIO()
            f.write('# '+data['pcode']+', '+str(int(time))+ ' days from infection\n')
            f.write('# '+fmt(data['alpha']))
            for acp in ac.T:
                f.write(fmt([str(acpn) for acpn in acp]))

            zf.writestr(str(int(time))+'_days.tsv', f.getvalue())



# Script
if __name__ == '__main__':

    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name

        # TODO: export whle trajectories

        # Sample by sample
        for i, sample in enumerate(patient.itersamples()):
            samplename = patient.code+'_sample_'+str(i+1)

            for protein in proteins:
                fn_in = sample.get_allele_counts_filename(protein, type='aa')
                fn_out = get_fn2_out(samplename, protein, 'npy', type='aa')
                if (not os.path.isfile(fn_out)) and (not os.path.islink(fn_out)):
                    os.symlink(fn_in, fn_out)

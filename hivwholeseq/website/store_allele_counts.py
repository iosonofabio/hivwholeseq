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
Value (TSV) format a.k.a. tab-delimited CSV. The comment character is # (hash).
Each file contains a separate time point. The estimated number of days from
infection is in the filename and again in the first line of the file.

Each row is a position in the genome, in nucleotides. Each column is a possible
nucleotide. The - (dash) indicates deletions, N an ambiguous nucleotide. The
latter should be rare and can be ignored for most purposes.

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

        ## Allele count trajectories
        #fn_in = get_allele_count_trajectories_filename(pname, 'genomewide')
        #npz = np.load(fn_in)
        #ind = npz['ind']
        #act = npz['act']

        #times = patient.times[ind]

        ## Write to file
        ## 1. Binary
        #fn_out = get_fn_out(patient.code, 'genomewide', 'npz')
        #np.savez(fn_out, times=times, act=act, alpha=alphal)

        ## 2. For download
        #fn_out = get_fn_out(patient.code, 'genomewide', 'zip')
        #save_for_download(fn_out, pcode=patient.code,
        #                  times=times, act=act, alpha=alphal)


        ## Coverage trajectories (binary format only)
        #cov = act.sum(axis=1)

        ## Write output
        #fn_out = get_coverage_filename(patient.code, 'genomewide')
        #np.savez(fn_out, times=times, cov=cov)


        # Sample by sample
        for i, sample in enumerate(patient.itersamples()):
            samplename = patient.code+'_sample_'+str(i+1)

            #for region in ['F'+str(j) for j in xrange(1, 7)] + ['genomewide']:
            #    fn_in = sample.get_allele_counts_filename(region)
            #    fn_out = get_fn2_out(samplename, region, 'npy')
            #    if (not os.path.isfile(fn_out)) and (not os.path.islink(fn_out)):
            #        os.symlink(fn_in, fn_out)

            for protein in proteins:
                fn_in = sample.get_allele_counts_filename(protein, type='aa')
                fn_out = get_fn2_out(samplename, protein, 'npy', type='aa')
                if (not os.path.isfile(fn_out)) and (not os.path.islink(fn_out)):
                    os.symlink(fn_in, fn_out)

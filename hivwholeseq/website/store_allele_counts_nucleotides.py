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
import shutil
import numpy as np

from hivwholeseq.utils.miseq import alphal
from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.website.filenames import get_allele_count_trajectories_filename as get_fn_out_traj
from hivwholeseq.website.filenames import get_allele_counts_filename as get_fn_out_sample
from hivwholeseq.website.filenames import get_coverage_filename


# Globals
s = \
'''The data is stored in the following format. Each file is in Tab Separated
Value (TSV) format. Each file contains a separate time point. The estimated
number of days from infection is in the filename and again in the first line
of the file.

Each row is a position in the genome, in nucleotides. The coordinates refer to
the reference sequence for each patient, which can be downloaded in FASTA or
Genbank (GB) format from the website.

Each column is a possible nucleotide, according to the alphabet ACGT-N. The -
(dash) indicates deletions (with respect to the patient reference sequence),
N an ambiguous nucleotide. The latter indicates ambiguous nucleotides, should
be rare and can be ignored for most purposes.

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

        print 'Trajectory'
        for ac, time in izip(data['act'], data['times']):
            f = StringIO.StringIO()
            f.write('# '+data['pcode']+', '+str(int(time))+ ' days from infection\n')
            f.write('# '+fmt(data['alpha']))
            for acp in ac.T:
                f.write(fmt([str(acpn) for acpn in acp]))

            print 'Sample time:', int(time)
            zf.writestr(str(int(time))+'_days.tsv', f.getvalue())



# Script
if __name__ == '__main__':

    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name

        # Allele count trajectories
        (act, ind) = patient.get_allele_count_trajectories('genomewide')
        times = patient.times[ind]

        # Write to file
        # 1. Binary
        fn_out = get_fn_out_traj(patient.code, 'genomewide', 'npz')
        np.savez(fn_out, times=times, act=act, alpha=alphal)

        # 2. For download
        fn_out = get_fn_out_traj(patient.code, 'genomewide', 'zip')
        save_for_download(fn_out, pcode=patient.code,
                          times=times, act=act, alpha=alphal)


        # Coverage trajectories (binary format only)
        cov = act.sum(axis=1)

        # Write output
        fn_out = get_coverage_filename(patient.code, 'genomewide')
        np.savez(fn_out, times=times, cov=cov)

        # Sample by sample
        print 'Sample by sample'
        for i, sample in enumerate(patient.itersamples()):
            samplename = patient.code+'_sample_'+str(i+1)
            print "Sample time:", sample['days since infection'], sample.name

            for region in ['F'+str(j) for j in xrange(1, 7)] + ['genomewide']:
                fn_in = sample.get_allele_counts_filename(region)
                fn_out = get_fn_out_sample(samplename, region, 'npy')
                if os.path.isfile(fn_out):
                    os.remove(fn_out)
                elif os.path.islink(fn_out):
                    os.unlink(fn_out)

                if os.path.isfile(fn_in):
                    shutil.copy(fn_in, fn_out)

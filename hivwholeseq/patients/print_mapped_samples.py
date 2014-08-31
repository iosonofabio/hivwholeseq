# vim: fdm=marker
'''
author:     Fabio Zanini
date:       29/06/14
content:    Print a table with what samples contain mapped and filtered reads
            from what PCR, and the same for allele counts.
'''
# Modules
import os
import argparse
import numpy as np


from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_filtered_filename, get_allele_counts_filename



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele frequency trajectories')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--stage', default='filtered',
                        help='Analysis type to print (filtered, counts)')

    args = parser.parse_args()
    pnames = args.patients
    VERBOSE = args.verbose
    stage = args.stage

    patients = load_patients() 
    if pnames is not None:
        patients = patients.loc[patients.index.isin(pnames)]

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        
        n_samples_tot = len(patient.samples)

        patient.discard_nonsequenced_samples()
        n_samples = len(patient.samples)
        n_samples_nonseq = n_samples_tot - n_samples

        # Time x fragment x PCR
        has_reads = np.zeros((n_samples, 6, 2), bool)
        for i, (samplename_pat, sample) in enumerate(patient.samples.iterrows()):
            for j, fragment in enumerate(['F'+str(n+1) for n in xrange(6)]):
                for k, PCR in enumerate((1, 2)):
                    if stage == 'filtered':
                        fn = get_mapped_filtered_filename(pname, samplename_pat, fragment, PCR=PCR)
                    elif stage == 'counts':
                        fn = get_allele_counts_filename(pname, samplename_pat, fragment, PCR=PCR)
                    if os.path.isfile(fn):
                        has_reads[i, j, k] = True

        # Print table
        print pname
        print '# samples:', n_samples_tot, '(seq: '+str(n_samples)+', nonseq:'+str(n_samples_nonseq)+')'
        print '-' * (12 + 1 + 1 + (1 + 3 + 1 + 1) * 6)
        for i, row in enumerate(has_reads):
            line = '{:12s}'.format(patient.samples.iloc[i].name)
            line = line+' |'
            for j in xrange(6):
                cell = row[j]
                if cell[0] and cell[1]:
                    line = line+' x x |'
                elif cell[0] and not cell[1]:
                    line = line+' x   |'
                elif not cell[0] and cell[1]:
                    line = line+'   x |'
                else:
                    line = line+'     |'
            print line
        print '-' * (12 + 1 + 1 + (1 + 3 + 1 + 1) * 6)
        print ''
                    


        





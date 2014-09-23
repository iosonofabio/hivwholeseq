# vim: fdm=marker
'''
author:     Fabio Zanini
date:       31/08/14
content:    Check status of patients: initial reference, genomewide reference,
            mapped reads, filtered reads, allele counts, allele frequency
            trajectories, linkage data structures.
'''
# Modules
import os
import argparse
from hivwholeseq.patients.patients import load_patients, load_patient, Patient, \
        SamplePat
from hivwholeseq.patients.filenames import get_decontaminate_summary_filename


# Globals
title_len = 15
cell_len = 7



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    VERBOSE = args.verbose

    patients = load_patients()
    for pname, p in patients.iterrows():
        p = Patient(p)

        n_samples = len(p.samples)
        p.discard_nonsequenced_samples()
        n_samples_seq = len(p.samples)
        print p.name, '# samples:', str(n_samples)+ ' ('+str(n_samples_seq)+' sequenced)'
        title = 'Samples'
        line = ('{:<'+str(title_len)+'}').format(title+':')
        line = line+' '.join(p.samples.index.tolist())
        print line

        fn = p.folder
        title = 'Folder'
        line = ('{:<'+str(title_len)+'}').format(title+':')
        if os.path.isdir(fn):
            status = 'OK'
        else:
            status = 'MISS'
        line = line + ('{:<'+str(cell_len)+'}').format(status)
        print line

        if status != 'OK':
            print ''
            continue

        title = 'References'
        line = ('{:<'+str(title_len)+'}').format(title+':')
        stati = []
        for fragment in ('F'+str(i+1) for i in xrange(6)):
            fn = p.get_reference_filename(fragment)
            if os.path.isfile(fn):
                status = 'OK'
            else:
                status = 'MISS'
            stati.append(status)
            line = line + fragment + ': ' + ('{:>'+str(cell_len - len(fragment) - 1)+'}').format(status) + '  '
        print line

        if frozenset(stati) != frozenset(['OK']):
            print ''
            continue

        title = 'Genome ref'
        line = ('{:<'+str(title_len)+'}').format(title+':')
        fn = p.get_reference_filename('genomewide')
        if os.path.isfile(fn):
            status = 'OK'
        else:
            status = 'MISS'
        line = line + ('{:<'+str(cell_len)+'}').format(status)
        print line

        if status != 'OK':
            print ''
            continue

        title = 'Map + filter'
        line = ('{:<'+str(title_len)+'}').format(title+':')
        print line
        for samplename, sample in p.samples.iterrows():
            sample = SamplePat(sample)
            title = sample.name
            line = ('{:<'+str(title_len)+'}').format(title+':')
            
            stati = []
            for fragment in ('F'+str(i+1) for i in xrange(6)):
                fn = sample.get_mapped_filtered_filename(fragment, PCR=1)
                if os.path.isfile(fn):
                    status = 'OK'
                else:
                    fn = sample.get_mapped_filtered_filename(fragment, PCR=2)
                    if os.path.isfile(fn):
                        status = 'PCR2'
                    else:
                        status = 'MISS'
                stati.append(status)
                line = line + fragment + ': ' + ('{:>'+str(cell_len - len(fragment) - 1)+'}').format(status) + '  '
            print line

        title = 'Decontaminate'
        line = ('{:<'+str(title_len)+'}').format(title+':')
        print line
        for samplename, sample in p.samples.iterrows():
            sample = SamplePat(sample)
            title = sample.name
            line = ('{:<'+str(title_len)+'}').format(title+':')
            
            stati = []
            for fragment in ('F'+str(i+1) for i in xrange(6)):
                fn = get_decontaminate_summary_filename(sample.patient, samplename,
                                                        fragment, PCR=1)
                if os.path.isfile(fn):
                    status = 'OK'
                else:
                    fn = get_decontaminate_summary_filename(sample.patient, samplename,
                                                            fragment, PCR=1)
                    if os.path.isfile(fn):
                        status = 'PCR2'
                    else:
                        status = 'MISS'
                stati.append(status)
                line = line + fragment + ': ' + ('{:>'+str(cell_len - len(fragment) - 1)+'}').format(status) + '  '
            print line

        title = 'Consensus'
        line = ('{:<'+str(title_len)+'}').format(title+':')
        print line
        for samplename, sample in p.samples.iterrows():
            sample = SamplePat(sample)
            title = sample.name
            line = ('{:<'+str(title_len)+'}').format(title+':')
            
            stati = []
            for fragment in ('F'+str(i+1) for i in xrange(6)):
                fn = sample.get_consensus_filename(fragment, PCR=1)
                if os.path.isfile(fn):
                    status = 'OK'
                else:
                    fn = sample.get_consensus_filename(fragment, PCR=2)
                    if os.path.isfile(fn):
                        status = 'PCR2'
                    else:
                        status = 'MISS'
                stati.append(status)
                line = line + fragment + ': ' + ('{:>'+str(cell_len - len(fragment) - 1)+'}').format(status) + '  '
            print line

        title = 'Allele counts'
        line = ('{:<'+str(title_len)+'}').format(title+':')
        print line
        for samplename, sample in p.samples.iterrows():
            sample = SamplePat(sample)
            title = sample.name
            line = ('{:<'+str(title_len)+'}').format(title+':')
            
            stati = []
            for fragment in ('F'+str(i+1) for i in xrange(6)):
                fn = sample.get_allele_counts_filename(fragment, PCR=1)
                if os.path.isfile(fn):
                    status = 'OK'
                else:
                    fn = sample.get_allele_counts_filename(fragment, PCR=2)
                    if os.path.isfile(fn):
                        status = 'PCR2'
                    else:
                        status = 'MISS'
                stati.append(status)
                line = line + fragment + ': ' + ('{:>'+str(cell_len - len(fragment) - 1)+'}').format(status) + '  '
            print line


        print ''

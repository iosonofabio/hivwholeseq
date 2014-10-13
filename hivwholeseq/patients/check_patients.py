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
import datetime
from hivwholeseq.patients.patients import load_patients, load_patient, Patient, \
        SamplePat
from hivwholeseq.patients.filenames import get_decontaminate_summary_filename
from hivwholeseq.generic_utils import modification_date


# Globals
title_len = 15
cell_len = 7


# Functions
def pretty_print_info(p, title, name, method, name_requisite=None, mod_dates={},
                      fragments=['F'+str(i+1) for i in xrange(6)]):
    '''Pretty printer for patient pipeline info'''
    import os, sys
    from hivwholeseq.patients.samples import SamplePat

    stati = set()
    
    line = ('{:<'+str(title_len)+'}').format(title+':')
    print line
    for samplename, sample in p.samples.iterrows():
        sample = SamplePat(sample)
        title = sample.name
        line = ('{:<'+str(title_len)+'}').format(title+':')
        
        for fragment in fragments:
            if isinstance(method, basestring) and hasattr(sample, method):
                fun = getattr(sample, method)
                fn = fun(fragment)
            else:
                fn = method(sample.patient, samplename, fragment)
            if os.path.isfile(fn):
                md = modification_date(fn)

                if name_requisite is None:
                    status = 'OK'
                    mod_dates[(name, fragment, samplename)] = md

                elif ((name_requisite, fragment, samplename) in mod_dates) and \
                   (md > mod_dates[(name_requisite, fragment, samplename)]):
                    status = 'OK'
                    mod_dates[(name, fragment, samplename)] = md

                elif ((name_requisite, fragment) in mod_dates) and \
                   (md > mod_dates[(name_requisite, fragment)]):
                    status = 'OK'
                    mod_dates[(name, fragment, samplename)] = md

                else:
                    status = 'OLD'

            else:
                status = 'MISS'
            stati.add(status)
            line = line+fragment+': '+\
                ('{:>'+str(cell_len - len(fragment) - 1)+'}').format(status)+'  '
        print line


    if 'OLD' in stati:
        sys.exit()



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    VERBOSE = args.verbose
    pnames = args.patients

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for pname, p in patients.iterrows():
        p = Patient(p)

        mod_dates = {}

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
                mod_dates[('reference', fragment)] = modification_date(fn)
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
            mod_dates[('reference', 'genomewide')] = modification_date(fn)
        else:
            status = 'MISS'
        line = line + ('{:<'+str(cell_len)+'}').format(status)
        print line

        if status != 'OK':
            print ''
            continue

        pretty_print_info(p, 'Map + filter', 'filter',
                          'get_mapped_filtered_filename',
                          'reference', mod_dates)

        pretty_print_info(p, 'Decontaminate', 'decontaminate',
                          get_decontaminate_summary_filename,
                          'filter', mod_dates)

        pretty_print_info(p, 'Consensus', 'consensus',
                          'get_consensus_filename',
                          'decontaminate', mod_dates)

        pretty_print_info(p, 'Allele counts', 'allele counts',
                          'get_allele_counts_filename',
                          'decontaminate', mod_dates)

        pretty_print_info(p, 'Allele cocounts', 'allele cocounts',
                          'get_allele_cocounts_filename',
                          'decontaminate', mod_dates)

        print ''

        break

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
def check_reference_overlap(p, VERBOSE=0):
    '''Check whether the reference from the various fragments overlap correctly'''
    from seqanpy import align_ladder
    from hivwholeseq.sequence_utils import pretty_print_pairwise_ali

    fragments = ['F'+str(i+1) for i in xrange(6)]
    title = 'Overlaps'
    line = ('{:<'+str(title_len)+'}').format(title+':')
    stati = []
    for i in xrange(len(fragments) - 1):
        ref1 = p.get_reference(fragments[i])
        ref2 = p.get_reference(fragments[i+1])
        (score, ali1, ali2) = align_ladder(ref1, ref2,
                                           score_gapopen=-10,
                                           score_gapext=-1)

        start2 = len(ali2) - len(ali2.lstrip('-'))
        end1 = len(ali1.rstrip('-'))

        if VERBOSE >= 4:
            pretty_print_pairwise_ali((ali1[start2: end1], ali2[start2: end1]),
                                      name1=fragments[i],
                                      name2=fragments[i+1],
                                      width=100)
        
        if ali1[start2: end1].count('-') == ali2[start2: end1].count('-'):
            status = 'OK'
        else:
            status = 'GAPS'
            import ipdb; ipdb.set_trace()

        line = line+fragments[i]+': '+\
            ('{:>'+str(cell_len - len(fragments[i]) - 1)+'}').format(status)+'  '
        stati.append(status)

    print line

    if 'GAPS' in stati:
        raise ValueError('GAPS status found') 


def pretty_print_info_patient(p, title, name, method, name_requisite=None, mod_dates={},
                              VERBOSE=0):
    '''Pretty printer for whole-patient info, fragment by fragment'''
    import os, sys

    line = ('{:<'+str(title_len)+'}').format(title+':')
    stati = []
    for fragment in ('F'+str(i+1) for i in xrange(6)):
        if isinstance(method, basestring):
            fun = getattr(p, method)
            fn = fun(fragment)
        else:
            fn = method(p.name, fragment)

        if os.path.isfile(fn):
            md = modification_date(fn)
            mod_dates[('refmap', fragment)] = md

            if name_requisite is None:
                status = 'OK'

            elif ((name_requisite, fragment) in mod_dates):
                if md > mod_dates[(name_requisite, fragment)]:
                    status = 'OK'
                else:
                    status = 'OLD'

        else:
            status = 'MISS'

        stati.append(status)
        line = line + fragment + ': ' + ('{:>'+str(cell_len - len(fragment) - 1)+'}').format(status) + '  '
    print line


def pretty_print_info(p, title, name, method, name_requisite=None, mod_dates={},
                      VERBOSE=0):
    '''Pretty printer for patient pipeline info'''
    import os, sys
    from hivwholeseq.patients.samples import SamplePat
    from hivwholeseq.mapping_utils import get_number_reads

    # NOTE: this function is used to check both entire patients and single samples
    if isinstance(p, SamplePat):
        sample_iter = [(p.name, p)]
    else:
        sample_iter = p.samples.iterrows()

    fragments=['F'+str(i+1) for i in xrange(6)]

    stati = set()
    line = ('{:<'+str(title_len)+'}').format(title+':')
    print line
    for samplename, sample in sample_iter:
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
                mod_dates[(name, fragment, samplename)] = md

                if name_requisite is None:
                    status = 'OK'

                elif ((name_requisite, fragment, samplename) in mod_dates):
                    if md > mod_dates[(name_requisite, fragment, samplename)]:
                        status = 'OK'
                    else:
                        status = 'OLD'

                elif ((name_requisite, fragment) in mod_dates):
                    if md > mod_dates[(name_requisite, fragment)]:
                        status = 'OK'
                    else:
                        status = 'OLD'

                elif 'contaminated' in sample[fragment]:
                    status = 'CONT'
                
                else:
                    status = 'ERROR'

            else:
                status = 'MISS'

            # Check the number of reads if requested
            if (status == 'OK') and (fn[-3:] == 'bam') and (VERBOSE >= 3):
                status = str(get_number_reads(fn))

            stati.add(status)
            line = line+fragment+': '+\
                ('{:>'+str(cell_len - len(fragment) - 1)+'}').format(status)+'  '
        print line


    if 'OLD' in stati:
        raise ValueError('OLD status found') 


def pretty_print_info_genomewide(p, title, name, method, mod_dates={}, VERBOSE=0,
                                 require_all=True):
    '''Pretty printer for patient pipeline info'''

    def check_requisite_genomewide(md, name_requisite, samplename, mod_dates,
                                   require_all=True):
        '''Check requisites for genomewide observables'''
        stati = []
        fragments=['F'+str(i+1) for i in xrange(6)]
        for fragment in fragments:
            if (name_requisite, fragment, samplename) not in mod_dates:
                stati.append('MISS')
            elif md < mod_dates[(name_requisite, fragment, samplename)]:
                stati.append('OLD')
            else:
                stati.append('OK')

        if 'OLD' in stati:
            return 'OLD'
        else:
            if require_all:
                if 'MISS' in stati:
                    return 'MISS'
                else:
                    return 'OK'
            else:
                if 'OK' in stati:
                    return 'OK'
                else:
                    return 'MISS'

    def check_contamination_genomewide(sample):
        '''Check whether any of the fragment samples is contaminated'''
        fragments=['F'+str(i+1) for i in xrange(6)]
        for fragment in fragments:
            if 'contaminated' in sample[fragment]:
                return True
        return False

    import os, sys
    from hivwholeseq.patients.samples import SamplePat

    # NOTE: this function is used to check both entire patients and single samples
    if isinstance(p, SamplePat):
        sample_iter = [(p.name, p)]
    else:
        sample_iter = p.samples.iterrows()

    stati = set()    
    line = ('{:<'+str(title_len)+'}').format(title+':')
    print line
    for samplename, sample in sample_iter:
        sample = SamplePat(sample)
        title = sample.name
        line = ('{:<'+str(title_len)+'}').format(title+':')
        
        if isinstance(method, basestring) and hasattr(sample, method):
            fun = getattr(sample, method)
            fn = fun('genomewide')
        else:
            fn = method(sample.patient, samplename, 'genomewide')
        if os.path.isfile(fn):
            md = modification_date(fn)
            mod_dates[(name, 'genomewide', samplename)] = md

            if name is None:
                status = 'OK'

            elif check_contamination_genomewide(sample):
                status = 'CONT'

            else:
                status = check_requisite_genomewide(md, name, samplename, mod_dates,
                                                    require_all=require_all)

        else:
            status = 'MISS'

        # Check the number of reads if requested
        if (status == 'OK') and (fn[-3:] == 'bam') and (VERBOSE >= 3):
            status = str(get_number_reads(fn))

        stati.add(status)
        line = line + ('{:<'+str(cell_len)+'}').format(status)
        print line

    if 'OLD' in stati:
        raise ValueError('OLD status found') 



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

        check_reference_overlap(p)

        # NOTE: we modify reference annotations all the time, but not the sequence
        # FIXME: find a better way to do this?
        pretty_print_info(p, 'Map + filter', 'filter',
                          'get_mapped_filtered_filename',
                          None,#'reference',
                          mod_dates)

        from hivwholeseq.patients.filenames import get_mapped_filtered_filename
        pretty_print_info(p, 'Decontaminate', 'decontaminate',
                          lambda pn, sn, fr: get_mapped_filtered_filename(pn, sn, fr, decontaminated=True),
                          'filter', mod_dates)

        pretty_print_info(p, 'Consensus', 'consensus',
                          'get_consensus_filename',
                          'decontaminate', mod_dates)

        pretty_print_info_genomewide(p, 'Cons genomewide', 'consensus',
                                     'get_consensus_filename',
                                     mod_dates)

        pretty_print_info(p, 'Allele counts', 'allele counts',
                          'get_allele_counts_filename',
                          'decontaminate', mod_dates)

        pretty_print_info(p, 'Allele cocounts', 'allele cocounts',
                          'get_allele_cocounts_filename',
                          'decontaminate', mod_dates)

        pretty_print_info_genomewide(p, 'Allele counts genomewide', 'allele counts',
                                     'get_allele_counts_filename',
                                     mod_dates, require_all=False)

        pretty_print_info_patient(p, 'Maps to HXB2', 'reference',
                                  'get_map_coordinates_reference_filename',
                                  'reference', mod_dates)


        print ''


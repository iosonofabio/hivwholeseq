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
from hivwholeseq.patients.patients import (load_patients, load_patient,
                                           iterpatient, SamplePat)
from hivwholeseq.patients.filenames import get_decontaminate_summary_filename
from hivwholeseq.utils.generic import modification_date
from hivwholeseq.utils.argparse import PatientsAction



# Globals
title_len = 15
cell_len = 7



# Functions
def check_reference_overlap(p, VERBOSE=0):
    '''Check whether the reference from the various fragments overlap correctly'''
    from seqanpy import align_ladder
    from hivwholeseq.utils.sequence import pretty_print_pairwise_ali

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


def print_info_patient(p, title, name, method, name_requisite=None,
                       VERBOSE=0):
    '''Pretty printer for whole-patient info, fragment by fragment'''
    import os, sys

    mod_dates = p.mod_dates

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


def print_info(p, title, name, method, name_requisite=None, VERBOSE=0):
    '''Pretty printer for patient pipeline info'''
    import os, sys
    from hivwholeseq.patients.samples import SamplePat
    from hivwholeseq.utils.mapping import get_number_reads

    mod_dates = p.mod_dates

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
                        print fn, md, mod_dates[(name_requisite, fragment, samplename)]

                elif ((name_requisite, fragment) in mod_dates):
                    if md > mod_dates[(name_requisite, fragment)]:
                        status = 'OK'
                    else:
                        status = 'OLD'

                        # NOTE: on Nov 13, 2014 I updated the mod dates of all
                        # references by mistake, without actually changing the
                        # sequences (ironically, probably testing a backup system
                        # for the refs themselves). So if the requisite is a ref
                        # seq and the date is this one, it's OK
                        if ((name_requisite == 'reference') and
                            mod_dates[(name_requisite, fragment)].date() == \
                            datetime.date(2014, 11, 13)):
                            status = 'OK'


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


def print_info_genomewide(p, title, name, method, VERBOSE=0, require_all=True):
    '''Pretty printer for patient pipeline info'''

    mod_dates = p.mod_dates

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


def check_pipeline_patient(p, VERBOSE=0):
    '''Check patient pipeline'''
    from hivwholeseq.utils.exceptions import PipelineError

    def print_info_summary(p):
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
            raise PipelineError('Missing patient folder')

        return status

    def print_info_references(p):
        '''Print info on references'''
        title = 'References'
        line = ('{:<'+str(title_len)+'}').format(title+':')
        stati = []
        for fragment in ('F'+str(i+1) for i in xrange(6)):
            fn = p.get_reference_filename(fragment)
            if os.path.isfile(fn):
                status = 'OK'
                p.mod_dates[('reference', fragment)] = modification_date(fn)
            else:
                status = 'MISS'
            stati.append(status)
            line = line + fragment + ': ' + ('{:>'+str(cell_len - len(fragment) - 1)+'}').format(status) + '  '
        print line
    
        if frozenset(stati) != frozenset(['OK']):
            print ''
            raise PipelineError('Amplicon reference failed!')
    
        title = 'Genome ref'
        line = ('{:<'+str(title_len)+'}').format(title+':')
        fn = p.get_reference_filename('genomewide', 'fasta')
        if os.path.isfile(fn):
            status = 'OK'
            p.mod_dates[('reference', 'genomewide')] = modification_date(fn)
        else:
            status = 'MISS'
        line = line + ('{:<'+str(cell_len)+'}').format(status)
        print line
    
        if status != 'OK':
            print ''
            raise PipelineError('Genomewide reference failed!')
    
        check_reference_overlap(p)

        title = 'Annotated'
        line = ('{:<'+str(title_len)+'}').format(title+':')
        fn = p.get_reference_filename('genomewide', 'gb')
        if os.path.isfile(fn):
            md = modification_date(fn)
            if md >= p.mod_dates[('reference', 'genomewide')]:
                status = 'OK'
            else:
                status = 'OLD'
        else:
            status = 'MISS'
        line = line + ('{:<'+str(cell_len)+'}').format(status)
        print line
        if status != 'OK':
            print ''
            raise PipelineError('Annotated reference failed!')


    p.mod_dates = {}
    
    print_info_summary(p)

    print_info_references(p)

    from hivwholeseq.patients.filenames import get_mapped_filtered_filename
    print_info(p, 'Map + filter', 'filter',
               lambda pn, sn, fr: get_mapped_filtered_filename(pn, sn, fr, decontaminated=False),
               'reference')

    print_info(p, 'Decontaminate', 'decontaminate',
               lambda pn, sn, fr: get_mapped_filtered_filename(pn, sn, fr, decontaminated=True),
               'filter')

    print_info(p, 'Consensus', 'consensus',
               'get_consensus_filename',
               'decontaminate')

    print_info_genomewide(p, 'Cons genomewide', 'consensus',
                          'get_consensus_filename')

    print_info(p, 'Allele counts', 'allele counts',
               'get_allele_counts_filename',
               'decontaminate')

    print_info(p, 'Allele cocounts', 'allele cocounts',
               'get_allele_cocounts_filename',
               'decontaminate')

    print_info_genomewide(p, 'Allele counts genomewide', 'allele counts',
                          'get_allele_counts_filename',
                          require_all=False)

    print_info_patient(p, 'Maps to HXB2', 'reference',
                       'get_map_coordinates_reference_filename',
                       'reference')


    print ''




# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--patients', action=PatientsAction,
                        help='Patient to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    VERBOSE = args.verbose
    pnames = args.patients

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for pname, p in iterpatient(patients):
        check_pipeline_patient(p, VERBOSE=VERBOSE)



# vim: fdm=marker
'''
author:     Fabio Zanini
date:       15/02/15
content:    Copy the data folder.
'''
# Modules
import os
import sys
import shutil

from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.samples import load_samples_sequenced, SamplePat
from hivwholeseq.utils.generic import mkdirs



# Globals
fragments = ['F'+str(i) for i in xrange(1, 7)]



# Functions
def copy_reference(dst_fn):
    '''Copy reference sequence(s) with annotations'''
    from hivwholeseq.reference import load_custom_reference_filename
    fn_src = load_custom_reference_filename('HXB2', format='gb')
    fn_dst = dst_fn+os.path.basename(fn_src)
    shutil.copy(fn_src, fn_dst)


def copy_initial_reference(patient, dst_fn):
    '''Copy initial patient mapping reference'''
    ref_fn = dst_fn+'reference/'
    mkdirs(ref_fn)

    for fragment in fragments:
        fn_src = patient.get_reference_filename(fragment)
        fn_dst = ref_fn+os.path.basename(fn_src)
        shutil.copy(fn_src, fn_dst)

    fn_src = patient.get_reference_filename('genomewide', format='gb')
    fn_dst = ref_fn+os.path.basename(fn_src)
    shutil.copy(fn_src, fn_dst)


def copy_folder(patient, dst_fn, foldername):
    '''Copy a whole folder'''
    src_fn = patient.folder+foldername+os.sep
    map_fn = dst_fn+foldername+os.sep

    mkdirs(map_fn)
    for fn_src in src_fn.listdir():
        copy(src_fn+fn_src, map_fn+fn_src)


def copy_glob(sample, dst_fn, pattern):
    '''Copy a glob pattern'''
    import glob
    src_fn = sample.get_foldername()
    fns_src = glob.glob(src_fn+pattern)
    for fn_src in fns_src:
        fn_dst = dst_fn+os.path.basename(fn_src)
        shutil.copy(fn_src, fn_dst)


def copy_reads(sample, dst_fn):
    '''Copy mapped reads'''
    dst_fn = dst_fn+'mapped_to_initial'+os.sep
    for fragment in fragments:
        fn_src = sample.get_mapped_filtered_filename(fragment, decontaminated=True)
        if not os.path.isfile(fn_src):
            continue

        fn_dst = dst_fn+os.path.basename(fn_src)
        shutil.copy(fn_src, fn_dst)



# Script
if __name__ == '__main__'

    parser = argparse.ArgumentParser(description='Copy data folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('destination',
                        help='Destination folder')
    parser.add_argument('--strip-PCR1', action='store_true',
                        help='Strip the ../PCR1/ part of the file tree')

    args = parser.parse_args()
    dst_fn = args.destination.lstrip(os.sep)+os.sep
    stirp_PCR1 = args.strip_PCR1

    patients_fn = dst_fn+'patients/'
    ref_fn = dst_fn+'reference/'

    print 'Make root folders'
    mkdirs(patients_fn)
    mkdirs(ref_fn)

    print 'Reference sequences'
    copy_reference(ref_fn)

    
    patients = load_patients()
    for pname, patient in patients.iterrows():
        print pname
        patient = Patient(patient)

        print 'Make folder'
        pat_fn = patients_fn+pname+os.sep
        mkdirs(pat_fn)

        print 'Mapping reference'
        copy_initial_reference(patient, pat_fn)

        
        print 'Coordinate maps'
        copy_folder(patient, pat_fn, 'coordinate_maps')


        print 'Alignments'
        copy_folder(patient, pat_fn, 'alignments')


        print 'Trees'
        copy_folder(patient, pat_fn, 'trees')


        print 'Haplotypes'
        copy_folder(patient, pat_fn, 'haplotypes')


        print 'Samples'
        for samplename, sample in patient.samples.iterrows():
            print samplename
            sample = SamplePat(sample)

            print 'Make folder'
            sm_fn = pat_fn+samplename+os.sep
            if not strip_PCR1:
                sm_fn += 'PCR1'+os.sep
            mkdirs(sm_fn)


            print 'Consensus'
            copy_glob(sample, sm_fn, 'consensus')
            

            print 'Allele counts'
            copy_glob(sample, sm_fn, 'allele_counts')
            

            print 'Allele cocounts'
            copy_glob(sample, sm_fn, 'allele_cocounts')


            print 'Reads'
            copy_reads(sample, sm_fn)

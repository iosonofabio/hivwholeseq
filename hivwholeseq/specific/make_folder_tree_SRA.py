# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/05/14
content:    Make a folder tree for the Short Read Archive (SRA) from an old
            sequencing run.
'''
# Modules
import os
import sys
from itertools import izip
import argparse
import subprocess as sp

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.datasets import MiSeq_runs



# Globals
root_folder = '/ebio/ag-neher/share/public_data/SRA/'
miseq_name = 'M01346'



# Functions
def convert_date(date):
    '''Convert date to SRA format'''
    (y, m, d) = date.split('-')
    y = y[-2:]
    datenew = y+m+d
    return datenew


def convert_run(run):
    '''Convert run to SRA format'''
    runnew = '00'+run[-2:]
    return runnew


def find_flowcell_code(dataset):
    '''Find flowcell code'''
    data_folder = dataset['folder']
    subdirs = filter(lambda x: 'adapterID' in x, os.listdir(data_folder))
    fn = data_folder+subdirs[0]+'/read1.fastq.gz'
    import gzip
    from Bio import SeqIO
    with gzip.open(fn, 'rb') as f:
        read_iter = SeqIO.parse(f, 'fastq')
        read = read_iter.next()
        rname = read.id
    flowcell = rname.split(':')[2]

    # Get it from the folders in abt6_ga2
    if '000000' not in flowcell:
        import re
        seq_run = int(re.findall(r'\d+$', dataset['name'])[0]) 
        root_folder_ga2 = '/ebio/abt6_ga2/images/'
        fns = [fn for fn in os.listdir(root_folder_ga2)
               if '{:04d}'.format(seq_run)+'_00000' in fn]
        if len(fns) == 1:
            rname = fns[0]
            flowcell = rname.split('_')[3]
        else:
            raise ValueError('Could not find flowcell')

    return flowcell


def make_root_folder(date, seqrun, fc):
    foldn = root_folder+date+'_'+miseq_name+'_'+seqrun+'_'+fc
    if os.path.isdir(foldn):
        print foldn+' already exists! Skipping...'
    else:
        mkdirs(foldn)
        print foldn+' created.'
    return foldn+'/'


def make_subdir_tree(folder, projectname, samplenames):
    folder = folder.rstrip('/')+'/'
    root_subfolder = folder+'Data/Intensities/FASTQ_deplex_1/'
    mkdirs(root_subfolder+'Undetermined_indices/Sample_lane1')
    for samplename in samplenames:
        mkdirs(root_subfolder+'Project_'+projectname+'/Sample_'+samplename)


def make_project_name(description):
    if 'Nextera' in description:
        if 'BluePippin' in description:
            projectname = 'HIVSwedenNexteraXTBluePippinAGNeher'
        else:
            projectname = 'HIVSwedenNexteraXTAGNeher'
    elif 'TruSeq' in description:
        projectname = 'HIVSwedenTruSeqNanoAGNeher'
    else:
        projectname = None
    return projectname


def convert_date_to_samplesheet(date):
    y = '20'+date[:2]
    m = date[2:4].lstrip('0')
    d = date[4:].lstrip('0')
    return '/'.join((m, d, y))


def convert_adapter(adaID, libtype):
    if libtype == 'NexteraXT':
        from hivwholeseq.sequencing.adapter_info import Nextera_XT_i7, Nextera_XT_i5
        (ada7, ada5) = tuple([int(x[1:]) for x in adaID.split('-')])
        adaI7 = 'N7'+'{:02d}'.format(ada7)
        adaI5 = 'S5'+'{:02d}'.format(ada5)
        seqI7 = Nextera_XT_i7[ada7]
        seqI5 = Nextera_XT_i5[ada5]
        return (adaI7, seqI7, adaI5, seqI5)

    elif libtype == 'TruSeq':
        from hivwholeseq.sequencing.adapter_info import TrueSeq_LT
        ada = int(adaID[2:])
        adaI7 = 'A'+'{:03d}'.format(ada)
        seqI7 = TrueSeq_LT[ada]
        return (adaI7, seqI7)

    else:
        raise ValueError('Other than NexteraXT and TruSeq not implemented')


def make_samplesheet_general(folder, projectname, date, seqrun, fc, n_cycles,
                             samplenames, adaIDs):
    ss_filename = folder+'SampleSheet.csv'
    if 'Nextera' in projectname:
        with open(ss_filename, 'w') as f:
            f.write('[Header]\n')
            f.write('IEMFileVersion,4\n')
            f.write('Investigator Name,Fabio Zanini\n')
            f.write('Experiment Name,run'+seqrun[1:]+'\n')
            f.write('Date,'+convert_date_to_samplesheet(date)+'\n')
            f.write('Workflow,GenerateFASTQ\n')
            f.write('Application,FASTQ Only\n')
            f.write('Assay,Nextera XT\n')
            f.write('Description,'+projectname+'\n')
            f.write('Chemistry,Amplicon\n')
            f.write('\n')

            f.write('[Reads]\n')
            f.write(str(n_cycles // 2)+'\n')
            f.write(str(n_cycles // 2)+'\n')
            f.write('\n')

            f.write('[Settings]\n')
            f.write('Adapter,CTGTCTCTTATACACATCT\n')
            f.write('\n')

            f.write('[Data]\n')
            f.write('Sample_ID,Sample_Name,Sample_Plate,Sample_Well,'+\
                    'I7_Index_ID,index,I5_Index_ID,index2,'+\
                    'Sample_Project,Description\n')
            for samplename, adaID in izip(samplenames, adaIDs):
                (adaI7, seqI7, adaI5, seqI5) = convert_adapter(adaID, 'NexteraXT')
                f.write(','.join((samplename, '', '', '',
                                  adaI7, seqI7, adaI5, seqI5,
                                  projectname, 'FabioZanini'))+'\n')

    elif 'TruSeq' in projectname:
        with open(ss_filename, 'w') as f:
            f.write('[Header]\n')
            f.write('IEMFileVersion,4\n')
            f.write('Investigator Name,Fabio Zanini\n')
            f.write('Experiment Name,run'+seqrun[1:]+'\n')
            f.write('Date,'+convert_date_to_samplesheet(date)+'\n')
            f.write('Workflow,GenerateFASTQ\n')
            f.write('Application,FASTQ Only\n')
            f.write('Assay,TruSeq LT\n')
            f.write('Description,'+projectname+'\n')
            f.write('Chemistry,Default\n')
            f.write('\n')

            f.write('[Reads]\n')
            f.write(str(n_cycles // 2)+'\n')
            f.write(str(n_cycles // 2)+'\n')
            f.write('\n')

            f.write('[Settings]\n')
            f.write('ReverseComplement,0\n')
            f.write('Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n')
            f.write('AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n')
            f.write('\n')

            f.write('[Data]\n')
            f.write('Sample_ID,Sample_Name,Sample_Plate,Sample_Well,'+\
                    'I7_Index_ID,index,'+\
                    'Sample_Project,Description\n')
            for samplename, adaID in izip(samplenames, adaIDs):
                (adaI7, seqI7) = convert_adapter(adaID, 'TruSeq')
                f.write(','.join((samplename, '', '', '',
                                  adaI7, seqI7,
                                  projectname, 'FabioZanini'))+'\n')


    else:
        raise ValueError('Other than NexteraXT and TruSeq not implemented')

    print 'General samplesheet created'


def make_runparameters(seqrun, fc, target_folder):
    '''Make runParameters XML file'''
    import re
    import hivwholeseq
    src_fn = os.path.dirname(hivwholeseq.__file__)+'/specific/runParameters_'+miseq_name+'.xml'
    dst_fn = target_folder+'runParameters.xml'

    with open(src_fn, 'r') as fin, open(dst_fn, 'w') as fout:
        for line in fin:

            # Flow cell name
            if '000000000-A8GYG' in line:
                line = line.replace('000000000-A8GYG', fc)

            # Run number
            if '<RunNumber>' in line:
                line = re.sub(r'<RunNumber>\d+', r'<RunNumber>'+str(int(seqrun)), line)

            fout.write(line)

    print 'runParameters.xml copied and corrected'


def copy_fastq_undetermined(data_folder, target_folder, symlink=False):
    '''Copy the fastq.gz or make links'''

    target_folder = target_folder.rstrip('/')+'/'
    root_subfolder = target_folder+'Data/Intensities/FASTQ_deplex_1/Undetermined_indices/Sample_lane1/'

    for i in ['1', '2']:
        src_fn = data_folder+'unclassified_reads/read'+i+'.fastq.gz'
        dst_fn = root_subfolder+'lane1_Undetermined_L001_R'+i+'_001.fastq.gz'
        if not symlink:
            out = sp.check_output(['ln', '-P', src_fn, dst_fn])
            if len(out):
                print out
        else:        
            os.symlink(src_fn, dst_fn)

    print 'Fastq files for undetermined indices linked'


def copy_fastq_samples(data_folder, target_folder, projectname, samplenames, adaIDs, symlink=False):
    '''Copy the fastq.gz or make links'''

    target_folder = target_folder.rstrip('/')+'/'
    root_subfolder = target_folder+'Data/Intensities/FASTQ_deplex_1/Project_'+projectname+'/'
    for samplename, adaID in izip(samplenames, adaIDs):
        src_fold = data_folder+'adapterID_'+adaID+'/'
        dst_fold = root_subfolder+'Sample_'+samplename+'/'

        if 'Nextera' in projectname:
            (adaI7, seqI7, adaI5, seqI5) = convert_adapter(adaID, 'NexteraXT')
            adastr = seqI7+'-'+seqI5
        elif 'TruSeq' in projectname:
            (adaI7, seqI7) = convert_adapter(adaID, 'TruSeq')
            adastr = seqI7
        else:
            raise ValueError('Other than NexteraXT and TruSeq not implemented')

        for i in ['1', '2']:
            src_fn = src_fold+'read'+i+'.fastq.gz'
            dst_fn = dst_fold+samplename+'_'+adastr+'_L001_R'+i+'_001.fastq.gz'


            if not symlink:
                out = sp.check_output(['ln', '-P', src_fn, dst_fn])
                if len(out):
                    print out
            else:
                os.symlink(src_fn, dst_fn)

        print 'Fastq files for '+samplename+' linked'



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Demultiplex HIV reads')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28, test_tiny)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--symlink', action='store_true',
                        help='Use symlinks instead of hardlinks')

    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose
    symlink = args.symlink

    dataset = MiSeq_runs[seq_run]

    # Get details
    data_folder = dataset['folder']
    n_cycles = dataset['n_cycles']
    date = convert_date(dataset['date'])
    seqrun = convert_run(seq_run)
    fc = find_flowcell_code(dataset)
    projectname = make_project_name(dataset['description'])
    if projectname is None:
        print 'Library is only a test. Skipping...'
        sys.exit()

    samplenames = map(lambda x: x.replace('_', '-'), dataset['samples'])
    adaIDs = dataset['adapters']

    folder = make_root_folder(date, seqrun, fc)
    
    make_subdir_tree(folder, projectname, samplenames)
    
    make_samplesheet_general(folder, projectname, date, seqrun, fc, n_cycles,
                             samplenames, adaIDs)

    make_runparameters(seqrun, fc, folder)

    copy_fastq_undetermined(data_folder, folder, symlink=symlink)

    copy_fastq_samples(data_folder, folder, projectname,
                       samplenames, adaIDs, symlink=symlink)


# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/12/13
content:    Module with all submit functions for the cluster. With this we can
            keep all cluster-specific code in one place.
'''
# Globals
import subprocess as sp

from . import JOBDIR, JOBLOGERR, JOBLOGOUT



# Functions
def nothing():
    '''Test function'''
    print 'JOBDIR', JOBDIR
    print 'JOBLOGERR', JOBLOGERR
    print 'JOBLOGOUT', JOBLOGOUT


def empty_log_folders():
    '''Empty log folders of old files'''
    import os, glob
    for folder in (JOBLOGOUT, JOBLOGERR):
        for fn in glob.glob(folder+'*'):
            # Keep .gitignore for the folders themselves
            if '.gitignore' not in fn:
                os.remove(fn)


# SEQUENCING
def fork_check_pipeline(seq_runs, adaIDs=None, pats=False,
                        detail=1, VERBOSE=0):
    '''Submit checking status of the pipeline to the cluster'''
    if VERBOSE:
        print 'Forking to the cluster'

    JOBSCRIPT = JOBDIR+'sequencing/check_pipeline.py'
    cluster_time = '00:59:59'
    vmem = '1G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'pipe',
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--runs', seq_runs,
                 '--adaIDs', adaIDs,
                 '--detail', detail,
                ]
    if not pats:
        call_list.append('--nopatients')

    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_quality_along_read(seq_run, VERBOSE=0, maxreads=-1, savefig=True):
    '''Submit quality check along read to the cluster'''
    if VERBOSE:
        print 'Forking to the cluster'

    JOBSCRIPT = JOBDIR+'sequencing/check_quality_along_read.py'
    cluster_time = '00:59:59'
    vmem = '1G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'quaalo',
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--verbose', VERBOSE,
                 '--maxreads', maxreads,
                ]
    if not savefig:
        call_list.append('--no-savefig')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_demultiplex(seq_run, VERBOSE=0, maxreads=-1, summary=True):
    '''Submit demultiplex script to the cluster'''
    if VERBOSE:
        print 'Forking to the cluster'

    JOBSCRIPT = JOBDIR+'demultiplex.py'
    cluster_times = ['0:59:59', '23:59:59']
    vmem = '8G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'demux',
                 '-l', 'h_rt='+cluster_times[not (0 < maxreads < 1e6)],
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--verbose', VERBOSE,
                 '--maxreads', maxreads,
                ]
    if not summary:
        call_list.append('--no-summary')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_trim(seq_run, adaID, VERBOSE=0, summary=True):
    '''Submit trim and divide script to the cluster for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID

    JOBSCRIPT = JOBDIR+'sequencing/trim_reads_lowq.py'
    cluster_time = '3:59:59'
    vmem = '1G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'pm '+adaID,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--adaIDs', adaID,
                 '--verbose', VERBOSE,
                ]
    if not summary:
        call_list.append('--no-summary')
    call_list = map(str, call_list)
    if VERBOSE >= 2:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_premap(seq_run, adaID, VERBOSE=0, threads=1, maxreads=-1,
                subsrate=0.05, gapopen=40, gapextend=3,
                reference='HXB2', summary=True, trimmed=False):
    '''Submit premap script to the cluster for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID

    JOBSCRIPT = JOBDIR+'sequencing/premap_to_reference.py'
    # It is hard to tell whether 1h is sufficient, because the final sorting takes
    # quite some time. So for now give up and require 2h.
    cluster_time = ['71:59:59', '3:59:59']
    vmem = '8G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'pm '+adaID,
                 '-l', 'h_rt='+cluster_time[threads >= 30],
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--adaIDs', adaID,
                 '--verbose', VERBOSE,
                 '--threads', threads,
                 '--reference', reference,
                 '--subsrate', subsrate,
                 '--gapopen', gapopen,
                 '--gapextend', gapextend,
                ]
    if not summary:
        call_list.append('--no-summary')
    if maxreads != -1:
        call_list.extend(['--maxreads', maxreads])
    if trimmed:
        call_list.append('--trimmed')
        
    call_list = map(str, call_list)
    if VERBOSE >= 2:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_premapped_coverage(samplename, VERBOSE=0, maxreads=-1):
    '''Submit premap coverage to the cluster'''
    if VERBOSE:
        print 'Forking to the cluster'

    JOBSCRIPT = JOBDIR+'sequencing/check_premapped_coverage.py'
    cluster_time = '00:59:59'
    vmem = '1G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'pcov '+samplename,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--samples', samplename,
                 '--verbose', VERBOSE,
                 '--maxreads', maxreads,
                 '--persist',
                ]
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_trim_and_divide(seq_run, adaID, VERBOSE=0, maxreads=-1, minisize=100,
                         summary=True):
    '''Submit trim and divide script to the cluster for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID

    JOBSCRIPT = JOBDIR+'sequencing/trim_and_divide.py'
    cluster_time = '2:59:59'
    vmem = '1G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'trdv '+adaID,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--adaIDs', adaID,
                 '--verbose', VERBOSE,
                 '--maxreads', maxreads,
                 '--minisize', minisize,
                ]
    if not summary:
        call_list.append('--no-summary')
    call_list = map(str, call_list)
    if VERBOSE >= 2:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_build_consensus_iterative(seq_run, adaID, fragment, n_reads=1000,
                                   iterations_max=0, VERBOSE=0, summary=True):
    '''Submit build consensus script to the cluster for each adapter ID and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'build_consensus_iterative.py'
    cluster_time = '0:59:59'
    vmem = '2G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'cb '+adaID+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '-n', n_reads,
                 '--iterations', iterations_max,
                 '--verbose', VERBOSE,
                ]
    if not summary:
        call_list.append('--no-summary')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_build_consensus(seq_run, adaID, fragment,
                         block_len_initial=100, n_reads_per_ali=31,
                         store_allele_counts=False, VERBOSE=0,
                         summary=True):
    '''Submit build consensus script to the cluster for each adapter ID and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'sequencing/build_consensus.py'
    if fragment == 'genomewide':
        cluster_time = '23:59:59'
    else:
        cluster_time = '0:59:59'

    vmem = '2G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'c '+adaID+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--block-length', block_len_initial,
                 '--reads-per-alignment', n_reads_per_ali,
                 '--verbose', VERBOSE,
                ]
    if not summary:
        call_list.append('--no-summary')
    if store_allele_counts:
        call_list.append('--allele-counts')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_map_to_consensus(seq_run, adaID, fragment, VERBOSE=3,
                          threads=1, maxreads=-1, filter_reads=False,
                          summary=True, rescue=False):
    '''Submit map script for each adapter ID and fragment
    
    Note on cluster runtime: we require less than 1 hr ONLY for tests, i.e. if
    maxreads is less or equal 10k. There is too much stochasticity in the cluster
    occupancy for anything else: the master thread is going to be cut down prematurely.
    '''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'sequencing/map_to_consensus.py'
    cluster_time = ['23:59:59', '0:59:59']
    vmem = '8G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'm '+adaID+' '+fragment,
                 '-l', 'h_rt='+cluster_time[0 < maxreads <= 10000],
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                 '--threads', threads,
                 '--maxreads', maxreads
                ]
    if filter_reads:
        call_list.append('--filter')
    if not summary:
        call_list.append('--no-summary')
    if rescue:
        call_list.append('--rescue')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_filter_mapped(seq_run, adaID, fragment, VERBOSE=0, maxreads=-1,
                       max_mismatches=30, susp_mismatches=25, summary=True):
    '''Submit filter script for each adapter ID and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'sequencing/filter_mapped_reads.py'
    cluster_time = '71:59:59'
    vmem = '2G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'f '+adaID+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                 '--max-mismatches', max_mismatches,
                 '--suspicious-mismatches', susp_mismatches,
                ]
    if not summary:
        call_list.append('--no-summary')
    if maxreads > 0:
        call_list.extend(['--maxreads', maxreads])
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_get_allele_counts(seq_run, adaID, fragment, VERBOSE=3):
    '''Submit get allele counts script to the cluster for each adapter ID and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'sequencing/get_allele_counts.py'
    cluster_time = '0:59:59'
    vmem = '4G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'acn '+adaID+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                ]
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_filter_allele_frequencies(seq_run, adaID, fragment, VERBOSE=3, summary=True):
    '''Submit filter allele frequencies to the cluster for each adapter ID and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'sequencing/filter_allele_frequencies.py'
    cluster_time = '0:59:59'
    vmem = '2G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'faf '+adaID+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                ]
    if not summary:
        call_list.append('--no-summary')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_extract_mutations(seq_run, adaID, VERBOSE=0, summary=True):
    '''Submit extract mutation to the cluster for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID

    JOBSCRIPT = JOBDIR+'sequencing/extract_mutations.py'
    cluster_time = '0:59:59'
    vmem = '8G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'exm_'+adaID,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--adaID', adaID,
                ]
    if not summary:
        call_list.append('--no-summary')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_get_coallele_counts(data_folder, adaID, fragment, VERBOSE=3, summary=True):
    '''Fork coallele counts for each adapter ID and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'sequencing/get_coallele_counts.py'
    cluster_time = '0:59:59'
    vmem = '8G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'ca '+adaID+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                ]
    if not summary:
        call_list.append('--no-summary')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_split_for_mapping(seq_run, adaID, fragment, VERBOSE=0, maxreads=-1, chunk_size=10000):
    '''Fork split reads for mapping for each adapter ID and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'sequencing/split_reads_for_mapping.py'
    cluster_time = '0:59:59'
    vmem = '8G'
    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'sfm '+adaID+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                 '--maxreads', maxreads,
                 '--chunksize', chunk_size,
                ]
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


# PHIX
def fork_get_allele_counts_phix(seq_run, maxreads=-1, VERBOSE=0, qual_min=None):
    '''Fork self for each adapter ID'''
    JOBSCRIPT = JOBDIR+'phix/get_allele_counts_phiX.py'
    cluster_time = '0:59:59'
    vmem = '8G'

    if VERBOSE:
        print 'Forking to the cluster'

    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'acpX'+seq_run,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--maxreads', maxreads,
                 '--verbose', VERBOSE,
                ]
    if qual_min is not None:
        call_list.extend(['--qual_min', qual_min])
    call_list = map(str, call_list)
    if VERBOSE >= 2:
        print ' '.join(call_list)
    sp.call(call_list)


# PATIENTS
def fork_map_to_initial_reference(samplename, fragment,
                                  VERBOSE=0, threads=1,
                                  n_pairs=-1, filter_reads=False,
                                  summary=True,
                                  only_chunks=[None],
                                  filtered=True):
    '''Fork to the cluster for each sample and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: sample '+samplename+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'store/map_to_initial_reference.py'
    cluster_time = ['23:59:59', '0:59:59']
    vmem = '8G'

    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'm '+samplename[:4]+' '+fragment,
                 '-l', 'h_rt='+cluster_time[(only_chunks != [None]) or (0 < n_pairs <= 10000)],
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--samples', samplename,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                 '--threads', threads,
                 '--maxreads', n_pairs,
                 '--skiphash',
                ]
    if filter_reads:
        call_list.append('--filter')
    if not summary:
        call_list.append('--no-summary')
    if only_chunks != [None]:
        call_list = call_list + ['--chunks'] + only_chunks
    if not filtered:
        call_list.append('--unfiltered')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_filter_mapped_init(samplename, fragment,
                            VERBOSE=0, n_pairs=-1,
                            PCR=1,
                            summary=True):
    '''Fork to the cluster for each sample and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: sample '+samplename+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'store/filter_mapped_reads.py'
    cluster_time = '23:59:59'
    vmem = '8G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'fmi '+samplename+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--samples', samplename,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                 '--maxreads', n_pairs,
                 '--PCR', PCR,
                ]
    if not summary:
        call_list.append('--no-summary')
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)


def fork_build_consensus_patient(samplename_pat, fragment, VERBOSE=0, PCR=1,
                                 block_len=100, n_reads_per_ali=31):
    '''Fork to the cluster for each sample and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: '+samplename_pat+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'store/store_consensus.py'
    cluster_time = '0:59:59'
    vmem = '2G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'ci'+fragment+samplename_pat,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--samples', samplename_pat,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                 '--PCR', PCR,
                 '--block-length', block_len,
                 '--reads-per-alignment', n_reads_per_ali,
                 '--save',
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)


def fork_get_allele_counts_patient(samplename, fragment, VERBOSE=0, PCR=1, qual_min=30):
    '''Fork to the cluster for each sample and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: '+samplename+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'store/store_allele_counts.py'
    cluster_time = '0:59:59'
    vmem = '2G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'ac'+fragment+samplename,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--samples', samplename,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                 '--save',
                 '--qualmin', qual_min,
                 '--PCR', PCR,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)


def fork_get_allele_counts_aa_patient(samplename, protein, VERBOSE=0, PCR=1, qual_min=30):
    '''Fork to the cluster for each sample and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: '+samplename+', protein '+protein

    JOBSCRIPT = JOBDIR+'store/store_allele_counts_aa.py'
    cluster_time = '0:59:59'
    vmem = '2G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'aac'+protein+samplename,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--samples', samplename,
                 '--proteins', protein,
                 '--verbose', VERBOSE,
                 '--save',
                 '--qualmin', qual_min,
                 '--PCR', PCR,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)


def fork_get_cocounts_patient(samplename, fragment, VERBOSE=0,
                              PCR=1, qual_min=30,
                              maxreads=-1, use_tests=False):
    '''Fork to the cluster for each patient, sample, and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: sample '+samplename+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'store/store_allele_cocounts.py'
    cluster_time = '23:59:59'
    vmem = '8G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'cc'+fragment+samplename,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--fragments', fragment,
                 '--samples', samplename,
                 '--verbose', VERBOSE,
                 '--maxreads', maxreads,
                 '--save',
                 '--qualmin', qual_min,
                 '--PCR', PCR,
                ]
    if use_tests:
        qsub_list.append('--tests')
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)


def fork_compress_cocounts_patient(samplename, fragment, VERBOSE=0,
                                   PCR=1, qual_min=30):
    '''Fork to the cluster for each patient, sample, and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: sample '+samplename+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'store/compress_allele_cocounts.py'
    cluster_time = '23:59:59'
    vmem = '8G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'ccc'+fragment+samplename,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--fragments', fragment,
                 '--samples', samplename,
                 '--verbose', VERBOSE,
                 '--qualmin', qual_min,
                 '--PCR', PCR,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)



def fork_decontaminate_reads_patient(samplename, fragment, VERBOSE=0, PCR=None,
                                     maxreads=-1, summary=True):
    '''Fork to the cluster the decontamination of reads'''
    if VERBOSE:
        print 'Fork to cluster: sample', samplename, fragment

    JOBSCRIPT = JOBDIR+'store/decontaminate_reads.py'
    cluster_time = '71:59:59'
    vmem = '2G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'de'+fragment+samplename,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--samples', samplename,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                ]
    if PCR is not None:
        qsub_list.extend(['--PCR', PCR])
    if maxreads != -1:
        qsub_list.extend(['--maxreads', maxreads])
    if not summary:
        qsub_list.append('--no-summary')
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)


def fork_get_allele_frequency_trajectory(pname, fragment, VERBOSE=0):
    '''Fork to the cluster for each sample and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: patient '+pname+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'patients/get_allele_frequency_trajectories.py'
    cluster_time = '0:59:59'
    vmem = '2G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'aft '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--patient', pname,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                 '--save',
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)


def fork_store_haplotypes_scan(pname, width, gap, start, end, VERBOSE=0,
                               freqmin=0.01, countmin=3):
    '''Fork to the cluster for each patient'''
    if VERBOSE:
        print 'Forking to the cluster: patient '+pname
    JOBSCRIPT = JOBDIR+'store/store_haplotypes_scan.py'
    cluster_time = '23:59:59'
    vmem = '2G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'scan '+pname,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--patients', pname,
                 '--width', width,
                 '--gap', gap,
                 '--start', start,
                 '--end', end,
                 '--freqmin', freqmin,
                 '--countmin', countmin,
                 '--verbose', VERBOSE,
                 '--save',
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)


# WEBSITE
def fork_store_haplotypes_website(pname, region, VERBOSE=0):
    '''Fork to the cluster for each patient and region'''
    if VERBOSE:
        print 'Forking to the cluster: patient '+pname+', region '+region

    JOBSCRIPT = JOBDIR+'website/store_haplotype_alignments_trees.py'
    cluster_time = '23:59:59'
    vmem = '8G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'spa'+pname+region,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--patients', pname,
                 '--regions', region,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)



# THEORY
def fork_calculate_beta_SFS(alpha, N, VERBOSE=0):
    '''Fork to the cluster in parallel for each alpha and N'''
    if VERBOSE:
        print 'Forking to the cluster: alpha '+str(alpha)+', N '+str(N)

    JOBSCRIPT = JOBDIR+'theory/build_betatree_sfs.py'
    cluster_time = '23:59:59'
    vmem = '4G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'SFS '+str(alpha)+' '+str(N),
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--alphas', alpha,
                 '--Ns', N,
                 '--verbose', VERBOSE,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)


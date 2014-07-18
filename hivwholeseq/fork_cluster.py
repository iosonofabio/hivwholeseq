# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/12/13
content:    Module with all submit functions for the cluster. With this we can
            keep all cluster-specific code in one place.
'''
# Globals
import subprocess as sp

import hivwholeseq
JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
JOBLOGOUT = JOBDIR+'logout/'
JOBLOGERR = JOBDIR+'logerr/'



# Functions
def empty_log_folders():
    '''Empty log folders of old files'''
    import shutil, os

    shutil.rmtree(JOBLOGOUT)
    os.mkdir(JOBLOGOUT)

    shutil.rmtree(JOBLOGERR)
    os.mkdir(JOBLOGERR)


def fork_quality_along_read(seq_run, VERBOSE=0, maxreads=-1, savefig=True):
    '''Submit quality check along read to the cluster'''
    if VERBOSE:
        print 'Forking to the cluster'

    JOBSCRIPT = JOBDIR+'quality_along_read.py'
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
    '''Submit premap script to the cluster for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID

    JOBSCRIPT = JOBDIR+'trim_reads_lowq.py'
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


def fork_premap(seq_run, adaID, VERBOSE=0, threads=1,
                reference='HXB2', summary=True, trimmed=False):
    '''Submit premap script to the cluster for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID

    JOBSCRIPT = JOBDIR+'premap_to_reference.py'
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
                ]
    if not summary:
        call_list.append('--no-summary')
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

    JOBSCRIPT = JOBDIR+'check_premapped_coverage.py'
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

    JOBSCRIPT = JOBDIR+'trim_and_divide.py'
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

    JOBSCRIPT = JOBDIR+'build_consensus.py'
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
                          summary=True):
    '''Submit map script for each adapter ID and fragment
    
    Note on cluster runtime: we require less than 1 hr ONLY for tests, i.e. if
    maxreads is less or equal 10k. There is too much stochasticity in the cluster
    occupancy for anything else: the master thread is going to be cut down prematurely.
    '''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'map_to_consensus.py'
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
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_filter_mapped(seq_run, adaID, fragment, VERBOSE=0, summary=True):
    '''Submit filter script for each adapter ID and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'filter_mapped_reads.py'
    cluster_time = '0:59:59'
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
                ]
    if not summary:
        call_list.append('--no-summary')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_get_allele_counts(seq_run, adaID, fragment, VERBOSE=3):
    '''Submit get allele counts script to the cluster for each adapter ID and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+adaID+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'get_allele_counts.py'
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

    JOBSCRIPT = JOBDIR+'filter_allele_frequencies.py'
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

    JOBSCRIPT = JOBDIR+'extract_mutations.py'
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

    JOBSCRIPT = JOBDIR+'get_coallele_counts.py'
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

    JOBSCRIPT = JOBDIR+'split_reads_for_mapping.py'
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



# PATIENTS
def fork_map_to_initial_consensus(samplename, fragment,
                                  VERBOSE=0, threads=1,
                                  n_pairs=-1, filter_reads=False,
                                  summary=True,
                                  only_chunks=[None]):
    '''Fork to the cluster for each sample and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: sample '+samplename+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'patients/map_to_initial_consensus.py'
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
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_paste_mapped_chunks_to_initial_consensus(pname, samplename, fragment,
                                                  VERBOSE=0, filter_reads=False,
                                                  summary=True):
    '''Fork to the cluster for each sample and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: patient '+pname+', sample '+\
                samplename+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'patients/paste_mapped_chunks.py'
    cluster_time = '0:59:59'
    vmem = '8G'

    call_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'pc '+samplename+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--patient', pname,
                 '--samples', samplename,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                ]
    if not summary:
        call_list.append('--no-summary')
    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    return sp.check_output(call_list)


def fork_filter_mapped_init(samplename, fragment,
                            VERBOSE=0, n_pairs=-1,
                            summary=True):
    '''Fork to the cluster for each sample and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: sample '+samplename+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'patients/filter_mapped_reads.py'
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
                ]
    if not summary:
        call_list.append('--no-summary')
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)


def fork_get_allele_counts_patient(samplename_pat, fragment, VERBOSE=0, qual_min=30):
    '''Fork to the cluster for each sample and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: '+samplename_pat+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'patients/get_allele_counts.py'
    cluster_time = '0:59:59'
    vmem = '2G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'ac '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--samples', samplename_pat,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                 '--save',
                 '--qualmin', qual_min,
                ]
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


def fork_get_cocounts_patient(pname, samplename, fragment, VERBOSE=0,
                              maxreads=-1, use_tests=False):
    '''Fork to the cluster for each patient, sample, and fragment'''
    if VERBOSE:
        print 'Forking to the cluster: patient '+pname+', sample'+samplename+', fragment '+fragment

    JOBSCRIPT = JOBDIR+'patients/get_allele_cocounts.py'
    cluster_time = '23:59:59'
    vmem = '8G'

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'coco '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--patient', pname,
                 '--fragments', fragment,
                 '--samples', samplename,
                 '--verbose', VERBOSE,
                 '--maxreads', maxreads,
                ]
    if use_tests:
        qsub_list.append('--tests')
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    return sp.check_output(qsub_list)


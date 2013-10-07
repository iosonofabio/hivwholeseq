# vim: fdm=marker
'''
author:     Fabio Zanini
date:       12/09/13
content:    Build consensus after map to HXB2 via samtools. Note: this requires
            the mapped reads to HXB2, and does not handle inserts.
'''
import argparse
import subprocess as sp

from mapping.datasets import MiSeq_runs
from mapping.adapter_info import load_adapter_table, foldername_adapter
from mapping.filenames import get_HXB2_fragmented


# Globals
# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'build_consensus_samtools.py'
cluster_time = '0:59:59'
vmem = '8G'


# Functions
def fork_self(miseq_run, adaID, fragment, VERBOSE=0):
    '''Fork self for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+'{:02d}'.format(adaID)

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'asb '+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', miseq_run,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE >= 2:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def sort_bam_file(data_folder, adaID, fragment, VERBOSE=0):
    '''Sort reads in BAM file'''
    mapped_reads_filename = data_folder+'subsample/'+foldername_adapter(adaID)+\
            'map_iter/mapped_to_HXB2_'+fragment+'.bam'
    sorted_reads_filename = data_folder+'subsample/'+foldername_adapter(adaID)+\
            'map_iter/mapped_to_HXB2_'+fragment+'_sorted'

    call_list = ['samtools', 'sort',
                 mapped_reads_filename, sorted_reads_filename]
    if VERBOSE >= 2:
        print ' '.join(call_list)
    sp.call(call_list)


def build_consensus(data_folder, adaID, fragment, VERBOSE=0):
    '''Assemble reads into a consensus'''
    # Input files
    sorted_reads_filename = data_folder+'subsample/'+foldername_adapter(adaID)+\
            'map_iter/mapped_to_HXB2_'+fragment+'_sorted.bam'

    # Output files
    output_filename = data_folder+'subsample/'+foldername_adapter(adaID)+\
            'map_iter/consensus_0_'+fragment+'_samtools.fastq'

    # Call samtools
    with open(output_filename, 'w') as fo:
        p1 = sp.Popen(['samtools',
                       'mpileup', '-u',
                       '-f', get_HXB2_fragmented(fragment),
                       sorted_reads_filename],               
                      stdout=sp.PIPE)
        p2 = sp.Popen(['bcftools',
                       'view', '-cg', '-'], stdin=p1.stdout, stdout=sp.PIPE)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        p3 = sp.Popen(['/ebio/ag-neher/share/programs/bundles/samtools-0.1.19/bcftools/vcfutils.pl',
                       'vcf2fq'], stdin=p2.stdout, stdout=fo)
        p2.stdout.close()  # Allow p2 to receive a SIGPIPE if p3 exits.

        # Wait until the end of the consensus building before regaining control
        p3.wait()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Assemble consensus de novo')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', type=int, nargs='+',
                        help='Adapter ID')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
 
    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    submit = args.submit

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over all requested samples
    for adaID in adaIDs:
        for fragment in fragments:

            # Submit to the cluster self if requested
            if submit:
                fork_self(miseq_run, adaID, fragment, VERBOSE=VERBOSE)
                continue

            # sort BAM file and build consensus
            sort_bam_file(data_folder, adaID, fragment, VERBOSE=VERBOSE)
            build_consensus(data_folder, adaID, fragment, VERBOSE=VERBOSE)

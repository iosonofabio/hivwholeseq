# vim: fdm=marker
'''
author:     Fabio Zanini
date:       15/06/14
content:    Check the status of the pipeline for one or more sequencing samples.
'''
# Modules
import os
import sys
import argparse

from hivwholeseq.generic_utils import getchar
from hivwholeseq.samples import SampleSeq, load_sequencing_run
from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.samples import load_samples_sequenced as lss


# Globals
len_fr = 8
len_msg = 5
spacing_fragments = 4



# Functions
def print_info(name, gfn, frs=None):
    '''Print info on these files'''
    print '{:<20s}'.format(name+':'),

    if frs is None:
        fn = gfn()
        if os.path.isfile(fn):
            print 'OK'
        else:
            print 'MISS'

    else:
        msg = []
        for fr in frs:
            ms = ('{:<'+str(len_fr)+'s}').format(fr+':')
            fn = gfn(fr)
            if os.path.isfile(fn):
                msg.append(ms+('{:>'+str(len_msg)+'}').format('OK'))
            else:
                msg.append(ms+('{:>'+str(len_msg)+'}').format('MISS'))
        print (' ' * spacing_fragments).join(msg)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Check sequencing run for missing parts of the analysis',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--runs', required=True, nargs='+',
                        help='Seq runs to analyze (e.g. Tue28, test_tiny)')
    parser.add_argument('--adaIDs', nargs='+',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--nopatients', action='store_false', dest='use_pats',
                        help='Include non-patient samples (e.g. reference strains)')
    parser.add_argument('--interactive', action='store_true',
                        help='Interactive mode')
    
    args = parser.parse_args()
    seq_runs = args.runs
    adaIDs = args.adaIDs
    use_pats = args.use_pats
    use_interactive = args.interactive

    samples_pat = lssp()
    samples = lss()

    samples = samples.loc[samples['seq run'].isin(seq_runs)]
    if adaIDs is not None:
        samples = samples.loc[samples.adapter.isin(adaIDs)]

    if len(seq_runs) >= 2:
        samples.sort(columns=['patient sample', 'seq run'], inplace=True)

    for samplename, sample in samples.iterrows():
        sample = SampleSeq(sample)
        print sample.name, 'seq:', sample['seq run'], sample.adapter,
        if sample['patient sample'] == 'nan':
            print 'not a patient sample',
            if use_pats:
                print '(skip)'
                continue
            else:
                print ''
        else:
            sample_pat = samples_pat.loc[sample['patient sample']]
            print 'patient: '+sample_pat.patient

        print_info('Premapped', sample.get_premapped_filename)
        print_info('Divided', sample.get_divided_filename, sample.regions_complete)
        print_info('Consensus', sample.get_consensus_filename, sample.regions_generic)
        print_info('Mapped', lambda x: sample.get_mapped_filename(x, filtered=False), sample.regions_generic)
        print_info('Filtered', lambda x: sample.get_mapped_filename(x, filtered=True), sample.regions_generic)

        print ''

        if use_interactive:
            print 'Press q to exit',
            sys.stdout.flush()
            ch = getchar()
            if ch.lower() in ['q']:
                print 'stopped'
                break
            else:
                sys.stdout.write("\x1b[1A")
                print ''

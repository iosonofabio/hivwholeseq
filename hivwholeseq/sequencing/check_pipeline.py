#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       15/06/14
content:    Check the status of the pipeline for one or more sequencing samples.
'''
# Modules
import os
import sys
from itertools import izip
import argparse
from Bio import SeqIO

from hivwholeseq.generic_utils import getchar
from hivwholeseq.sequencing.samples import SampleSeq, load_sequencing_run
from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.patients.patients import SamplePat
from hivwholeseq.sequencing.samples import load_samples_sequenced as lss
from hivwholeseq.mapping_utils import get_number_reads
from hivwholeseq.fork_cluster import fork_check_pipeline as fork_self


# Globals
len_fr = 8
len_msg = 6
spacing_fragments = 4



# Functions
def check_status(sample, step, detail=1):
    '''Check for a sample a certain step of the pipeline at a certain detail'''
    if detail == 1:
        if step == 'premapped':
            return [os.path.isfile(sample.get_premapped_filename())]
        elif step == 'divided':
            return [(fr, os.path.isfile(sample.get_divided_filename(fr)))
                    for fr in sample.regions_complete]
        elif step == 'consensus':
            return [(fr, os.path.isfile(sample.get_consensus_filename(fr)))
                    for fr in sample.regions_generic]
        elif step == 'mapped':
            return [(fr, os.path.isfile(sample.get_mapped_filename(fr, filtered=False)))
                    for fr in sample.regions_generic]
        elif step == 'filtered':
            return [(fr, os.path.isfile(sample.get_mapped_filename(fr, filtered=True)))
                    for fr in sample.regions_generic]
        elif step == 'mapped_initial':
            return [(fr, os.path.isfile(sample.get_mapped_to_initial_filename(fr)))
                    for fr in sample.regions_generic]
        elif step == 'mapped_filtered':
            # Check whether the mapped filtered is older than the mapped_initial
            from hivwholeseq.generic_utils import modification_date
            out = []
            for fr in sample.regions_generic:
                fn_mi = sample.get_mapped_to_initial_filename(fr)
                fn_mf = sample.get_mapped_filtered_filename(fr)
                if not os.path.isfile(fn_mf):
                    out.append((fr, False))
                    continue

                if not os.path.isfile(fn_mi):
                    out.append((fr, True))
                    continue

                md_mi = modification_date(fn_mi)
                md_mf = modification_date(fn_mf)
                if md_mf < md_mi:
                    out.append((fr, 'OLD'))
                else:
                    out.append((fr, True))
            return out

    elif detail == 2:
        if step in ('filtered', 'consensus'):
            return check_status(sample, step, detail=3)
        else:
            return check_status(sample, step, detail=1)

    elif detail == 3:
        if step == 'premapped':
            if os.path.isfile(sample.get_premapped_filename()):
                return [get_number_reads(sample.get_premapped_filename())]
            else:
                return [False]

        elif step == 'divided':
            stati = []
            for fr in sample.regions_complete:
                fn = sample.get_divided_filename(fr)
                if os.path.isfile(fn):
                    status = (fr, get_number_reads(fn))
                else:
                    status = (fr, False)
                stati.append(status)
            return stati

        elif step == 'consensus':
            stati = []
            for fr in sample.regions_generic:
                fn = sample.get_consensus_filename(fr)
                if os.path.isfile(fn):
                    status = (fr, len(SeqIO.read(fn, 'fasta')))
                else:
                    status = (fr, False)
                stati.append(status)
            return stati

        elif step == 'mapped':
            stati = []
            for fr in sample.regions_generic:
                fn = sample.get_mapped_filename(fr, filtered=False)
                if os.path.isfile(fn):
                    status = (fr, get_number_reads(fn))
                else:
                    status = (fr, False)
                stati.append(status)
            return stati

        elif step == 'filtered':
            stati = []
            for fr in sample.regions_generic:
                fn = sample.get_mapped_filename(fr, filtered=True)
                if os.path.isfile(fn):
                    status = (fr, get_number_reads(fn))
                else:
                    status = (fr, False)
                stati.append(status)
            return stati

        # TODO: add mapped_to_initial and downstream
        elif step in ('mapped_initial', 'mapped_filtered'):
            return check_status(sample, step, detail=1)


def print_info(name, status, detail=1):
    '''Print info on these files'''
    print '{:<20s}'.format(name+':'),

    if name.lower() in ['premapped']:
        status = status[0]
        if status == True:
            print 'OK'
        elif status == False:
            print 'MISS'
        else:
            print str(status)

    else:
        stati = list(status)
        msg = []
        for (fr, status) in stati:
            ms = ('{:<'+str(len_fr)+'s}').format(fr+':')
            if status == True:
                msg.append(ms+('{:>'+str(len_msg)+'}').format('OK'))
            elif status == False:
                msg.append(ms+('{:>'+str(len_msg)+'}').format('MISS'))
            else:
                msg.append(ms+('{:>'+str(len_msg)+'}').format(str(status)))
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
    parser.add_argument('--detail', type=int, default=1,
                        help='Include details on number of reads, length of consensus')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    
    args = parser.parse_args()
    seq_runs = args.runs
    adaIDs = args.adaIDs
    use_pats = args.use_pats
    use_interactive = args.interactive
    detail = args.detail
    submit = args.submit

    if submit:
        fork_self(seq_runs, adaIDs=adaIDs,
                  pats=use_pats,
                  detail=detail)
        sys.exit()

    samples_pat = lssp(include_wrong=True)
    samples = lss()

    samples = samples.loc[samples['seq run'].isin(seq_runs)]
    if adaIDs is not None:
        samples = samples.loc[samples.adapter.isin(adaIDs)]

    if len(seq_runs) >= 2:
        samples.sort(columns=['patient sample', 'seq run'], inplace=True)

    for isa, (samplename, sample) in enumerate(samples.iterrows()):
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

        steps = ['premapped', 'divided', 'consensus', 'mapped', 'filtered',
                 'mapped_initial', 'mapped_filtered']
        for step in steps:
            status = check_status(sample, step, detail=detail)
            print_info(step.capitalize(), status, detail=detail)

        if (isa != len(samples) - 1):
            print ''

        if use_interactive and (isa != len(samples) - 1):
            print 'Press q to exit',
            sys.stdout.flush()
            ch = getchar()
            if ch.lower() in ['q']:
                print 'stopped'
                break
            else:
                sys.stdout.write("\x1b[1A")
                print ''

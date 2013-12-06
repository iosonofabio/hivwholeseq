# vim: fdm=marker
'''
author:     Fabio Zanini
date:       14/10/13
content:    Prepare reference files for NL4-3 and similia: split fragments, etc.
'''
# Modules
import argparse
from itertools import izip
import numpy as np
from Bio import SeqIO
from Bio import pairwise2

from mapping.datasets import MiSeq_runs
from mapping.miseq import alpha
from mapping.filenames import get_NL43_entire, get_NL43_fragmented, \
        get_F10_entire, get_F10_fragmented
from mapping.primer_info import primers_inner as pr



# Globals
#             run adaID  sample  ref_function     fragmented_function
references = {(28, 2): ('NL4-3', get_NL43_entire, get_NL43_fragmented),
              (28, 7): ('F10', get_F10_entire, get_F10_fragmented)}



# Functions
def find_fragments(miseq_run, adaID, trim_primers=True):
    '''Find the fragments'''
    # Get the annotated sequence
    sample, get_reference_filename, _ = references[(miseq_run, adaID)]
    refseq = SeqIO.read(get_reference_filename(), 'fasta')

    # Look for the primers
    # Exact matches will not work because of ambiguous primers
    # Make 30 overlapping windows and search in those
    refseq_frag = {}
    n_windows = 30
    windows = [(max(0, pos - 20), refseq[max(0, pos - 20): \
                           min(len(refseq), pos + len(refseq) / n_windows + 20)])
               for pos in len(refseq) / n_windows * np.arange(n_windows)]

    for (fragment, prfs) in pr.iteritems():

        poss = {}
        for prf, key in izip(prfs, ('fwd', 'rev')):
            prfl = list(prf)
            score = 2 * sum(1 for c in prfl if c in alpha[:4]) \
                    -0.5 * sum(1 for c in prfl if c not in alpha[:4])
            
            for pos, win in windows:
                a = pairwise2.align.localms(str(win.seq), prf, 2, -0.5, -2.5, -0.1)[0]
                posp = pos + len(a[1]) - len(a[1].lstrip('-'))
                
                # The F5 primers are a bit fuzzy
                if 'F5' in fragment:
                    score_gap = 10
                else:
                    score_gap = 5
                # LTR crap: F6 rev is found at the beginning
                if (a[2] >= score - score_gap) and ((not poss) or (posp > poss['fwd'])):
                    poss[key] = posp
                    break

        # SF162 is only a partial reference, hence it's a mess
        if len(poss) == 2:
            if trim_primers:
                rsf = refseq[poss['fwd'] + len(prfs[0]): poss['rev']]
            else:
                rsf = refseq[poss['fwd']: poss['rev'] + len(prfs[1])]
            rsf.name = rsf.name+'_fragment_'+fragment
            refseq_frag[fragment] = rsf

        elif len(poss) == 1:
            # If we have only rev (e.g. F10 F1), take everything to rev
            if 'rev' in poss:
                if trim_primers:
                    rsf = refseq[:poss['rev']]
                else:
                    rsf = refseq[:poss['rev'] + len(prfs[1])]
                rsf.name = rsf.name+'_fragment_'+fragment
                refseq_frag[fragment] = rsf

            else:
                # Only fwd is unheard of
                import ipdb; ipdb.set_trace()
        else:
            # No primer match is unheard of
            import ipdb; ipdb.set_trace()

    return refseq_frag


def write_output_fragment(miseq_run, adaID, fragment, seq, trim_primers=True, VERBOSE=0):
    '''Write the output of a fragment'''
    get_fragmented_filename = references[(miseq_run, adaID)][2]
    filename = get_fragmented_filename(fragment, trim_primers=trim_primers)
    SeqIO.write(seq, filename, 'fasta')
    if VERBOSE:
        print 'File written:', filename



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Prepare reference sequences')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--trim_primers', action='store_true',
                        help='Trim the inner PCR primers')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    VERBOSE = args.verbose
    trim_primers = args.trim_primers

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Iterate over adapter IDs
    for adaID in adaIDs:
        refseq_frag = find_fragments(miseq_run, adaID, trim_primers=trim_primers) 

        sample = references[(miseq_run, adaID)][0]
        for fragment, seq in refseq_frag.iteritems():
            if 1000 < len(seq) < 2000:
                write_output_fragment(miseq_run, adaID, fragment, seq, VERBOSE=VERBOSE,
                                      trim_primers=trim_primers)
                if VERBOSE:
                    print fragment, len(seq)
            else:
                if VERBOSE:
                    print 'Fragment', fragment, 'not found!'



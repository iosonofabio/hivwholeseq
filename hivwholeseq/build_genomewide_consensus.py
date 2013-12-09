# vim: fdm=marker
'''
author:     Fabio Zanini
date:       24/10/13
content:    Build genome wide consensus from the 6 fragments.
'''
# Modules
import argparse
from itertools import izip
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as SeqIO

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.filenames import get_consensus_filename, get_merged_consensus_filename
from hivwholeseq.check_overlaps import get_overlapping_fragments, get_overlap, \
        check_overlap_consensus
from hivwholeseq.samples import samples



# Functions
def merge_consensi(data_folder, adaID, fragments, VERBOSE=0):
    '''Merge consensi at overlapping pairs'''
    import warnings

    consensi = {frag: SeqIO.read(get_consensus_filename(data_folder, adaID, frag,
                                             trim_primers=True), 'fasta')
                for frag in fragments}

    pairs = get_overlapping_fragments(fragments)
    overlaps = {}
    for (frag1, frag2) in pairs:
        overlap = get_overlap(data_folder, adaID, frag1, frag2, VERBOSE=VERBOSE)
        is_diff = check_overlap_consensus(data_folder, adaID, frag1, frag2,
                                          overlap, VERBOSE=VERBOSE)
        if is_diff:
            warnings.warn(frag1+' and '+frag2+' have different consensi.', RuntimeWarning)
        overlaps[(frag1, frag2)] = overlap

    consensus = []
    fragments = sorted(fragments)
    for i, frag in enumerate(fragments):
        # If the start is not an overlap, start a new consensus and copy all
        if (i == 0) or (fragments[i-1], frag) not in overlaps:
            cons = [[frag], str(consensi[frag].seq)]
            consensus.append(cons)

        # copy from the end of the overlap on
        else:
            cons = consensus[-1]
            (_, start, _) = overlaps[(fragments[i-1], frag)]
            cons[0].append(frag)
            cons[1] = cons[1]+str(consensi[frag][start:].seq)

    # Make SeqRecords out of consensi
    for i, (frags, cons) in enumerate(consensus):
        name = 'adaID_'+str(adaID)+'_'+'-'.join(frags)
        rec = SeqRecord(Seq(cons, IUPAC.ambiguous_dna),
                        id=name, name=name)
        consensus[i] = (frags, rec)

    return consensus



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Merge consensi from fragments')
    parser.add_argument('--run', required=True,
                        help='MiSeq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    
    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = dataset['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs


    # Iterate over adaIDs
    for adaID in adaIDs:

        samplename = dataset['samples'][dataset['adapters'].index(adaID)]
        fragments = [fr[:2] for fr in samples[samplename]['fragments']]

        consensus = merge_consensi(data_folder, adaID, fragments, VERBOSE=VERBOSE)

        for (frags, cons) in consensus:
            output_filename = get_merged_consensus_filename(data_folder, adaID, frags)
            SeqIO.write(cons, output_filename, 'fasta')


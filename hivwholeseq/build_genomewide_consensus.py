# vim: fdm=marker
'''
author:     Fabio Zanini
date:       24/10/13
content:    Build genome wide consensus and allele frequencies from the 6 fragments.
'''
# Modules
import argparse
from itertools import izip
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as SeqIO

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.filenames import get_consensus_filename, \
        get_allele_frequencies_filename,\
        get_merged_consensus_filename, get_merged_allele_frequencies_filename
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
            cons[0].append(frag)
            tmp = overlaps[(fragments[i-1], frag)]
            if tmp is not None:
                (_, start, _) = tmp
                cons[1] = cons[1]+str(consensi[frag][start:].seq)
            else:
                cons[1] = cons[1]+('N' * 10)+str(consensi[frag].seq)

    # Make SeqRecords out of consensi
    for i, (frags, cons) in enumerate(consensus):
        name = 'adaID_'+str(adaID)+'_'+'-'.join(frags)
        rec = SeqRecord(Seq(cons, IUPAC.ambiguous_dna),
                        id=name, name=name)
        consensus[i] = (frags, rec)

    return consensus


def merge_allele_frequencies(data_folder, adaID, fragments, VERBOSE=0):
    '''Merge allele frequencies at overlapping pairs'''
    import warnings
    import numpy as np

    consensi = {frag: SeqIO.read(get_consensus_filename(data_folder, adaID, frag,
                                             trim_primers=True), 'fasta')
                for frag in fragments}
    nus = {frag: np.load(get_allele_frequencies_filename(data_folder, adaID, frag))
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

    nu = []
    fragments = sorted(fragments)
    for i, frag in enumerate(fragments):
        # If the start is not an overlap, start a new chunk and copy all
        if (i == 0) or (fragments[i-1], frag) not in overlaps:
            nuf = [[frag], nus[frag]]
            nu.append(nuf)

        # else, copy from the end of the overlap on
        # FIXME: we could average the consensus zone out of indels...
        else:
            nuf = nu[-1]
            nuf[0].append(frag)
            tmp = overlaps[(fragments[i-1], frag)]
            if tmp is not None:
                (_, start, _) = tmp
                #(recursion is not the most efficient but -- oh, well)
                nuf[1] = np.concatenate([nuf[1], nus[frag][:, start:]], axis=1)
            else:
                tmp = np.zeros((nuf[1].shape[0], 10), float)
                tmp[-1] = 1
                nuf[1] = np.concatenate([nuf[1], tmp, nus[frag][:, start:]], axis=1)

    return nu



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Merge consensi and allele frequencies')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--nus', action='store_true',
                        help='Also create genome-wide allele frequencies')
    
    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    do_nus = args.nus
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = dataset['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    for adaID in adaIDs:

        samplename = dataset['samples'][dataset['adapters'].index(adaID)]
        fragments = [fr[:2] for fr in samples[samplename]['fragments']]

        # Write one or more merged consensi
        consensus = merge_consensi(data_folder, adaID, fragments, VERBOSE=VERBOSE)
        for (frags, cons) in consensus:
            output_filename = get_merged_consensus_filename(data_folder, adaID, frags)
            SeqIO.write(cons, output_filename, 'fasta')

        # Write allele frequencies
        if do_nus:
            nu = merge_allele_frequencies(data_folder, adaID, fragments, VERBOSE=VERBOSE)
            for (frags, nuf) in nu:
                output_filename = get_merged_allele_frequencies_filename(data_folder, adaID, frags)
                nuf.dump(output_filename)

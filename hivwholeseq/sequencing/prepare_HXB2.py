# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/08/13
content:    Chop HXB2 into the 6 fragments, which are used as "chromosomes" for
            mapping.
'''
# Modules
import os
import re
import Bio.SeqIO as SeqIO

from hivwholeseq.reference import load_HXB2
from hivwholeseq.sequencing.primer_info import primers_coordinates_HXB2_inner as pci
from hivwholeseq.sequencing.primer_info import primers_coordinates_HXB2_outer as pco
from hivwholeseq.sequencing.filenames import get_HXB2_fragmented, get_HXB2_entire



# Script
if __name__ == '__main__':

    # Make output folder if necessary
    dirname = os.path.dirname(get_HXB2_entire())
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # Get the annotated sequence
    seq = load_HXB2()

    # 1. Copy the entire reference verbatim
    SeqIO.write(seq, get_HXB2_entire(), 'fasta')

    # 2. Make a cropped sequence from F1o to F6o (outer primers), to reduce LTR
    # degeneracy problems during premapping
    start = pco['F1'][0][0]
    end = pco['F6'][1][1]
    seq_cropped = seq[start: end]
    seq_cropped.id = seq_cropped.name = seq.name+'_cropped_F1o_F6o'
    seq_cropped.description = seq.description+' (cropped to fragments F1-F6, outer primers)'
    SeqIO.write(seq_cropped, get_HXB2_entire(cropped=True), 'fasta')

    # 3. Extract fragments coordinates and sequences (this is sample specific)
    frag_coos = {f: [co[0][0], co[1][1]] for f, co in pci.iteritems()}
    frag_seqs = {f: seq[co[0]: co[1]] for f, co in frag_coos.iteritems()}
    
    # Extract fragments without the inner primers
    frag_coos_trimmed = {f: [co[0][1], co[1][0]] for f, co in pci.iteritems()}
    frag_seqs_trimmed = {f: seq[co[0]: co[1]] for f, co in frag_coos_trimmed.iteritems()}

    # Relabel sequences
    for f in frag_seqs:
        # Relabel fragments with primers
        fs = frag_seqs[f]
        fs.name = 'HXB2_reference_'+f
        fs.id = fs.name
        fs.description = re.sub('_', ' ', fs.id)

        # Relabel fragments without primers
        fs_trimmed = frag_seqs_trimmed[f]
        fs_trimmed.name = 'HXB2_reference_'+f+'_trim_primers'
        fs_trimmed.id = fs_trimmed.name
        fs_trimmed.description = re.sub('_', ' ', fs_trimmed.id)

    # Write to file both the fragmented (in separated files)
    for f in frag_seqs:
        SeqIO.write(frag_seqs[f], get_HXB2_fragmented(f), 'fasta')
        SeqIO.write(frag_seqs_trimmed[f],
                    get_HXB2_fragmented(f, trim_primers=True), 'fasta')
        print f

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

from mapping.reference import load_HXB2
from mapping.primer_info import primers_coordinates_HXB2_inner as pci
from mapping.filenames import get_HXB2_fragmented, get_HXB2_entire



# Script
if __name__ == '__main__':

    # Get the annotated sequence
    seq = load_HXB2()

    # Extract fragments coordinates and sequences (this is sample specific)
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

    # Make a cropped sequence from F1 to F6, to reduce LTR degeneracy problems
    start_F1 = frag_coos['F1'][0]
    end_F6 = frag_coos['F6'][1]
    seq_cropped = seq[start_F1: end_F6]
    seq_cropped.id = seq_cropped.name = seq.name+'_cropped_F1_F6'
    seq_cropped.description = seq.description+' (cropped to fragments F1-F6)'

    # Make output folder if necessary
    dirname = os.path.dirname(get_HXB2_entire())
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # Write to file both the entire, the cropped, and the fragmented (in separated files)
    SeqIO.write(seq, get_HXB2_entire(), 'fasta')
    SeqIO.write(seq_cropped, get_HXB2_entire(cropped=True), 'fasta')
    for f in frag_seqs:
        SeqIO.write(frag_seqs[f], get_HXB2_fragmented(f), 'fasta')
        SeqIO.write(frag_seqs_trimmed[f],
                    get_HXB2_fragmented(f, trim_primers=True), 'fasta')
        print f

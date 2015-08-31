# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/11/14
content:    Test reading SHAPE data for the NL4-3 sequence.
'''
# Globals
_seqs = None
_shape = None


# Functions
def load_SHAPE():
    '''Load the SHAPE data of a reference sequence (for now only NL4-3)'''
    global _seqs
    global _shape

    if _seqs is None:

        import numpy as np
        import pandas as pd
        from hivwholeseq.filenames import reference_folder

        fn = reference_folder+'SHAPE/nmeth.3029-S2.xlsx'
        table = pd.read_excel(fn, index_col=0, parse_cols='A:D')

        seqm = table.loc[:, 'Nt identity']
        seqs = str(''.join(seqm))
        
        shape = table.loc[:, '1M7 SHAPE MaP']
        shape = np.ma.masked_where(shape < -990, shape)

        _seqs = str(seqs)
        _shape = shape.copy()

    else:
        seqs = str(_seqs)
        shape = _shape.copy()
    
    return (seqs, shape)


def add_SHAPE_to_seqrecord(seqrecord, VERBOSE=0):
    '''Add SHAPE info to an HIV sequence (no need to be genomewide)'''
    import numpy as np
    from seqanpy import align_overlap

    global _seqs
    global _shape

    if _seqs is None:
        seqs, shape = load_SHAPE()
    else:
        # No need to copy, they are nodified within this function and no pointers
        # are transmitted outside
        seqs = _seqs
        shape = _shape

    (score, ali1, ali2) = align_overlap(seqs, seqrecord)
    if VERBOSE >= 3:
        from hivwholeseq.utils.sequence import pretty_print_pairwise_ali
        pretty_print_pairwise_ali([ali1, ali2], width=100,
                                  name1='SHAPE',
                                  name2=seqrecord.name)

    alim1 = np.fromstring(ali1, 'S1')
    alim2 = np.fromstring(ali2, 'S1')
    indshape = list((alim1 != '-').nonzero()[0])
    indout = (alim2 != '-').nonzero()[0]
    shapeout = np.ma.masked_all(len(seqrecord))
    for posout, posali in enumerate(indout):
        if alim1[posali] == alim2[posali]:
            posshape = indshape.index(posali)
            shapeout[posout] = shape[posshape]

    seqrecord.letter_annotations['SHAPE'] = shapeout



# Script
if __name__ == '__main__':

    from hivwholeseq.reference import load_custom_reference
    ref = load_custom_reference('NL4-3')
    refs = ''.join(ref)

    seqs, shape = load_SHAPE()

    from seqanpy import align_overlap
    from hivwholeseq.utils.sequence import pretty_print_pairwise_ali
    (score, ali1, ali2) = align_overlap(refs, seqs)
    pretty_print_pairwise_ali([ali1, ali2], name1='NL4-3', name2='SHAPE', width=100)
    if score == len(ali2.replace('-', '')) * 3:
        print 'The sequence is NL32, cut at both edges'

    from hivwholeseq.patients.patients import load_patient
    patient = load_patient('p1')
    refpat = patient.get_reference('genomewide')

    add_SHAPE_to_seqrecord(refpat)


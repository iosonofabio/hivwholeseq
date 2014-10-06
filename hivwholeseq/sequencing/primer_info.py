# vim: fdm=marker
'''
author:     Fabio Zanini
date:       24/09/13
content:    Information about the primers used for PCR in Sweden.

            Note: the old coordinates were wrong.
'''
# Modules
from Bio.Seq import reverse_complement as rc


# Globals
fragments_genes = {'gag': ['F1', 'F2'],
                   'pol': ['F1', 'F2', 'F3', 'F4'],
                   'env': ['F4', 'F5', 'F6'],
                   'vif': ['F4'],
                   'vpr': ['F4'],
                   'vpu': ['F4'],
                   'tat': [['F4'], ['F6']],
                   'rev': [['F4'], ['F6']],
                   #'nef': ['F6']
                   }

fragments_exons = {g: frs for (g, frs) in fragments_genes.iteritems()
                   if g not in ('tat', 'rev')}
fragments_exons['tat1'] = fragments_genes['tat'][0]
fragments_exons['tat2'] = fragments_genes['tat'][1]
fragments_exons['rev1'] = fragments_genes['rev'][0]
fragments_exons['rev2'] = fragments_genes['rev'][1]

fragments_RNA_structures = {'RRE': ['F5', 'F6'],
                            "LTR5'": ['F1'],
                            "LTR3'": ['F6']}

fragments_other = {'env peptide': ['F4', 'F5'],
                   'psi': ['F1']}

# Note: the reverse primers get reverse complemented (so everything is positive sense)
primers_inner = {'F1': ['AAGTAGTGTGTGCCCGTCTGT', rc('TGCCAAAGAGTGATYTGAGGG')],
                 'F2': ['GGGCTGTTGGARATGTGG', rc('ACAAACTCCCAYTCAGGAATCCA')],
                 'F3a': ['GAAAGCATAGTRATATGGGGAAA', rc('CACCTGCCATCTGTTTTCCATA')],
                 'F4': ['TGGAAAGGTGAAGGGGCAG', rc('GTACACAGGCATGTGTRGCCCA')],
                 'F5a': ['TAAGAGAAAGAGCAGAAGACAGTGG', rc('CCAAATYCCYAGGAGCTGTTGATC')],
                 'F5b': ['TCTATTATGGRGTACCTGTRTGG', rc('CCAAATYCCYAGGAGCTGTTG')],
                 'F6': ['CAGGAAGCACTATGGGCGC', rc('CCAGAGAGCTCCCAGG')],
                }

primers_outer = {'F1': ['CTCAATAAAGCTTGCCTTGAGTGC', rc('ACTGTATCATCTGCTCCTGTRTCT')],
                 'F2': ['AAATTGCAGGGCYCCTAG', rc('CTRTTAGCTGCCCCATCTACATAG')],
                 'F3a': ['CACACTAATGATGTAARACARTTAACAG', rc('TTCCATGTTYTAATCCTCATCCTGTCTAC')],
                 # NOTE: F3b rev might be degenerate, ask Lina (she's on holiday now)!!
                 'F3b': ['CACACTAATGATGTAARACARTTAACAG', 'CAAGAATAAGTTCAGAAGTACACATCCC'],
                 'F4': ['CGGGTTTATTWCAGRGACAGCAGA', rc('GGGGTTAAYTTTACACATGGYTTTA')],
                 'F5a': ['GGCATYTCCTATGGCAGGAAGAAG', rc('GTGGTGCARATGAGTTTTCCAGAGCA')],
                 'F5b': ['AGAATAAGAGAAAGAGCAGAAGA', rc('ATGAGTTTTCCAGAGCANCCCCA')],
                 'F6': ['GGGTTCTTRGGARCAGCAGGAAG', rc('ATTGAGGCTTAAGCAGTGGGTTC')],
                }

# All primers together (e.g. F5ai is the inner primer for F5, version a)
primers_PCR = dict([(fi+'i', si) for (fi, si) in primers_inner.iteritems()] + \
                   [(fo+'o', so) for (fo, so) in primers_outer.iteritems()])

# Note: the reverse are sorted already, and coordinates start from 0 and end at
# the last nucleotide + 1 (a la Python).
primers_coordinates_HXB2_inner = {'F1': [[550, 571], [2251, 2272]],
                                  'F2': [[2021, 2039], [3776, 3799]],
                                  'F3a': [[3680, 3703], [5039, 5061]],
                                  'F4': [[4955, 4974], [6428, 6450]],
                                  'F5a': [[6197, 6222], [7988, 8012]],
                                  'F5b': [[6336, 6359], [7991, 8012]],
                                  'F6': [[7800, 7819], [9566, 9582]],
                                 }

primers_coordinates_HXB2_outer = {'F1': [[523, 547], [2323, 2347]],
                                  'F2': [[1997, 2015], [3868, 3892]],
                                  'F3a': [[3629, 3657], [5076, 5105]],
                                  'F3b': [[3629, 3657], [5185, 5213]],
                                  'F4': [[4898, 4922], [6570, 6595]],
                                  'F5a': [[5959, 5982], [8015, 8041]],
                                  'F5b': [[6193, 6216], [8009, 8032]],
                                  'F6': [[7784, 7807], [9591, 9614]],
                                 }

primers_coordinates_HXB2 = dict([(fi+'i', si) for (fi, si) in primers_coordinates_HXB2_inner.iteritems()] + \
                                [(fo+'o', so) for (fo, so) in primers_coordinates_HXB2_outer.iteritems()])

# Functions
def find_fragment(refseq, fragment, threshold=0.7):
    '''Find the coordinate of one fragment in the refseq'''
    import numpy as np

    refs = ''.join(refseq)
    refm = np.fromstring(refs, 'S1')
    rl = len(refm)

    (prfwd, prrev) = primers_PCR[fragment]

    seed = np.ma.array(np.fromstring(prfwd, 'S1'))
    seed[seed == 'Y'] = np.ma.masked
    seed[seed == 'R'] = np.ma.masked
    seed[seed == 'W'] = np.ma.masked
    seed[seed == 'N'] = np.ma.masked
    sl = len(seed)
    n_matches = np.array([(seed == refm[i: i + sl]).sum() for i in xrange(len(refm) - sl)], int)
    poss_start = np.argsort(n_matches)[::-1][:5]
    if n_matches[poss_start[0]] < threshold * (-seed.mask).sum():
        raise ValueError('Start of fragment not found')

    seed = np.ma.array(np.fromstring(prrev, 'S1'))
    seed[seed == 'Y'] = np.ma.masked
    seed[seed == 'R'] = np.ma.masked
    seed[seed == 'W'] = np.ma.masked
    seed[seed == 'N'] = np.ma.masked
    sl = len(seed)
    n_matches = np.array([(seed == refm[i: i + sl]).sum() for i in xrange(len(refm) - sl)], int)
    poss_end = np.argsort(n_matches)[::-1][:5]
    if n_matches[poss_end[0]] < threshold * (-seed.mask).sum():
        raise ValueError('End of fragment not found')

    found = False
    for pos_start in poss_start:
        for pos_end in poss_end:
            if 1000 < pos_end - pos_start < 2300:
                found = True
                break
        if found:
            break

    if not found:
        raise ValueError('No suitable combination of fwd/rev primers found')

    return (pos_start, pos_end)


def find_primer_seq(seq, primer, from_right=False, threshold=0.7):
    '''Find an ambiguous primer in a sequence'''
    import numpy as np

    refs = ''.join(seq)
    refm = np.fromstring(refs, 'S1')
    rl = len(refm)

    seed = np.ma.array(np.fromstring(primer, 'S1'))
    seed[seed == 'Y'] = np.ma.masked
    seed[seed == 'R'] = np.ma.masked
    seed[seed == 'W'] = np.ma.masked
    seed[seed == 'N'] = np.ma.masked
    sl = len(seed)
    n_matches = np.array([(seed == refm[i: i + sl]).sum() for i in xrange(len(refm) - sl)], int)
    if not from_right:
        pos_start = np.argmax(n_matches)
    else:
        pos_start = len(n_matches) - 1 - np.argmax(n_matches[::-1])

    if n_matches[pos_start] < threshold * (-seed.mask).sum():
        raise ValueError('Start of fragment not found')

    return pos_start


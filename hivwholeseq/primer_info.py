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
# Note: the reverse primers get reverse complemented
primers_inner = {'F1': ['AAGTAGTGTGTGCCCGTCTGT', rc('TGCCAAAGAGTGATYTGAGGG')],
                 'F2': ['GGGCTGTTGGARATGTGG', rc('ACAAACTCCCAYTCAGGAATCCA')],
                 'F3': ['GAAAGCATAGTRATATGGGGAAA', rc('CACCTGCCATCTGTTTTCCATA')],
                 'F4': ['TGGAAAGGTGAAGGGGCAG', rc('GTACACAGGCATGTGTRGCCCA')],
                 'F5a': ['TAAGAGAAAGAGCAGAAGACAGTGG', rc('CCAAATYCCYAGGAGCTGTTGATC')],
                 'F5b': ['TCTATTATGGRGTACCTGTRTGG', rc('CCAAATYCCYAGGAGCTGTTG')],
                 'F6': ['CAGGAAGCACTATGGGCGC', rc('CCAGAGAGCTCCCAGG')],
                }

primers_outer = {'F1': ['CTCAATAAAGCTTGCCTTGAGTGC', rc('ACTGTATCATCTGCTCCTGTRTCT')],
                 'F2': ['AAATTGCAGGGCYCCTAG', rc('CTRTTAGCTGCCCCATCTACATAG')],
                 'F3': ['CACACTAATGATGTAARACARTTAACAG', rc('TTCCATGTTYTAATCCTCATCCTGTCTAC')],
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
                                  'F3': [[3680, 3703], [5039, 5061]],
                                  'F4': [[4955, 4974], [6428, 6450]],
                                  'F5a': [[6197, 6222], [7988, 8012]],
                                  'F5b': [[6336, 6359], [7991, 8012]],
                                  'F6': [[7800, 7819], [9566, 9582]],
                                 }

primers_coordinates_HXB2_outer = {'F1': [[523, 547], [2323, 2347]],
                                  'F2': [[1997, 2015], [3868, 3892]],
                                  'F3': [[3629, 3657], [5076, 5105]],
                                  'F4': [[4898, 4922], [6570, 6595]],
                                  'F5a': [[5959, 5982], [8015, 8041]],
                                  'F5b': [[6193, 6216], [8009, 8032]],
                                  'F6': [[7784, 7807], [9591, 9614]],
                                 }


# Functions
def get_fragment_positions(fragments):
    '''Get a dictionary of positions for the selected fragments'''
    #TODO
    pass

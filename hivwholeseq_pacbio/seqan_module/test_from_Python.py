# vim: fdm=marker
'''
author:     Fabio Zanini
date:       30/01/14
content:    Test of the library from Python
'''

# Script
if __name__ == '__main__':

    import seqanpy as sap

    # Test find seed
    refseq = 'AAAGGTAGCTTACAGGCT'
    seed = 'TAGCTAAC'

    output = sap.find_seed(refseq, seed)
    print output
    print refseq[output[0]: output[0] + len(seed)]
    print seed


# vim: fdm=indent
'''
author:     Fabio Zanini
date:       31/08/15
content:    Module that contains fast conversion factories between private and
            anonymyzed versions of data.

            Such a module is useful during research development, when consensus
            on anonymized patient/sample names has not been reached yey.

            NOTE: remember to always anonymyze data before publishing them!
'''
# pdict is a dictionary with the anonymyzed patient names as keys and the actual
# names as values. It must be invertible (do not reuse values).
pdict = {'p1': 'WRITE_ME',
        }

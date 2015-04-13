# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/10/14
content:    Support module for parsing input arguments
'''
# Modules
from __future__ import absolute_import
from argparse import Action



# Globals
_pdict = {'p1': '20097',
          'p2': '15363',
          'p3': '15823',
          'p4': '15313',
          'p5': '15376',
          'p6': '20529',
          'p7': '15107',
          'p8': '9669',
          'p9': '15241',
          'p10': '15034',
          'p11': '15319',
         }



# Classes
class RoiAction(Action):
    '''Argparse action for a genomic region of interest'''

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(RoiAction, self).__init__(option_strings, dest, nargs=3, **kwargs)


    def __call__(self, parser, namespace, values, option_string=None):
        from argparse import ArgumentTypeError

        if len(values) == 0:
            setattr(namespace, self.dest, None)
            return

        try:
            fragment = values[0]
            start = int(values[1])
            if not (values[2] == '+oo'):
                end = int(values[2])
            else:
                end = values[2]

        except ValueError:
            raise ArgumentTypeError('Arguments not good for a ROI.')

        roi = (fragment, start, end)
        setattr(namespace, self.dest, roi)


class RoisAction(Action):
    '''Argparse action for genomic regions of interest'''

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(RoisAction, self).__init__(option_strings, dest, nargs='+', **kwargs)


    def __call__(self, parser, namespace, values, option_string=None):
        from argparse import ArgumentTypeError

        if len(values) == 0:
            setattr(namespace, self.dest, [])
            return

        if len(values) % 3:
            raise ArgumentTypeError('nargs for ROIs must be a multiple of 3.')

        rois = []
        for i in xrange(len(values) // 3):
            try:
                fragment = values[3 * i]
                start = int(values[3 * i + 1])
                end = int(values[3 * i + 2])

            except ValueError:
                raise ArgumentTypeError('Arguments not good for a ROI.')

            rois.append((fragment, start, end))

        setattr(namespace, self.dest, rois)


class PatientAction(Action):
    '''Argparse action for a patient'''

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(PatientAction, self).__init__(option_strings, dest, nargs=1, **kwargs)


    def __call__(self, parser, namespace, values, option_string=None):
        from argparse import ArgumentTypeError

        if len(values) == 0:
            setattr(namespace, self.dest, None)
            return

        pnames = [_pdict.get(v, _pdict.get('p'+str(v), v)) for v in values]

        pv = _pdict.values()
        for pname in pnames:
            if pname not in pv:
                raise ArgumentTypeError('Patient not recognized: '+str(pname))

        setattr(namespace, self.dest, pnames[0])


class PatientsAction(Action):
    '''Action for patients'''

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(PatientsAction, self).__init__(option_strings, dest, nargs='+', **kwargs)


    def __call__(self, parser, namespace, values, option_string=None):
        from argparse import ArgumentTypeError

        if len(values) == 0:
            setattr(namespace, self.dest, [])
            return

        pnames = [_pdict.get(v, _pdict.get('p'+str(v), v)) for v in values]

        pv = _pdict.values()
        for pname in pnames:
            if pname not in pv:
                raise ArgumentTypeError('Patient not recognized: '+str(pname))

        setattr(namespace, self.dest, pnames)

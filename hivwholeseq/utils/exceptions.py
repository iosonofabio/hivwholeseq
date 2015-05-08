# vim: fdm=indent
'''
author:     Fabio Zanini
date:       06/05/15
content:    Custom exceptions.
'''
class RoiError(ValueError):
    '''Error related to a Region Of Interest'''
    pass


class NoDataError(ValueError):
    '''Error in case we have no data covering this situation'''
    pass

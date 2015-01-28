# vim: fdm=marker
'''
author:     Fabio Zanini
date:       31/10/14
content:    Support module to add useful plot functions, e.g. logit scales.
'''
# Modules
import numpy as np
from matplotlib import transforms as mtransforms
from matplotlib.ticker import Formatter, FixedLocator
from matplotlib.ticker import NullFormatter, ScalarFormatter, FuncFormatter
from matplotlib import scale as mscale



# Classes
class LogitTransform(mtransforms.Transform):
    """
    Logit transform: x -> log10(x / (1 - x)).
    """
    input_dims = 1
    output_dims = 1
    is_separable = True

 
    def transform_non_affine(self, a):
        """
        This transform takes an Nx1 ``numpy`` array and returns a
        transformed copy. Importantly, the ``transform`` method
        *must* return an array that is the same shape as the input
        array, since these values need to remain synchronized with
        values in the other dimension.
        """
        masked = np.ma.masked_where((a <= 0) | (a >= 1), a)
        if masked.mask.any():
            return np.ma.log(masked / (1.0 - masked))
        else:
            return np.log(a / (1.0 - a))


    def inverted(self):
        """
        Override this method so matplotlib knows how to get the
        inverse transform for this transform.
        """
        return LogitTransformInverted()


class LogitTransformInverted(mtransforms.Transform):
    """
    Inverted logit transform: x -> 1 / (1 + exp(-y)).
    """
    input_dims = 1
    output_dims = 1
    is_separable = True
 
    def transform_non_affine(self, a):
        """
        This transform takes an Nx1 ``numpy`` array and returns a
        transformed copy.
        """
        return 1.0 / (1.0 + np.exp(-a))


    def inverted(self):
        """
        Override this method so matplotlib knows how to get the
        inverse transform for this transform.
        """
        return LogitTransform()


class LogitScale(mscale.ScaleBase):
    """
    A standard logit scale. Care is taken so values beyond ]0, 1[ are not
    plotted.

    Only base 10 ticks are implemented (e.g. [0.01, 0.1, 0.5, 0.9, 0.99]).
    """
    name = 'logit'
    thresholds = [-10, 10]


    def __init__(self, axis, **kwargs):
        pass


    def get_transform(self):
        """
        Logit transform: x -> log10(x / (1 - x)).
        """
        return LogitTransform()


    def set_default_locators_and_formatters(self, axis):
        """Set default ticks and labels"""

        # FIXME: pretty rough, ouch!
        logrange = self.thresholds
        tickleft = np.logspace(logrange[0], -1, -logrange[0])
        tickcenter = [0.5]
        tickright = 1 - np.logspace(-logrange[1], -1, logrange[1])[::-1]
        tickloc = np.concatenate([tickleft, tickcenter, tickright])

        tickleft = [[10**po * x for x in xrange(2 , 10)]
                    for po in xrange(logrange[0] - 1, -1)]
        tickcenter = [[0.1 * x for x in xrange(2 , 9)]]
        tickright = [[1 - 10**po * (10 - x) for x in xrange(2, 10)]
                     for po in xrange(-2, logrange[1] - 2, -1)]
        ticklocminor = np.concatenate(tickleft + tickcenter + tickright)

        axis.set_major_locator(FixedLocator(tickloc))
        axis.set_minor_locator(FixedLocator(ticklocminor))

        axis.set_major_formatter(FuncFormatter(logit_formatter))
        axis.set_minor_formatter(NullFormatter())

    
    def limit_range_for_scale(self, vmin, vmax, minpos):
        """
        Limit the domain to values in ]0, 1[.
        """
        #TODO: this is not very easy with the stupid minpos argument only
        return (max(vmin, 10**(self.thresholds[0])),
                min(vmax, 1 - 10**(-self.thresholds[1])))


# Functions
def logit_formatter(x, pos):
    '''Format as 10^x for < 0.5, 1 - 10^x for > 0.5'''
    if x < 0.4:
        s = r'10^{'+r'{:1.0f}'.format(np.round(np.log10(x), 0))+r'}'

    elif x > 1.0 - 0.4:
        s = r'1 - 10^{'+r'{:1.0f}'.format(np.round(np.log10(1 - x), 0))+r'}'

    else:
        s = str(np.round(x, 1))

    return r'$\mathdefault{'+s+r'}$'


def finalize():
    '''Finalize the module by registering its components import matplotlib'''
    mscale.register_scale(LogitScale)



# Finalize
finalize()



# Script (for debugging)
if __name__ == '__main__':


    x = [0.00001, 0.01, 0.02, 0.1, 0.2, 0.5, 0.9, 0.98, 0.99]
    y = np.arange(len(x))

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.set_xscale('logit')

    ax.grid(True)

    plt.ion()
    plt.show()

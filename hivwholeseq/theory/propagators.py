# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/10/14
content:    Get and plot the propagator according to various versions of Katya's
            paper, i.e. Kosheleva and Desai 2013.
'''
# Modules
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt



# Functions
def propagator_neutral(xi, t, xf=np.linspace(0, 1, 100)[1:-1], n=10000):
    '''Neutral process propagator (Kimura solution with Gegenbauer polynomials)
    
    Note: for speed reasons, 4-5 float matrices of size
                    
                    len(xi) * len(xf) * n

          are stored in memory. If that's too big, call this function with less
          aggressive vectorization, e.g. within one or more for cycles.
    '''
    from scipy import special as sp
    r = 1 - 2 * xi
    z = 1 - 2 * xf
    i_array = np.arange(1, n + 1)

    if np.isscalar(xi) and np.isscalar(xf):
        i = i_array
        res = (2 * i+1) * (1 - r**2) / i / (i+1) * \
               sp.eval_gegenbauer(i-1, 1, r) * \
               sp.eval_gegenbauer(i-1, 1, z) * \
               np.exp(-i * (i+1) * t / 4)
        return (xf, res)

    elif np.isscalar(xi):
        zs = np.repeat(z, len(i_array))
        i = np.tile(i_array, len(z))
        res = (2 * i+1) * (1 - r**2) / i / (i+1) * \
              sp.eval_gegenbauer(i-1, 1, r) * \
              sp.eval_gegenbauer(i-1, 1, zs) * \
              np.exp(-i * (i+1) * t / 4)
        res = res.reshape((len(z), len(i_array))).sum(axis=1)
        return (xf, res)

    elif np.isscalar(xf):
        rs = np.repeat(r, len(i_array))
        i = np.tile(i_array, len(r))
        res = (2 * i+1) * (1 - rs**2) / i / (i+1) * \
              sp.eval_gegenbauer(i-1, 1, rs) * \
              sp.eval_gegenbauer(i-1, 1, z) * \
              np.exp(-i * (i+1) * t / 4)
        res = res.reshape((len(r), len(i_array))).sum(axis=1)
        return (xf, res)

    else:
        # Note: this function minimizes Gegenbauer polynomial evaluations, but
        # the cost is a rather opaque (vectorized) code.
        i = np.tile(i_array, len(r) * len(z))
        i = i.reshape((len(r), len(z), len(i_array)))

        # r enters appears both alone and via the Gegenbauer
        rs = np.tile(r, len(i_array) * len(z))
        rs = rs.reshape((len(i_array), len(z), len(r))).swapaxes(0, 2)
        rs_gb = sp.eval_gegenbauer(i[:, 0].ravel() - 1, 1, rs[:, 0].ravel())
        rs_gb = np.tile(rs_gb, len(z))
        rs_gb = rs_gb.reshape((len(z), len(r), len(i_array))).swapaxes(0, 1)

        # z only appears via the Gegenbauer polymonial
        zs_gb = sp.eval_gegenbauer(i[0].ravel() - 1, 1, np.repeat(z, len(i_array)))
        zs_gb = np.tile(zs_gb, len(r))
        zs_gb = zs_gb.reshape((len(r), len(z), len(i_array)))

        res = (2 * i+1) * (1 - rs**2) / i / (i+1) * \
              rs_gb * \
              zs_gb * \
              np.exp(-i * (i+1) * t / 4)
        res = res.sum(axis=2)
        return (xf, res)


    #res = np.zeros_like(xf)
    #for i in i_array:
    #    res += (2 * i+1) * (1 - r**2) / i / (i+1) * \
    #           sp.eval_gegenbauer(i-1, 1, r) * \
    #           sp.eval_gegenbauer(i-1, 1, z) * \
    #           np.exp(-i * (i+1) * t / 4)

    #return (xf, res)


def propagator_BSC(xi, t, xf=np.linspace(0, 1, 1000)[1:-1]):
    '''BSC propagator
    
    See eq.15 of the Kosheleva and Desai (2013) paper.
    
                                          sin(pi a) xi (1 - xi)    
    rho(xf | xi) = -------------------------------------------------------------
                   xf (1 - xf) pi [ (1 - xi)^2   (xf / (1 - xf))^a + 
                                      xi^2       ((1 - xf) / xf)^a +
                                   2 xi (1 - xi)     cos(pi a)       ]
    
    '''
    a = np.exp(-t)
    xii = 1.0 - xi
    xfi = 1.0 - xf
    pi = np.pi
    num = np.sin(pi * a) * xi * xii
    den = xf * xfi * pi * (xii**2 * (xf / xfi)**a + \
                           xi**2 * (xfi / xf)**a + \
                           2 * xi * xii * np.cos(pi * a))

    return (xf, num / den)


def rho16(xi, q, xf=np.linspace(0, 1, 1000)):
    '''rho(x_k | x_k-1) = x_k-1 (1 - x_k-1) / (q * dx^2) -- q = 8 in HIV'''
    dx2 = (xf - xi)**2
    rho = xi * (1.0 - xi) / (q * dx2)
    return (xf, rho)



# Script
if __name__ == '__main__':

        ## Eq. 15 has the long formula
        #for t in np.logspace(-3, 0, 5):
        #    fig, ax = plt.subplots()
        #    xis = [0.13, 0.29, 0.5]
        #    for i, xi in enumerate(xis):
        #        xf, rho = propagator_BSC(xi, t)
        #        ax.plot(xf, rho, color=cm.jet(1.0 * i / len(xis)), lw=2,
        #                label='$x_i = '+str(xi)+'$')
        #        ax.axvline(xi, lw=2, color=cm.jet(1.0 * i / len(xis)))
        #    ax.set_ylabel('$\\rho (xf | xi)$')
        #    ax.set_xlabel('Final frequency')
        #    ax.set_title('rho15, $t = 10^{'+str(np.log10(t))+'}$')
        #    ax.set_yscale('log')
        #    ax.grid(True)

        # Check vectorization of gegenbauer polynomials
        xis = np.array([1e-2, 1e-1, 5e-1, 9e-1, 9.9e-1])
        xfs = np.array([0.03, 0.1, 0.3, 0.7, 0.9, 0.97])
        t = 1
        (_, res) = propagator_neutral(xis, t, xf=xfs, n=10000)
        res2 = []
        for i, xi in enumerate(xis):
            (_, r) = propagator_neutral(xi, t, xf=xfs, n=10000)
            res2.append(r)
            if (r != res[i]).any():
                print 'Error'
                break


        ## Eq. 16 contains the q --> +oo asymptote
        #fig, ax = plt.subplots()
        #q = 8
        #xis = [0.13, 0.29, 0.5]
        #for i, xi in enumerate(xis):
        #    xf, rho = rho16(xi, q)
        #    ax.plot(xf, rho, color=cm.jet(1.0 * i / len(xis)), lw=2,
        #            label='$x_i = '+str(xi)+'$')
        #ax.set_ylabel('$\\rho (xf | xi)$')
        #ax.set_xlabel('Final frequency')
        #ax.set_title('rho16')
        #ax.grid(True)

        plt.ion()
        plt.show()

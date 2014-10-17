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
def propagator_neutral(x,p,t,n):
    '''Neutral process propagator (Kimura solution with Gegenbauer polynomials)'''
    #TODO: vectorize this function
    from scipy import special as sp
    r = 1 - 2 * p
    z = 1 - 2 * x
    i_array = np.arange(1, n + 1)
    res = np.zeros_like(x)
    for i in i_array:
        res += (2 * i+1) * (1 - r**2) / i / (i+1) * \
               sp.eval_gegenbauer(i-1, 1, r) * \
               sp.eval_gegenbauer(i-1, 1, z) * \
               np.exp(-i * (i+1) * t / 4)

    return res


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

        # Eq. 15 has the long formula
        for t in np.logspace(-3, 0, 5):
            fig, ax = plt.subplots()
            xis = [0.13, 0.29, 0.5]
            for i, xi in enumerate(xis):
                xf, rho = propagator_BSC(xi, t)
                ax.plot(xf, rho, color=cm.jet(1.0 * i / len(xis)), lw=2,
                        label='$x_i = '+str(xi)+'$')
                ax.axvline(xi, lw=2, color=cm.jet(1.0 * i / len(xis)))
            ax.set_ylabel('$\\rho (xf | xi)$')
            ax.set_xlabel('Final frequency')
            ax.set_title('rho15, $t = 10^{'+str(np.log10(t))+'}$')
            ax.set_yscale('log')
            ax.grid(True)


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

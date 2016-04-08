import numpy as np
from gpaw.xc.kernel import XCKernel


class LB94(XCKernel):
    """Correction to LDA to resemble asymptotic -1/r potential.

    See:

        van Leeuwen and Baerends, Phys.Rev.A vol 49 (1994) 2421
    """

    def __init__(self, beta=0.05):
        XCKernel.__init__(self, 'LDA')
        self.name = 'LB94'
        self.type = 'GGA'
        self.beta = beta

    def calculate(self, e_g, n_sg, dedn_sg,
                  sigma_xg=None, dedsigma_xg=None,
                  tau_sg=None, dedtau_sg=None):
        XCKernel.calculate(self, e_g, n_sg, dedn_sg, sigma_xg, dedsigma_xg)
        for s, n_g in enumerate(n_sg):
            n_g = n_g * len(n_sg)
            n_g[n_g < 1e-10] = 1e-10
            y_g = n_g**(1 / 3.0)
            x_g = sigma_xg[2 * s]**0.5 / (y_g * n_g) * len(n_sg)
            x_g[x_g > 500] = 0.0
            dedn_sg[s] -= (self.beta * x_g**2 * y_g /
                           (1 + 3 * self.beta * x_g * np.arcsinh(x_g)))
        dedsigma_xg[:] = 0.0

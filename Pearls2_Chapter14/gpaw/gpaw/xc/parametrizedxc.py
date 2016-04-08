import numpy as np

from gpaw.xc.libxc import LibXC
from gpaw.xc.kernel import XCKernel


class ParametrizedKernel(XCKernel):
    def __init__(self, name):
        self.name = name
        self.kernels = []
        self.coefs = []
        self.type = 'LDA'
        for x in name.split('+'):
            c, n = x.split('_', 1)
            self.coefs.append(float(c))
            kernel = LibXC(n)
            self.kernels.append(kernel)
            if kernel.type == 'GGA':
                self.type = 'GGA'

    def calculate(self, e_g, n_sg, dedn_sg,
                  sigma_xg=None, dedsigma_xg=None,
                  tau_sg=None, dedtau_sg=None):
        e_g[:] = 0.0
        e_g_tmp = np.empty_like(e_g)
        dedn_sg_tmp = np.empty_like(dedn_sg)

        if self.type == 'GGA':
            dedsigma_xg[:] = 0.0
            dedsigma_xg_tmp = np.empty_like(dedsigma_xg)
        else:
            dedsigma_xg_tmp = None
            
        for kernel, c in zip(self.kernels, self.coefs):
            dedn_sg_tmp[:] = 0.0
            if self.type == 'GGA':
                dedsigma_xg_tmp[:] = 0.0
            kernel.calculate(e_g_tmp, n_sg, dedn_sg_tmp,
                             sigma_xg, dedsigma_xg_tmp,
                             tau_sg, dedtau_sg)
            e_g += c * e_g_tmp
            dedn_sg += c * dedn_sg_tmp
            if self.type == 'GGA':
                dedsigma_xg += c * dedsigma_xg_tmp

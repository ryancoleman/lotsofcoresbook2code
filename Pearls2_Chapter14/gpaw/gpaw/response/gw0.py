from __future__ import print_function

import numpy as np

from gpaw.response.g0w0 import G0W0


class GW0(G0W0):
    def __init__(self, calc, filename='gw0', **kwargs):
        G0W0.__init__(self, calc, filename, savew=True, **kwargs)
        try:
            self.qp_xsin = np.load(filename + '-gw0.npy')
        except IOError:
            self.qp_xsin = np.empty((1,) + self.shape)
            
        self.iteration = len(self.qp_xsin)
    
    def get_k_point(self, s, K, n1, n2):
        kpt = G0W0.get_k_point(self, s, K, n1, n2)
        if self.iteration > 1:
            b1, b2 = self.bands
            m1 = max(b1, n1)
            m2 = min(b2, n2)
            if m2 > m1:
                k = self.calc.wfs.kd.bz2ibz_k[K]
                i = self.kpts.index(k)
                qp_n = self.qp_xsin[-1, s, i, m1 - b1:m2 - b1]
                kpt.eps_n[m1 - n1:m2 - n1] = qp_n
                
        return kpt
    
    def update_qp_energies(self):
        print('Iteration:', self.iteration, file=self.fd)
        
        if self.iteration == 1:
            self.qp_xsin[0] = self.eps_sin
            
        self.Z_sin = 1 / (1 - self.dsigma_sin)
        self.qp_sin = self.eps_sin + self.Z_sin * (
            self.sigma_sin + self.exx_sin -
            self.vxc_sin)
        
        self.iteration += 1
        qp_xsin = np.empty((self.iteration,) + self.shape)
        qp_xsin[:-1] = self.qp_xsin
        qp_xsin[-1] = self.qp_sin
        self.qp_xsin = qp_xsin
        
        if self.world.rank == 0:
            with open(self.filename + '.qp.npy', 'w') as fd:
                np.save(fd, self.qp_xsin)

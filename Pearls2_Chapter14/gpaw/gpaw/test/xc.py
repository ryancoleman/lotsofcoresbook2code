import numpy as np
from gpaw.xc.libxc import LibXC, short_names
from gpaw.xc.kernel import XCKernel, codes
from gpaw.test import equal

funcs = []
modes = []
for name in short_names:
    try:
        LibXC(name)
    except NameError:
        continue
    funcs.append(name)
    modes.append(0)
for name in codes:
    funcs.append(name)
    modes.append(1)


def create_xc(func, mode):
    isinstance(func, str)
    isinstance(mode, int)
    if mode == 0:
        xc = LibXC(func)
    else:
        xc = XCKernel(func)
    return xc

    
def f1(n_xg, xc):
    e_g = np.empty_like(n_xg[0])
    n_sg = n_xg[:1]
    sigma_xg = n_xg[1:2]
    tau_sg = n_xg[2:]
    dedn_sg = np.zeros_like(n_sg)
    dedsigma_xg = np.zeros_like(sigma_xg)
    dedtau_sg = np.zeros_like(tau_sg)
    xc.calculate(e_g, n_sg, dedn_sg, sigma_xg, dedsigma_xg, tau_sg, dedtau_sg)
    return e_g, np.concatenate((dedn_sg, dedsigma_xg, dedtau_sg))

    
def f2(n_xg, xc):
    e_g = np.empty_like(n_xg[0])
    n_sg = n_xg[:2]
    sigma_xg = n_xg[2:5]
    tau_sg = n_xg[5:]
    dedn_sg = np.zeros_like(n_sg)
    dedsigma_xg = np.zeros_like(sigma_xg)
    dedtau_sg = np.zeros_like(tau_sg)
    xc.calculate(e_g, n_sg, dedn_sg, sigma_xg, dedsigma_xg, tau_sg, dedtau_sg)
    return e_g, np.concatenate((dedn_sg, dedsigma_xg, dedtau_sg))

eps = 1.0e-6

n_xg = np.array(
    [[0.2, 0.01, 0.4],
     [0.2, 0.1, 0.5],
     [0.01, 0.01, 0.2],
     [0.1, 0.3, 0.5]]).T.copy()

for i, func in enumerate(funcs):
    xc = create_xc(funcs[i], modes[i])
    e0_g, d0_xg = f1(n_xg, xc)
    d_xg = np.empty_like(d0_xg)
    for x, n_g in enumerate(n_xg):
        m_xg = n_xg.copy()
        m_xg[x] += eps
        d_xg[x] = 0.5 * f1(m_xg, xc)[0] / eps
        m_xg[x] -= 2 * eps
        d_xg[x] -= 0.5 * f1(m_xg, xc)[0] / eps
    ns_xg = np.empty((7, len(n_g)))
    ns_xg[:2] = n_xg[0] / 2
    ns_xg[2:5] = n_xg[1] / 4
    ns_xg[5:] = n_xg[2] / 2
    es_g, ds_xg = f2(ns_xg, xc)
    error = (abs(d0_xg - d_xg).max() +
             abs(es_g - e0_g).max() +
             abs(ds_xg[:2] - d0_xg[0]).max() +
             abs(ds_xg[2:5].sum(0) / 4 - d0_xg[1]).max() +
             abs(ds_xg[5:] - d0_xg[2]).max())
    equal(error, 0, 6e-9)
    del xc

# Numbers from old lxc_xc.py test:
na = 2.0
nb = 1.0
sigma0 = 2.0  # (0.0, 1.0, 1.0)
sigma1 = 2.0
sigma2 = 5.0  # (1.0, 2.0, 0.0)
taua = (3 * np.pi**2)**(2. / 3.) * na**(5. / 3.) / 2 * sigma0
taub = (3 * np.pi**2)**(2. / 3.) * nb**(5. / 3.) / 2 * sigma2

n_xg = np.array(
    [[na, nb, sigma0, sigma1, sigma2, taua, taub],
     [0.1, 0.1, 0.025, 0.025, 0.025, 0.25, 0.25],
     [0.1, 0.1, 0.125, 0.125, 0.125, 0.0025, 0.025],
     [0.1, 0.1, 0.01, 0.01, 0.01, 0.2, 0.2],
     [0.1, 0.2, 0.1, -0.08, 0.10, 0.01, 0.05],
     [0.1, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01],
     [0.1, 0.1, 0.1, 0.15, 0.20, 0.01, 0.05]]).T.copy()


for i, func in enumerate(funcs):
    xc = create_xc(funcs[i], modes[i])
    if xc.type == 'MGGA':
        N_xg = n_xg[:, :1].copy()
    else:
        N_xg = n_xg
    e0_g, d0_xg = f2(N_xg, xc)
    d_xg = np.empty_like(d0_xg)
    for x, n_g in enumerate(N_xg):
        m_xg = N_xg.copy()
        m_xg[x] += eps
        d_xg[x] = 0.5 * f2(m_xg, xc)[0] / eps
        m_xg[x] -= 2 * eps
        d_xg[x] -= 0.5 * f2(m_xg, xc)[0] / eps
    equal(abs(d0_xg - d_xg).max(), 0, 2e-8)
    del xc

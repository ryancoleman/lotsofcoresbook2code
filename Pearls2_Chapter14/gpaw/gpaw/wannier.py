from math import pi
from pickle import load, dump

import numpy as np
from ase.units import Bohr

from _gpaw import localize
from gpaw.utilities.tools import dagger, lowdin


class Wannier:
    def __init__(self, calc=None, spin=0, nbands=None):
        self.spin = spin
        self.Z_nnc = None
        if calc is not None:
            if not calc.wfs.gd.orthogonal:
                raise NotImplementedError('Wannier function analysis ' +
                                          'requires an orthogonal cell.')
            self.cell_c = calc.wfs.gd.cell_cv.diagonal() * Bohr
            if nbands is None:
                nbands = calc.get_number_of_bands()
            self.Z_nnc = np.empty((nbands, nbands, 3), complex)
            print("calculating Z_nnc")
            for c in range(3):
                G_c = np.zeros(3)
                G_c[c] = 1
                self.Z_nnc[:, :, c] = calc.get_wannier_integrals(
                    spin, 0, 0, G_c, nbands)
            self.value = 0.0
            self.U_nn = np.identity(nbands)

    def load(self, filename):
        self.cell_c, self.Z_nnc, self.value, self.U_nn = load(open(filename))

    def dump(self, filename):
        dump((self.cell_c, self.Z_nnc, self.value, self.U_nn), filename)
        
    def localize(self, eps=1e-5, iterations=-1):
        i = 0
        while i != iterations:
            value = localize(self.Z_nnc, self.U_nn)
            print(i, value)
            if value - self.value < eps:
                break
            i += 1
            self.value = value
        return value # / Bohr**6

    def get_centers(self):
        scaled_c = -np.angle(self.Z_nnc.diagonal()).T / (2 * pi)
        return (scaled_c % 1.0) * self.cell_c

    def get_function(self, calc, n, pad=True):
        if pad:
            return calc.wfs.gd.zero_pad(self.get_function(calc, n, False))
        psit_nG = calc.wfs.kpt_u[self.spin].psit_nG[:]
        psit_nG = psit_nG.reshape((calc.wfs.bd.nbands, -1))
        return np.dot(self.U_nn[:, n],
                       psit_nG).reshape(calc.wfs.gd.n_c) / Bohr**1.5

    def get_hamiltonian(self, calc):
        # U^T diag(eps_n) U
        eps_n = calc.get_eigenvalues(kpt=0, spin=self.spin)
        return np.dot(dagger(self.U_nn) * eps_n, self.U_nn)


class LocFun(Wannier):
    def localize(self, calc, M=None, T=0, projections=None, ortho=True,
                 verbose=False):
        # M is size of Hilbert space to fix. Default is ~ number of occ. bands.
        if M is None:
            M = 0
            f_n = calc.get_occupation_numbers(0, self.spin)
            while f_n[M] > .01:
                M += 1

        if projections is None:
            projections = single_zeta(calc, self.spin, verbose=verbose)
        
        self.U_nn, self.S_jj = get_locfun_rotation(projections, M, T, ortho)
        if self.Z_nnc is None:
            self.value = 1
            return self.value
        
        self.Z_jjc = np.empty(self.S_jj.shape + (3,))
        for c in range(3):
            self.Z_jjc[:, :, c] = np.dot(dagger(self.U_nn),
                                     np.dot(self.Z_nnc[:, :, c], self.U_nn))
        
        self.value = np.sum(np.abs(self.Z_jjc.diagonal())**2)
        
        return self.value # / Bohr**6

    def get_centers(self):
        z_jjc = np.empty(self.S_jj.shape+(3,))
        for c in range(3):
            z_jjc = np.dot(dagger(self.U_nn),
                            np.dot(self.Z_nnc[:,:,c], self.U_nn))

        scaled_c = -np.angle(z_jjc.diagonal()).T / (2 * pi)
        return (scaled_c % 1.0) * self.cell_c

    def get_eigenstate_centers(self):
        scaled_c = -np.angle(self.Z_nnc.diagonal()).T / (2 * pi)
        return (scaled_c % 1.0) * self.cell_c
    
    def get_proj_norm(self, calc):
        return np.array([np.linalg.norm(U_j) for U_j in self.U_nn])
            

def get_locfun_rotation(projections_nj, M=None, T=0, ortho=False):
    """Mikkel Strange's localized functions.
    
    projections_nj = <psi_n|p_j>
    psi_n: eigenstates
    p_j: localized function
    Nw =  number of localized functions
    M = Number of fixed states
    T = Number of virtual states to exclude (from above) 
    """

    Nbands, Nw = projections_nj.shape
    if M is None:
        M = Nw
    L = Nw - M # Extra degrees of freedom
    V = Nbands - M - T# Virtual states
    a0_nj = projections_nj[:M, :]
    a0_vj = projections_nj[M:M + V, :]

    if V == 0:
        D_jj = np.dot(dagger(projections_nj), projections_nj)
        U_nj = 1.0 / np.sqrt(D_jj.diagonal()) * projections_nj
        S_jj = np.dot(dagger(U_nj), U_nj)
        assert np.diagonal(np.linalg.cholesky(S_jj)).min() > .01, \
               'Close to linear dependence.'
        if ortho:
            lowdin(U_nj, S_jj)
            S_jj = np.identity(len(S_jj))
        return U_nj, S_jj

    #b_v, b_vv = np.linalg.eigh(np.dot(a0_vj, dagger(a0_vj)))
    #T_vp = b_vv[:, np.argsort(-b_v)[:L]]
    b_j, b_jj = np.linalg.eigh(np.dot(dagger(a0_vj), a0_vj))
    T_vp = np.dot(a0_vj, b_jj[:, np.argsort(-b_j.real)[:L]])

    R_vv = np.dot(T_vp, dagger(T_vp))
    D_jj = np.dot(dagger(a0_nj), a0_nj) + np.dot(dagger(a0_vj),
                                                   np.dot(R_vv, a0_vj))
    D2_j = 1.0 / np.sqrt(D_jj.diagonal())
    ap_nj = D2_j * a0_nj
    ap_vj = D2_j * np.dot(R_vv, a0_vj)
    S_jj = np.dot(dagger(ap_nj), ap_nj) + np.dot(dagger(ap_vj), ap_vj)

    # Check for linear dependencies
    Scd = np.diagonal(np.linalg.cholesky(S_jj)).min()
    if Scd < 0.01:
        print(('Warning: possibly near linear dependence.\n'
               'Minimum eigenvalue of cholesky decomposition is %s' % Scd))

    if ortho:
        lowdin(ap_nj, S_jj)
        lowdin(ap_vj, S_jj)
        S_jj = np.identity(len(S_jj))

    U_nj = np.concatenate([ap_nj.flat, ap_vj.flat]).reshape(M+V, Nw)
    return U_nj, S_jj


def single_zeta(paw, spin, verbose=False):
    angular = [['1'],
               ['y', 'z', 'x'],
               ['xy', 'yz', '3z^2-r^2', 'xz', 'x^2-y^2'],
               ['3x^2y-y^3', 'xyz', '5yz^2-yr^2', '5z^3-3zr^2',
                '5xz^2-xr^2', 'x^2z-y^2z', 'x^3-3xy^2'],
               ]
    if verbose:
        print('index atom orbital')
    p_jn = []
    for a, P_ni in paw.wfs.kpt_u[spin].P_ani.items():
        setup = paw.wfs.setups[a]
        i = 0
        for l, n in zip(setup.l_j, setup.n_j):
            if n < 0:
                break
            for j in range(i, i + 2 * l + 1):
                p_jn.append(P_ni[:, j])
                if verbose:
                    print('%5i %4i %s_%s' % (len(p_jn), a,
                                             'spdf'[l], angular[l][j - i]))
            i += 2 * l + 1
    projections_nj = dagger(np.array(p_jn))
    assert projections_nj.shape[0] >= projections_nj.shape[1]
    return projections_nj

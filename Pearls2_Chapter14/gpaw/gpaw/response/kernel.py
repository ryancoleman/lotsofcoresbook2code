# -*- coding: utf-8
# In an attempt to appease epydoc and still have readable docstrings,
# the vertical bar | is represented by u'\u2758' in this module.
"""This module defines Coulomb and XC kernels for the response model.
"""

import numpy as np
#from math import sqrt, pi, sin, cos, exp
from gpaw.utilities.blas import gemmdot
from gpaw.xc import XC
from gpaw.sphere.lebedev import weight_n, R_nv
from gpaw.mpi import world, rank, size
from ase.dft.kpoints import monkhorst_pack


def v3D_Coulomb(qG):
    """Coulomb Potential in the 3D Periodic Case.
    
    Periodic calculation, no cutoff.

    ::
    
      v3D = 4 pi / |q+G|^2
    """
    
    v_q = 1. / (qG**2).sum(axis=1)
    return v_q


def v2D_Coulomb(qG, N_p, N_np, R):
    """Coulomb Potential in the 2D Periodic Case.

    Slab/Surface/Layer calculation.
    Cutoff in non-periodic N_np direction.
    No cutoff in periodic N_p directions.

    ::

      v2D = 4 pi/|q+G|^2 * [1 + exp(-|G_p|R)*[(G_n/|G_p|)*sin(G_n R)
                                                        - cos(G_n R)]
    """

    G_nR = qG[:, N_np[0]] * R
    G_pR = (qG[:, N_p[0]]**2 + qG[:, N_p[1]]**2)**0.5 * R

    v_q = 1. / (qG**2).sum(axis=1)
    v_q *= 1. + np.exp(-G_pR) * ((G_nR / G_pR) * np.sin(G_nR)
                                 - np.cos(G_nR))
    return v_q.astype(complex)


def v1D_Coulomb(qG, N_p, N_np, R):
    """Coulomb Potential in the 1D Periodic Case.
    
    Nanotube/Nanowire/Atomic Chain calculation.
    Cutoff in non-periodic N_np directions.
    No cutoff in periodic N_p direction.

    ::
    
      v1D = 4 pi/|q+G|^2 * [1 + |G_n|R J_1(|G_n|R) K_0(G_p R)
                              - G_p R J_0(|G_n| R) K_1(G_p R)]
    """

    from scipy.special import j1, k0, j0, k1

    G_nR = (qG[:, N_np[0]]**2 + qG[:, N_np[1]]**2)**0.5 * R
    G_pR = abs(qG[:, N_p[0]]) * R
    v_q = 1. / (qG**2).sum(axis=1)
    v_q *= (1. + G_nR * j1(G_nR) * k0(G_pR)
            - G_pR * j0(G_nR) * k1(G_pR))
    return v_q


def v0D_Coulomb(qG, R):
    """Coulomb Potential in the 0D Non-Periodic Case

    Isolated System/Molecule calculation.
    Spherical cutoff in all directions.

    ::

      v0D = 4 pi/|G|^2 * [1 - cos(|G|R)]
    """

    qGR = ((qG.sum(axis=1))**2)**0.5 * R
    v_q = 1. / (qG**2).sum(axis=1)
    v_q *= 1. - np.cos(qGR)
    return v_q


class BaseCoulombKernel:
    """This is an abstract class for calculating the Coulomb kernel.

    """
    def __init__(self, vcut):
        self.vcut = vcut

    def v_Coulomb(self, qG):
        """Returns the Coulomb potential in k-space.
        """
        pass

    def calculate_Kc(self,
                     q_c,
                     Gvec_Gc,
                     bcell_cv,
                     integrate_gamma=False,
                     N_k=None,
                     symmetric=True):
        """Symmetric Coulomb kernel.
        """
        pass

    def calculate_Kc_q(self,
                       bcell_cv,
                       N_k,
                       q_qc=None,
                       Gvec_c=None,
                       N=100):
        """What does this do? XXX
        """
        pass


class CoulombKernel3D(BaseCoulombKernel):
    """Class for calculating the Coulomb kernel for a 3D calculation.

    vcut = [False, False, False].
    """
    def v_Coulomb(self, qG):
        """Coulomb Potential in the 3D Periodic Case.
        Periodic calculation, no cutoff.
        """
        return v3D_Coulomb(qG)


class CoulombKernel2D(BaseCoulombKernel):
    """Class for calculating the Coulomb kernel for a 2D calculation.
    
    Slab/Surface/Layer calculation.
    vcut = [False, False, True].
    Cutoff in non-periodic N_np direction.
    No cutoff in periodic N_p directions.
    R = 1/2 max cell[N_np].
    """
    def __init__(self, vcut, cell):
        # Periodic Directions
        self.N_p = vcut.argsort()[:2]
        # Non-Periodic Direction
        self.N_np = vcut.argsort()[-1:]
        # 2D Radial Cutoff is one half the cell dimension in Bohr
        self.R = cell[self.N_np, self.N_np] / 2.
        BaseCoulombKernel.__init__(self, vcut)

    def v_Coulomb(self, qG):
        """Coulomb Potential in the 2D Periodic Case.

        Slab/Surface/Layer calculation.
        Cutoff in non-periodic N_np direction.
        No cutoff in periodic N_p directions.
        """
        
        return v2D_Coulomb(qG, self.N_p, self.N_np, self.R)


class CoulombKernel1D(BaseCoulombKernel):
    """Class for calculating the Coulomb kernel for a 1D calculation.
 
    Nanotube/Nanowire/Atomic Chain calculation.
    vcut = [False, True, True].
    Cutoff in non-periodic N_np directions.
    No cutoff in periodic N_p direction.
    R = 1/2 max cell[N_np].
    """
    def __init__(self, vcut, cell):
        from scipy.special import j1, k0, j0, k1
        # Periodic Direction
        self.N_p = vcut.argsort()[:1]
        # Non-Periodic Directions
        self.N_np = vcut.argsort()[-2:]
        # 2D Radial Cutoff is one half the cell dimension in Bohr
        self.R = max(cell[self.N_np[0], self.N_np[0]],
                     cell[self.N_np[1], self.N_np[1]]) / 2.
        self.j1 = j1
        self.k0 = k0
        self.j0 = j0
        self.k1 = k1
        BaseCoulombKernel.__init__(self, vcut)

    def v_Coulomb(self, qG):
        """Coulomb Potential in the 1D Periodic Case.

        Nanotube/Nanowire/Atomic Chain calculation.
        Cutoff in non-periodic N_np directions.
        No cutoff in periodic N_p direction.
        """

        return v1D_Coulomb(qG, self.N_p, self.N_np, self.R)


class CoulombKernel0D(BaseCoulombKernel):
    """Class for calculating the Coulomb kernel for a 0D calculation.
    
    Isolated System/Molecule calculation.
    vcut = [True, True, True].
    Spherical cutoff in all directions.
    R = 1/2 max cell.    
    """
    def __init__(self, vcut, cell):
        self.R = max(np.diag(cell)) / 2.
        BaseCoulombKernel.__init__(self, vcut)

    def v_Coulomb(self, qG):
        """Coulomb Potential in the 0D Non-Periodic Case

        Isolated System/Molecule calculation.
        Spherical cutoff in all directions.
        """

        return v0D_Coulomb(qG, self.R)


class CoulombKernel(BaseCoulombKernel):
    """Class for calculating the Coulomb kernel.

    vcut = None (Default).
    Applies Coulomb cutoff in non-periodic directions.
    vcut = '3D' => [False, False, False].
    vcut = '2D' => [False, False, True].
    vcut = '1D' => [False, True, True].
    vcut = '0D' => [True, True, True].
    vcut is a numpy arry of type bool of length 3.
    
    """
    def __init__(self, vcut, pbc, cell):
        if vcut is None:
            vcut = pbc - True
        elif isinstance(vcut, str):
            # Assume vcut = 'nD', where n is number of periodic directions
            # Should be deprecated
            n_p = int(vcut[0])
            vcut = np.ones(3, bool)
            vcut[:n_p] = False

        # 3D periodic case
        if vcut.sum() == 0:
            self.impl = CoulombKernel3D(vcut)
        # 2D periodic case
        elif vcut.sum() == 1:
            self.impl = CoulombKernel2D(vcut, cell)
        # 1D periodic case
        elif vcut.sum() == 2:
            self.impl = CoulombKernel1D(vcut, cell)
        # 0D periodic case
        elif vcut.sum() == 3:
            self.impl = CoulombKernel0D(vcut, cell)
        
        BaseCoulombKernel.__init__(self, vcut)

    def v_Coulomb(self, qG):
        """Returns vcut dependent Coulomb potential.
        
        Uses v_Coulomb as implemented in self.impl
        """
        
        return self.impl.v_Coulomb(qG)

    def calculate_Kc(self,
                     q_c,
                     Gvec_Gc,
                     bcell_cv,
                     integrate_gamma=False,
                     N_k=None,
                     symmetric=True):

        """Symmetric Coulomb kernel"""

        if integrate_gamma:
            assert N_k is not None
            if (q_c == 0).all():
                # Only to avoid dividing by zero below!
                # The term is handled separately later
                q_c = [1.e-12, 0., 0.]

        qG_c = np.zeros((len(Gvec_Gc), 3))
        qG_c[:, 0] = Gvec_Gc[:, 0] + q_c[0]
        qG_c[:, 1] = Gvec_Gc[:, 1] + q_c[1]
        qG_c[:, 2] = Gvec_Gc[:, 2] + q_c[2]

        qG = np.dot(qG_c, bcell_cv)

        # Calculate coulomb kernel
        Kc_G = self.v_Coulomb(qG)

        if integrate_gamma and np.abs(q_c).sum() < 1e-8:
            bvol = np.abs(np.linalg.det(bcell_cv)) / (N_k[0] * N_k[1] * N_k[2])
            R_q0 = (3. * bvol / (4. * np.pi))**(1. / 3.)
            Kc_G[0] = 4. * np.pi * R_q0 / bvol

        if symmetric:
            Kc_G = np.array(Kc_G, np.complex)**0.5
            Kc_G = np.outer(Kc_G.conj(), Kc_G)
            return 4. * np.pi * Kc_G
        else:
            return 4. * np.pi * (Kc_G * np.ones([len(Kc_G), len(Kc_G)])).T

    def calculate_Kc_q(self,
                       bcell_cv,
                       N_k,
                       q_qc=None,
                       Gvec_c=None,
                       N=100):

        if q_qc is None:
            q_qc = np.array([[1.e-12, 0., 0.]])
        if Gvec_c is None:
            Gvec_c = np.array([0, 0, 0])

        bcell = bcell_cv.copy()
        bcell[0] /= N_k[0]
        bcell[1] /= N_k[1]
        bcell[2] /= N_k[2]
        bvol = np.abs(np.linalg.det(bcell))

        gamma = np.where(np.sum(abs(q_qc), axis=1) < 1e-5)[0][0]
        if (Gvec_c == 0).all():
            q_qc = q_qc.copy()
            # Just to avoid zero division
            q_qc[gamma] = [0.001, 0., 0.]

        q_qc[:, 0] += Gvec_c[0]
        q_qc[:, 1] += Gvec_c[1]
        q_qc[:, 2] += Gvec_c[2]
        qG = np.dot(q_qc, bcell_cv)

        K0_q = self.v_Coulomb(qG)
        if (Gvec_c == 0).all():
            R_q0 = (3. * bvol / (4. * np.pi))**(1. / 3.)
            K0_q[gamma] = 4. * np.pi * R_q0 / bvol

        bvol = np.abs(np.linalg.det(bcell))
        # Integrated Coulomb kernel
        qt_qc = monkhorst_pack((N, N, N))
        qt_qcv = np.dot(qt_qc, bcell)
        dv = bvol / len(qt_qcv)
        K_q = np.zeros(len(q_qc), complex)

        for i in range(len(q_qc)):
            q = qt_qcv.copy() + qG[i]
            v_q = v3D_Coulomb(q)
            K_q[i] = np.sum(v_q) * dv / bvol

        return 4. * np.pi * K_q, 4. * np.pi * K0_q


def calculate_Kc(q_c,
                 Gvec_Gc,
                 acell_cv,
                 bcell_cv,
                 pbc,
                 vcut=None,
                 integrate_gamma=False,
                 N_k=None,
                 symmetric=True):

    """Symmetric Coulomb kernel"""

    if integrate_gamma:
        assert N_k is not None
        if (q_c == 0).all():
            # Only to avoid dividing by zero below!
            # The term is handled separately later
            q_c = [1.e-12, 0., 0.]

    #G_p = np.array(pbc, float)
    #G_n = np.array([1,1,1]) - G_p
    if vcut is None:
        vcut = '3D'
    G_p = np.array(pbc, bool).argsort()[-int(vcut[0]):]
    G_n = np.array(pbc, bool).argsort()[:int(vcut[0])]

    if vcut is None or vcut == '3D':
        pass
    elif vcut == '2D':
        # R is half the length of cell in non-periodic direction
        if pbc.all():  # G_n.sum() < 1e-8: # default dir is z
            G_n = np.array([2])  # np.array([0,0,1])
            G_p = np.array([0, 1])  # np.array([1,1,0])
        acell_n = acell_cv[G_n, G_n]
        R = max(acell_n) / 2.
    elif vcut == '1D':
        # R is the radius of the cylinder containing the cell.
        if pbc.all():  # G_n.sum() < 1e-8:
            raise ValueError('Check boundary condition ! ')
        acell_n = acell_cv
        R = max(acell_n[G_n[0], G_n[0]],
                acell_n[G_n[1], G_n[1]]) / 2.
    elif vcut == '0D':
        # R is the minimum radius of a sphere containing the cell.
        acell_n = acell_cv
        R = min(acell_n[0, 0],
                acell_n[1, 1],
                acell_n[2, 2]) / 2.
    else:
        NotImplemented

    qG_c = np.zeros((len(Gvec_Gc), 3))
    qG_c[:, 0] = Gvec_Gc[:, 0] + q_c[0]
    qG_c[:, 1] = Gvec_Gc[:, 1] + q_c[1]
    qG_c[:, 2] = Gvec_Gc[:, 2] + q_c[2]

    qG = np.dot(qG_c, bcell_cv)

    # Calculate Coulomb kernel
    if vcut is None or vcut == '3D':
        Kc_G = v3D_Coulomb(qG)
    elif vcut == '2D':
        Kc_G = v2D_Coulomb(qG, G_p, G_n, R)
    elif vcut == '1D':
        Kc_G = v1D_Coulomb(qG, G_p, G_n, R)
    elif vcut == '0D':
        Kc_G = v0D_Coulomb(qG, R)
    else:
        NotImplemented

    if integrate_gamma and np.abs(q_c).sum() < 1e-8:
        bvol = np.abs(np.linalg.det(bcell_cv)) / (N_k[0] * N_k[1] * N_k[2])
        R_q0 = (3. * bvol / (4. * np.pi))**(1. / 3.)
        Kc_G[0] = 4. * np.pi * R_q0 / bvol

    if symmetric:
        #print "Kc_G"
        #print Kc_G
        Kc_G = np.array(Kc_G, np.complex)**0.5
        #print "Kc_G^0.5"
        #print Kc_GC
        #print "Outer Product"
        Kc_G = np.outer(Kc_G.conj(), Kc_G)
        #print Kc_G
        return 4. * np.pi * Kc_G  # np.outer(Kc_G.conj(), Kc_G)
    else:
        return 4. * np.pi * (Kc_G * np.ones([len(Kc_G), len(Kc_G)])).T


def calculate_Kxc(gd,
                  nt_sG,
                  npw,
                  Gvec_Gc,
                  nG,
                  vol,
                  bcell_cv,
                  R_av,
                  setups,
                  D_asp,
                  functional='ALDA',
                  density_cut=None):
    """ALDA kernel"""

    # The soft part
    #assert np.abs(nt_sG[0].shape - nG).sum() == 0
    if functional == 'ALDA_X':
        x_only = True
        A_x = -3. / 4. * (3. / np.pi)**(1. / 3.)
        nspins = len(nt_sG)
        assert nspins in [1, 2]
        fxc_sg = nspins**(1. / 3.) * 4. / 9. * A_x * nt_sG**(-2. / 3.)
    else:
        assert len(nt_sG) == 1
        x_only = False
        fxc_sg = np.zeros_like(nt_sG)
        xc = XC(functional[1:])
        xc.calculate_fxc(gd, nt_sG, fxc_sg)

    if density_cut is not None:
        fxc_sg[np.where(nt_sG * len(nt_sG) < density_cut)] = 0.0

    # FFT fxc(r)
    nG0 = nG[0] * nG[1] * nG[2]
    tmp_sg = [np.fft.fftn(fxc_sg[s]) * vol / nG0 for s in range(len(nt_sG))]

    r_vg = gd.get_grid_point_coordinates()
    Kxc_sGG = np.zeros((len(fxc_sg), npw, npw), dtype=complex)
    for s in range(len(fxc_sg)):
        for iG in range(npw):
            for jG in range(npw):
                dG_c = Gvec_Gc[iG] - Gvec_Gc[jG]
                if (nG / 2 - np.abs(dG_c) > 0).all():
                    index = (dG_c + nG) % nG
                    Kxc_sGG[s, iG, jG] = tmp_sg[s][index[0],
                                                   index[1],
                                                   index[2]]
                else:  # not in the fft index
                    dG_v = np.dot(dG_c, bcell_cv)
                    dGr_g = gemmdot(dG_v, r_vg, beta=0.0)
                    Kxc_sGG[s, iG, jG] = gd.integrate(np.exp(-1j * dGr_g)
                                                      * fxc_sg[s])

    # The PAW part
    KxcPAW_sGG = np.zeros_like(Kxc_sGG)
    dG_GGv = np.zeros((npw, npw, 3))
    for iG in range(npw):
        for jG in range(npw):
            dG_c = Gvec_Gc[iG] - Gvec_Gc[jG]
            dG_GGv[iG, jG] = np.dot(dG_c, bcell_cv)

    for a, setup in enumerate(setups):
        if rank == a % size:
            rgd = setup.xc_correction.rgd
            n_qg = setup.xc_correction.n_qg
            nt_qg = setup.xc_correction.nt_qg
            nc_g = setup.xc_correction.nc_g
            nct_g = setup.xc_correction.nct_g
            Y_nL = setup.xc_correction.Y_nL
            dv_g = rgd.dv_g

            D_sp = D_asp[a]
            B_pqL = setup.xc_correction.B_pqL
            D_sLq = np.inner(D_sp, B_pqL.T)
            nspins = len(D_sp)

            f_sg = rgd.empty(nspins)
            ft_sg = rgd.empty(nspins)

            n_sLg = np.dot(D_sLq, n_qg)
            nt_sLg = np.dot(D_sLq, nt_qg)

            # Add core density
            n_sLg[:, 0] += np.sqrt(4. * np.pi) / nspins * nc_g
            nt_sLg[:, 0] += np.sqrt(4. * np.pi) / nspins * nct_g

            coefatoms_GG = np.exp(-1j * np.inner(dG_GGv, R_av[a]))
            for n, Y_L in enumerate(Y_nL):
                w = weight_n[n]
                f_sg[:] = 0.0
                n_sg = np.dot(Y_L, n_sLg)
                if x_only:
                    f_sg = nspins * (4 / 9.) * A_x * (nspins * n_sg)**(-2 / 3.)
                else:
                    xc.calculate_fxc(rgd, n_sg, f_sg)

                ft_sg[:] = 0.0
                nt_sg = np.dot(Y_L, nt_sLg)
                if x_only:
                    ft_sg = nspins * (4 / 9.) * (A_x
                                                 * (nspins * nt_sg)**(-2 / 3.))
                else:
                    xc.calculate_fxc(rgd, nt_sg, ft_sg)
                for i in range(len(rgd.r_g)):
                    coef_GG = np.exp(-1j * np.inner(dG_GGv, R_nv[n])
                                     * rgd.r_g[i])
                    for s in range(len(f_sg)):
                        KxcPAW_sGG[s] += w * np.dot(coef_GG,
                                                    (f_sg[s, i] -
                                                     ft_sg[s, i])
                                                    * dv_g[i]) * coefatoms_GG

    world.sum(KxcPAW_sGG)
    Kxc_sGG += KxcPAW_sGG

    return Kxc_sGG / vol


def calculate_Kc_q(acell_cv,
                   bcell_cv,
                   pbc,
                   N_k,
                   vcut=None,
                   integrated=True,
                   q_qc=np.array([[1.e-12, 0., 0.]]),
                   Gvec_c=np.array([0, 0, 0]),
                   N=100):
    # get cutoff parameters
    #G_p = np.array(pbc, int)
    # Normal Direction
    #G_n = np.array([1,1,1]) - G_p
    if vcut is None:
        vcut = '3D'
    G_p = np.array(pbc, bool).argsort()[-int(vcut[0]):]
    G_n = np.array(pbc, bool).argsort()[:int(vcut[0])]

    if vcut is None or vcut == '3D':
        pass
    elif vcut == '2D':
        if pbc.all():  # G_n.sum() < 1e-8: # default dir is z
            G_n = np.array([2])
            G_p = np.array([0, 1])
        acell_n = acell_cv[G_n, G_n]
        R = max(acell_n) / 2.
    elif vcut == '1D':
        # R is the radius of the cylinder containing the cell.
        if pbc.all():  # G_n.sum() < 1e-8:
            raise ValueError('Check boundary condition ! ')
        acell_n = acell_cv
        R = max(acell_n[G_n[0], G_n[0]], acell_n[G_n[1], G_n[1]]) / 2.
    elif vcut == '0D':
        # R is the minimum radius of a sphere containing the cell.
        acell_n = acell_cv
        R = min(acell_n[0, 0], acell_n[1, 1], acell_n[2, 2]) / 2.
    else:
        NotImplemented

    bcell = bcell_cv.copy()
    bcell[0] /= N_k[0]
    bcell[1] /= N_k[1]
    bcell[2] /= N_k[2]
    bvol = np.abs(np.linalg.det(bcell))

    gamma = np.where(np.sum(abs(q_qc), axis=1) < 1e-5)[0][0]

    if (Gvec_c[G_p] == 0).all():
        q_qc[gamma][G_p[0]] = 1.e-12  # Just to avoid zero division

    q_qc[:, 0] += Gvec_c[0]
    q_qc[:, 1] += Gvec_c[1]
    q_qc[:, 2] += Gvec_c[2]
    qG = np.dot(q_qc, bcell_cv)

    if vcut is None or vcut == '3D':
        K0_q = v3D_Coulomb(qG)
        if (Gvec_c == 0).all():
            R_q0 = (3. * bvol / (4. * np.pi))**(1. / 3.)
            K0_q[gamma] = 4. * np.pi * R_q0 / bvol
    elif vcut == '2D':
        K0_q = v2D_Coulomb(qG, G_p, G_n, R)
        if (Gvec_c == 0).all():
            R_q0 = (3 * bvol / (4. * np.pi))**(1. / 3.)
            K0_q[gamma] = (4. * np.pi * R_q0 / bvol)
    elif vcut == '1D':
        K0_q = v1D_Coulomb(qG, G_p, G_n, R)
    elif vcut == '0D':
        K0_q = v0D_Coulomb(qG, R)
    else:
        NotImplemented

    bvol = np.abs(np.linalg.det(bcell))
    # Integrated Coulomb kernel
    qt_qc = monkhorst_pack((N, N, N))
    qt_qcv = np.dot(qt_qc, bcell)
    dv = bvol / len(qt_qcv)
    K_q = np.zeros(len(q_qc), complex)

    for i in range(len(q_qc)):
        q = qt_qcv.copy() + qG[i]
        if vcut is None or vcut == '3D':
            v_q = v3D_Coulomb(q)
        elif vcut == '2D':
            v_q = v2D_Coulomb(q, G_p, G_n, R)
        elif vcut == '1D':
            v_q = v1D_Coulomb(q, G_p, G_n, R)
        elif vcut == '0D':
            v_q = v0D_Coulomb(q, R)
        else:
            NotImplemented
        K_q[i] = np.sum(v_q) * dv / bvol

    return 4. * np.pi * K_q, 4. * np.pi * K0_q

""" This module implements absorbing boundary conditions to minimize errors
caused by boundary reflections. The first classes use negative imaginary
potentials of different shapes which are known to absorb waves. The last
class PML uses the perfectly matched layer approach.

    For information about imaginary potentials please see:

    Takashi Nakatsukasa, Kazuhiro Yabana, J. Chem. Phys. 114, 2550 (2001);
    doi:10.1063/1.1338527
    Daniel Neuhasuer and Michael Baer,  J. Chem. Phys. 90, 4351 (1989);
    doi:10.1063/1.456646

    About PML:
    Chunxiong Zheng, Journal of Computational  Physics,
    volume 227, Issue 1(Nov 2007) pages 537-556;
    doi:10.1016/j.jcp.2007.08.004
"""

import numpy as np

from ase.units import Bohr

class DummyAbsorbingBoundary:
    """ Virtual and not usable class for abosrbing boundaries."""

    def __init__(self):
        self.v_imag = None
        self.type = None

    def set_up(self):
        pass

    def get_potential_matrix(self):
        return self.v_imag



class LinearAbsorbingBoundary(DummyAbsorbingBoundary):

    """
    This class uses linear negative imaginary potential for absorption.
    For distances larger than abc_range, the negative imaginary potential
    increases linearly. Positive float abc_strength describes the steepness
    of the slope. High values of abc_strength cause more reflection from
    the potential and small values do not absorb enough.

    Parameters:
      abc_range: Positive float number. Absorbing potential is zero
                 where distance from the middle point or chosen positions
                 in the grid is smaller and starts to increase linearly
                 after this point.
      abc_strength: Positive float number. A good value for this is usually
                    about 0.01

      positions: An array of positions in the grid. Each of these points
                 has abc_range of free space around them. If this is left
                 as None, the middle point of the grid is used. You can use
                 the position array of an ASE Atoms object here.
    """

    def __init__(self, abc_range, abc_strength, positions = None):
        DummyAbsorbingBoundary.__init__(self)
        self.abc_strength = abc_strength
        self.abc_r = abc_range / Bohr
        self.positions = positions / Bohr
        self.type = 'IPOT'


    def set_up(self, gd):
        """
        Creates the potential matrix self.v_imag.

        Parameters:

        gd: grid descriptor
        """

        assert gd.orthogonal

        #self.v_imag = np.zeros((gd.n_c[0],gd.n_c[1],gd.n_c[2]),dtype=complex)
        self.v_imag = gd.zeros(dtype=complex)

        # If positions array wasn't given, uses the middle point of the grid.
        if self.positions is None:
            self.positions=[np.array([gd.N_c[0] *gd.h_cv[0, 0] * 0.5, # middle
                                      gd.N_c[1] *gd.h_cv[1, 1] * 0.5,
                                      gd.N_c[2] *gd.h_cv[2, 2] * 0.5])]

        for i in range(gd.n_c[0]):
            x = (i + gd.beg_c[0]) * gd.h_cv[0, 0]
            for j in range(gd.n_c[1]):
                y = (j + gd.beg_c[1]) * gd.h_cv[1, 1]
                for k in range(gd.n_c[2]):
                    z = (k + gd.beg_c[2]) * gd.h_cv[2, 2]
                    position=np.array([x,y,z])
                    # Calculates the distance from the nearest chosen point
                    # in the grid
                    position_vectors = self.positions-position

                    r = np.linalg.norm(position_vectors[0])
                    for vector in position_vectors:
                        if np.linalg.norm(vector)<r:
                            r = np.linalg.norm(vector)

                    if r > self.abc_r:
                        self.v_imag[i][j][k] = \
                            -(0+1j) * self.abc_strength*(r-self.abc_r)




class P4AbsorbingBoundary(DummyAbsorbingBoundary):
    """
    The negative imaginary potential used by this class are 4th degree
    polynomials which are constructed so that the value is zero at the
    beginning and abc_strength at the end. The derivative is
    zero at the begininning and zero at the end.

    Parameters:
      abc_range: Positive float number. Absorbing potential is zero where
                 distance from the middle point or chosen positions in the
                 grid is smaller and starts to increase after this point.

      abc_strength: Positive float number.

      positions: An array of positions in the grid. Each of these points has
                 abc_range of free space around them. If this is left as None,
                 the middle point of the grid is used. You can use the position
                 array of an ASE Atoms object here. You have to define the
                 width parameter if you use this.

      width:   The width of the absorbing layer. If you don't define positions
               array a value for this is automatically generated so that the
               width is from abc_range to the end of the grid (works best for
               cube).
               If you use the atom-centered model you have to define the
               width yourself.
    """

    def __init__(self, abc_range, abc_strength, positions = None, width=None):
        DummyAbsorbingBoundary.__init__(self)
        self.abc_r = abc_range / Bohr
        self.abc_strength = abc_strength
        self.positions = positions / Bohr
        self.width = width / Bohr
        self.type = 'IPOT'

    def set_up(self,gd):

        #self.v_imag = np.zeros((gd.n_c[0],gd.n_c[1],gd.n_c[2]),dtype=complex)
        self.v_imag = gd.zeros(dtype=complex)
        vo = self.abc_strength

        # If positions array wasn't given, uses the middle point of the
        # grid as the center.
        if self.positions is None:
            self.positions=[np.array([gd.N_c[0] *gd.h_cv[0, 0] *0.5, # middle
                                      gd.N_c[1] *gd.h_cv[1, 1] *0.5,
                                      gd.N_c[2] *gd.h_cv[2, 2] * 0.5])]

        if self.width is None:
            self.width = np.linalg.norm(self.positions[0]) / np.sqrt(3) - self.abc_r

        for i in range(gd.n_c[0]):
            x = (i + gd.beg_c[0]) * gd.h_cv[0, 0]
            for j in range(gd.n_c[1]):
                y = (j + gd.beg_c[1]) * gd.h_cv[1, 1]
                for k in range(gd.n_c[2]):
                    z = (k + gd.beg_c[2]) * gd.h_cv[2, 2]
                    position=np.array([x,y,z])

                    position_vectors = self.positions-position

                    r = np.linalg.norm(position_vectors[0])
                    for vector in position_vectors:
                        if np.linalg.norm(vector)<r:
                            r = np.linalg.norm(vector)

                    if r > self.abc_r:
                        if r < self.abc_r+self.width:
                            self.v_imag[i][j][k] = (0+1j) \
                             * ((np.sqrt(self.abc_strength) -
                                 np.sqrt(self.abc_strength)/self.width**2
                                 * (r-self.abc_r)**2)**2 - self.abc_strength)
                        else:
                            self.v_imag[i][j][k] = -self.abc_strength*(0+1j)
                    #print i, j, k, self.v_imag[i][j][k]


class PML:
    """
    Important! You must use the BiCGStab solver or your time propagation
    will likely crash. Give the parameter solver='BiCGStab' as you create
    the TDDFT object.

    Using PML makes the time progation slower so reserve twice more time.
    As imaginary potential is usually almost equally good, this is mostly
    for testing (until improved?).

    And now for something completely different, a Perfectly matched layer.
    Note that you can't give positions array as parameter for this class.

    Parameters:
      abc_range:     Starting from the center of the grid, the amount of
                     free space before PML starts to affect.

      abc_strength:   Positive float, good amount should be around 0.1
    """

    def __init__(self,abc_range,abc_strength):

        self.abc_range = abc_range / Bohr
        self.abc_strength = abc_strength
        self.type = 'PML'


    def set_up(self,gd):
        r"""Set up matrices for PML.

        Creates the matrixes needed when the PML is applied in tdopers.py
        when applying time-dependent hamiltonian.

        Perfectly matched layer is applied as potential Vpml = Tpml-T,
        Where  T = -.5*\nabla^{2}\psi  (Use latex for equations)

        Tpml we get from

        'A perfectly matched layer approach to the nonlinear Schrodinger wave
        equations',
        Journal of Computational  Physics,
        volume 227, Issue 1(Nov 2007) pages 537-556,
        Author: Chunxiong Zheng

        T_{pml} = -0.5*(G\nabla G\nabla \psi + G^{2}\nabla^{2}\psi)

        where G = \frac{1}{1+R\sigma}

        V_{pml} = -0.5 * (G\nabla G\nabla \psi + (G^{2}-1)\nabla^{2}\psi)

        This is probably not the most optimal approach and slows the
        propagation.

        The matrixes created here are G and the gradients of G in all
        directions.

        """
        R = (0+1j) # Complex number, has to be in the first quadrant.
        self.G = gd.zeros(dtype=complex)
        self.G[:] = 1.0

        self.dG = gd.zeros(n=3, dtype=complex)


        r0=np.array([gd.N_c[0] * gd.h_cv[0, 0] * 0.5,
                     gd.N_c[1] * gd.h_cv[1, 1] * 0.5,
                     gd.N_c[2] * gd.h_cv[2, 2] * 0.5]) # middle point
        for i in range(gd.n_c[0]):
            x = (i + gd.beg_c[0]) * gd.h_cv[0, 0]
            for j in range(gd.n_c[1]):
                y = (j + gd.beg_c[1]) * gd.h_cv[1, 1]
                for k in range(gd.n_c[2]):
                    z = (k + gd.beg_c[2]) * gd.h_cv[2, 2]
                    position=np.array([x,y,z])
                    r = np.linalg.norm(position - r0)


                    if r > self.abc_range:
                        self.G[i][j][k] = (1.0 + R * self.abc_strength * (r-self.abc_range)**2)**-1.0
                        self.dG[0][i][j][k] = \
                          -(1.0+R*self.abc_strength*(r-self.abc_range)**2.0)**-2.0*2.0*R*self.abc_strength*(r-self.abc_range)*(x-r0[0])/r

                        self.dG[1][i][j][k] = -(1.0+R*self.abc_strength*(r-self.abc_range)**2.0)**-2.0*2.0*R*self.abc_strength*(r-self.abc_range)*(y-r0[1])/r

                        self.dG[2][i][j][k] = -(1.0+R*self.abc_strength*(r-self.abc_range)**2.0)**-2.0*2.0*R*self.abc_strength*(r-self.abc_range)*(z-r0[2])/r

    def get_G(self):
        return self.G

    def get_dG(self):
        return self.dG

    def isPML(self):
        return True



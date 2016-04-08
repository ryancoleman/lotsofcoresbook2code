import numpy as np

from ase.units import Bohr

"""This module defines different external potentials to be used in
time-independent and time-dependent calculations."""


class ExternalPotential:
    """ External potential

    """
    def __init__(self, vext_g=None, gd=None):
        """Initialize with a grid and the corresponding grid descriptor.

        Grid values should be in Hartree.
        """
        self.vext_g = vext_g
        self.gd = gd
        if self.gd is not None:
            assert gd.orthogonal
            if np.alltrue(vext_g.shape == gd.get_size_of_global_array()):
                # this is a global array and has to be distributed
                self.vext_g = self.gd.zeros()
                self.gd.distribute(vext_g, self.vext_g)

    def get_potential(self, gd=None):
        if self.gd is None:
            self.gd = gd
        else:
            if gd is not None:
                # make sure we are talking about the same grid
                assert gd == self.gd
        return self.vext_g

    def get_taylor(self, position=None, spos_c=None):
        """Get the Taylor expansion around a point

        position [Angstrom]"""
        # use only 0 order term, i.e. the value
        return [[self.get_value(position, spos_c)]]

    def get_value(self, position=None, spos_c=None):
        """The potential value (as seen by an electron)
        at a certain grid point.

        position [Angstrom]
        spos_c scaled position on the grid"""
        if spos_c is None:
            spos_c = self.gd.scale_position(position / Bohr)
        g_c = self.gd.get_nearest_grid_point(spos_c)
        g_c -= (g_c == self.gd.n_c)  # force point to this domain
        return self.vext_g[tuple(g_c)]

    def get_nuclear_energy(self, nucleus):
        """Return the energy contribution of the bare nucleus."""
        return 0.0  # don't assume anything about the nucleus

    def add_linear_field(self, wfs, spos_ac, a_nG, b_nG, strength, kpt):
        """Adds (does NOT apply) linear field
        f(x,y,z) = str_x * x + str_y * y + str_z * z to wavefunctions.

        Parameters
        ----------
        pt_nuclei: List of ?LocalizedFunctions?
            Projectors (paw.pt_nuclei)
        a_nG:
            the wavefunctions
        b_nG:
            the result
        strength: float[3]
            strength of the linear field
        kpt: KPoint
            K-point
        """

        gd = wfs.gd

        # apply local part of x to smooth wavefunctions psit_n
        for i in range(gd.n_c[0]):
            x = (i + gd.beg_c[0]) * gd.h_cv[0, 0]
            b_nG[:, i, :, :] += (strength[0] * x) * a_nG[:, i, :, :]

        # FIXME: combine y and z to one vectorized operation,
        # i.e., make yz-array and take its product with a_nG

        # apply local part of y to smooth wavefunctions psit_n
        for i in range(gd.n_c[1]):
            y = (i + gd.beg_c[1]) * gd.h_cv[1, 1]
            b_nG[:, :, i, :] += (strength[1] * y) * a_nG[:, :, i, :]

        # apply local part of z to smooth wavefunctions psit_n
        for i in range(gd.n_c[2]):
            z = (i + gd.beg_c[2]) * gd.h_cv[2, 2]
            b_nG[:, :, :, i] += (strength[2] * z) * a_nG[:, :, :, i]

        # apply the non-local part for each nucleus

        # number of wavefunctions, psit_nG
        n = len(a_nG)
        P_ani = wfs.pt.dict(n)
        wfs.pt.integrate(a_nG, P_ani, kpt.q)

        coef_ani = {}
        for a, P_ni in P_ani.items():
            c0 = np.dot(spos_ac[a] * gd.cell_cv.diagonal(), strength)
            cxyz = strength
            # calculate coefficient
            # ---------------------
            #
            # coeffs_ni =
            #   P_nj * c0 * 1_ij
            #   + P_nj * cx * x_ij
            #
            # where (see spherical_harmonics.py)
            #
            #   1_ij = sqrt(4pi) Delta_0ij
            #   y_ij = sqrt(4pi/3) Delta_1ij
            #   z_ij = sqrt(4pi/3) Delta_2ij
            #   x_ij = sqrt(4pi/3) Delta_3ij
            # ...

            Delta_iiL = wfs.setups[a].Delta_iiL

            #   1_ij = sqrt(4pi) Delta_0ij
            #   y_ij = sqrt(4pi/3) Delta_1ij
            #   z_ij = sqrt(4pi/3) Delta_2ij
            #   x_ij = sqrt(4pi/3) Delta_3ij
            oneij = np.sqrt(4 * np.pi) \
                * np.dot(P_ni, Delta_iiL[:, :, 0])
            yij = np.sqrt(4 * np.pi / 3) \
                * np.dot(P_ni, Delta_iiL[:, :, 1])
            zij = np.sqrt(4 * np.pi / 3) \
                * np.dot(P_ni, Delta_iiL[:, :, 2])
            xij = np.sqrt(4 * np.pi / 3) \
                * np.dot(P_ni, Delta_iiL[:, :, 3])

            # coefficients
            # coefs_ni = sum_j ( <phi_i| f(x,y,z) | phi_j>
            #                    - <phit_i| f(x,y,z) | phit_j> ) P_nj
            coef_ani[a] = (c0 * oneij +
                           cxyz[0] * xij + cxyz[1] * yij + cxyz[2] * zij)

        # add partial wave pt_nG to psit_nG with proper coefficient
        wfs.pt.add(b_nG, coef_ani, kpt.q)


class ConstantPotential(ExternalPotential):
    """Constant potential for tests."""
    def __init__(self, constant=1.):
        self.constant = constant
        ExternalPotential.__init__(self)

    def get_potential(self, gd):
        if self.vext_g is None:
            self.gd = gd
            self.vext_g = gd.zeros() + self.constant
        return self.vext_g

    def get_ion_energy_and_forces(self, atoms):
        """Return the ionic energy and force contribution."""
        forces = np.zeros((len(atoms), 3))
        energy = 0
        return energy, forces


class ElectrostaticPotential(ExternalPotential):
    """External electrostatic potential

    The action of the external potential on the nucleus is defined in the
    electrostatic case.
    """
    def get_ion_energy_and_forces(self, atoms):
        """Return the ionic energy and force contribution."""
        forces = np.zeros((len(atoms), 3))
        energy = 0
        for i, atom in enumerate(atoms):
            taylor = self.get_taylor(atom.position)
#            print "pos, taylor=", atom.position, taylor
            Z = atom.number
            energy -= Z * taylor[0][0]
            if len(taylor) > 1:
                # see spherical_harmonics.py for the assignment
                forces[i] += Z * np.array([taylor[1][2],   # x
                                           taylor[1][0],   # y
                                           taylor[1][1]])  # z
        return energy, forces


class ConstantElectricField(ElectrostaticPotential):
    """External constant electric field"""
    def __init__(self, strength, direction=[0, 0, 1], center=None):
        """
        strength: field strength [atomic units]
        direction: polarisation direction
        center: the center of zero field [Angstrom]
        """
        self.strength = strength
        if center is None:
            self.center = None
        else:
            self.center = np.array(center) / Bohr

        # normalise the direction
        dir = np.array(direction, float)
        dir /= np.sqrt(np.dot(dir, dir))
        self.direction = dir

    def get_potential(self, gd=None):
        """Get the potential on the grid."""

        if hasattr(self, 'potential'):
            if gd == self.gd or gd is None:
                # nothing changed
                return self.potential

        for c in range(3):
            if self.direction[c] != 0 and gd.pbc_c[c]:
                raise NotImplementedError('ConstantElectricField is ' +
                                          'not suitable for periodic ' +
                                          'boundary conditions.')

        self.gd = gd

        if self.center is None:
            # use the center of the grid as default
            self.center = 0.5 * gd.cell_cv.sum(0)

        r_Rv = gd.get_grid_point_coordinates().transpose((1, 2, 3, 0))
        self.potential = self.strength * np.dot(r_Rv - self.center,
                                                self.direction)
        return self.potential

    def get_taylor(self, position=None, spos_c=None):
        """Get the Taylor expansion around a point

        position [Angstrom]"""
        if position is None:
            gd = self.gd
            pos = np.dot(gd.N_c * spos_c, gd.h_cv) * Bohr
        else:
            pos = position
        # see spherical_harmonics.py for the assignment
        return [[self.get_value(position=pos)],
                [self.strength * self.direction[1],   # y
                 self.strength * self.direction[2],   # z
                 self.strength * self.direction[0]]]  # x

    def get_value(self, position=None, spos_c=None):
        """The potential value (as seen by an electron)
        at a certain grid point.

        position [Angstrom]
        spos_c scaled position on the grid"""
        gd = self.gd
        if position is None:
            vr = np.dot(gd.N_c * spos_c, gd.h_cv) - self.center
        else:
            vr = position / Bohr - self.center
        return self.strength * np.dot(vr, self.direction)

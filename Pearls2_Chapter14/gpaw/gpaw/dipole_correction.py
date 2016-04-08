import numpy as np

from gpaw.utilities import erf


class DipoleCorrection:
    """Dipole-correcting wrapper around another PoissonSolver."""
    def __init__(self, poissonsolver, direction):
        """Construct dipole correction object.

        poissonsolver is a GPAW Poisson solver.
        direction is 0, 1 or 2 and specifies x, y or z.
        """
        self.corrector = DipoleCorrector(direction)
        self.poissonsolver = poissonsolver

    def get_stencil(self):
        return self.poissonsolver.get_stencil()

    def set_grid_descriptor(self, gd):
        if gd.pbc_c[self.corrector.c]:
            raise ValueError('System must be non-periodic along '
                             'dipole correction axis')
        self.poissonsolver.set_grid_descriptor(gd)

    def get_description(self):
        poissondesc = self.poissonsolver.get_description()
        desc = 'Dipole correction along %s-axis' % 'xyz'[self.corrector.c]
        return '\n'.join([poissondesc, desc])

    def initialize(self):
        self.poissonsolver.initialize()

    def solve(self, phi, rho, **kwargs):
        gd = self.poissonsolver.gd
        drho, dphi = self.corrector.get_dipole_correction(gd, rho)
        phi -= dphi
        iters = self.poissonsolver.solve(phi, rho + drho, **kwargs)
        phi += dphi
        return iters

    def estimate_memory(self, mem):
        self.poissonsolver.estimate_memory(mem)


class DipoleCorrector:
    def __init__(self, direction):
        self.c = direction

    def get_dipole_correction(self, gd, rhot_g):
        """Get dipole corrections to charge and potential.

        Returns arrays drhot_g and dphit_g such that if rhot_g has the
        potential phit_g, then rhot_g + drhot_g has the potential
        phit_g + dphit_g, where dphit_g is an error function.

        The error function is chosen so as to be largely constant at the
        cell boundaries and beyond.
        """
        # This implementation is not particularly economical memory-wise

        c = self.c

        # Right now the dipole correction must be along one coordinate
        # axis and orthogonal to the two others.  The two others need not
        # be orthogonal to each other.
        for c1 in range(3):
            if c1 != c:
                if np.vdot(gd.cell_cv[c], gd.cell_cv[c1]) > 1e-12:
                    raise ValueError('Dipole correction axis must be '
                                     'orthogonal to the two other axes.')

        moment = gd.calculate_dipole_moment(rhot_g)[c]
        if abs(moment) < 1e-12:
            return gd.zeros(), gd.zeros()

        r_g = gd.get_grid_point_coordinates()[c]
        cellsize = abs(gd.cell_cv[c, c])
        sr_g = 2.0 / cellsize * r_g - 1.0  # sr ~ 'scaled r'
        alpha = 12.0  # should perhaps be variable
        drho_g = sr_g * np.exp(-alpha * sr_g**2)
        moment2 = gd.calculate_dipole_moment(drho_g)[c]
        factor = -moment / moment2
        drho_g *= factor
        phifactor = factor * (np.pi / alpha)**1.5 * cellsize**2 / 4.0
        dphi_g = -phifactor * erf(sr_g * np.sqrt(alpha))
        return drho_g, dphi_g

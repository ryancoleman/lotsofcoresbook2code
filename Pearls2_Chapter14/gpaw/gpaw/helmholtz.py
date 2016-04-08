from math import pi, sqrt, fabs
import numpy as np
from numpy import exp

from gpaw.poisson import PoissonSolver
from gpaw.utilities import cerf, erf
from gpaw.utilities.gauss import Gaussian
from gpaw.fd_operators import FDOperator, laplace
from gpaw.transformers import Transformer


class HelmholtzGaussian(Gaussian):
    def get_phi(self, k2):
        """Get the solution of the Helmholtz equation for a Gaussian."""
#   This should lead to very big errors
        r = np.sqrt(self.r2)
        k = sqrt(k2)
        sigma = 1. / sqrt(2 * self.a)
        p = sigma * k / sqrt(2)
        i = 1j
        rhop = r / sqrt(2) / sigma + i * p
        rhom = r / sqrt(2) / sigma - i * p
        h = np.sin(k * r)

        # the gpaw-Gaussian is sqrt(4 * pi) times a 3D normalized Gaussian
        return sqrt(4 * pi) * exp(-p**2) / r / 2 * (
            np.cos(k * r) * (cerf(rhop) + cerf(rhom)) +
            i * np.sin(k * r) * (2 + cerf(rhop) - cerf(rhom)))
        

class ScreenedPoissonGaussian(Gaussian):
    def get_phi(self, mu2):
        """Get the solution of the screened Poisson equation for a Gaussian.
        
           The Gaussian is centered to middle of grid-descriptor."""
        r = np.sqrt(self.r2)
        mu = sqrt(mu2)
        sigma = 1. / sqrt(2 * self.a)
        sig2 = sigma**2
        mrho = (sig2 * mu - r) / (sqrt(2) * sigma)
        prho = (sig2 * mu + r) / (sqrt(2) * sigma)

        def erfc(values):
            return 1. - erf(values)

        # the gpaw-Gaussian is sqrt(4 * pi) times a 3D normalized Gaussian
        return sqrt(4 * pi) * exp(sig2 * mu2 / 2.0) / (2 * r) * (
                exp(-mu * r) * erfc(mrho) - exp(mu * r) * erfc(prho))
        

class HelmholtzOperator(FDOperator):
    def __init__(self, gd, scale=1.0, n=1, dtype=float, k2=0.0):
        """Helmholtz for general non orthorhombic grid.

        gd: GridDescriptor
            Descriptor for grid.
        scale: float
            Scaling factor.  Use scale=-0.5 for a kinetic energy operator.
        n: int
            Range of stencil.  Stencil has O(h^(2n)) error.
        dtype: float or complex
            Datatype to work on.
        """

        # Order the 13 neighbor grid points:
        M_ic = np.indices((3, 3, 3)).reshape((3, -3)).T[-13:] - 1
        u_cv = gd.h_cv / (gd.h_cv**2).sum(1)[:, np.newaxis]**0.5
        u2_i = (np.dot(M_ic, u_cv)**2).sum(1)
        i_d = u2_i.argsort()

        m_mv = np.array([(2, 0, 0), (0, 2, 0), (0, 0, 2),
                         (0, 1, 1), (1, 0, 1), (1, 1, 0)])
        # Try 3, 4, 5 and 6 directions:
        for D in range(3, 7):
            h_dv = np.dot(M_ic[i_d[:D]], gd.h_cv)
            A_md = (h_dv**m_mv[:, np.newaxis, :]).prod(2)
            a_d, residual, rank, s = np.linalg.lstsq(A_md, [1, 1, 1, 0, 0, 0])
            if residual.sum() < 1e-14:
                assert rank == D, 'You have a weird unit cell!'
                # D directions was OK
                break

        a_d *= scale
        offsets = [(0, 0, 0)]
        coefs = [laplace[n][0] * a_d.sum()]
        coefs[0] += k2 * scale
        for d in range(D):
            M_c = M_ic[i_d[d]]
            offsets.extend(np.arange(1, n + 1)[:, np.newaxis] * M_c)
            coefs.extend(a_d[d] * np.array(laplace[n][1:]))
            offsets.extend(np.arange(-1, -n - 1, -1)[:, np.newaxis] * M_c)
            coefs.extend(a_d[d] * np.array(laplace[n][1:]))

        FDOperator.__init__(self, coefs, offsets, gd, dtype)
        
        self.description = (
            '%d*%d+1=%d point O(h^%d) finite-difference Helmholtz' %
            ((self.npoints - 1) // n, n, self.npoints, 2 * n))


class HelmholtzSolver(PoissonSolver):
    """Solve the Helmholtz or screened Poisson equations.

       The difference between the Helmholtz equation:

           (Laplace + k^2) phi = n

       and the screened Poisson equation:

           (Laplace - mu^2) phi = n

       is only the sign of the added inhomogenity. Because of
       this we can use one class to solve both. So if k2 is
       greater zero we'll try to solve the Helmhlotz equation,
       otherwise we'll try to solve the screened Poisson equation.
    """

    def __init__(self, k2=0.0, nn='M', relax='GS', eps=2e-10):
        assert k2 <= 0, 'Currently only defined for k^2<=0'
        PoissonSolver.__init__(self, nn, relax, eps)
        self.k2 = k2
        
    def set_grid_descriptor(self, gd):
        # Should probably be renamed initialize
        self.gd = gd
        self.dv = gd.dv

        gd = self.gd
        scale = -0.25 / pi

        if self.nn == 'M':
            raise ValueError(
                'Helmholtz not defined for Mehrstellen stencil')
        self.operators = [HelmholtzOperator(gd, scale, self.nn, k2=self.k2)]
        self.B = None

        self.interpolators = []
        self.restrictors = []

        level = 0
        self.presmooths = [2]
        self.postsmooths = [1]

        # Weights for the relaxation,
        # only used if 'J' (Jacobi) is chosen as method
        self.weights = [2.0 / 3.0]

        while level < 4:
            try:
                gd2 = gd.coarsen()
            except ValueError:
                break
            self.operators.append(HelmholtzOperator(gd2, scale, 1,
                                                    k2=self.k2))
            self.interpolators.append(Transformer(gd2, gd))
            self.restrictors.append(Transformer(gd, gd2))
            self.presmooths.append(4)
            self.postsmooths.append(4)
            self.weights.append(1.0)
            level += 1
            gd = gd2

        self.levels = level

        if self.relax_method == 1:
            self.description = 'Gauss-Seidel'
        else:
            self.description = 'Jacobi'
        self.description += ' solver with %d multi-grid levels' % (level + 1)
        self.description += '\nStencil: ' + self.operators[0].description

    def load_gauss(self):
        """Load the gaussians."""
        if not hasattr(self, 'rho_gauss'):
            if self.k2 > 0:
                gauss = HelmholtzGaussian(self.gd)
            else:
                gauss = ScreenedPoissonGaussian(self.gd)
            self.rho_gauss = gauss.get_gauss(0)
            self.phi_gauss = gauss.get_phi(abs(self.k2))

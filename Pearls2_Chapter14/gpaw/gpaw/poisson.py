# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

from math import pi

import numpy as np
from numpy.fft import fftn, ifftn, fft2, ifft2

from gpaw.transformers import Transformer
from gpaw.fd_operators import Laplace, LaplaceA, LaplaceB
from gpaw import PoissonConvergenceError
from gpaw.utilities.blas import axpy
from gpaw.utilities.gauss import Gaussian
from gpaw.utilities.ewald import madelung
from gpaw.utilities.tools import construct_reciprocal
import _gpaw


class PoissonSolver:
    def __init__(self, nn=3, relax='J', eps=2e-10, maxiter=1000,
                 remove_moment=None):
        self.relax = relax
        self.nn = nn
        self.eps = eps
        self.charged_periodic_correction = None
        self.maxiter = maxiter
        self.remove_moment = remove_moment

        # Relaxation method
        if relax == 'GS':
            # Gauss-Seidel
            self.relax_method = 1
        elif relax == 'J':
            # Jacobi
            self.relax_method = 2
        else:
            raise NotImplementedError('Relaxation method %s' % relax)
        
        self.description = None

    def get_stencil(self):
        return self.nn

    def set_grid_descriptor(self, gd):
        # Should probably be renamed initialize
        self.gd = gd
        self.dv = gd.dv

        gd = self.gd
        scale = -0.25 / pi

        if self.nn == 'M':
            if not gd.orthogonal:
                raise RuntimeError('Cannot use Mehrstellen stencil with '
                                   'non orthogonal cell.')

            self.operators = [LaplaceA(gd, -scale)]
            self.B = LaplaceB(gd)
        else:
            self.operators = [Laplace(gd, scale, self.nn)]
            self.B = None

        self.interpolators = []
        self.restrictors = []

        level = 0
        self.presmooths = [2]
        self.postsmooths = [1]

        # Weights for the relaxation,
        # only used if 'J' (Jacobi) is chosen as method
        self.weights = [2.0 / 3.0]

        while level < 8:
            try:
                gd2 = gd.coarsen()
            except ValueError:
                break
            self.operators.append(Laplace(gd2, scale, 1))
            self.interpolators.append(Transformer(gd2, gd))
            self.restrictors.append(Transformer(gd, gd2))
            self.presmooths.append(4)
            self.postsmooths.append(4)
            self.weights.append(1.0)
            level += 1
            gd = gd2

        self.levels = level

    def get_description(self):
        name = {1: 'Gauss-Seidel', 2: 'Jacobi'}[self.relax_method]
        coarsest_grid = self.operators[-1].gd.N_c
        coarsest_grid_string = ' x '.join([str(N) for N in coarsest_grid])
        assert self.levels + 1 == len(self.operators)
        lines = ['%s solver with %d multi-grid levels'
                 % (name, self.levels + 1),
                 '    Coarsest grid: %s points' % coarsest_grid_string]
        if coarsest_grid.max() > 24:
            lines.extend(['    Warning: Coarse grid has more than 24 points.',
                          '             More multi-grid levels recommended.'])
        lines.extend(['    Stencil: %s' % self.operators[0].description,
                      '    Tolerance: %e' % self.eps,
                      '    Max iterations: %d' % self.maxiter])
        return '\n'.join(lines)

    def initialize(self, load_gauss=False):
        # Should probably be renamed allocate
        gd = self.gd
        self.rhos = [gd.empty()]
        self.phis = [None]
        self.residuals = [gd.empty()]
        for level in range(self.levels):
            gd2 = gd.coarsen()
            self.phis.append(gd2.empty())
            self.rhos.append(gd2.empty())
            self.residuals.append(gd2.empty())
            gd = gd2
        assert len(self.phis) == len(self.rhos)
        level += 1
        assert level == self.levels

        self.step = 0.66666666 / self.operators[0].get_diagonal_element()
        self.presmooths[level] = 8
        self.postsmooths[level] = 8

        if load_gauss:
            self.load_gauss()

    def load_gauss(self):
        if not hasattr(self, 'rho_gauss'):
            gauss = Gaussian(self.gd)
            self.rho_gauss = gauss.get_gauss(0)
            self.phi_gauss = gauss.get_gauss_pot(0)

    def solve(self, phi, rho, charge=None, eps=None, maxcharge=1e-6,
              zero_initial_phi=False):

        if eps is None:
            eps = self.eps
        actual_charge = self.gd.integrate(rho)
        background = (actual_charge / self.gd.dv /
                      self.gd.get_size_of_global_array().prod())

        if self.remove_moment:
            assert not self.gd.pbc_c.any()
            if not hasattr(self, 'gauss'):
                self.gauss = Gaussian(self.gd)
            rho_neutral = rho.copy()
            phi_cor_L = []
            for L in range(self.remove_moment):
                phi_cor_L.append(self.gauss.remove_moment(rho_neutral, L))
            # Remove multipoles for better initial guess
            for phi_cor in phi_cor_L:
                phi -= phi_cor

            niter = self.solve_neutral(phi, rho_neutral, eps=eps)
            # correct error introduced by removing multipoles
            for phi_cor in phi_cor_L:
                phi += phi_cor

            return niter
        if charge is None:
            charge = actual_charge
        if abs(charge) <= maxcharge:
            # System is charge neutral. Use standard solver
            return self.solve_neutral(phi, rho - background, eps=eps)

        elif abs(charge) > maxcharge and self.gd.pbc_c.all():
            # System is charged and periodic. Subtract a homogeneous
            # background charge

            # Set initial guess for potential
            if zero_initial_phi:
                phi[:] = 0.0

            iters = self.solve_neutral(phi, rho - background, eps=eps)
            return iters

        elif abs(charge) > maxcharge and not self.gd.pbc_c.any():
            # The system is charged and in a non-periodic unit cell.
            # Determine the potential by 1) subtract a gaussian from the
            # density, 2) determine potential from the neutralized density
            # and 3) add the potential from the gaussian density.

            # Load necessary attributes
            self.load_gauss()

            # Remove monopole moment
            q = actual_charge / np.sqrt(4 * pi)  # Monopole moment
            rho_neutral = rho - q * self.rho_gauss  # neutralized density

            # Set initial guess for potential
            if zero_initial_phi:
                phi[:] = 0.0
            else:
                axpy(-q, self.phi_gauss, phi)  # phi -= q * self.phi_gauss

            # Determine potential from neutral density using standard solver
            niter = self.solve_neutral(phi, rho_neutral, eps=eps)

            # correct error introduced by removing monopole
            axpy(q, self.phi_gauss, phi)  # phi += q * self.phi_gauss

            return niter
        else:
            # System is charged with mixed boundaryconditions
            msg = 'Charged systems with mixed periodic/zero'
            msg += ' boundary conditions'
            raise NotImplementedError(msg)

    def solve_neutral(self, phi, rho, eps=2e-10):
        self.phis[0] = phi

        if self.B is None:
            self.rhos[0][:] = rho
        else:
            self.B.apply(rho, self.rhos[0])

        niter = 1
        maxiter = self.maxiter
        while self.iterate2(self.step) > eps and niter < maxiter:
            niter += 1
        if niter == maxiter:
            msg = 'Poisson solver did not converge in %d iterations!' % maxiter
            raise PoissonConvergenceError(msg)

        # Set the average potential to zero in periodic systems
        if np.alltrue(self.gd.pbc_c):
            phi_ave = self.gd.comm.sum(np.sum(phi.ravel()))
            N_c = self.gd.get_size_of_global_array()
            phi_ave /= np.product(N_c)
            phi -= phi_ave

        return niter

    def iterate(self, step, level=0):
        residual = self.residuals[level]
        niter = 0
        while True:
            niter += 1
            if level > 0 and niter == 1:
                residual[:] = -self.rhos[level]
            else:
                self.operators[level].apply(self.phis[level], residual)
                residual -= self.rhos[level]
            error = self.gd.comm.sum(np.vdot(residual, residual))
            if niter == 1 and level < self.levels:
                self.restrictors[level].apply(residual, self.rhos[level + 1])
                self.phis[level + 1][:] = 0.0
                self.iterate(4.0 * step, level + 1)
                self.interpolators[level].apply(self.phis[level + 1], residual)
                self.phis[level] -= residual
                continue
            residual *= step
            self.phis[level] -= residual
            if niter == 2:
                break

        return error

    def iterate2(self, step, level=0):
        """Smooths the solution in every multigrid level"""

        residual = self.residuals[level]

        if level < self.levels:
            self.operators[level].relax(self.relax_method,
                                        self.phis[level],
                                        self.rhos[level],
                                        self.presmooths[level],
                                        self.weights[level])

            self.operators[level].apply(self.phis[level], residual)
            residual -= self.rhos[level]
            self.restrictors[level].apply(residual,
                                          self.rhos[level + 1])
            self.phis[level + 1][:] = 0.0
            self.iterate2(4.0 * step, level + 1)
            self.interpolators[level].apply(self.phis[level + 1], residual)
            self.phis[level] -= residual

        self.operators[level].relax(self.relax_method,
                                    self.phis[level],
                                    self.rhos[level],
                                    self.postsmooths[level],
                                    self.weights[level])
        if level == 0:
            self.operators[level].apply(self.phis[level], residual)
            residual -= self.rhos[level]
            error = self.gd.comm.sum(np.dot(residual.ravel(),
                                            residual.ravel())) * self.dv
            return error

    def estimate_memory(self, mem):
        # XXX Memory estimate works only for J and GS, not FFT solver
        # Poisson solver appears to use same amount of memory regardless
        # of whether it's J or GS, which is a bit strange

        gdbytes = self.gd.bytecount()
        nbytes = -gdbytes  # No phi on finest grid, compensate ahead
        for level in range(self.levels):
            nbytes += 3 * gdbytes  # Arrays: rho, phi, residual
            gdbytes //= 8
        mem.subnode('rho, phi, residual [%d levels]' % self.levels, nbytes)

    def __repr__(self):
        template = 'PoissonSolver(relax=\'%s\', nn=%s, eps=%e)'
        representation = template % (self.relax, repr(self.nn), self.eps)
        return representation


class NoInteractionPoissonSolver:
    relax_method = 0
    nn = 1

    def get_description(self):
        return 'No interaction'
    
    def get_stencil(self):
        return 1

    def solve(self, phi, rho, charge):
        return 0

    def set_grid_descriptor(self, gd):
        pass

    def initialize(self):
        pass


class FFTPoissonSolver(PoissonSolver):
    """FFT Poisson solver for general unit cells."""

    relax_method = 0
    nn = 999

    def __init__(self, eps=2e-10):
        self.charged_periodic_correction = None
        self.remove_moment = None
        self.eps = eps

    def get_description(self):
        return 'FFT'

    def set_grid_descriptor(self, gd):
        assert gd.pbc_c.all()
        self.gd = gd

    def initialize(self):
        if self.gd.comm.rank == 0:
            self.k2_Q, self.N3 = construct_reciprocal(self.gd)

    def solve_neutral(self, phi_g, rho_g, eps=None):
        if self.gd.comm.size == 1:
            phi_g[:] = ifftn(fftn(rho_g) * 4.0 * pi / self.k2_Q).real
        else:
            rho_g = self.gd.collect(rho_g)
            if self.gd.comm.rank == 0:
                globalphi_g = ifftn(fftn(rho_g) * 4.0 * pi / self.k2_Q).real
            else:
                globalphi_g = None
            self.gd.distribute(globalphi_g, phi_g)
        return 1


class FixedBoundaryPoissonSolver(PoissonSolver):
    """Solve the Poisson equation with FFT in two directions,
    and with central differential method in the third direction."""

    def __init__(self, nn=1):
        # XXX How can this work when it does not call __init___
        # on PoissonSolver? -askhl
        self.nn = nn
        self.charged_periodic_correction = None
        assert self.nn == 1

    def set_grid_descriptor(self, gd):
        assert gd.pbc_c.all()
        assert gd.orthogonal
        self.gd = gd

    def initialize(self, b_phi1, b_phi2):
        distribution = np.zeros([self.gd.comm.size], int)
        if self.gd.comm.rank == 0:
            gd = self.gd
            N_c1 = gd.N_c[:2, np.newaxis]
            i_cq = np.indices(gd.N_c[:2]).reshape((2, -1))
            i_cq += N_c1 // 2
            i_cq %= N_c1
            i_cq -= N_c1 // 2
            B_vc = 2.0 * np.pi * gd.icell_cv.T[:2, :2]
            k_vq = np.dot(B_vc, i_cq)
            k_vq *= k_vq
            k_vq2 = np.sum(k_vq, axis=0)
            k_vq2 = k_vq2.reshape(-1)

            b_phi1 = fft2(b_phi1, None, (0, 1))
            b_phi2 = fft2(b_phi2, None, (0, 1))

            b_phi1 = b_phi1[:, :, -1].reshape(-1)
            b_phi2 = b_phi2[:, :, 0].reshape(-1)

            loc_b_phi1 = np.array_split(b_phi1, self.gd.comm.size)
            loc_b_phi2 = np.array_split(b_phi2, self.gd.comm.size)
            loc_k_vq2 = np.array_split(k_vq2, self.gd.comm.size)

            self.loc_b_phi1 = loc_b_phi1[0]
            self.loc_b_phi2 = loc_b_phi2[0]
            self.k_vq2 = loc_k_vq2[0]

            for i in range(self.gd.comm.size):
                distribution[i] = len(loc_b_phi1[i])
            self.gd.comm.broadcast(distribution, 0)

            for i in range(1, self.gd.comm.size):
                self.gd.comm.ssend(loc_b_phi1[i], i, 135)
                self.gd.comm.ssend(loc_b_phi2[i], i, 246)
                self.gd.comm.ssend(loc_k_vq2[i], i, 169)
        else:
            self.gd.comm.broadcast(distribution, 0)
            self.loc_b_phi1 = np.zeros([distribution[self.gd.comm.rank]],
                                       dtype=complex)
            self.loc_b_phi2 = np.zeros([distribution[self.gd.comm.rank]],
                                       dtype=complex)
            self.k_vq2 = np.zeros([distribution[self.gd.comm.rank]])
            self.gd.comm.receive(self.loc_b_phi1, 0, 135)
            self.gd.comm.receive(self.loc_b_phi2, 0, 246)
            self.gd.comm.receive(self.k_vq2, 0, 169)

        k_distribution = np.arange(np.sum(distribution))
        self.k_distribution = np.array_split(k_distribution,
                                             self.gd.comm.size)

        self.d1, self.d2, self.d3 = self.gd.N_c
        self.r_distribution = np.array_split(np.arange(self.d3),
                                             self.gd.comm.size)
        self.comm_reshape = not (self.gd.parsize_c[0] == 1
                                 and self.gd.parsize_c[1] == 1)

    def solve(self, phi_g, rho_g, charge=None):
        if charge is None:
            actual_charge = self.gd.integrate(rho_g)
        else:
            actual_charge = charge

        if self.charged_periodic_correction is None:
            self.charged_periodic_correction = madelung(self.gd.cell_cv)

        background = (actual_charge / self.gd.dv /
                      self.gd.get_size_of_global_array().prod())

        self.solve_neutral(phi_g, rho_g - background)
        phi_g += actual_charge * self.charged_periodic_correction

    def scatter_r_distribution(self, global_rho_g, dtype=float):
        d1 = self.d1
        d2 = self.d2
        comm = self.gd.comm
        index = self.r_distribution[comm.rank]
        if comm.rank == 0:
            rho_g1 = global_rho_g[:, :, index]
            for i in range(1, comm.size):
                ind = self.r_distribution[i]
                comm.ssend(global_rho_g[:, :, ind].copy(), i, 178)
        else:
            rho_g1 = np.zeros([d1, d2, len(index)], dtype=dtype)
            comm.receive(rho_g1, 0, 178)
        return rho_g1

    def gather_r_distribution(self, rho_g, dtype=complex):
        comm = self.gd.comm
        index = self.r_distribution[comm.rank]
        d1, d2, d3 = self.d1, self.d2, self.d3
        if comm.rank == 0:
            global_rho_g = np.zeros([d1, d2, d3], dtype)
            global_rho_g[:, :, index] = rho_g
            for i in range(1, comm.size):
                ind = self.r_distribution[i]
                rho_gi = np.zeros([d1, d2, len(ind)], dtype)
                comm.receive(rho_gi, i, 368)
                global_rho_g[:, :, ind] = rho_gi
        else:
            comm.ssend(rho_g, 0, 368)
            global_rho_g = None
        return global_rho_g

    def scatter_k_distribution(self, global_rho_g):
        comm = self.gd.comm
        index = self.k_distribution[comm.rank]
        if comm.rank == 0:
            rho_g = global_rho_g[index]
            for i in range(1, comm.size):
                ind = self.k_distribution[i]
                comm.ssend(global_rho_g[ind], i, 370)
        else:
            rho_g = np.zeros([len(index), self.d3], dtype=complex)
            comm.receive(rho_g, 0, 370)
        return rho_g

    def gather_k_distribution(self, phi_g):
        comm = self.gd.comm
        index = self.k_distribution[comm.rank]
        d12 = self.d1 * self.d2
        if comm.rank == 0:
            global_phi_g = np.zeros([d12, self.d3], dtype=complex)
            global_phi_g[index] = phi_g
            for i in range(1, comm.size):
                ind = self.k_distribution[i]
                phi_gi = np.zeros([len(ind), self.d3], dtype=complex)
                comm.receive(phi_gi, i, 569)
                global_phi_g[ind] = phi_gi
        else:
            comm.ssend(phi_g, 0, 569)
            global_phi_g = None
        return global_phi_g

    def solve_neutral(self, phi_g, rho_g):
        # b_phi1 and b_phi2 are the boundary Hartree potential values
        # of left and right sides

        if self.comm_reshape:
            global_rho_g0 = self.gd.collect(rho_g)
            rho_g1 = self.scatter_r_distribution(global_rho_g0)
        else:
            rho_g1 = rho_g

        # use copy() to avoid the C_contiguous=False
        rho_g2 = fft2(rho_g1, None, (0, 1)).copy()

        global_rho_g = self.gather_r_distribution(rho_g2)
        if self.gd.comm.rank == 0:
            global_rho_g.shape = (self.d1 * self.d2, self.d3)
        rho_g3 = self.scatter_k_distribution(global_rho_g)

        du0 = np.zeros(self.d3 - 1, dtype=complex)
        du20 = np.zeros(self.d3 - 2, dtype=complex)
        h2 = self.gd.h_cv[2, 2] ** 2

        phi_g1 = np.zeros(rho_g3.shape, dtype=complex)
        index = self.k_distribution[self.gd.comm.rank]
        for phi, rho, rv2, bp1, bp2, i in zip(phi_g1, rho_g3,
                                              self.k_vq2,
                                              self.loc_b_phi1,
                                              self.loc_b_phi2,
                                              range(len(index))):
            A = np.zeros(self.d3, dtype=complex) + 2 + h2 * rv2
            phi = rho * np.pi * 4 * h2
            phi[0] += bp1
            phi[-1] += bp2
            du = du0 - 1
            dl = du0 - 1
            du2 = du20 - 1
            _gpaw.linear_solve_tridiag(self.d3, A, du, dl, du2, phi)
            phi_g1[i] = phi

        global_phi_g = self.gather_k_distribution(phi_g1)
        if self.gd.comm.rank == 0:
            global_phi_g.shape = (self.d1, self.d2, self.d3)
        phi_g2 = self.scatter_r_distribution(global_phi_g, dtype=complex)
        # use copy() to avoid the C_contiguous=False
        phi_g3 = ifft2(phi_g2, None, (0, 1)).real.copy()
        if self.comm_reshape:
            global_phi_g = self.gather_r_distribution(phi_g3, dtype=float)
            self.gd.distribute(global_phi_g, phi_g)
        else:
            phi_g[:] = phi_g3

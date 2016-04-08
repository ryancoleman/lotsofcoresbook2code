# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""Grid-descriptors

This module contains classes defining two kinds of grids:

* Uniform 3D grids.
* Radial grids.
"""

from math import pi

import numpy as np

import _gpaw
import gpaw.mpi as mpi
from gpaw.domain import Domain
from gpaw.utilities import mlsqr
from gpaw.utilities.blas import rk, r2k, gemv, gemm
from gpaw.mic.micblas import gemm as mic_gemm
from gpaw.mic.micblas import rk as mic_rk
from gpaw.mic.micblas import r2k as mic_r2k
from gpaw.mic import stream

import pymic as mic

# Remove this:  XXX
assert (-1) % 3 == 2
assert (np.array([-1]) % 3)[0] == 2

NONBLOCKING = False


class GridBoundsError(ValueError):
    pass


class GridDescriptor(Domain):
    """Descriptor-class for uniform 3D grid

    A ``GridDescriptor`` object holds information on how functions, such
    as wave functions and electron densities, are discreticed in a
    certain domain in space.  The main information here is how many
    grid points are used in each direction of the unit cell.

    There are methods for tasks such as allocating arrays, performing
    symmetry operations and integrating functions over space.  All
    methods work correctly also when the domain is parallelized via
    domain decomposition.

    This is how a 2x2x2 3D array is laid out in memory::

        3-----7
        |\    |\
        | \   | \
        |  1-----5      z
        2--|--6  |   y  |
         \ |   \ |    \ |
          \|    \|     \|
           0-----4      +-----x

    Example:

     >>> a = np.zeros((2, 2, 2))
     >>> a.ravel()[:] = range(8)
     >>> a
     array([[[0, 1],
             [2, 3]],
            [[4, 5],
             [6, 7]]])
     """

    ndim = 3  # dimension of ndarrays

    def __init__(self, N_c, cell_cv=(1, 1, 1), pbc_c=True,
                 comm=None, parsize=None):
        """Construct grid-descriptor object.

        parameters:

        N_c: 3 ints
            Number of grid points along axes.
        cell_cv: 3 float's or 3x3 floats
            Unit cell.
        pbc_c: one or three bools
            Periodic boundary conditions flag(s).
        comm: MPI-communicator
            Communicator for domain-decomposition.
        parsize: tuple of 3 ints, a single int or None
            Number of domains.

        Note that if pbc_c[c] is True, then the actual number of gridpoints
        along axis c is one less than N_c[c].

        Attributes:

        ==========  ========================================================
        ``dv``      Volume per grid point.
        ``h_cv``    Array of the grid spacing along the three axes.
        ``N_c``     Array of the number of grid points along the three axes.
        ``n_c``     Number of grid points on this CPU.
        ``beg_c``   Beginning of grid-point indices (inclusive).
        ``end_c``   End of grid-point indices (exclusive).
        ``comm``    MPI-communicator for domain decomposition.
        ==========  ========================================================

        The length unit is Bohr.
        """

        if isinstance(pbc_c, int):
            pbc_c = (pbc_c,) * 3
        if comm is None:
            comm = mpi.world

        self.N_c = np.array(N_c, int)
        if (self.N_c != N_c).any():
            raise ValueError('Non-int number of grid points %s' % N_c)
        
        Domain.__init__(self, cell_cv, pbc_c, comm, parsize, self.N_c)
        self.rank = self.comm.rank

        parsize_c = self.parsize_c
        n_c, remainder_c = divmod(self.N_c, parsize_c)

        self.beg_c = np.empty(3, int)
        self.end_c = np.empty(3, int)

        self.n_cp = []
        for c in range(3):
            n_p = (np.arange(parsize_c[c] + 1) * float(self.N_c[c]) /
                   parsize_c[c])
            n_p = np.around(n_p + 0.4999).astype(int)
            
            if not self.pbc_c[c]:
                n_p[0] = 1

            if not np.alltrue(n_p[1:] - n_p[:-1]):
                raise ValueError('Grid too small!')
                    
            self.beg_c[c] = n_p[self.parpos_c[c]]
            self.end_c[c] = n_p[self.parpos_c[c] + 1]
            self.n_cp.append(n_p)
            
        self.n_c = self.end_c - self.beg_c

        self.h_cv = self.cell_cv / self.N_c[:, np.newaxis]
        self.volume = abs(np.linalg.det(self.cell_cv))
        self.dv = self.volume / self.N_c.prod()

        self.orthogonal = not (self.cell_cv -
                               np.diag(self.cell_cv.diagonal())).any()

        # Sanity check for grid spacings:
        h_c = self.get_grid_spacings()
        if max(h_c) / min(h_c) > 1.3:
            raise ValueError('Very anisotropic grid spacings: %s' % h_c)

    def new_descriptor(self, N_c=None, cell_cv=None, pbc_c=None,
                       comm=None, parsize=None):
        """Create new descriptor based on this one.

        The new descriptor will use the same class (possibly a subclass)
        and all arguments will be equal to those of this descriptor
        unless new arguments are provided."""
        if N_c is None:
            N_c = self.N_c
        if cell_cv is None:
            cell_cv = self.cell_cv
        if pbc_c is None:
            pbc_c = self.pbc_c
        if comm is None:
            comm = self.comm
        if parsize is None:
            parsize = self.parsize_c
        return self.__class__(N_c, cell_cv, pbc_c, comm, parsize)

    def get_grid_spacings(self):
        L_c = (np.linalg.inv(self.cell_cv)**2).sum(0)**-0.5
        return L_c / self.N_c

    def get_size_of_global_array(self, pad=False):
        if pad:
            return self.N_c
        else:
            return self.N_c - 1 + self.pbc_c

    def flat_index(self, G_c):
        g1, g2, g3 = G_c - self.beg_c
        return g3 + self.n_c[2] * (g2 + g1 * self.n_c[1])
    
    def get_slice(self):
        return [slice(b - 1 + p, e - 1 + p) for b, e, p in
                zip(self.beg_c, self.end_c, self.pbc_c)]

    def zeros(self, n=(), dtype=float, global_array=False, pad=False, usemic=False):
        """Return new zeroed 3D array for this domain.

        The type can be set with the ``dtype`` keyword (default:
        ``float``).  Extra dimensions can be added with ``n=dim``.  A
        global array spanning all domains can be allocated with
        ``global_array=True``."""

        array = self._new_array(n, dtype, True, global_array, pad)
        if usemic:
            oa = stream.bind(array)
            stream.sync()
            return oa
        else:
            return array
    
    def empty(self, n=(), dtype=float, global_array=False, pad=False, usemic=False):
        """Return new uninitialized 3D array for this domain.

        The type can be set with the ``dtype`` keyword (default:
        ``float``).  Extra dimensions can be added with ``n=dim``.  A
        global array spanning all domains can be allocated with
        ``global_array=True``."""

        array = self._new_array(n, dtype, False, global_array, pad)
        if usemic:
            oa = stream.bind(array)
            stream.sync()
            return oa
        else:
            return array
        
    def _new_array(self, n=(), dtype=float, zero=True,
                   global_array=False, pad=False):
        if global_array:
            shape = self.get_size_of_global_array(pad)
        else:
            shape = self.n_c
            
        if isinstance(n, int):
            n = (n,)

        shape = n + tuple(shape)

        if zero:
            return np.zeros(shape, dtype)
        else:
            return np.empty(shape, dtype)
        
    def integrate(self, a_xg, b_yg=None,
                  global_integral=True, hermitian=False,
                  _transposed_result=None):
        """Integrate function(s) over domain.

        a_xg: ndarray
            Function(s) to be integrated.
        b_yg: ndarray
            If present, integrate a_xg.conj() * b_yg.
        global_integral: bool
            If the array(s) are distributed over several domains, then the
            total sum will be returned.  To get the local contribution
            only, use global_integral=False.
        hermitian: bool
            Result is hermitian.
        _transposed_result: ndarray
            Long story.  Don't use this unless you are a method of the
            MatrixOperator class ..."""
        
        xshape = a_xg.shape[:-3]
        
        if b_yg is None:
            # Only one array:
            result = a_xg.reshape(xshape + (-1,)).sum(axis=-1) * self.dv
            if global_integral:
                if result.ndim == 0:
                    result = self.comm.sum(result)
                else:
                    self.comm.sum(result)
            return result

        if isinstance(a_xg, mic.OffloadArray):
            # offload arrays have to be contiguous in any case
            A_xg = a_xg
            B_yg = b_yg
        else:
            A_xg = np.ascontiguousarray(a_xg.reshape((-1,) + a_xg.shape[-3:]))
            B_yg = np.ascontiguousarray(b_yg.reshape((-1,) + b_yg.shape[-3:]))

        if _transposed_result is None:
            result_yx = np.zeros((len(B_yg), len(A_xg)), A_xg.dtype)
        else:
            result_yx = _transposed_result
            global_integral = False

        if isinstance(a_xg, mic.OffloadArray):
            result_yx_mic = stream.bind(result_yx)
            stream.sync()
            # result_yx_mic.fillfrom(result_yx)
            # result_yx_mic.array[:] = result_yx[:]
            # result_yx_mic.update_device()

        if a_xg is b_yg:
            if isinstance(a_xg, mic.OffloadArray):
                # dsyrk performs badly in MIC so use dgemm here
                # mic_rk(self.dv, A_xg, 0.0, result_yx_mic)
                mic_gemm(self.dv, A_xg, A_xg, 0.0, result_yx_mic, 'c')
            else:
                rk(self.dv, A_xg, 0.0, result_yx)
        elif hermitian:
            if isinstance(a_xg, mic.OffloadArray):
                mic_r2k(self.dv, A_xg, B_yg, 0.0, result_yx_mic)
            else:
                r2k(0.5 * self.dv, A_xg, B_yg, 0.0, result_yx)
        else:
            if isinstance(a_xg, mic.OffloadArray):
                mic_gemm(self.dv, A_xg, B_yg, 0.0, result_yx_mic, 'c')
            else:
                gemm(self.dv, A_xg, B_yg, 0.0, result_yx, 'c')
        
        if isinstance(a_xg, mic.OffloadArray):
            result_yx_mic.update_host()
            stream.sync()

        if global_integral:
            self.comm.sum(result_yx)

        yshape = b_yg.shape[:-3]
        result = result_yx.T.reshape(xshape + yshape)
        
        if result.ndim == 0:
            return result.item()
        else:
            return result

    def gemm(self, alpha, psit_nG, C_mn, beta, newpsit_mG):
        """Helper function for MatrixOperator class."""
        gemm(alpha, psit_nG, C_mn, beta, newpsit_mG)

    def gemv(self, alpha, psit_nG, C_n, beta, newpsit_G, trans='t'):
        """Helper function for CG eigensolver."""
        gemv(alpha, psit_nG, C_n, beta, newpsit_G, trans)

    def coarsen(self):
        """Return coarsened `GridDescriptor` object.

        Reurned descriptor has 2x2x2 fewer grid points."""
        
        if np.sometrue(self.N_c % 2):
            raise ValueError('Grid %s not divisible by 2!' % self.N_c)

        return self.new_descriptor(self.N_c // 2)

    def refine(self):
        """Return refined `GridDescriptor` object.

        Returned descriptor has 2x2x2 more grid points."""
        return self.new_descriptor(self.N_c * 2)
    
    def get_boxes(self, spos_c, rcut, cut=True):
        """Find boxes enclosing sphere."""
        N_c = self.N_c
        ncut = rcut * (self.icell_cv**2).sum(axis=1)**0.5 * self.N_c
        npos_c = spos_c * N_c
        beg_c = np.ceil(npos_c - ncut).astype(int)
        end_c = np.ceil(npos_c + ncut).astype(int)

        if cut:
            for c in range(3):
                if not self.pbc_c[c]:
                    if beg_c[c] < 0:
                        beg_c[c] = 0
                    if end_c[c] > N_c[c]:
                        end_c[c] = N_c[c]
        else:
            for c in range(3):
                if (not self.pbc_c[c] and
                    (beg_c[c] < 0 or end_c[c] > N_c[c])):
                    msg = ('Box at %.3f %.3f %.3f crosses boundary.  '
                           'Beg. of box %s, end of box %s, max box size %s' %
                           (tuple(spos_c) + (beg_c, end_c, self.N_c)))
                    raise GridBoundsError(msg)
                    
        range_c = ([], [], [])
        
        for c in range(3):
            b = beg_c[c]
            e = b
            
            while e < end_c[c]:
                b0 = b % N_c[c]
               
                e = min(end_c[c], b + N_c[c] - b0)

                if b0 < self.beg_c[c]:
                    b1 = b + self.beg_c[c] - b0
                else:
                    b1 = b
                    
                e0 = b0 - b + e
                              
                if e0 > self.end_c[c]:
                    e1 = e - (e0 - self.end_c[c])
                else:
                    e1 = e
                if e1 > b1:
                    range_c[c].append((b1, e1))
                b = e
        
        boxes = []

        for b0, e0 in range_c[0]:
            for b1, e1 in range_c[1]:
                for b2, e2 in range_c[2]:
                    b = np.array((b0, b1, b2))
                    e = np.array((e0, e1, e2))
                    beg_c = np.array((b0 % N_c[0], b1 % N_c[1], b2 % N_c[2]))
                    end_c = beg_c + e - b
                    disp = (b - beg_c) / N_c
                    beg_c = np.maximum(beg_c, self.beg_c)
                    end_c = np.minimum(end_c, self.end_c)
                    if (beg_c[0] < end_c[0] and
                        beg_c[1] < end_c[1] and
                        beg_c[2] < end_c[2]):
                        boxes.append((beg_c, end_c, disp))

        return boxes

    def get_nearest_grid_point(self, spos_c, force_to_this_domain=False):
        """Return index of nearest grid point.
        
        The nearest grid point can be on a different CPU than the one the
        nucleus belongs to (i.e. return can be negative, or larger than
        gd.end_c), in which case something clever should be done.
        The point can be forced to the grid descriptors domain to be
        consistent with self.get_rank_from_position(spos_c).
        """
        g_c = np.around(self.N_c * spos_c).astype(int)
        if force_to_this_domain:
            for c in range(3):
                g_c[c] = max(g_c[c], self.beg_c[c])
                g_c[c] = min(g_c[c], self.end_c[c] - 1)
        return g_c - self.beg_c

    def plane_wave(self, k_c):
        """Evaluate plane wave on grid.

        Returns::

               _ _
              ik.r
             e    ,

        where the wave vector is given by k_c (in units of reciprocal
        lattice vectors)."""

        N_c = self.N_c
        return np.exp(2j * pi * np.dot(np.indices(N_c).T, k_c / N_c).T)

    def symmetrize(self, a_g, op_scc, ft_sc=None):
        if len(op_scc) == 1:
            return
        
        if ft_sc is not None and not ft_sc.any():
            ft_sc = None
            
        A_g = self.collect(a_g)
        if self.comm.rank == 0:
            B_g = np.zeros_like(A_g)
            for s, op_cc in enumerate(op_scc):
                if ft_sc is None:
                    _gpaw.symmetrize(A_g, B_g, op_cc)
                else:
                    _gpaw.symmetrize_ft(A_g, B_g, op_cc, ft_sc[s])
        else:
            B_g = None
        self.distribute(B_g, a_g)
        a_g /= len(op_scc)
    
    def collect(self, a_xg, broadcast=False):
        """Collect distributed array to master-CPU or all CPU's."""
        if self.comm.size == 1:
            return a_xg

        xshape = a_xg.shape[:-3]

        # Collect all arrays on the master:
        if self.rank != 0:
            # There can be several sends before the corresponding receives
            # are posted, so use syncronous send here
            self.comm.ssend(a_xg, 0, 301)
            if broadcast:
                A_xg = self.empty(xshape, a_xg.dtype, global_array=True)
                self.comm.broadcast(A_xg, 0)
                return A_xg
            else:
                return None

        # Put the subdomains from the slaves into the big array
        # for the whole domain:
        A_xg = self.empty(xshape, a_xg.dtype, global_array=True)
        parsize_c = self.parsize_c
        r = 0
        for n0 in range(parsize_c[0]):
            b0, e0 = self.n_cp[0][n0:n0 + 2] - self.beg_c[0]
            for n1 in range(parsize_c[1]):
                b1, e1 = self.n_cp[1][n1:n1 + 2] - self.beg_c[1]
                for n2 in range(parsize_c[2]):
                    b2, e2 = self.n_cp[2][n2:n2 + 2] - self.beg_c[2]
                    if r != 0:
                        a_xg = np.empty(xshape +
                                        ((e0 - b0), (e1 - b1), (e2 - b2)),
                                        a_xg.dtype.char)
                        self.comm.receive(a_xg, r, 301)
                    A_xg[..., b0:e0, b1:e1, b2:e2] = a_xg
                    r += 1
        if broadcast:
            self.comm.broadcast(A_xg, 0)
        return A_xg

    def distribute(self, B_xg, b_xg):
        """Distribute full array B_xg to subdomains, result in b_xg.

        B_xg is not used by the slaves (i.e. it should be None on all slaves)
        b_xg must be allocated on all nodes and will be overwritten.
        """

        if self.comm.size == 1:
            b_xg[:] = B_xg
            return
        
        if self.rank != 0:
            self.comm.receive(b_xg, 0, 42)
            return
        else:
            parsize_c = self.parsize_c
            requests = []
            r = 0
            for n0 in range(parsize_c[0]):
                b0, e0 = self.n_cp[0][n0:n0 + 2] - self.beg_c[0]
                for n1 in range(parsize_c[1]):
                    b1, e1 = self.n_cp[1][n1:n1 + 2] - self.beg_c[1]
                    for n2 in range(parsize_c[2]):
                        b2, e2 = self.n_cp[2][n2:n2 + 2] - self.beg_c[2]
                        if r != 0:
                            a_xg = B_xg[..., b0:e0, b1:e1, b2:e2].copy()
                            request = self.comm.send(a_xg, r, 42, NONBLOCKING)
                            # Remember to store a reference to the
                            # send buffer (a_xg) so that is isn't
                            # deallocated:
                            requests.append((request, a_xg))
                        else:
                            b_xg[:] = B_xg[..., b0:e0, b1:e1, b2:e2]
                        r += 1
                        
            for request, a_xg in requests:
                self.comm.wait(request)
        
    def zero_pad(self, a_xg):
        """Pad array with zeros as first element along non-periodic directions.

        Should only be invoked on global arrays.
        """
        assert np.all(a_xg.shape[-3:] == (self.N_c + self.pbc_c - 1))
        if self.pbc_c.all():
            return a_xg

        npbx, npby, npbz = 1 - self.pbc_c
        b_xg = np.zeros(a_xg.shape[:-3] + tuple(self.N_c), dtype=a_xg.dtype)
        b_xg[..., npbx:, npby:, npbz:] = a_xg
        return b_xg

    def calculate_dipole_moment(self, rho_g):
        """Calculate dipole moment of density."""
        rho_01 = rho_g.sum(axis=2)
        rho_02 = rho_g.sum(axis=1)
        rho_cg = [rho_01.sum(axis=1), rho_01.sum(axis=0), rho_02.sum(axis=0)]
        rhog_c = [np.dot(np.arange(self.beg_c[c], self.end_c[c]), rho_cg[c])
                  for c in range(3)]
        d_c = -np.dot(rhog_c, self.h_cv) * self.dv
        self.comm.sum(d_c)
        return d_c

    def wannier_matrix(self, psit_nG, psit_nG1, G_c, nbands=None):
        """Wannier localization integrals

        The soft part of Z is given by (Eq. 27 ref1)::

            ~       ~     -i G.r   ~
            Z   = <psi | e      |psi >
             nm       n             m
                    
        psit_nG and psit_nG1 are the set of wave functions for the two
        different spin/kpoints in question.

        ref1: Thygesen et al, Phys. Rev. B 72, 125119 (2005)
        """

        if nbands is None:
            nbands = len(psit_nG)

        if nbands == 0:
            return np.zeros((0, 0), complex)

        e_G = np.exp(-2j * pi * np.dot(np.indices(self.n_c).T +
                                       self.beg_c, G_c / self.N_c).T)
        a_nG = (e_G * psit_nG[:nbands].conj()).reshape((nbands, -1))
        return np.inner(a_nG,
                        psit_nG1[:nbands].reshape((nbands, -1))) * self.dv

    def find_center(self, a_R):
        """Calculate center of positive function."""
        assert self.orthogonal
        r_vR = self.get_grid_point_coordinates()
        a_R = a_R.astype(complex)
        center = []
        for L, r_R in zip(self.cell_cv.diagonal(), r_vR):
            z = self.integrate(a_R, np.exp(2j * pi / L * r_R))
            center.append(np.angle(z) / (2 * pi) * L % L)
        return np.array(center)
        
    def bytecount(self, dtype=float):
        """Get the number of bytes used by a grid of specified dtype."""
        return long(np.prod(self.n_c)) * np.array(1, dtype).itemsize

    def get_grid_point_coordinates(self, dtype=float, global_array=False):
        """Construct cartesian coordinates of grid points in the domain."""
        r_vG = np.dot(np.indices(self.n_c, dtype).T + self.beg_c,
                      self.h_cv).T.copy()
        if global_array:
            return self.collect(r_vG, broadcast=True)  # XXX waste!
        else:
            return r_vG

    def get_grid_point_distance_vectors(self, r_v, mic=True, dtype=float):
        """Return distances to a given vector in the domain.

        mic: if true adopts the mininimum image convention
        procedure by W. Smith in 'The Minimum image convention in
        Non-Cubic MD cells' March 29, 1989
        """
        s_Gc = (np.indices(self.n_c, dtype).T + self.beg_c) / self.N_c
        cell_cv = self.N_c * self.h_cv
        s_Gc -= np.linalg.solve(cell_cv.T, r_v)
        
        if mic:
            # XXX do the correction twice works better
            s_Gc -= self.pbc_c * (2 * s_Gc).astype(int)
            s_Gc -= self.pbc_c * (2 * s_Gc).astype(int)
            # sanity check
            assert((s_Gc * self.pbc_c >= -0.5).all())
            assert((s_Gc * self.pbc_c <= 0.5).all())
                
        return np.dot(s_Gc, cell_cv).T.copy()

    def interpolate_grid_points(self, spos_nc, vt_g, target_n, use_mlsqr=True):
        """Return interpolated values.

        Calculate interpolated values from array vt_g based on the
        scaled coordinates on spos_c.

        Uses moving least squares algorithm by default, or otherwise
        trilinear interpolation.
        
        This doesn't work in parallel, since it would require
        communication between neighbouring grid.  """

        assert mpi.world.size == 1

        if use_mlsqr:
            mlsqr(3, 2.3, spos_nc, self.N_c, self.beg_c, vt_g, target_n)
        else:
            for n, spos_c in enumerate(spos_nc):
                g_c = self.N_c * spos_c - self.beg_c

                # The begin and end of the array slice
                bg_c = np.floor(g_c).astype(int)
                Bg_c = np.ceil(g_c).astype(int)

                # The coordinate within the box (bottom left = 0,
                # top right = h_c)
                dg_c = g_c - bg_c
                Bg_c %= self.N_c

                target_n[n] = (
                    vt_g[bg_c[0], bg_c[1], bg_c[2]] *
                    (1.0 - dg_c[0]) * (1.0 - dg_c[1]) * (1.0 - dg_c[2]) +
                    vt_g[Bg_c[0], bg_c[1], bg_c[2]] *
                    (0.0 + dg_c[0]) * (1.0 - dg_c[1]) * (1.0 - dg_c[2]) +
                    vt_g[bg_c[0], Bg_c[1], bg_c[2]] *
                    (1.0 - dg_c[0]) * (0.0 + dg_c[1]) * (1.0 - dg_c[2]) +
                    vt_g[Bg_c[0], Bg_c[1], bg_c[2]] *
                    (0.0 + dg_c[0]) * (0.0 + dg_c[1]) * (1.0 - dg_c[2]) +
                    vt_g[bg_c[0], bg_c[1], Bg_c[2]] *
                    (1.0 - dg_c[0]) * (1.0 - dg_c[1]) * (0.0 + dg_c[2]) +
                    vt_g[Bg_c[0], bg_c[1], Bg_c[2]] *
                    (0.0 + dg_c[0]) * (1.0 - dg_c[1]) * (0.0 + dg_c[2]) +
                    vt_g[bg_c[0], Bg_c[1], Bg_c[2]] *
                    (1.0 - dg_c[0]) * (0.0 + dg_c[1]) * (0.0 + dg_c[2]) +
                    vt_g[Bg_c[0], Bg_c[1], Bg_c[2]] *
                    (0.0 + dg_c[0]) * (0.0 + dg_c[1]) * (0.0 + dg_c[2]))

    def __eq__(self, other):
        return (self.dv == other.dv and
                (self.h_cv == other.h_cv).all() and
                (self.N_c == other.N_c).all() and
                (self.n_c == other.n_c).all() and
                (self.beg_c == other.beg_c).all() and
                (self.end_c == other.end_c).all())

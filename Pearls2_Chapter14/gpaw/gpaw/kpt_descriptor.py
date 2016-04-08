# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""K-point/spin combination-descriptors

This module contains classes for defining combinations of two indices:

* Index k for irreducible kpoints in the 1st Brillouin zone.
* Index s for spin up/down if spin-polarized (otherwise ignored).

"""

import numpy as np

from ase.dft.kpoints import monkhorst_pack, get_monkhorst_pack_size_and_offset

from gpaw.kpoint import KPoint
import gpaw.mpi as mpi
import _gpaw


def to1bz(bzk_kc, cell_cv):
    """Wrap k-points to 1. BZ.

    Return k-points wrapped to the 1. BZ.

    bzk_kc: (n,3) ndarray
        Array of k-points in units of the reciprocal lattice vectors.
    cell_cv: (3,3) ndarray
        Unit cell.
    """

    B_cv = 2.0 * np.pi * np.linalg.inv(cell_cv).T
    K_kv = np.dot(bzk_kc, B_cv)
    N_xc = np.indices((3, 3, 3)).reshape((3, 27)).T - 1
    G_xv = np.dot(N_xc, B_cv)

    bz1k_kc = bzk_kc.copy()

    # Find the closest reciprocal lattice vector:
    for k, K_v in enumerate(K_kv):
        # If a k-point has the same distance to several reciprocal
        # lattice vectors, we don't want to pick a random one on the
        # basis of numerical noise, so we round off the differences
        # between the shortest distances to 6 decimals and chose the
        # one with the lowest index.
        d = ((G_xv - K_v)**2).sum(1)
        x = (d - d.min()).round(6).argmin()
        bz1k_kc[k] -= N_xc[x]

    return bz1k_kc


class KPointDescriptor:
    """Descriptor-class for k-points."""

    def __init__(self, kpts, nspins=1, collinear=True):
        """Construct descriptor object for kpoint/spin combinations (ks-pair).

        Parameters
        ----------
        kpts: None, sequence of 3 ints, or (n,3)-shaped array
            Specification of the k-point grid. None=Gamma, list of
            ints=Monkhorst-Pack, ndarray=user specified.
        nspins: int
            Number of spins.

        Attributes
        ===================  =================================================
        ``N_c``               Number of k-points in the different directions.
        ``nspins``            Number of spins in total.
        ``mynspins``          Number of spins on this CPU.
        ``nibzkpts``          Number of irreducible kpoints in 1st BZ.
        ``nks``               Number of k-point/spin combinations in total.
        ``mynks``             Number of k-point/spin combinations on this CPU.
        ``gamma``             Boolean indicator for gamma point calculation.
        ``comm``              MPI-communicator for kpoint distribution.
        ``weight_k``          Weights of each k-point
        ``ibzk_kc``           Unknown
        ``ibzk_qc``           Unknown
        ``sym_k``             Unknown
        ``time_reversal_k``   Unknown
        ``bz2ibz_k``          Unknown
        ``ibz2bz_k``          Unknown
        ``bz2bz_ks``          Unknown
        ``symmetry``          Object representing symmetries
        ===================  =================================================
        """

        if kpts is None:
            self.bzk_kc = np.zeros((1, 3))
            self.N_c = np.array((1, 1, 1), dtype=int)
            self.offset_c = np.zeros(3)
        elif isinstance(kpts[0], int):
            self.bzk_kc = monkhorst_pack(kpts)
            self.N_c = np.array(kpts, dtype=int)
            self.offset_c = np.zeros(3)
        else:
            self.bzk_kc = np.array(kpts, float)
            try:
                self.N_c, self.offset_c = \
                    get_monkhorst_pack_size_and_offset(self.bzk_kc)
            except ValueError:
                self.N_c = None
                self.offset_c = None

        self.collinear = collinear
        self.nspins = nspins
        self.nbzkpts = len(self.bzk_kc)

        # Gamma-point calculation?
        self.gamma = self.nbzkpts == 1 and np.allclose(self.bzk_kc, 0)
            
        # Point group and time-reversal symmetry neglected:
        self.weight_k = np.ones(self.nbzkpts) / self.nbzkpts
        self.ibzk_kc = self.bzk_kc.copy()
        self.sym_k = np.zeros(self.nbzkpts, int)
        self.time_reversal_k = np.zeros(self.nbzkpts, bool)
        self.bz2ibz_k = np.arange(self.nbzkpts)
        self.ibz2bz_k = np.arange(self.nbzkpts)
        self.bz2bz_ks = np.arange(self.nbzkpts)[:, np.newaxis]
        self.nibzkpts = self.nbzkpts
        self.nks = self.nibzkpts * self.nspins

        self.set_communicator(mpi.serial_comm)

        if self.gamma:
            self.description = '1 k-point (Gamma)'
        else:
            self.description = '%d k-points' % self.nbzkpts
            if self.N_c is not None:
                self.description += (': %d x %d x %d Monkhorst-Pack grid' %
                                     tuple(self.N_c))
                if self.offset_c.any():
                    self.description += ' + ['
                    for x in self.offset_c:
                        if x != 0 and abs(round(1 / x) - 1 / x) < 1e-12:
                            self.description += '1/%d,' % round(1 / x)
                        else:
                            self.description += '%f,' % x
                    self.description = self.description[:-1] + ']'

    def __len__(self):
        """Return number of k-point/spin combinations of local CPU."""

        return self.mynks

    def set_symmetry(self, atoms, symmetry, comm=None):
        """Create symmetry object and construct irreducible Brillouin zone.

        atoms: Atoms object
            Defines atom positions and types and also unit cell and
            boundary conditions.
        symmetry: Symmetry object
            Symmetry object.
        """

        self.symmetry = symmetry

        for c, periodic in enumerate(atoms.pbc):
            if not periodic and not np.allclose(self.bzk_kc[:, c], 0.0):
                raise ValueError('K-points can only be used with PBCs!')

        # Find symmetry operations of atoms:
        symmetry.analyze(atoms.get_scaled_positions())

        if symmetry.time_reversal or symmetry.point_group:
            (self.ibzk_kc, self.weight_k,
             self.sym_k,
             self.time_reversal_k,
             self.bz2ibz_k,
             self.ibz2bz_k,
             self.bz2bz_ks) = symmetry.reduce(self.bzk_kc, comm)

        # Number of irreducible k-points and k-point/spin combinations.
        self.nibzkpts = len(self.ibzk_kc)
        if self.collinear:
            self.nks = self.nibzkpts * self.nspins
        else:
            self.nks = self.nibzkpts

    def set_communicator(self, comm):
        """Set k-point communicator."""

        # Ranks < self.rank0 have mynks0 k-point/spin combinations and
        # ranks >= self.rank0 have mynks0+1 k-point/spin combinations.
        mynks0, x = divmod(self.nks, comm.size)
        self.rank0 = comm.size - x
        self.comm = comm

        # My number and offset of k-point/spin combinations
        self.mynks = self.get_count()
        self.ks0 = self.get_offset()

        if self.nspins == 2 and comm.size == 1:  # NCXXXXXXXX
            # Avoid duplicating k-points in local list of k-points.
            self.ibzk_qc = self.ibzk_kc.copy()
            self.weight_q = self.weight_k
        else:
            self.ibzk_qc = np.vstack((self.ibzk_kc,
                                      self.ibzk_kc))[self.get_slice()]
            self.weight_q = np.hstack((self.weight_k,
                                       self.weight_k))[self.get_slice()]

    def copy(self, comm=mpi.serial_comm):
        """Create a copy with shared symmetry object."""
        kd = KPointDescriptor(self.bzk_kc, self.nspins)
        kd.weight_k = self.weight_k
        kd.ibzk_kc = self.ibzk_kc
        kd.sym_k = self.sym_k
        kd.time_reversal_k = self.time_reversal_k
        kd.bz2ibz_k = self.bz2ibz_k
        kd.ibz2bz_k = self.ibz2bz_k
        kd.bz2bz_ks = self.bz2bz_ks
        kd.symmetry = self.symmetry
        kd.nibzkpts = self.nibzkpts
        kd.nks = self.nks
        kd.set_communicator(comm)
        return kd

    def create_k_points(self, gd):
        """Return a list of KPoints."""

        sdisp_cd = gd.sdisp_cd

        kpt_u = []

        for ks in range(self.ks0, self.ks0 + self.mynks):
            s, k = divmod(ks, self.nibzkpts)
            q = (ks - self.ks0) % self.nibzkpts
            if self.collinear:
                weight = self.weight_k[k] * 2 / self.nspins
            else:
                weight = self.weight_k[k]
            if self.gamma:
                phase_cd = np.ones((3, 2), complex)
            else:
                phase_cd = np.exp(2j * np.pi *
                                  sdisp_cd * self.ibzk_kc[k, :, np.newaxis])
            kpt_u.append(KPoint(weight, s, k, q, phase_cd))

        return kpt_u

    def collect(self, a_ux, broadcast=True):
        """Collect distributed data to all."""

        if self.comm.rank == 0 or broadcast:
            xshape = a_ux.shape[1:]
            a_skx = np.empty((self.nspins, self.nibzkpts) + xshape, a_ux.dtype)
            a_Ux = a_skx.reshape((-1,) + xshape)
        else:
            a_skx = None

        if self.comm.rank > 0:
            self.comm.send(a_ux, 0)
        else:
            u1 = self.get_count(0)
            a_Ux[0:u1] = a_ux
            requests = []
            for rank in range(1, self.comm.size):
                u2 = u1 + self.get_count(rank)
                requests.append(self.comm.receive(a_Ux[u1:u2], rank,
                                                  block=False))
                u1 = u2
            assert u1 == len(a_Ux)
            self.comm.waitall(requests)

        if broadcast:
            self.comm.broadcast(a_Ux, 0)

        return a_skx

    def transform_wave_function(self, psit_G, k, index_G=None, phase_G=None):
        """Transform wave function from IBZ to BZ.

        k is the index of the desired k-point in the full BZ.
        """

        s = self.sym_k[k]
        time_reversal = self.time_reversal_k[k]
        op_cc = np.linalg.inv(self.symmetry.op_scc[s]).round().astype(int)

        # Identity
        if (np.abs(op_cc - np.eye(3, dtype=int)) < 1e-10).all():
            if time_reversal:
                return psit_G.conj()
            else:
                return psit_G
        # General point group symmetry
        else:
            ik = self.bz2ibz_k[k]
            kibz_c = self.ibzk_kc[ik]
            b_g = np.zeros_like(psit_G)
            kbz_c = np.dot(self.symmetry.op_scc[s], kibz_c)
            if index_G is not None:
                assert index_G.shape == psit_G.shape == phase_G.shape,\
                    'Shape mismatch %s vs %s vs %s' % (index_G.shape,
                                                       psit_G.shape,
                                                       phase_G.shape)
                _gpaw.symmetrize_with_index(psit_G, b_g, index_G, phase_G)
            else:
                _gpaw.symmetrize_wavefunction(psit_G, b_g, op_cc.copy(),
                                              np.ascontiguousarray(kibz_c),
                                              kbz_c)

            if time_reversal:
                return b_g.conj()
            else:
                return b_g

    def get_transform_wavefunction_index(self, nG, k):
        """Get the "wavefunction transform index".

        This is a permutation of the numbers 1, 2, .. N which
        associates k + q to some k, and where N is the total
        number of grid points as specified by nG which is a
        3D tuple.

        Returns index_G and phase_G which are one-dimensional
        arrays on the grid."""

        s = self.sym_k[k]
        op_cc = np.linalg.inv(self.symmetry.op_scc[s]).round().astype(int)

        # General point group symmetry
        if (np.abs(op_cc - np.eye(3, dtype=int)) < 1e-10).all():
            nG0 = np.prod(nG)
            index_G = np.arange(nG0).reshape(nG)
            phase_G = np.ones(nG)
        else:
            ik = self.bz2ibz_k[k]
            kibz_c = self.ibzk_kc[ik]
            index_G = np.zeros(nG, dtype=int)
            phase_G = np.zeros(nG, dtype=complex)

            kbz_c = np.dot(self.symmetry.op_scc[s], kibz_c)
            _gpaw.symmetrize_return_index(index_G, phase_G, op_cc.copy(),
                                          np.ascontiguousarray(kibz_c),
                                          kbz_c)
        return index_G, phase_G

    def find_k_plus_q(self, q_c, kpts_k=None):
        """Find the indices of k+q for all kpoints in the Brillouin zone.

        In case that k+q is outside the BZ, the k-point inside the BZ
        corresponding to k+q is given.

        Parameters
        ----------
        q_c: ndarray
            Coordinates for the q-vector in units of the reciprocal
            lattice vectors.
        kpts_k: list of ints
            Restrict search to specified k-points.

        """
        k_x = kpts_k
        if k_x is None:
            return self.find_k_plus_q(q_c, range(self.nbzkpts))

        i_x = []
        for k in k_x:
            kpt_c = self.bzk_kc[k] + q_c
            d_kc = kpt_c - self.bzk_kc
            d_k = abs(d_kc - d_kc.round()).sum(1)
            i = d_k.argmin()
            if d_k[i] > 1e-8:
                raise RuntimeError('Could not find k+q!')
            i_x.append(i)

        return i_x

    def get_bz_q_points(self, first=False):
        """Return the q=k1-k2. q-mesh is always Gamma-centered."""
        shift_c = 0.5 * ((self.N_c + 1) % 2) / self.N_c
        bzq_qc = monkhorst_pack(self.N_c) + shift_c
        if first:
            return to1bz(bzq_qc, self.symmetry.cell_cv)
        else:
            return bzq_qc

    def get_ibz_q_points(self, bzq_qc, op_scc):
        """Return ibz q points and the corresponding symmetry operations that
        work for k-mesh as well."""

        ibzq_qc_tmp = []
        ibzq_qc_tmp.append(bzq_qc[-1])
        weight_tmp = [0]

        for i, op_cc in enumerate(op_scc):
            if np.abs(op_cc - np.eye(3)).sum() < 1e-8:
                identity_iop = i
                break

        ibzq_q_tmp = {}
        iop_q = {}
        timerev_q = {}
        diff_qc = {}

        for i in range(len(bzq_qc) - 1, -1, -1):  # loop opposite to kpoint
            try:
                ibzk, iop, timerev, diff_c = self.find_ibzkpt(
                    op_scc, ibzq_qc_tmp, bzq_qc[i])
                find = False
                for ii, iop1 in enumerate(self.sym_k):
                    if iop1 == iop and self.time_reversal_k[ii] == timerev:
                        find = True
                        break
                if find is False:
                    raise ValueError('cant find k!')

                ibzq_q_tmp[i] = ibzk
                weight_tmp[ibzk] += 1.
                iop_q[i] = iop
                timerev_q[i] = timerev
                diff_qc[i] = diff_c
            except ValueError:
                ibzq_qc_tmp.append(bzq_qc[i])
                weight_tmp.append(1.)
                ibzq_q_tmp[i] = len(ibzq_qc_tmp) - 1
                iop_q[i] = identity_iop
                timerev_q[i] = False
                diff_qc[i] = np.zeros(3)

        # reverse the order.
        nq = len(ibzq_qc_tmp)
        ibzq_qc = np.zeros((nq, 3))
        ibzq_q = np.zeros(len(bzq_qc), dtype=int)
        for i in range(nq):
            ibzq_qc[i] = ibzq_qc_tmp[nq - i - 1]
        for i in range(len(bzq_qc)):
            ibzq_q[i] = nq - ibzq_q_tmp[i] - 1
        self.q_weights = np.array(weight_tmp[::-1]) / len(bzq_qc)
        return ibzq_qc, ibzq_q, iop_q, timerev_q, diff_qc

    def find_ibzkpt(self, symrel, ibzk_kc, bzk_c):
        """Find index in IBZ and related symmetry operations."""
        find = False
        ibzkpt = 0
        iop = 0
        timerev = False

        for sign in (1, -1):
            for ioptmp, op in enumerate(symrel):
                for i, ibzk in enumerate(ibzk_kc):
                    diff_c = bzk_c - sign * np.dot(op, ibzk)
                    if (np.abs(diff_c - diff_c.round()) < 1e-8).all():
                        ibzkpt = i
                        iop = ioptmp
                        find = True
                        if sign == -1:
                            timerev = True
                        break
                if find:
                    break
            if find:
                break

        if not find:
            raise ValueError('Cant find corresponding IBZ kpoint!')
        return ibzkpt, iop, timerev, diff_c.round()

    def where_is_q(self, q_c, bzq_qc):
        """Find the index of q points in BZ."""
        d_qc = q_c - bzq_qc
        d_q = abs(d_qc - d_qc.round()).sum(1)
        q = d_q.argmin()
        if d_q[q] > 1e-8:
            raise RuntimeError('Could not find q!')
        return q

    def get_count(self, rank=None):
        """Return the number of ks-pairs which belong to a given rank."""

        if rank is None:
            rank = self.comm.rank
        assert rank in xrange(self.comm.size)
        mynks0 = self.nks // self.comm.size
        mynks = mynks0
        if rank >= self.rank0:
            mynks += 1
        return mynks

    def get_offset(self, rank=None):
        """Return the offset of the first ks-pair on a given rank."""

        if rank is None:
            rank = self.comm.rank
        assert rank in xrange(self.comm.size)
        mynks0 = self.nks // self.comm.size
        ks0 = rank * mynks0
        if rank >= self.rank0:
            ks0 += rank - self.rank0
        return ks0

    def get_rank_and_index(self, s, k):
        """Find rank and local index of k-point/spin combination."""

        u = self.where_is(s, k)
        rank, myu = self.who_has(u)
        return rank, myu

    def get_slice(self, rank=None):
        """Return the slice of global ks-pairs which belong to a given rank."""

        if rank is None:
            rank = self.comm.rank
        assert rank in xrange(self.comm.size)
        mynks, ks0 = self.get_count(rank), self.get_offset(rank)
        uslice = slice(ks0, ks0 + mynks)
        return uslice

    def get_indices(self, rank=None):
        """Return the global ks-pair indices which belong to a given rank."""

        uslice = self.get_slice(rank)
        return np.arange(*uslice.indices(self.nks))

    def get_ranks(self):
        """Return array of ranks as a function of global ks-pair indices."""

        ranks = np.empty(self.nks, dtype=int)
        for rank in range(self.comm.size):
            uslice = self.get_slice(rank)
            ranks[uslice] = rank
        assert (ranks >= 0).all() and (ranks < self.comm.size).all()
        return ranks

    def who_has(self, u):
        """Convert global index to rank information and local index."""

        mynks0 = self.nks // self.comm.size
        if u < mynks0 * self.rank0:
            rank, myu = divmod(u, mynks0)
        else:
            rank, myu = divmod(u - mynks0 * self.rank0, mynks0 + 1)
            rank += self.rank0
        return rank, myu

    def global_index(self, myu, rank=None):
        """Convert rank information and local index to global index."""

        if rank is None:
            rank = self.comm.rank
        assert rank in xrange(self.comm.size)
        ks0 = self.get_offset(rank)
        u = ks0 + myu
        return u

    def what_is(self, u):
        """Convert global index to corresponding kpoint/spin combination."""

        s, k = divmod(u, self.nibzkpts)
        return s, k

    def where_is(self, s, k):
        """Convert kpoint/spin combination to the global index thereof."""

        u = k + self.nibzkpts * s
        return u

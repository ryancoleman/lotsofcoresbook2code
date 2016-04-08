# Copyright (C) 2003  CAMP
# Copyright (C) 2014 R. Warmbier Materials for Energy Research Group,
# Wits University
# Please see the accompanying LICENSE file for further information.
from __future__ import print_function, division

import functools

import numpy as np

import _gpaw
import gpaw.mpi as mpi
from gpaw.utilities import gcd


def frac(f, n=2 * 3 * 4 * 5):
    if not isinstance(f, (int, float)):
        return np.array([frac(a, n) for a in f]).T
    if f == 0:
        return 0, 1
    x = n * f
    if abs(x - round(x)) > 1e-6:
        raise ValueError
    x = int(round(x))
    d = gcd(x, n)
    return x // d, n // d


def sfrac(f):
    if f == 0:
        return '0'
    return '%d/%d' % frac(f)


class Symmetry:
    """Interface class for determination of symmetry, point and space groups.

    It also provides to apply symmetry operations to kpoint grids,
    wavefunctions and forces.
    """
    def __init__(self, id_a, cell_cv, pbc_c=np.ones(3, bool), tolerance=1e-7,
                 point_group=True, time_reversal=True, symmorphic=True):
        """Construct symmetry object.

        Parameters:

        id_a: list of int
            Numbered atomic types
        cell_cv: array(3,3), float
            Cartesian lattice vectors
        pbc_c: array(3), bool
            Periodic boundary conditions.
        tolerance: float
            Tolerance for symmetry determination.
        symmorphic: bool
            Switch for the use of non-symmorphic symmetries aka: symmetries
            with fractional translations.  Default is to use only symmorphic
            symmetries.
        point_group: bool
            Use point-group symmetries.
        time_reversal: bool
            Use time-reversal symmetry.
        tolerance: float
            Relative tolerance.

        Attributes:

        op_scc:
            Array of rotation matrices
        ft_sc:
            Array of fractional translation vectors
        a_sa:
            Array of atomic indices after symmetry operation
        has_inversion:
            (bool) Have inversion
        """

        self.id_a = id_a
        self.cell_cv = np.array(cell_cv, float)
        assert self.cell_cv.shape == (3, 3)
        self.pbc_c = np.array(pbc_c, bool)
        self.tol = tolerance
        self.symmorphic = symmorphic
        self.point_group = point_group
        self.time_reversal = time_reversal

        # Disable fractional translations for non-periodic boundary conditions:
        if not self.pbc_c.all():
            self.symmorphic = True

        self.op_scc = np.identity(3, int).reshape((1, 3, 3))
        self.ft_sc = np.zeros((1, 3))
        self.a_sa = np.arange(len(id_a)).reshape((1, -1))
        self.has_inversion = False
        self.gcd_c = np.ones(3, int)

    def analyze(self, spos_ac):
        """Determine list of symmetry operations.

        First determine all symmetry operations of the cell. Then call
        ``prune_symmetries`` to remove those symmetries that are not satisfied
        by the atoms.
        """
        if self.point_group:
            self.find_lattice_symmetry()
            self.prune_symmetries_atoms(spos_ac)

    def find_lattice_symmetry(self):
        """Determine list of symmetry operations."""

        # Symmetry operations as matrices in 123 basis
        self.op_scc = []

        # Metric tensor
        metric_cc = np.dot(self.cell_cv, self.cell_cv.T)

        # Generate all possible 3x3 symmetry matrices using base-3 integers
        power = (6561, 2187, 729, 243, 81, 27, 9, 3, 1)

        # operation is a 3x3 matrix, with possible elements -1, 0, 1, thus
        # there are 3**9 = 19683 possible matrices
        for base3id in xrange(19683):
            op_cc = np.empty((3, 3), dtype=int)
            m = base3id
            for ip, p in enumerate(power):
                d, m = divmod(m, p)
                op_cc[ip // 3, ip % 3] = 1 - d

            # The metric of the cell should be conserved after applying
            # the operation
            opmetric_cc = np.dot(np.dot(op_cc, metric_cc), op_cc.T)

            if np.abs(metric_cc - opmetric_cc).sum() > self.tol:
                continue

            # Operation must not swap axes that are not both periodic
            pbc_cc = np.logical_and.outer(self.pbc_c, self.pbc_c)
            if op_cc[~(pbc_cc | np.identity(3, bool))].any():
                continue

            # Operation must not invert axes that are not periodic
            pbc_cc = np.logical_and.outer(self.pbc_c, self.pbc_c)
            if not (op_cc[np.diag(~self.pbc_c)] == 1).all():
                continue

            # operation is a valid symmetry of the unit cell
            self.op_scc.append(op_cc)

        self.op_scc = np.array(self.op_scc)
        self.ft_sc = np.zeros((len(self.op_scc), 3))

    def prune_symmetries_atoms(self, spos_ac):
        """Remove symmetries that are not satisfied by the atoms."""

        if len(spos_ac) == 0:
            return

        # Build lists of atom numbers for each type of atom - one
        # list for each combination of atomic number, setup type,
        # magnetic moment and basis set:
        a_ij = {}
        for a, id in enumerate(self.id_a):
            if id in a_ij:
                a_ij[id].append(a)
            else:
                a_ij[id] = [a]

        a_j = a_ij[self.id_a[0]]  # just pick the first species

        # if supercell disable fractional translations:
        if not self.symmorphic:
            op_cc = np.identity(3, int)
            ftrans_sc = spos_ac[a_j[1:]] - spos_ac[a_j[0]]
            ftrans_sc -= np.rint(ftrans_sc)
            for ft_c in ftrans_sc:
                a_a = self.check_one_symmetry(spos_ac, op_cc, ft_c, a_ij)
                if a_a is not None:
                    self.symmorphic = True
                    break

        symmetries = []
        ftsymmetries = []

        # go through all possible symmetry operations
        for op_cc in self.op_scc:
            # first ignore fractional translations
            a_a = self.check_one_symmetry(spos_ac, op_cc, [0, 0, 0], a_ij)
            if a_a is not None:
                symmetries.append((op_cc, [0, 0, 0], a_a))
            elif not self.symmorphic:
                # check fractional translations
                sposrot_ac = np.dot(spos_ac, op_cc)
                ftrans_jc = sposrot_ac[a_j] - spos_ac[a_j[0]]
                ftrans_jc -= np.rint(ftrans_jc)
                for ft_c in ftrans_jc:
                    try:
                        nom_c, denom_c = frac(ft_c)
                    except ValueError:
                        continue
                    ft_c = nom_c / denom_c
                    a_a = self.check_one_symmetry(spos_ac, op_cc, ft_c, a_ij)
                    if a_a is not None:
                        ftsymmetries.append((op_cc, ft_c, a_a))
                        for c, d in enumerate(denom_c):
                            self.gcd_c[c] = gcd(self.gcd_c[c] * d, d)

        # Add symmetry operations with fractional translations at the end:
        symmetries.extend(ftsymmetries)
        self.op_scc = np.array([sym[0] for sym in symmetries])
        self.ft_sc = np.array([sym[1] for sym in symmetries])
        self.a_sa = np.array([sym[2] for sym in symmetries])

        inv_cc = -np.eye(3, dtype=int)
        self.has_inversion = (self.op_scc == inv_cc).all(2).all(1).any()

    def check_one_symmetry(self, spos_ac, op_cc, ft_c, a_ij):
        """Checks whether atoms satisfy one given symmetry operation."""

        a_a = np.zeros(len(spos_ac), int)
        for a_j in a_ij.values():
            spos_jc = spos_ac[a_j]
            for a in a_j:
                spos_c = np.dot(spos_ac[a], op_cc)
                sdiff_jc = spos_c - spos_jc - ft_c
                sdiff_jc -= sdiff_jc.round()
                indices = np.where(abs(sdiff_jc).max(1) < self.tol)[0]
                if len(indices) == 1:
                    j = indices[0]
                    a_a[a] = a_j[j]
                else:
                    assert len(indices) == 0
                    return

        return a_a

    def check(self, spos_ac):
        """Check if positions satisfy symmetry operations."""

        nsymold = len(self.op_scc)
        self.prune_symmetries_atoms(spos_ac)
        if len(self.op_scc) < nsymold:
            raise RuntimeError('Broken symmetry!')

    def reduce(self, bzk_kc, comm=None):
        """Reduce k-points to irreducible part of the BZ.

        Returns the irreducible k-points and the weights and other stuff.

        """
        nbzkpts = len(bzk_kc)
        U_scc = self.op_scc
        nsym = len(U_scc)

        time_reversal = self.time_reversal and not self.has_inversion
        bz2bz_ks = map_k_points(bzk_kc, U_scc, time_reversal, comm, self.tol)

        bz2bz_k = -np.ones(nbzkpts + 1, int)
        ibz2bz_k = []
        for k in range(nbzkpts - 1, -1, -1):
            # Reverse order looks more natural
            if bz2bz_k[k] == -1:
                bz2bz_k[bz2bz_ks[k]] = k
                ibz2bz_k.append(k)
        ibz2bz_k = np.array(ibz2bz_k[::-1])
        bz2bz_k = bz2bz_k[:-1].copy()

        bz2ibz_k = np.empty(nbzkpts, int)
        bz2ibz_k[ibz2bz_k] = np.arange(len(ibz2bz_k))
        bz2ibz_k = bz2ibz_k[bz2bz_k]

        weight_k = np.bincount(bz2ibz_k) * (1.0 / nbzkpts)

        # Symmetry operation mapping IBZ to BZ:
        sym_k = np.empty(nbzkpts, int)
        for k in range(nbzkpts):
            # We pick the first one found:
            try:
                sym_k[k] = np.where(bz2bz_ks[bz2bz_k[k]] == k)[0][0]
            except IndexError:
                print(nbzkpts)
                print(k)
                print(bz2bz_k)
                print(bz2bz_ks[bz2bz_k[k]])
                print(np.shape(np.where(bz2bz_ks[bz2bz_k[k]] == k)))
                print(bz2bz_k[k])
                print(bz2bz_ks[bz2bz_k[k]] == k)
                raise

        # Time-reversal symmetry used on top of the point group operation:
        if time_reversal:
            time_reversal_k = sym_k >= nsym
            sym_k %= nsym
        else:
            time_reversal_k = np.zeros(nbzkpts, bool)

        assert (ibz2bz_k[bz2ibz_k] == bz2bz_k).all()
        for k in range(nbzkpts):
            sign = 1 - 2 * time_reversal_k[k]
            dq_c = (np.dot(U_scc[sym_k[k]], bzk_kc[bz2bz_k[k]]) -
                    sign * bzk_kc[k])
            dq_c -= dq_c.round()
            assert abs(dq_c).max() < 1e-10

        return (bzk_kc[ibz2bz_k], weight_k,
                sym_k, time_reversal_k, bz2ibz_k, ibz2bz_k, bz2bz_ks)

    def check_grid(self, N_c):
        """Check that symmetries are comensurate with grid."""
        for U_cc, ft_c in zip(self.op_scc, self.ft_sc):
            assert not (U_cc * N_c - (U_cc.T * N_c).T).any()
            t_c = ft_c * N_c
            assert np.allclose(t_c, t_c.round())

    def symmetrize(self, a, gd):
        """Symmetrize array."""
        gd.symmetrize(a, self.op_scc, self.ft_sc)

    def symmetrize_positions(self, spos_ac):
        """Symmetrizes the atomic positions."""
        spos_tmp_ac = np.zeros_like(spos_ac)
        spos_new_ac = np.zeros_like(spos_ac)
        for i, op_cc in enumerate(self.op_scc):
            spos_tmp_ac[:] = 0.
            for a in range(len(spos_ac)):
                spos_c = np.dot(spos_ac[a], op_cc) - self.ft_sc[i]
                # Bring back the negative ones:
                spos_c = spos_c - np.floor(spos_c + 1e-5)
                spos_tmp_ac[self.a_sa[i][a]] += spos_c
            spos_new_ac += spos_tmp_ac

        spos_new_ac /= len(self.op_scc)
        return spos_new_ac

    def symmetrize_wavefunction(self, a_g, kibz_c, kbz_c, op_cc,
                                time_reversal):
        """Generate Bloch function from symmetry related function in the IBZ.

        a_g: ndarray
            Array with Bloch function from the irreducible BZ.
        kibz_c: ndarray
            Corresponing k-point coordinates.
        kbz_c: ndarray
            K-point coordinates of the symmetry related k-point.
        op_cc: ndarray
            Point group operation connecting the two k-points.
        time-reversal: bool
            Time-reversal symmetry required in addition to the point group
            symmetry to connect the two k-points.
        """

        # Identity
        if (np.abs(op_cc - np.eye(3, dtype=int)) < 1e-10).all():
            if time_reversal:
                return a_g.conj()
            else:
                return a_g
        # Inversion symmetry
        elif (np.abs(op_cc + np.eye(3, dtype=int)) < 1e-10).all():
            return a_g.conj()
        # General point group symmetry
        else:
            import _gpaw
            b_g = np.zeros_like(a_g)
            if time_reversal:
                # assert abs(np.dot(op_cc, kibz_c) - -kbz_c) < tol
                _gpaw.symmetrize_wavefunction(a_g, b_g, op_cc.T.copy(),
                                              kibz_c, -kbz_c)
                return b_g.conj()
            else:
                # assert abs(np.dot(op_cc, kibz_c) - kbz_c) < tol
                _gpaw.symmetrize_wavefunction(a_g, b_g, op_cc.T.copy(),
                                              kibz_c, kbz_c)
                return b_g

    def symmetrize_forces(self, F0_av):
        """Symmetrize forces."""
        F_ac = np.zeros_like(F0_av)
        for map_a, op_cc in zip(self.a_sa, self.op_scc):
            op_vv = np.dot(np.linalg.inv(self.cell_cv),
                           np.dot(op_cc, self.cell_cv))
            for a1, a2 in enumerate(map_a):
                F_ac[a2] += np.dot(F0_av[a1], op_vv)
        return F_ac / len(self.op_scc)

    def print_symmetries(self, fd):
        """Print symmetry information."""

        p = functools.partial(print, file=fd)

        n = len(self.op_scc)
        nft = self.ft_sc.any(1).sum()
        p('Symmetries present (total):', n)
        if not self.symmorphic:
            p('Symmetries with fractional translations:', nft)

        # X-Y grid of symmetry matrices:
        p()
        nx = 6 if self.symmorphic else 3
        ns = len(self.op_scc)
        y = 0
        for y in range((ns + nx - 1) // nx):
            for c in range(3):
                for x in range(nx):
                    s = x + y * nx
                    if s == ns:
                        break
                    op_c = self.op_scc[s, c]
                    ft = self.ft_sc[s, c]
                    p('  (%2d %2d %2d)' % tuple(op_c), end='')
                    if not self.symmorphic:
                        p(' + (%4s)' % sfrac(ft), end='')
                p()
            p()


def map_k_points(bzk_kc, U_scc, time_reversal, comm=None, tol=1e-11):
    """Find symmetry relations between k-points.

    This is a Python-wrapper for a C-function that does the hard work
    which is distributed over comm.

    The map bz2bz_ks is returned.  If there is a k2 for which::

      = _    _    _
      U q  = q  + N,
       s k1   k2

    where N is a vector of integers, then bz2bz_ks[k1, s] = k2, otherwise
    if there is a k2 for which::

      = _     _    _
      U q  = -q  + N,
       s k1    k2

    then bz2bz_ks[k1, s + nsym] = k2, where nsym = len(U_scc).  Otherwise
    bz2bz_ks[k1, s] = -1.
    """

    if comm is None or isinstance(comm, mpi.DryRunCommunicator):
        comm = mpi.serial_comm

    nbzkpts = len(bzk_kc)
    ka = nbzkpts * comm.rank // comm.size
    kb = nbzkpts * (comm.rank + 1) // comm.size
    assert comm.sum(kb - ka) == nbzkpts

    if time_reversal:
        U_scc = np.concatenate([U_scc, -U_scc])

    bz2bz_ks = np.zeros((nbzkpts, len(U_scc)), int)
    bz2bz_ks[ka:kb] = -1
    _gpaw.map_k_points(np.ascontiguousarray(bzk_kc),
                       np.ascontiguousarray(U_scc), tol, bz2bz_ks, ka, kb)
    comm.sum(bz2bz_ks)
    return bz2bz_ks


def atoms2symmetry(atoms, id_a=None):
    """Create symmetry object from atoms object."""
    if id_a is None:
        id_a = atoms.get_atomic_numbers()
    symmetry = Symmetry(id_a, atoms.cell, atoms.pbc,
                        symmorphic=False,
                        time_reversal=False)
    symmetry.analyze(atoms.get_scaled_positions())
    return symmetry

    
def analyze_atoms(filename):
    """Analyse symmetry.
    
    filename: str
        filename containing atomic positions and unit cell."""
        
    import sys
    from ase.io import read
    atoms = read(filename)
    symmetry = atoms2symmetry(atoms)
    symmetry.print_symmetries(sys.stdout)

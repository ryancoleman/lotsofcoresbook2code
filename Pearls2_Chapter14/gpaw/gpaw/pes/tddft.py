"""PES using approximate LrTDDFT scheme.

"""
import numpy as np

from ase.units import Hartree

from gpaw.pes import BasePES
from gpaw.pes.state import State
from gpaw.utilities import packed_index

from numpy import sqrt, pi


class TDDFTPES(BasePES):

    def __init__(self, mother, excited_daughter, daughter=None,
                 shift=True, tolerance={}):
        self.tolerance = {
            'occupation': 1.e-10,
            'magnetic': 1.e-10,
            'grid': 0,
        }
        for key in tolerance.keys():
            if not key in self.tolerance:
                raise RuntimeError("Tolerance key '%s' not known."
                                   % key)
            self.tolerance[key] = tolerance[key]

        if excited_daughter.calculator is not None:
            self.c_d = excited_daughter.calculator
        else:
            self.c_d = daughter

        self.c_m = mother
        self.gd = self.c_m.wfs.gd
        self.lr_d = excited_daughter

        self.c_m.converge_wave_functions()
        self.c_d.converge_wave_functions()
        self.lr_d.diagonalize()

        self.check_systems()
        self.lr_d.jend = self.lr_d.kss[-1].j

        # Make good way for initialising these

        kmax = 0
        lmax = 0
        for kss in self.lr_d.kss:
            kmax = max(kmax, kss.i)
            lmax = max(lmax, kss.j)
        self.kmax = kmax + 1
        self.lmax = lmax + 1

        self.f = None
        self.be = None
        self.shift = shift

        def gs_orbitals(calc):
            indicees = []
            nbands = calc.get_number_of_bands()
            spin = (calc.get_number_of_spins() == 2)
            f_tolerance = (1 + spin) * self.tolerance['occupation']
            for kpt in calc.wfs.kpt_u:
                for i in range(nbands):
                    if kpt.f_n[i] > f_tolerance:
                        indicees.append(i + kpt.s * nbands)
                        if not spin:
                            indicees.append(i + nbands)
            return indicees
        self.gs_m = gs_orbitals(self.c_m)
        self.imax = len(self.gs_m)
        self.gs_d = gs_orbitals(self.c_d)

        if (len(self.gs_m) != len(self.gs_d) + 1):
            raise RuntimeError(('Mother valence %d does not correspond ' +
                                'to daughter valence %d. ' +
                                'Modify tolerance["occupation"] ?') %
                               (len(self.gs_m), len(self.gs_d)))

    def _calculate(self):

        self.ks_overlaps()
        self.single_overlaps()
        self.full_overlap_matrix()

        self._create_f()

    def ks_overlaps(self):
        """Evaluate KS overlaps of mother and daughter."""
        bands_m = self.c_m.get_number_of_bands()
        spin_m = self.c_m.get_number_of_spins() == 2
        bands_d = self.c_d.get_number_of_bands()
        spin_d = self.c_d.get_number_of_spins() == 2

        self.overlap = np.zeros((2 * bands_m, 2 * bands_d))
        for i_m in range(bands_m):
            for s_m in range(2):
                k_m = spin_m * s_m
                wf_m = self.c_m.wfs.kpt_u[k_m].psit_nG[i_m]

                for j_d in range(bands_d):
                    k_d = spin_d * s_m

                    wf_d = self.c_d.wfs.kpt_u[k_d].psit_nG[j_d]
                    me = self.gd.integrate(wf_m * wf_d)

                    i = s_m * bands_m + i_m
                    j = s_m * bands_d + j_d
                    self.overlap[i, j] = me + self._nuc_corr(i_m, j_d,
                                                             k_m, k_d)

    def single_overlaps(self):
        self.singles = np.zeros((self.imax, len(self.lr_d)))
        nbands_d = self.c_d.get_number_of_bands()

        for i, i_m in enumerate(self.gs_m):
            for kl, kss in enumerate(self.lr_d.kss):
                if kss.fij > self.tolerance['occupation']:
                    spin = kss.pspin

                    keep_row = list(self.gs_m)
                    keep_row.remove(i_m)

                    k_d = kss.i + spin * nbands_d
                    l_d = kss.j + spin * nbands_d
                    keep_col = list(self.gs_d)
                    keep_col.remove(k_d)
                    keep_col.append(l_d)

                    d_ikl = np.zeros((len(keep_row), len(keep_col)))

                    for col in range(len(keep_col)):
                        for row in range(len(keep_row)):
                            d_ikl[row, col] = self.overlap[keep_row[row],
                                                           keep_col[col]]

                    self.singles[i, kl] = np.linalg.det(d_ikl)

    def gs_gs_overlaps(self):
        """Evaluate overlap matrix of mother and daughter ground states.

        """
        g0 = np.zeros((self.imax))
        for i, i_m in enumerate(self.gs_m):
            keep_row = list(self.gs_m)
            keep_row.remove(i_m)

            keep_col = list(self.gs_d)
            d_i00 = np.zeros((len(keep_row), len(keep_col)))

            for col in range(len(keep_col)):
                for row in range(len(keep_row)):
                    d_i00[row, col] = self.overlap[keep_row[row],
                                                   keep_col[col]]

            g0[i] = (-1) ** (self.imax + i) * np.linalg.det(d_i00)
        return g0

    def full_overlap_matrix(self):
        """Full overlap matrix of mother and daughter many particle states.

        """
        self.g_Ii = np.zeros((len(self.lr_d) + 1, self.imax))
        self.g_Ii[0, :] = self.gs_gs_overlaps()

        for I in range(len(self.lr_d)):
            for i in range(self.imax):
                gi = 0
                for kl in range(len(self.lr_d)):
                    gi += self.lr_d[I].f[kl] * self.singles[i, kl]
                self.g_Ii[1 + I, i] = (-1.) ** (self.imax + i) * gi

    def _create_f(self):
        self.f = (self.g_Ii * self.g_Ii).sum(axis=1)

        if self.shift is True:
            shift = (self.c_d.get_potential_energy() -
                     self.c_m.get_potential_energy())
        else:
            shift = float(self.shift)

        self.be = (np.array([0] + list(self.lr_d.get_energies() * Hartree)) +
                   shift)

    def _nuc_corr(self, i_m, j_d, k_m, k_d):
        ma = 0

        for a, P_ni_m in self.c_m.wfs.kpt_u[k_m].P_ani.items():
            P_ni_d = self.c_d.wfs.kpt_u[k_d].P_ani[a]
            Pi_i = P_ni_m[i_m]
            Pj_i = P_ni_d[j_d]
            Delta_pL = self.c_m.wfs.setups[a].Delta_pL

            for i in range(len(Pi_i)):
                for j in range(len(Pj_i)):
                    pij = Pi_i[i] * Pj_i[j]
                    ij = packed_index(i, j, len(Pi_i))
                    ma += Delta_pL[ij, 0] * pij

        self.gd.comm.sum(ma)
        return sqrt(4 * pi) * ma

    def check_systems(self):
        """Check that mother and daughter systems correspond to each other.

        """
        gtol = self.tolerance['grid']
        mtol = self.tolerance['magnetic']
        if (np.abs(self.c_m.wfs.gd.cell_cv -
                   self.c_d.wfs.gd.cell_cv) > gtol).any():
            raise RuntimeError('Not the same grid:' +
                               str(self.c_m.wfs.gd.cell_cv) + ' !=' +
                               str(self.c_d.wfs.gd.cell_cv))
        if (np.abs(self.c_m.wfs.gd.h_cv -
                   self.c_d.wfs.gd.h_cv) > gtol).any():
            raise RuntimeError('Not the same grid')
        if (self.c_m.atoms.positions != self.c_m.atoms.positions).any():
            raise RuntimeError('Not the same atomic positions')
        if np.abs(np.abs(self.c_m.get_magnetic_moment() -
                         self.c_d.get_magnetic_moment()) - 1) > mtol:
            raise RuntimeError(('Mother (%g) ' %
                                self.c_m.get_magnetic_moment()) +
                               ('and daughter spin (%g) ' %
                                self.c_d.get_magnetic_moment()) +
                               'are not compatible')

    def Dyson_orbital(self, I):
        """Return the Dyson orbital corresponding to excition I."""
        if not hasattr(self, 'g'):
            self._calculate()
        if not hasattr(self, 'morbitals'):
            nbands = self.c_m.get_number_of_bands()
            spin = self.c_m.get_number_of_spins() == 2
            morbitals_ig = []
            for i in self.gs_m:
                k = int(i >= nbands) * spin
                i -= nbands * int(i >= nbands)
                morbitals_ig.append(self.c_m.wfs.kpt_u[k].psit_nG[i])
            self.morbitals_ig = np.array(morbitals_ig)

        dyson = State(self.gd)
        gridnn = self.gd.zeros()
        for i, g in enumerate(self.g_Ii[I]):
            gridnn += g * self.morbitals_ig[i]
        dyson.set_grid(gridnn)
        dyson.set_energy(-self.be[I])
        return dyson

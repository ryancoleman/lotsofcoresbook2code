# Copyright (C) 2008  CAMd
# Please see the accompanying LICENSE file for further information.

"""This module is used in delta self-consistent field (dSCF) calculations

dSCF is a simple 'ad hoc' method to estimating excitation energies within
DFT. The only difference to ordinary DFT is that one or more electrons(s)
are forced to occupy one or more predefined orbitals. The only restriction
on these orbitals is that they must be linear combinations of available
Kohn-Sham orbitals.

"""

import sys
import copy
import numpy as np

from ase.parallel import paropen
from ase.units import Hartree
from ase.utils import devnull

from gpaw import mpi
from gpaw.occupations import OccupationNumbers, FermiDirac


def dscf_calculation(paw, orbitals, atoms):
    """Helper function to prepare a calculator for a dSCF calculation

    Parameters
    ==========
    orbitals: list of lists
        Orbitals which one wants to occupy. The format is
        orbitals = [[1.0,orb1,0],[1.0,orb2,1],...], where 1.0 is the no.
        of electrons, orb1 and orb2 are the orbitals (see MolecularOrbitals
        below for an example of an orbital class). 0 and 1 represents the
        spin (up and down). This number is ignored in a spin-paired
        calculation.

    Example
    =======

    >>> atoms.set_calculator(calc)
    >>> e_gs = atoms.get_potential_energy() #ground state energy
    >>> sigma_star=MolecularOrbitals(calc, molecule=[0,1],
    >>>                              w=[[1.,0.,0.,0.],[-1.,0.,0.,0.]])
    >>> dscf_calculation(calc, [[1.0,sigma_star,1]], atoms)
    >>> e_exc = atoms.get_potential_energy() #excitation energy

    """

    # If the calculator has not been initialized the occupation object
    # is None
    if paw.occupations is None:
        paw.initialize(atoms)
    occ = paw.occupations
    if occ.width == 0:
        occ.width = 1e-6
    if isinstance(occ, OccupationsDSCF):
        paw.occupations.orbitals = orbitals
    else:
        new_occ = OccupationsDSCF(occ.width * Hartree, orbitals)
        paw.occupations = new_occ
    # If the calculator has already converged (for the ground state),
    # reset self-consistency and let the density be updated right away
    if paw.scf.converged:
        paw.scf.niter_fixdensity = 0
        paw.scf.reset()


class OccupationsDSCF(FermiDirac):
    """Occupation class.

    Corresponds to the ordinary FermiDirac class in occupation.py. Only
    difference is that it forces some electrons in the supplied orbitals
    in stead of placing all the electrons by a Fermi-Dirac distribution.
    """

    def __init__(self, width, orbitals):
        FermiDirac.__init__(self, width)
        self.orbitals = orbitals
        self.norbitals = len(self.orbitals)

        self.cnoe = 0.0
        for orb in self.orbitals:
            self.cnoe += orb[0]

    def set_number_of_electrons(self, wfs):
        self.nvalence = wfs.nvalence - self.cnoe

    def calculate(self, wfs):
        FermiDirac.calculate(self, wfs)

        # Get the expansion coefficients c_un for each dscf-orbital
        # and incorporate their respective occupations into kpt.ne_o
        c_oun = []
        for orb in self.orbitals:
            ef = self.fermilevel
            if self.fixmagmom:
                fermilevels = [ef + 0.5 * self.split, ef - 0.5 * self.split]
            else:
                fermilevels = ef
            c_oun.append(orb[1].expand(fermilevels, wfs))

        for u, kpt in enumerate(wfs.kpt_u):
            kpt.ne_o = np.zeros(self.norbitals, dtype=float)
            kpt.c_on = np.zeros((self.norbitals, len(kpt.f_n)), dtype=complex)

            for o, orb in enumerate(self.orbitals):
                # TODO XXX false if orb[0]<0 since abs(c_n)**2>0
                #kpt.c_on[o,:] = abs(orb[0])**0.5 * c_oun[o][u]

                kpt.ne_o[o] = orb[0]
                kpt.c_on[o, :] = c_oun[o][u]

                if wfs.nspins == 2:
                    assert orb[2] in range(2), 'Invalid spin index'

                    if orb[2] == kpt.s:
                        kpt.ne_o[o] *= kpt.weight
                    else:
                        kpt.ne_o[o] = 0.0
                else:
                    kpt.ne_o[o] *= 0.5 * kpt.weight

        # Correct the magnetic moment
        for orb in self.orbitals:
            if orb[2] == 0:
                self.magmom += orb[0]
            elif orb[2] == 1:
                self.magmom -= orb[0]

    def calculate_band_energy(self, wfs):
        FermiDirac.calculate_band_energy(self, wfs)

        de_band = 0.0
        for kpt in wfs.kpt_u:
            if hasattr(kpt, 'c_on'):
                for ne, c_n in zip(kpt.ne_o, kpt.c_on):
                    de_band += ne * np.dot(np.abs(c_n)**2, kpt.eps_n)
        self.e_band += wfs.band_comm.sum(wfs.kd.comm.sum(de_band))


class MolecularOrbital:
    """Class defining the orbitals that should be filled in a dSCF calculation.

    An orbital is defined through a linear combination of the atomic
    partial waves. In each self-consistent cycle the method expand
    is called. This method take the Kohn-Sham orbitals fulfilling the
    criteria given by Estart, Eend and nos and return the best
    possible expansion of the orbital in this basis. The integral
    of the Kohn-Sham all-electron wavefunction ``|u,n>`` (u being local spin
    and kpoint index) and the partial wave ``|\phi_i^a>`` is approximated
    by::

      wfs.kpt_u[u].P_ani = <\tilde p_i^a|\tilde\psi_{un}>.

    Parameters
    ----------
    paw: gpaw calculator instance
        The calculator used in the dSCF calculation.
    Estart: float
        Kohn-Sham orbitals with an energy above Efermi+Estart are used
        in the linear expansion.
    Eend: float
        Kohn-Sham orbitals with an energy below Efermi+Eend are used
        in the linear expansion.
    nos: int
        The maximum Number Of States used in the linear expansion.
    weights: dictionary
        The keys represent atoms and have values which are lists of weights
        of the contributing partial waves. The default value thuis corresponds
        to an antibonding 2\sigma orbital of atoms 0 and 1.
        Format::

          {1. atom: [weight of 1. projector function of the 1. atom,
                     weight of 2. projector function of the 1. atom, ...],
           2. atom: [weight of 1. projector function of the 2. atom,
                     weight of 2. projector function of the 2. atom, ...],
           ...}
    """

    def __init__(self, paw, Estart=0.0, Eend=1.e6,
                 nos=None, weights={0: [1], 1: [-1]}):

        self.fixmom = paw.input_parameters.fixmom
        self.w = weights
        self.Estart = Estart
        self.Eend = Eend
        self.nos = nos

    def expand(self, epsF, wfs):

        if wfs.nspins == 1:
            epsF = [epsF]
        elif not self.fixmom:
            epsF = [epsF, epsF]

        if self.nos is None:
            self.nos = wfs.bd.nbands

        c_un = []
        for u, kpt in enumerate(wfs.kpt_u):
            Porb_n = np.zeros(wfs.bd.nbands, dtype=complex)
            for a, P_ni in kpt.P_ani.items():
                if a in self.w.keys():
                    for i in range(len(self.w[a])):
                        Porb_n += self.w[a][i] * P_ni[:, i]
            wfs.gd.comm.sum(Porb_n)

            # Starting from KS orbitals with largest overlap,
            # fill in the expansion coeffients between Estart and Eend
            c_n = np.zeros(wfs.bd.nbands, dtype=complex)
            nos = 0
            bandpriority = np.argsort(abs(Porb_n)**2)[::-1]

            for n in bandpriority:
                if (kpt.eps_n[n] > epsF[kpt.s] + self.Estart and
                    kpt.eps_n[n] < epsF[kpt.s] + self.Eend):
                    c_n[n] = Porb_n[n].conj()
                    nos += 1
                if nos == self.nos:
                    break

            # Normalize expansion coefficients
            c_n /= np.sqrt(sum(abs(c_n)**2))
            c_un.append(c_n)

        return c_un


class AEOrbital:
    """Class defining the orbitals that should be filled in a dSCF calculation.

    An orbital is defined through a linear combination of KS orbitals
    which is determined by this class as follows: For each kpoint and spin
    we calculate the quantity ``c_n = <n|a>`` where ``|n>`` is the
    all-electron KS states in the calculation and ``|a>`` is the
    all-electron resonant state to be kept occupied. We can then
    write ``|a> = Sum(c_n|n>)`` and in each self-consistent cycle the
    method expand is called. This method take the Kohn-Sham
    orbitals fulfilling the criteria given by Estart, Eend and
    nos (Number Of States) and return the best possible expansion of
    the orbital in this basis.

    Parameters
    ----------
    paw: gpaw calculator instance
        The calculator used in the dSCF calculation.
    molecule: list of integers
        The atoms, which are a part of the molecule.
    Estart: float
        Kohn-Sham orbitals with an energy above Efermi+Estart are used
        in the linear expansion.
    Eend: float
        Kohn-Sham orbitals with an energy below Efermi+Eend are used
        in the linear expansion.
    nos: int
        The maximum Number Of States used in the linear expansion.
    wf_u: list of wavefunction arrays
        Wavefunction to be occupied on the kpts on this processor:

        wf_u = [kpt.psit_nG[n] for kpt in calc_mol.wfs.kpt_u]

    p_uai: list of dictionaries
        Projector overlaps with the wavefunction to be occupied for each
        kpoint. These are used when correcting to all-electron wavefunction
        overlaps:

        p_uai = [dict([(mol[a], P_ni[n]) for a, P_ni in kpt.P_ani.items()])
                 for kpt in paw.wfs.kpt_u]

        where mol is a list of atoms contributing to the state and n is the
        band. See examples in the dscf documentation on the gpaw webpage.
    """

    def __init__(self, paw, wf_u, p_uai, Estart=0.0, Eend=1.e6, nos=None,
                 txt='-'):

        self.fixmom = paw.input_parameters.fixmom
        self.wf_u = wf_u
        self.p_uai = p_uai
        self.Estart = Estart
        self.Eend = Eend
        self.nos = nos

        if txt is None:
            self.txt = devnull
        elif txt == '-':
            self.txt = sys.stdout
        elif isinstance(txt, str):
            self.txt = paropen(txt, 'w')
        else:
            self.txt = txt

    def expand(self, epsF, wfs):

        if wfs.nspins == 1:
            epsF = [epsF]
        elif not self.fixmom:
            epsF = [epsF, epsF]

        if self.nos is None:
            self.nos = wfs.bd.nbands

        # Check dimension of lists
        if len(self.wf_u) == len(wfs.kpt_u):
            wf_u = self.wf_u
            p_uai = self.p_uai
        else:
            raise RuntimeError('List of wavefunctions has wrong size')

        c_un = []
        p_u = []
        for u, kpt in enumerate(wfs.kpt_u):

            # Inner product of pseudowavefunctions
            wf = np.reshape(wf_u[u], -1)
            Wf_n = kpt.psit_nG
            Wf_n = np.reshape(Wf_n, (len(kpt.f_n), -1))
            Porb_n = np.dot(Wf_n.conj(), wf) * wfs.gd.dv

            # Correction to obtain inner product of AE wavefunctions
            for a, p_i in p_uai[u].items():
                for n in range(wfs.bd.nbands):
                    for i in range(len(p_i)):
                        for j in range(len(p_i)):
                            Porb_n[n] += (kpt.P_ani[a][n][i].conj() *
                                          wfs.setups[a].dO_ii[i][j] *
                                          p_i[j])
            wfs.gd.comm.sum(Porb_n)

            #print 'Kpt:', kpt.k, ' Spin:', kpt.s, \
            #      ' Sum_n|<orb|nks>|^2:', sum(abs(Porb_n)**2)
            p_u.append(np.array([sum(abs(Porb_n)**2)], dtype=float))

            # Starting from KS orbitals with largest overlap,
            # fill in the expansion coeffients
            c_n = np.zeros(wfs.bd.nbands, dtype=complex)
            nos = 0
            bandpriority = np.argsort(abs(Porb_n)**2)[::-1]

            for n in bandpriority:
                if (kpt.eps_n[n] > epsF[kpt.s] + self.Estart and
                    kpt.eps_n[n] < epsF[kpt.s] + self.Eend):
                    c_n[n] = Porb_n[n]
                    nos += 1
                if nos == self.nos:
                    break

            # Normalize expansion coefficients
            c_n /= np.sqrt(sum(abs(c_n)**2))
            c_un.append(c_n)

        for s in range(wfs.nspins):
            for k in range(wfs.kd.nibzkpts):
                p = wfs.collect_auxiliary(p_u, k, s)
                if wfs.world.rank == 0:
                    self.txt.write('Kpt: %d, Spin: %d, '
                                   'Sum_n|<orb|nks>|^2: %f\n' % (k, s, p))

        return c_un

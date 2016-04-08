# -*- coding: utf-8 -*-

"""This module defines an ELF class."""

from numpy import pi

from gpaw.fd_operators import Gradient
from gpaw.lfc import LocalizedFunctionsCollection as LFC


def _elf(nt_sg, nt_grad2_sg, taut_sg, ncut, spinpol):
    """Pseudo electron localisation function (ELF) as defined in
    Becke and Edgecombe, J. Chem. Phys., vol 92 (1990) 5397

    More comprehensive definition in
    M. Kohout and A. Savin, Int. J. Quantum Chem., vol 60 (1996) 875-882

    Arguments:
     =============== =====================================================
     ``nt_sg``       Pseudo valence density.
     ``nt_grad2_sg`` Squared norm of the density gradient.
     ``tau_sg``      Kinetic energy density.
     ``ncut``        Minimum density cutoff parameter.
     ``spinpol``     Boolean indicator for spin polarization.
     =============== =====================================================
    """

    # Fermi constant
    cF = 3.0 / 10 * (3 * pi**2)**(2.0 / 3.0)

    if spinpol:
        # Kouhut eq. (9)
        D0 = 2**(2.0/3.0) * cF * (nt_sg[0]**(5.0/3.0) + nt_sg[1]**(5.0/3.0))
        taut = taut_sg.sum(axis=0)
        D = taut - (nt_grad2_sg[0] / nt_sg[0] + nt_grad2_sg[1] / nt_sg[1]) / 8
    else:
        # Kouhut eq. (7)
        D0 = cF * nt_sg[0]**(5.0/3.0)
        taut = taut_sg[0]
        D = taut - nt_grad2_sg[0] / nt_sg[0] / 8

    elf_g = 1.0 / (1.0 + (D / D0)**2)

    if ncut is not None:
        nt = nt_sg.sum(axis=0)
        elf_g[nt < ncut] = 0.0

    return elf_g


class ELF:
    """ELF object for calculating the electronic localization function.

    Arguments:
     =============== =====================================================
     ``paw``         Instance of ``GPAW`` class.
     ``ncut``        Density cutoff below which the ELF is zero.
     =============== =====================================================
    """

    def __init__(self, paw=None, ncut=1e-6):
        """Create the ELF object."""

        self.gd = paw.wfs.gd
        self.paw = paw
        self.finegd = paw.density.finegd
        self.nspins = paw.density.nspins
        self.density = paw.density

        self.ncut = ncut
        self.spinpol = (self.nspins == 2)

        self.initialize(paw)

    def initialize(self, paw):

        if not paw.initialized:
            raise RuntimeError('PAW instance is not initialized')
        paw.converge_wave_functions()

        self.tauct = LFC(self.gd,
                         [[setup.tauct] for setup in self.density.setups],
                         forces=True, cut=True)
        spos_ac = paw.atoms.get_scaled_positions() % 1.0
        self.tauct.set_positions(spos_ac)

        self.taut_sg = None
        self.nt_grad2_sG = self.gd.empty(self.nspins)
        self.nt_grad2_sg = None

    def interpolate(self):

        self.density.interpolate_pseudo_density()

        if self.taut_sg is None:
            self.taut_sg = self.finegd.empty(self.nspins)
            self.nt_grad2_sg = self.finegd.empty(self.nspins)

        ddr_v = [Gradient(self.finegd, v, n=3).apply for v in range(3)]
        self.nt_grad2_sg[:] = 0.0
        d_g = self.finegd.empty()

        # Transfer the densities from the coarse to the fine grid
        for s in range(self.nspins):
            self.density.interpolator.apply(self.taut_sG[s],
                                            self.taut_sg[s])
            #self.density.interpolator.apply(self.nt_grad2_sG[s],
            #                                self.nt_grad2_sg[s])
            for v in range(3):
                ddr_v[v](self.density.nt_sg[s], d_g)
                self.nt_grad2_sg[s] += d_g**2.0

    def update(self):
        self.taut_sG = self.paw.wfs.calculate_kinetic_energy_density()

        # Add the pseudo core kinetic array
        for taut_G in self.taut_sG:
            self.tauct.add(taut_G, 1.0 / self.paw.wfs.nspins)

        # For periodic boundary conditions
        if self.paw.wfs.kd.symmetry is not None:
            self.paw.wfs.kd.symmetry.symmetrize(self.taut_sG[0], self.paw.wfs.gd)

        self.nt_grad2_sG[:] = 0.0

        d_G = self.gd.empty()

        for s in range(self.nspins):
            for v in range(3):
                self.paw.wfs.taugrad_v[v](self.density.nt_sG[s], d_G)
                self.nt_grad2_sG[s] += d_G**2.0

        # TODO are nct from setups usable for nt_grad2_sG ?

    def get_electronic_localization_function(self, gridrefinement=1,
                                             pad=True, broadcast=True):

        # Returns dimensionless electronic localization function
        if gridrefinement == 1:
            elf_G = _elf(self.density.nt_sG, self.nt_grad2_sG,
                         self.taut_sG, self.ncut, self.spinpol)
            elf_G = self.gd.collect(elf_G, broadcast)
            if pad:
                elf_G = self.gd.zero_pad(elf_G)
            return elf_G
        elif gridrefinement == 2:
            if self.nt_grad2_sg is None:
                self.interpolate()

            elf_g = _elf(self.density.nt_sg, self.nt_grad2_sg,
                         self.taut_sg, self.ncut, self.spinpol)
            elf_g = self.finegd.collect(elf_g, broadcast)
            if pad:
                elf_g = self.finegd.zero_pad(elf_g)
            return elf_g
        else:
            raise NotImplementedError('Arbitrary refinement not implemented')

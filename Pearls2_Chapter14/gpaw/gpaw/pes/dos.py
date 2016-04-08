"""Photoelectron spectra from the (shifted) DOS approach."""

import numpy as np
from ase.units import Hartree

from gpaw.pes import BasePES


class DOSPES(BasePES):

    """PES derived from density of states with shifted KS-energies.

    """

    def __init__(self, mother, daughter=None, shift=False, vde=None, f_min=0.1):
        self.c_m = mother
        self.c_d = daughter
        self.f = None
        self.be = None
        self.shift = shift
        self.vde = vde
        self.f_min = f_min

    def _calculate(self):
        """Evaluate energies and spectroscopic factors."""

        self.be = []
        self.f = []
        ex_m = []
        for spin in range(self.c_m.get_number_of_spins()):
            # use only the Gamma point
            eps_n = self.c_m.get_eigenvalues(0, spin)
            f_n = self.c_m.get_occupation_numbers(0, spin)
            for e, f in zip(eps_n, f_n):
                if f > self.f_min:
                    self.be.append(-e)
                    self.f.append(f)
                    ex_m.append(-e)
        self.be = np.array(self.be)

        # find HOMO energy
        ex_m.sort()
        e_HOMO = ex_m[0]

        if self.vde is not None:
            assert(self.shift is False)
            energy_shift = -e_HOMO + self.vde
        else:
            if self.shift is True:
                e_m = self.c_m.get_potential_energy()
                try:
                    energy_shift = float(self.c_d) - e_HOMO
                except AttributeError:
                    e_d = self.c_d.get_potential_energy()
                    energy_shift = e_d - e_m - e_HOMO
            else:
                energy_shift = float(self.shift)

        self.be += energy_shift

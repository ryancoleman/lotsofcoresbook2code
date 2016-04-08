from ase.units import Bohr

import gpaw.mpi as mpi
from gpaw.lrtddft.spectrum import Writer

ds_prefactor = {
    'Ang': Bohr ** 2,
    'a.u.': 1.0,
    'Mb': Bohr ** 2 * 100}


class PESpectrum(Writer):

    def __init__(self,
                 enlist,
                 folding='Gauss',
                 width=0.08):  # Gauss/Lorentz width
        Writer.__init__(self, folding, width)
        self.title = 'Photo emission spectrum'
        self.fields = 'Binding energy [eV]   '
        if folding is None:
            self.fields += 'Spectroscopic factor'
        else:
            self.fields += 'Folded spectroscopic factor'

        self.energies = enlist[0]
        self.values = []
        for val in enlist[1]:
            self.values.append([val])


class BasePES:

    def save_folded_pes(self,
                        filename=None,
                        width=0.08,  # Gauss/Lorentz width
                        emin=None,
                        emax=None,
                        de=None,
                        folding='Gauss',
                        comment=None):

        ew = self.get_energies_and_weights()
        if mpi.rank == mpi.MASTER:
            sp = PESpectrum(ew, folding, width)
            sp.write(filename, emin, emax, de, comment)

    def get_energies_and_weights(self):
        if self.be is None or self.f is None:
            self._calculate()

        return self.be, self.f

    def set_first_peak_energy(self, energy):
        self.first_peak_energy = energy

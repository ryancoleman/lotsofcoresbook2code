from ase import Atoms
from gpaw import GPAW
from gpaw.test import equal

for mode in ['fd', 'pw']:
    print(mode)
    hydrogen = Atoms('H',
                     cell=(2.5, 3, 3.5),
                     pbc=1,
                     calculator=GPAW(txt=None, mode=mode))
    hydrogen.get_potential_energy()
    dens = hydrogen.calc.density
    ham = hydrogen.calc.hamiltonian
    ham.poisson.eps = 1e-20
    dens.interpolate_pseudo_density()
    dens.calculate_pseudo_charge()
    ham.update(dens)
    ham.get_energy(hydrogen.calc.occupations)
    y = (ham.vt_sG[0, 0, 0, 0] - ham.vt_sG[0, 0, 0, 1]) * ham.gd.dv
    x = 0.0001
    dens.nt_sG[0, 0, 0, 0] += x
    dens.nt_sG[0, 0, 0, 1] -= x
    dens.interpolate_pseudo_density()
    dens.calculate_pseudo_charge()
    ham.update(dens)
    e1 = ham.get_energy(hydrogen.calc.occupations) - ham.Ekin
    dens.nt_sG[0, 0, 0, 0] -= 2 * x
    dens.nt_sG[0, 0, 0, 1] += 2 * x
    dens.interpolate_pseudo_density()
    dens.calculate_pseudo_charge()
    ham.update(dens)
    e2 = ham.get_energy(hydrogen.calc.occupations) - ham.Ekin
    equal(y, (e1 - e2) / (2 * x), 2e-8)

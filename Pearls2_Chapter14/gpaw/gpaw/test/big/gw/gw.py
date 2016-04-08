import ase.db
from ase.lattice import bulk

from gpaw import GPAW, PW, FermiDirac
from gpaw.response.g0w0 import G0W0


data = {
    'C': ['diamond', 3.553],
    'Si': ['diamond', 5.421],
    'Ge': ['diamond', 5.644],
    'SiC': ['zincblende', 4.346],
    'AlN': ['zincblende', 4.368],
    'AlP': ['zincblende', 5.451],
    'AlAs': ['zincblende', 5.649],
    'GaN': ['zincblende', 4.520],
    'GaP': ['zincblende', 5.439],
    'GaAs': ['zincblende', 5.640],
    'InP': ['zincblende', 5.858],
    'InAs': ['zincblende', 6.047],
    'InSb': ['zincblende', 6.468],
    'BN': ['zincblende', 3.615],
    'GaSb': ['zincblende', 6.136],
    'MgO': ['rocksalt', 4.213],
    'ZnO': ['zincblende', 4.580],
    'ZnS': ['zincblende', 5.420],
    'ZnSe': ['zincblende', 5.674],
    'ZnTe': ['zincblende', 6.079],
    'CdO': ['rocksalt', 4.695],
    'CdS': ['zincblende', 5.832],
    'CdSe': ['zincblende', 6.077],
    'CdTe': ['zincblende', 6.477]}


c = ase.db.connect('gw.db')

for name in data:
    id = c.reserve(name=name)
    if id is None:
        continue
        
    x, a = data[name]
    atoms = bulk(name, x, a=a)
    atoms.calc = GPAW(mode=PW(600),
                      xc='LDA',
                      parallel={'band': 1},
                      occupations=FermiDirac(0.02),
                      kpts={'size': (6, 6, 6), 'gamma': True},
                      txt='%s.txt' % name)
    atoms.get_potential_energy()
    atoms.calc.diagonalize_full_hamiltonian(nbands=400)
    atoms.calc.write(name, mode='all')
    n = int(atoms.calc.get_number_of_electrons()) // 2
    gw = G0W0(name, 'gw-' + name,
              nbands=400,
              kpts=[(0, 0, 0), (0.5, 0.5, 0.5), (0.5, 0.5, 0)],
              ecut=200,
              hilbert=True,
              fast=True,
              domega0=0.1,
              eta=0.2,
              bands=(0, n + 2))
    results = gw.calculate()
    c.write(atoms, name=name, data=results)
    del c[id]

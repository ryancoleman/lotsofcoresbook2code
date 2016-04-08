import sys

from ase import Atoms
from ase.data.s22 import data, s22
from ase.calculators.vdwcorrection import vdWTkatchenko09prl
from ase.parallel import parprint

from gpaw import GPAW, FermiDirac
from gpaw.cluster import Cluster
from gpaw.analyse.hirshfeld import HirshfeldDensity, HirshfeldPartitioning
from gpaw.analyse.vdwradii import vdWradii

h = 0.25
box = 3.0

molecule = 'Adenine-thymine_complex_stack'

Energy = {
    'PBE': [],
    'vdW-DF': [],
    'TS09': []}

for molecule in ['Adenine-thymine_complex_stack']:
    ss = Cluster(Atoms(data[molecule]['symbols'], 
                       data[molecule]['positions']))
    
    # split the structures
    s1 = ss.find_connected(0)
    s2 = ss.find_connected(-1)
    assert len(ss) == len(s1) + len(s2)
    
    c = GPAW(xc='PBE', h=h, nbands=-6,
             occupations=FermiDirac(width=0.1), txt=None)
    cdf = GPAW(xc='vdW-DF', h=h, nbands=-6, occupations=FermiDirac(width=0.1),
               txt=None)
    
    for s in [s1, s2, ss]:
        s.set_calculator(c)
        s.minimal_box(box, h=h)
        Energy['PBE'].append(s.get_potential_energy())
        cc = vdWTkatchenko09prl(HirshfeldPartitioning(c),
                                vdWradii(s.get_chemical_symbols(), 'PBE'))
        s.set_calculator(cc)
        Energy['TS09'].append(s.get_potential_energy())

        s.set_calculator(cdf)
        Energy['vdW-DF'].append(s.get_potential_energy())

    parprint('Coupled cluster binding energy',
             -data[molecule]['interaction energy CC'] * 1000, 'meV')
    for xc in ['PBE', 'vdW-DF', 'TS09']:
        ene = Energy[xc]
#        print xc, 'energies', ene
        parprint(xc, 'binding energy',
                 (ene[0] + ene[1] - ene[2]) * 1000, 'meV')

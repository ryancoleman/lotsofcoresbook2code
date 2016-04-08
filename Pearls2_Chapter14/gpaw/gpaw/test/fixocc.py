from __future__ import print_function
from ase.structure import molecule
from ase.parallel import parprint
from gpaw import GPAW
from gpaw.cluster import Cluster
from gpaw.occupations import FixedOccupations, ZeroKelvin
from gpaw.test import equal

h = 0.4
box = 2
nbands = 2
txt = '-'
txt = None

H2 = Cluster(molecule('H2'))
H2.minimal_box(box, h)
convergence = {'energy':0.01, 'eigenstates':1.e-3, 'density':1.e-2}

if 1:
    # test ZeroKelvin vs FixedOccupations
    c = GPAW(h=h, nbands=nbands, 
             occupations=ZeroKelvin(True), 
             convergence=convergence,
             txt=txt)
    H2.set_calculator(c)
    E_zk = H2.get_potential_energy()

    c = GPAW(h=h, nbands=nbands, 
             occupations=FixedOccupations([[2, 0]]), 
             convergence=convergence,
             txt=txt)
    H2.set_calculator(c)
    E_fo = H2.get_potential_energy()
    parprint(E_zk, E_fo)
    equal(E_zk, E_fo, 1.e-10)

if 1:
    # test spin-paired vs spin-polarized
    c = GPAW(h=h, nbands=nbands, 
             occupations=FixedOccupations([[1,1]]), 
             convergence=convergence,
             txt=txt)
    H2.set_calculator(c)
    E_ns = H2.get_potential_energy()
if 1:
    c = GPAW(h=h, nbands=nbands, spinpol=True,
             occupations=FixedOccupations([[0.5, 0.5]] * 2), 
             convergence=convergence,
             txt=txt)
    H2.set_calculator(c)
    E_sp = H2.get_potential_energy()
    parprint(E_ns, E_sp)
    equal(E_ns, E_sp, 1.e-6)

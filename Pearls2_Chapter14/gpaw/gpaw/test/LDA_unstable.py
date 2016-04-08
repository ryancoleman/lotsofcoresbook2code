# http://listserv.fysik.dtu.dk/pipermail/gpaw-developers/2014-February/004374.html
from __future__ import print_function
from ase import Atom, Atoms
from gpaw import GPAW
from gpaw.test import equal


for i in range(40):
    a = 6.
    b = a / 2
    mol = Atoms([Atom('O',(b, b, 0.1219 + b)),
                 Atom('H',(b, 0.7633 + b, -0.4876 + b)),
                 Atom('H',(b, -0.7633 + b, -0.4876 + b))],
                pbc=False, cell=[a, a, a])
    calc = GPAW(gpts=(24, 24, 24), nbands=4, mode='lcao', txt=None,
                xc='LDA')
    def stop():
        calc.scf.converged = True
    calc.attach(stop, 1)
    mol.set_calculator(calc)
    e = mol.get_potential_energy()
    if i == 0:
        eref = e
    if calc.wfs.world.rank == 0:
        print(repr(e))
    equal(e - eref, 0, 1.e-12)

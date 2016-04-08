from __future__ import print_function
from gpaw import GPAW, restart
from ase.structure import molecule
from gpaw.test import equal
Eini0 = -17.7436014436
Iini0 = 12
esolvers = ['cg', 'rmm-diis', 'dav']
E0 = {'cg':       -17.7439210822,
      'rmm-diis': -17.7439094546,
      'dav':      -17.7439760796}
I0 = {'cg': 6, 'rmm-diis': 6, 'dav': 8}


calc = GPAW(xc='LDA',
            eigensolver='cg',
            convergence={'eigenstates': 3.5e-5},
            #txt=None,
            dtype=complex)



mol = molecule('N2')
mol.center(vacuum=3.0)
mol.set_calculator(calc)

Eini = mol.get_potential_energy()
Iini = calc.get_number_of_iterations()
print(('%10s: %12.6f eV in %3d iterations' %
       ('init(cg)', Eini, Iini)))
equal(Eini, Eini0, 1E-8)

calc.write('N2.gpw', mode='all')
del calc, mol

E = {}
I = {}
for esolver in esolvers:

    mol, calc = restart('N2.gpw', txt=None)

    if (calc.wfs.dtype!=complex or
        calc.wfs.kpt_u[0].psit_nG.dtype!=complex):
        raise AssertionError('ERROR: restart failed to read complex WFS')
    
    calc.scf.reset()
    calc.set(convergence={'eigenstates': 3.5e-9})
    calc.set(eigensolver=esolver)

    E[esolver]=mol.get_potential_energy()
    I[esolver]=calc.get_number_of_iterations()
    
    print(('%10s: %12.6f eV in %3d iterations' %
           (esolver, E[esolver], I[esolver])))

for esolver in esolvers:
    print(esolver)
    equal(E[esolver], E0[esolver], 8E-5)
    

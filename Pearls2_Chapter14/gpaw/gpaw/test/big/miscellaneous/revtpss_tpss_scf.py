from __future__ import print_function
from ase import Atom, Atoms
from gpaw import GPAW

L = 6
name='N2'
a = Atoms([Atom('N', (L/2.+1.098/2.,L/2.,L/2.)),
           Atom('N', (L/2.-1.098/2.,L/2.,L/2.))], cell=(L, L, L), pbc=False)

calc=GPAW(h=0.22, xc='PBE',convergence={'eigenstates':1.0e-7}, stencils=(3,3), txt = name+'.txt',
          eigensolver='rmm-diis')
a.set_calculator(calc)
e_n2 = a.get_potential_energy()
n2t = calc.get_xc_difference('TPSS')
n2rt = calc.get_xc_difference('revTPSS')
#print >> open('NSCF.txt','a'), name+ "TPSS Energy", e_n2 + n2t, name+ "revTPSS Energy", e_n2 + n2rt

a.calc.set(xc='TPSS')
e_n2t = a.get_potential_energy()

a.calc.set(xc='revTPSS')
e_n2rt = a.get_potential_energy()

name='N'
b = Atoms([Atom('N', (L/2.,L/2.,L/2.),magmom=3)], cell=(L, L, L), pbc=False)

calc=GPAW(h=0.22, xc='PBE',convergence={'eigenstates':1.0e-7}, stencils=(3,3), txt = name+'.txt',
          eigensolver='rmm-diis',
          fixmom=True, hund=True)
b.set_calculator(calc)
e_n = b.get_potential_energy()
nt = calc.get_xc_difference('TPSS')
nrt = calc.get_xc_difference('revTPSS')
#print >> open('NSCF.txt','a'), name+ "TPSS Energy", e_n + nt, name+ "revTPSS Energy", e_n + nrt

b.calc.set(xc='TPSS')
e_nt = b.get_potential_energy()

b.calc.set(xc='revTPSS')
e_nrt = b.get_potential_energy()

print('Atm. Experiment  ', -228.5)
print('Atm. PBE         ', (e_n2-2*e_n)*23.06)
print('Atm. TPSS(nsc)   ', ((e_n2+n2t)-2*(e_n+nt))*23.06)
print('Atm. TPSS        ', (e_n2t-2*e_nt)*23.06)
print('Atm. revTPSS(nsc)', ((e_n2+n2rt)-2*(e_n+nrt))*23.06)
print('Atm. revTPSS     ', (e_n2rt-2*e_nrt)*23.06)


from ase import Atom, Atoms
from gpaw import GPAW, mpi
from gpaw.test import equal
from gpaw.lrtddft import LrTDDFT

from gpaw.pes.dos import DOSPES
from gpaw.pes.tddft import TDDFTPES

txt = None
R=0.7 # approx. experimental bond length
a = 3.0
c = 3.0
h = .3
H2 = Atoms([Atom('H', (a / 2, a / 2, (c - R) / 2)),
                Atom('H', (a / 2, a / 2, (c + R) / 2))],
               cell=(a, a, c))

H2_plus = Atoms([Atom('H', (a / 2, a / 2, (c - R) / 2)),
                Atom('H', (a / 2, a / 2, (c + R) / 2))],
               cell=(a, a, c))

xc='LDA'

calc = GPAW(gpts=(12, 12, 12), xc=xc, nbands=1,
            parallel={'domain': mpi.world.size},
            spinpol=True, txt=txt)
H2.set_calculator(calc)
e_H2 = H2.get_potential_energy()
niter_H2 = calc.get_number_of_iterations()


calc_plus = GPAW(gpts=(12, 12, 12), xc=xc, nbands=2, 
                 parallel={'domain': mpi.world.size},
                 spinpol=True, txt=txt)
calc_plus.set(charge=+1)
H2_plus.set_calculator(calc_plus)
e_H2_plus = H2_plus.get_potential_energy()
niter_H2_plus = calc.get_number_of_iterations()

out = 'dospes.dat'
pes = DOSPES(calc, calc_plus, shift=True)
pes.save_folded_pes(filename=out, folding=None)
pes.save_folded_pes(filename=None, folding=None)

# check for correct shift
VDE = calc_plus.get_potential_energy() - calc.get_potential_energy()
BE_HOMO = 1.e23
be_n, f_n = pes.get_energies_and_weights()
for be, f in zip(be_n, f_n):
    if f > 0.1 and be < BE_HOMO:
        BE_HOMO = be
equal(BE_HOMO, VDE)

lr = LrTDDFT(calc_plus, xc=xc)

out = 'lrpes.dat'
pes = TDDFTPES(calc, lr)
pes.save_folded_pes(filename=out, folding='Gauss')
pes.save_folded_pes(filename=None, folding=None)

energy_tolerance = 0.0001
niter_tolerance = 1
equal(e_H2, -3.90059, energy_tolerance)
equal(e_H2_plus, 10.5659703, energy_tolerance)

# io
out = 'lrpes.dat.gz'
lr.write(out)
lr = LrTDDFT(out)
lr.set_calculator(calc_plus)
pes = TDDFTPES(calc, lr)
pes.save_folded_pes(filename=None, folding=None)


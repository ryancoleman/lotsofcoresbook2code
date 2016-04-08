from ase.lattice.spacegroup import crystal
from gpaw import GPAW
from gpaw import PW
from gpaw.test import equal

name = 'quartz'
# no. 152 - trigonal

a = 5.032090
c = a * 1.0968533
p0 = (0.4778763, 0.0, 1. / 3.)
p1 = (0.4153076, 0.2531340, 0.2029893)

atoms = crystal(['Si', 'O'], basis=[p0, p1],
                spacegroup=152, cellpar=[a, a, c, 90, 90, 120])


## with fractional translations
calc = GPAW(mode=PW(),
            xc='LDA',
            kpts=(3, 3, 3),
            nbands=42,
            symmetry={'symmorphic': False},
            gpts=(20, 20, 24),
            eigensolver='rmm-diis')

atoms.set_calculator(calc)
energy_fractrans = atoms.get_potential_energy()

assert(len(calc.wfs.kd.ibzk_kc) == 7)
assert(len(calc.wfs.kd.symmetry.op_scc) == 6)

## without fractional translations
calc = GPAW(mode=PW(),
            xc='LDA',
            kpts=(3, 3, 3),
            nbands=42,
            gpts=(20, 20, 24),
            eigensolver='rmm-diis')

atoms.set_calculator(calc)
energy_no_fractrans = atoms.get_potential_energy()

assert(len(calc.wfs.kd.ibzk_kc) == 10)
assert(len(calc.wfs.kd.symmetry.op_scc) == 2)

equal(energy_fractrans, energy_no_fractrans, 1e-3)

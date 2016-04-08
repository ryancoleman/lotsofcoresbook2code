from ase.lattice.spacegroup import crystal
from gpaw import GPAW
from gpaw import PW
from gpaw.test import equal

name = 'cristobalite'
# no. 92 - tetragonal

a = 5.0833674
c = 7.0984738
p0 = (0.2939118, 0.2939118, 0.0)
p1 = (0.2412656, 0.0931314, 0.1739217)

atoms = crystal(['Si', 'O'], basis=[p0, p1],
                spacegroup=92, cellpar=[a, a, c, 90, 90, 90])


## with fractional translations
calc = GPAW(mode=PW(),
            xc='LDA',
            kpts=(3, 3, 2),
            nbands=40,
            symmetry={'symmorphic': False},
            gpts=(24, 24, 32),
            eigensolver='rmm-diis')

atoms.set_calculator(calc)
energy_fractrans = atoms.get_potential_energy()

assert(len(calc.wfs.kd.ibzk_kc) == 3)
assert(len(calc.wfs.kd.symmetry.op_scc) == 8)

## without fractional translations
calc = GPAW(mode=PW(),
            xc='LDA',
            kpts=(3, 3, 2),
            nbands=40,
            gpts=(24, 24, 32),
            eigensolver='rmm-diis')

atoms.set_calculator(calc)
energy_no_fractrans = atoms.get_potential_energy()

assert(len(calc.wfs.kd.ibzk_kc) == 6)
assert(len(calc.wfs.kd.symmetry.op_scc) == 2)

equal(energy_fractrans, energy_no_fractrans, 1e-7)

from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.calculators.emt import EMT
from ase.lattice.surface import fcc100, add_adsorbate

from gpaw import GPAW, PW


def aual100(site, height, calc=None):
    slab = fcc100('Al', size=(2, 2, 2))
    slab.center(axis=2, vacuum=3.0)
    add_adsorbate(slab, 'Au', height, site)
    mask = [atom.symbol == 'Al' for atom in slab]
    fixlayer = FixAtoms(mask=mask)
    slab.set_constraint(fixlayer)

    if calc is None:
        calc = GPAW(mode=PW(200), kpts=(2, 2, 1), xc='PBE', txt=site + '.txt',
                    eigensolver='rmm-diis', nbands=40)

    slab.set_calculator(calc)
    qn = QuasiNewton(slab, trajectory=site + calc.name + '.traj')
    qn.run(fmax=0.05)

    if isinstance(calc, GPAW):
        calc.write(site + '.gpw')

    return slab.get_potential_energy()

e_hollow = aual100('hollow', 1.6)
e_bridge = aual100('bridge', 2.0)
e_ontop = aual100('ontop', 2.4)
assert abs(e_bridge - e_hollow - 0.352) < 0.01
assert abs(e_ontop - e_hollow - 0.711) < 0.01

calc = EMT()
e_hollow = aual100('hollow', 1.6, calc)
e_bridge = aual100('bridge', 2.0, calc)
e_ontop = aual100('ontop', 2.4, calc)
assert abs(e_bridge - e_hollow - 0.401) < 0.01
assert abs(e_ontop - e_hollow - 0.745) < 0.01

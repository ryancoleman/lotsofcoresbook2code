from math import sqrt

from ase import Atoms

from gpaw import GPAW
from gpaw.test import equal
from gpaw.xc.vdw import VDWFunctional

energy_tolerance = 0.002


def test():
    vdw = VDWFunctional('vdW-DF', verbose=1)
    d = 3.9
    x = d / sqrt(3)
    L = 3.0 + 2 * 4.0
    dimer = Atoms('Ar2', [(0, 0, 0), (x, x, x)], cell=(L, L, L))
    dimer.center()
    calc = GPAW(h=0.2, xc='revPBE')
    dimer.set_calculator(calc)
    e2 = dimer.get_potential_energy()
    calc.write('Ar2.gpw')
    e2vdw = calc.get_xc_difference(vdw)
    e2vdwb = GPAW('Ar2.gpw').get_xc_difference(vdw)
    print(e2vdwb - e2vdw)
    assert abs(e2vdwb - e2vdw) < 1e-9
    del dimer[1]
    e = dimer.get_potential_energy()
    evdw = calc.get_xc_difference(vdw)

    E = 2 * e - e2
    Evdw = E + 2 * evdw - e2vdw
    print(E, Evdw)
    assert abs(E - -0.0048) < 6e-3, abs(E)
    assert abs(Evdw - +0.0223) < 3e-2, abs(Evdw)

    print(e2, e)
    equal(e2, -0.001581923, energy_tolerance)
    equal(e, -0.003224226, energy_tolerance)

test()
